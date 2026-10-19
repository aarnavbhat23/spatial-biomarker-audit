"""Physical cell visibility and logical reuse for residual pair sampling."""
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from shapely import box, union_all

from study_region_scale import frame
from study_residual_acquisition import nested_frames, residual_interval, residual_trial

ROOT=Path(__file__).resolve().parent


def grouped_indices(owner, number):
    order=np.argsort(owner,kind='stable')
    offsets=np.r_[0,np.cumsum(np.bincount(owner,minlength=number))]
    return [order[offsets[i]:offsets[i+1]] for i in range(number)]


def prepare_visibility(data,core):
    f=frame(data,core)
    edges=data['edges'];coords=data['coords']
    owned=grouped_indices(f['owners'],len(f['counts']))
    # Each node's incident edges; each edge is visited at most twice per trial.
    node_order=np.argsort(np.r_[edges[:,0],edges[:,1]],kind='stable')
    edge_ids=np.tile(np.arange(len(edges)),2)[node_order]
    offsets=np.r_[0,np.cumsum(np.bincount(edges.ravel(),minlength=len(coords)))]
    incident=[edge_ids[offsets[i]:offsets[i+1]] for i in range(len(coords))]
    tree=cKDTree(coords)
    visible=[]
    for rectangle in f['rectangles']:
        center=(rectangle[:2]+rectangle[2:])/2
        ids=np.asarray(tree.query_ball_point(center,(core+160)/2,p=np.inf),dtype=int)
        inside=np.all((coords[ids]>=rectangle[:2])&(coords[ids]<rectangle[2:]),axis=1)
        visible.append(ids[inside])
    return {'frame':f,'owned':owned,'incident':incident,'visible':visible}


def resolved_edges(edges,labels,seen):
    a,b=edges.T
    return ((seen[a]&~labels[a,0]) | (seen[b]&~labels[b,1]) | (seen[a]&seen[b]))


def cell_trial(data,prepared,seed,width,max_draws=2048):
    rng=np.random.default_rng(seed);f=prepared['frame']
    edges=data['edges'];labels=data['labels'];total=len(edges)
    outcomes=labels[edges[:,0],0]&labels[edges[:,1],1]
    truth=float(outcomes.mean())
    seen=np.zeros(len(labels),dtype=bool);unresolved=np.ones(total,dtype=bool)
    remaining=f['counts'].copy();known_positive=0
    values=[];lower=[];upper=[];rectangles=[];checkpoints=[]
    covered_all=True;stopped=False;lo,hi=0.,1.
    zero_shortcuts=0
    for step in range(1,max_draws+1):
        unknown=int(remaining.sum())
        if unknown==0:break
        selected=int(rng.choice(len(remaining),p=remaining/unknown))
        owned=prepared['owned'][selected]
        query=owned[unresolved[owned]]
        newly_seen=prepared['visible'][selected][~seen[prepared['visible'][selected]]]
        seen[newly_seen]=True
        # The selected halo must reveal both endpoints of every owned edge.
        if not np.all(seen[edges[query]]):
            raise RuntimeError('Acquisition did not cover selected edge endpoints')
        acquired_positive=int(outcomes[query].sum())
        lower.append(known_positive/total);upper.append((known_positive+unknown)/total)
        values.append(known_positive/total+unknown/total*acquired_positive/len(query))
        affected=np.unique(np.concatenate([prepared['incident'][i] for i in newly_seen]))
        affected=affected[unresolved[affected]]
        identified=affected[resolved_edges(edges[affected],labels,seen)]
        zero_shortcuts+=int(np.sum(~np.all(seen[edges[identified]],axis=1)))
        known_positive+=int(outcomes[identified].sum())
        unresolved[identified]=False
        remaining-=np.bincount(f['owners'][identified],minlength=len(remaining))
        if 'physical_rectangles' in prepared:
            rectangles.extend(prepared['physical_rectangles'][selected])
        else:
            rectangles.append(f['rectangles'][selected])
        current_low=known_positive/total
        current_high=(known_positive+int(remaining.sum()))/total
        if step in (8,16,32,64,128,256,512,1024,2048) or step==max_draws or current_high-current_low<=width:
            lo,hi=residual_interval(values,lower,upper,current_low,current_high)
            covered=lo-1e-12<=truth<=hi+1e-12
            covered_all=covered_all and covered
            checkpoints.append({'draws':step,'interval':[lo,hi],'covered':bool(covered)})
            if 0<=hi-lo<=width:
                stopped=True;break
    rects=np.asarray(rectangles)
    area=float(union_all(box(*rects.T)).area)
    summed_area=float(np.sum((rects[:,2]-rects[:,0])*(rects[:,3]-rects[:,1])))
    return {'draws':len(values),'interval':[lo,hi],'stopped':stopped,
            'covered_all_checkpoints':bool(covered_all),'empty_confidence_set':lo>hi,
            'area_native_pixels':area,'summed_area_native_pixels':summed_area,
            'observed_cells':int(seen.sum()),'resolved_edges':int((~unresolved).sum()),
            'zero_shortcut_edges':zero_shortcuts,'checkpoints':checkpoints}


def main():
    result={'seed':20261031,'scope':'Exploratory cell reuse on three development patients, 10 repeats; no fresh confirmation.',
            'patients':{}}
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        counts,positive,frames=nested_frames(data)
        entries=[]
        for core in (128,2048):
            prepared=prepare_visibility(data,core)
            for width in (.04,.005):
                for reuse in ('nested_core','cell_visibility'):
                    rows=[]
                    for repeat in range(10):
                        seed=20261031+1000*int(patient[3:])+repeat
                        row=(cell_trial(data,prepared,seed,width) if reuse=='cell_visibility' else
                             residual_trial(counts,positive,frames,seed,str(core),width,0))
                        rows.append(row)
                    entry={'core':core,'target_width':width,'reuse':reuse,'rows':rows,
                           'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                           'coverage_all_checkpoints':float(np.mean([r['covered_all_checkpoints'] for r in rows])),
                           'mean_draws':float(np.mean([r['draws'] for r in rows])),
                           'mean_area_native_pixels':float(np.mean([r['area_native_pixels'] for r in rows])),
                           'mean_summed_area_native_pixels':float(np.mean([r['summed_area_native_pixels'] for r in rows]))}
                    entries.append(entry)
                    print(patient,core,width,reuse,'area',round(entry['mean_area_native_pixels']),
                          'draws',entry['mean_draws'],'attain',entry['stop_fraction'],flush=True)
        result['patients'][patient]=entries
        (ROOT/'results/cell_reuse_study.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
