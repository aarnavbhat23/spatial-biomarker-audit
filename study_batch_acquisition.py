"""Batch-boundary betting and a deterministic graph-certificate control."""
import json
from pathlib import Path

import numpy as np

from study_cell_reuse import grouped_indices,resolved_edges
from study_raster_acquisition import prepare_frame_pairs
from study_residual_acquisition import residual_interval

ROOT=Path(__file__).resolve().parent


def physical_index(data,prepared,field):
    tiles,node_field=np.unique(data['coords']//field,axis=0,return_inverse=True)
    lookup={tuple(tile):i for i,tile in enumerate(tiles)}
    group_fields=[np.array([lookup[tuple(rect[:2]/field)] for rect in rectangles])
                  for rectangles in prepared['physical_rectangles']]
    return node_field,grouped_indices(node_field,len(tiles)),group_fields


def batch_values(outcomes,queries,known_positive,unknown,total):
    """Queries are frozen pre-batch; duplicate draws must remain duplicate values."""
    return [known_positive/total+unknown/total*float(outcomes[q].mean()) for q in queries]


def batch_trial(data,prepared,field,seed,batch_size,width=.005,max_draws=2048,greedy_batch=1):
    """batch_size=None selects deterministic greedy acquisition, not a sampler."""
    if batch_size is not None and batch_size<1:raise ValueError('Positive batch size required')
    if greedy_batch<1:raise ValueError('Positive greedy batch size required')
    rng=np.random.default_rng(seed);f=prepared['frame']
    edges=data['edges'];labels=data['labels'];total=len(edges)
    outcomes=labels[edges[:,0],0]&labels[edges[:,1],1]
    truth=float(outcomes.mean())
    _,field_nodes,group_fields=physical_index(data,prepared,field)
    group_ids=np.concatenate([np.full(len(ids),i) for i,ids in enumerate(group_fields)])
    field_ids=np.concatenate(group_fields)
    seen=np.zeros(len(labels),bool);acquired=np.zeros(len(field_nodes),bool)
    unresolved=np.ones(total,bool);remaining=f['counts'].copy();known_positive=0
    values=[];lower=[];upper=[];checkpoints=[];rounds=0;stopped=False;covered_all=True
    lo,hi=0.,1.
    while remaining.sum()>0 and (batch_size is None or len(values)<max_draws):
        unknown=int(remaining.sum());old_positive=known_positive
        if batch_size is None:
            scores=np.bincount(field_ids,weights=remaining[group_ids],minlength=len(field_nodes))
            scores[acquired]=-1
            fields=np.argsort(-scores,kind='stable')[:min(greedy_batch,int((scores>0).sum()))]
            queries=[]
        else:
            selected=rng.choice(len(remaining),size=min(batch_size,max_draws-len(values)),
                                replace=True,p=remaining/unknown)
            queries=[owned[unresolved[owned]] for owned in (prepared['owned'][i] for i in selected)]
            fields=np.unique(np.concatenate([group_fields[i] for i in selected]))
        new_fields=fields[~acquired[fields]]
        newly_seen=np.concatenate([field_nodes[i] for i in new_fields])
        acquired[fields]=True;seen[newly_seen]=True;rounds+=1
        if batch_size is not None:
            if not all(np.all(seen[edges[q]]) for q in queries):
                raise RuntimeError('Missing selected edge endpoints')
            values.extend(batch_values(outcomes,queries,old_positive,unknown,total))
            lower.extend([old_positive/total]*len(queries))
            upper.extend([(old_positive+unknown)/total]*len(queries))
        affected=np.unique(np.concatenate([prepared['incident'][i] for i in newly_seen]))
        affected=affected[unresolved[affected]]
        identified=affected[resolved_edges(edges[affected],labels,seen)]
        known_positive+=int(outcomes[identified].sum());unresolved[identified]=False
        remaining-=np.bincount(f['owners'][identified],minlength=len(remaining))
        lo,hi=known_positive/total,(known_positive+int(remaining.sum()))/total
        if batch_size is not None:lo,hi=residual_interval(values,lower,upper,lo,hi)
        covered=bool(lo-1e-12<=truth<=hi+1e-12);covered_all&=covered
        checkpoints.append({'round':rounds,'queries':len(values),'fields':int(acquired.sum()),
                            'interval':[lo,hi],'covered':covered})
        if 0<=hi-lo<=width:stopped=True;break
    return {'physical_fields':int(acquired.sum()),'rounds':rounds,'queries':len(values),
            'stopped':stopped,'covered_all_checkpoints':covered_all,'empty_confidence_set':lo>hi,
            'interval':[lo,hi],'checkpoints':checkpoints}


def main():
    report={'scope':'Development-only batch and deterministic certificate screen; hypothetical field costs.',
            'seed':20261105,'target_width':.005,'repetitions':5,'patients':{}}
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));entries=[]
        for field in (512,1024,2048):
            prepared=prepare_frame_pairs(data,field)
            methods=[(b,1) for b in (1,4,16,64)]+[(None,b) for b in (1,4,16,64)]
            for batch,greedy_batch in methods:
                repeats=1 if batch is None else 5
                rows=[batch_trial(data,prepared,field,20261105+1000*int(patient[3:])+repeat,batch,greedy_batch=greedy_batch)
                      for repeat in range(repeats)]
                entry={'field_native_pixels':field,'batch_size':batch,'greedy_batch':greedy_batch if batch is None else None,
                       'policy':'greedy_certificate' if batch is None else 'batch_pps','rows':rows,
                       'mean_fields':float(np.mean([r['physical_fields'] for r in rows])),
                       'mean_rounds':float(np.mean([r['rounds'] for r in rows])),
                       'mean_queries':float(np.mean([r['queries'] for r in rows])),
                       'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                       'coverage_all_checkpoints':float(np.mean([r['covered_all_checkpoints'] for r in rows]))}
                entry['cost_scenarios']={str(h):entry['mean_fields']+h*entry['mean_rounds'] for h in (0,1,10,100)}
                entries.append(entry)
                print(patient,field,batch,greedy_batch,'fields',entry['mean_fields'],'rounds',entry['mean_rounds'],
                      'attain',entry['stop_fraction'],flush=True)
        report['patients'][patient]=entries
        (ROOT/'results/batch_acquisition_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
