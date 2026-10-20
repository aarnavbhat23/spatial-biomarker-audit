"""Joint spatial certificates with shared observations and competing proposals."""
import json
from pathlib import Path
import numpy as np
from study_batch_acquisition import physical_index
from study_cell_reuse import resolved_edges
from study_raster_acquisition import prepare_frame_pairs
from study_residual_acquisition import residual_interval

ROOT=Path(__file__).resolve().parent


def envelope_distribution(remaining,totals,active):
    envelope=np.max(remaining[active]/totals[active,None],axis=0)
    mass=float(envelope.sum())
    if mass<=0:raise ValueError('No unresolved active claims')
    return envelope/mass,mass


def choose_distribution(remaining,totals,active,widths,policy):
    update=active.copy()
    if policy=='envelope':p,_=envelope_distribution(remaining,totals,active)
    elif policy=='mixture':
        rows=remaining[active];p=np.mean(rows/rows.sum(axis=1,keepdims=True),axis=0)
    elif policy=='uniform':
        p=np.any(remaining[active]>0,axis=0).astype(float);p/=p.sum()
    elif policy in ('focused','focused_reuse'):
        candidates=np.flatnonzero(active);chosen=candidates[np.argmax(widths[candidates])]
        p=remaining[chosen]/remaining[chosen].sum();update[:]=False;update[chosen]=True
        if policy=='focused_reuse':update=active&~np.any((remaining>0)&(p[None,:]==0),axis=1)
    else:raise ValueError(policy)
    caps=np.zeros(len(totals))
    for j in np.flatnonzero(update):
        if np.any((remaining[j]>0)&(p==0)):raise ValueError('Proposal misses unresolved claim support')
        caps[j]=np.max(remaining[j,p>0]/(totals[j]*p[p>0]))
    return p,update,caps


def claim_masks(data,radii=(40,60,80)):
    difference=data['coords'][data['edges'][:,0]]-data['coords'][data['edges'][:,1]]
    squared=np.sum(difference**2,axis=1)
    return np.stack([squared<=radius**2+1e-9 for radius in radii])


def multiclaim_trial(data,prepared,masks,field,seed,policy,width=.005,batch=16,max_queries=4096):
    rng=np.random.default_rng(seed);f=prepared['frame'];edges=data['edges'];labels=data['labels']
    outcomes=labels[edges[:,0],0]&labels[edges[:,1],1];totals=masks.sum(axis=1)
    if np.any(totals==0):raise ValueError('Each claim must contain edges')
    truth=np.sum(masks&outcomes[None,:],axis=1)/totals;q=len(totals);groups=len(f['counts'])
    remaining=np.stack([np.bincount(f['owners'][mask],minlength=groups) for mask in masks])
    _,field_nodes,group_fields=physical_index(data,prepared,field)
    group_ids=np.concatenate([np.full(len(ids),i) for i,ids in enumerate(group_fields)])
    field_ids=np.concatenate(group_fields)
    seen=np.zeros(len(labels),bool);acquired=np.zeros(len(field_nodes),bool)
    unresolved=np.ones(len(edges),bool);known=np.zeros(q,dtype=int);done=np.zeros(q,bool)
    values=[[] for _ in range(q)];lower=[[] for _ in range(q)];upper=[[] for _ in range(q)]
    intervals=np.tile([0.,1.],(q,1));queries_used=0;rounds=0;checkpoints=[];covered_all=True;empty_any=False
    while not done.all() and (policy=='greedy' or queries_used<max_queries):
        active=~done
        # Deterministically complete claims need no further randomized updates.
        empty=active&(remaining.sum(axis=1)==0)
        intervals[empty]=np.stack([known[empty]/totals[empty]]*2,axis=1);done[empty]=True
        if done.all():break
        active=~done
        old_known=known.copy()
        if policy=='greedy':
            group_score=np.sum(remaining[active]/totals[active,None],axis=0)
            score=np.bincount(field_ids,weights=group_score[group_ids],minlength=len(field_nodes));score[acquired]=-1
            fields=np.argsort(-score,kind='stable')[:min(batch,int((score>0).sum()))]
            selected=None
        else:
            p,update,caps=choose_distribution(remaining,totals,active,np.diff(intervals,axis=1).ravel(),policy)
            selected=rng.choice(groups,size=min(batch,max_queries-queries_used),p=p)
            frozen=[owned[unresolved[owned]] for owned in (prepared['owned'][g] for g in selected)]
            fields=np.unique(np.concatenate([group_fields[g] for g in selected]))
        new_fields=fields[~acquired[fields]];newly_seen=np.concatenate([field_nodes[i] for i in new_fields])
        acquired[fields]=True;seen[newly_seen]=True;rounds+=1
        if selected is not None:
            if not all(np.all(seen[edges[ids]]) for ids in frozen):raise RuntimeError('Unseen query endpoints')
            for j in np.flatnonzero(update):
                positives=np.array([np.sum(outcomes[ids]&masks[j,ids]) for ids in frozen])
                values[j].extend((old_known[j]/totals[j]+positives/(totals[j]*p[selected])).tolist())
                lower[j].extend([old_known[j]/totals[j]]*len(selected))
                upper[j].extend([old_known[j]/totals[j]+caps[j]]*len(selected))
            queries_used+=len(selected)
        affected=np.unique(np.concatenate([prepared['incident'][i] for i in newly_seen]))
        affected=affected[unresolved[affected]]
        identified=affected[resolved_edges(edges[affected],labels,seen)];unresolved[identified]=False
        for j in range(q):
            ids=identified[masks[j,identified]];known[j]+=int(outcomes[ids].sum())
            remaining[j]-=np.bincount(f['owners'][ids],minlength=groups)
            low=known[j]/totals[j];high=(known[j]+remaining[j].sum())/totals[j]
            if done[j]:intervals[j]=max(intervals[j,0],low),min(intervals[j,1],high)
            elif policy=='greedy':intervals[j]=low,high
            else:intervals[j]=residual_interval(values[j],lower[j],upper[j],low,high,alpha=.05/q)
        covered=bool(np.all((intervals[:,0]-1e-12<=truth)&(truth<=intervals[:,1]+1e-12)))
        covered_all&=covered;empty_any|=bool(np.any(intervals[:,0]>intervals[:,1]))
        lengths=np.diff(intervals,axis=1).ravel();done|=(lengths>=0)&(lengths<=width)
        checkpoints.append({'round':rounds,'queries':queries_used,'fields':int(acquired.sum()),
                            'intervals':intervals.tolist(),'joint_covered':covered,'done':done.tolist()})
    return {'fields':int(acquired.sum()),'rounds':rounds,'queries':queries_used,'stopped':bool(done.all()),
            'joint_covered_all_checkpoints':covered_all,'empty_set_any_checkpoint':empty_any,
            'intervals':intervals.tolist(),'truth':truth.tolist(),'checkpoints':checkpoints}


def main():
    report={'scope':'Exploratory shared nested distance claims, fixed graph and exact calls; no new patients.',
            'seed':20261109,'radii_native_pixels':[40,60,80],'joint_alpha':.05,'target_width':.005,
            'field_native_pixels':1024,'batch':16,'patients':{}}
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));prepared=prepare_frame_pairs(data,1024)
        masks=claim_masks(data);entries=[]
        for policy in ('envelope','mixture','uniform','focused','focused_reuse','greedy'):
            rows=[multiclaim_trial(data,prepared,masks,1024,20261109+1000*int(patient[3:])+repeat,policy)
                  for repeat in range(1 if policy=='greedy' else 5)]
            entry={'policy':policy,'rows':rows,'mean_fields':float(np.mean([r['fields'] for r in rows])),
                   'mean_rounds':float(np.mean([r['rounds'] for r in rows])),
                   'attainment':float(np.mean([r['stopped'] for r in rows])),
                   'joint_coverage_all_checkpoints':float(np.mean([r['joint_covered_all_checkpoints'] for r in rows]))}
            entries.append(entry);print(patient,policy,'fields',entry['mean_fields'],'rounds',entry['mean_rounds'],
                                        'attain',entry['attainment'],'coverage',entry['joint_coverage_all_checkpoints'],flush=True)
        report['patients'][patient]=entries
        (ROOT/'results/multiclaim_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
