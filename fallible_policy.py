"""Field review policies see only recorded calls and queried fallible responses."""
import numpy as np
from fallible_bounds import linear_optimum,review_bounds
from general_marker_bounds import fast_bounds


def fields_from_coordinates(coordinates):
    xy=np.asarray(coordinates,float);span=np.ptp(xy,axis=0)
    bins=np.minimum(3,((xy-xy.min(axis=0))/np.maximum(span,1e-12)*4).astype(int))
    return bins[:,0]+4*bins[:,1]


def endpoint_labels(raw,bound):
    worlds=[]
    for detail in bound['details']:
        world=raw.copy();world[detail['a_flips'],0]^=1;world[detail['b_flips'],1]^=1;worlds.append(world)
    return worlds


def choose_field(edges,raw,reference,budgets,errors,fields,bound,history,policy,rng):
    n=len(raw);u,v=edges.T;pending=reference[:,0]<0;available=np.unique(fields[pending])
    if not len(available):return None
    if policy=='uniform':return int(rng.choice(available))
    degree=np.bincount(u,minlength=n)+np.bincount(v,minlength=n)
    influence=np.bincount(u,weights=raw[v,1],minlength=n)+np.bincount(v,weights=raw[u,0],minlength=n)
    if policy in ('lookahead','witness'):
        low,high=endpoint_labels(raw,bound)
        changed=low[u,0]*low[v,1]!=high[u,0]*high[v,1]
        disagreement=np.bincount(u,weights=changed*(low[u,0]!=high[u,0]),minlength=n)
        disagreement+=np.bincount(v,weights=changed*(low[v,1]!=high[v,1]),minlength=n)
    if policy=='lookahead':
        anchor=np.stack([linear_optimum(np.zeros(n),raw[:,k],reference[:,k],budgets[k],errors[k])[1] for k in (0,1)],axis=1)
    candidates=[];current=history[1]-history[0]
    for field in available:
        mask=pending&(fields==field);cost=1+.02*mask.sum()
        if policy=='degree':score=degree[mask].sum()/cost
        elif policy=='influence':score=influence[mask].sum()/cost
        elif policy=='witness':score=disagreement[mask].sum()/cost
        elif policy=='lookahead':
            widths=[]
            for proxy in (anchor,low,high):
                future=reference.copy();future[mask]=proxy[mask]
                result=review_bounds(edges,raw,future,budgets,errors,witnesses=False)
                if result['empty']:raise AssertionError('Feasible proxy produced an empty set')
                lo=max(history[0],result['outer_bounds'][0]);hi=min(history[1],result['outer_bounds'][1])
                widths.append(hi-lo)
            score=(current-max(widths))/cost
        else:raise ValueError('Unknown policy')
        candidates.append((float(score),-int(field)))
    return -max(candidates)[1]


def run_policy(edges,raw,budgets,errors,fields,review,policy,seed,target_width=.005):
    reference=np.full_like(raw,-1);rng=np.random.default_rng(seed);trace=[];cells=0;rounds=0
    history=fast_bounds(edges,raw,*budgets)['outer_bounds']
    while True:
        bound=review_bounds(edges,raw,reference,budgets,errors)
        if bound['empty']:return {'empty':True,'trace':trace}
        history=[max(history[0],bound['outer_bounds'][0]),min(history[1],bound['outer_bounds'][1])]
        step={'round':rounds,'reviewed_cells':cells,'outer_bounds':history.copy(),'attained_bounds':bound['attained_bounds'],
              'endpoint_gaps':[bound['attained_bounds'][0]-history[0],history[1]-bound['attained_bounds'][1]]}
        trace.append(step)
        if history[1]-history[0]<=target_width+1e-12 or np.all(reference>=0):break
        field=choose_field(edges,raw,reference,budgets,errors,fields,bound,history,policy,rng)
        ids=np.flatnonzero((fields==field)&(reference[:,0]<0));response=np.asarray(review(ids),int)
        if response.shape!=(len(ids),2) or not np.isin(response,[0,1]).all():raise ValueError('Malformed reference response')
        reference[ids]=response;cells+=len(ids);rounds+=1;step['next_field']=field;step['new_cells']=len(ids)
    return {'empty':False,'attained':history[1]-history[0]<=target_width+1e-12,'fields':rounds,'cells':cells,'rounds':rounds,
            'quoted_cost':{str(latency):rounds*(1+latency)+.02*cells for latency in (0,1,8)},'trace':trace}
