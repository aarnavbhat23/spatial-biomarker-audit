"""Feasible categorical witnesses and conditional pairwise review lower bounds."""
import numpy as np
from contrast_count_slices import exact_numerator


def categorical_witness(edges,raw,budget_a,budget_b,maximize,start=None,max_moves=50):
    raw=np.asarray(raw);edges=np.asarray(edges)
    if raw.ndim!=2 or raw.shape[1]!=2 or not np.isin(raw,[0,1]).all() or np.any(raw.sum(axis=1)>1):
        raise ValueError('Need mutually exclusive binary identities')
    if edges.ndim!=2 or edges.shape[1]!=2 or not np.issubdtype(edges.dtype,np.integer):raise ValueError('Need directed integer edges')
    n=len(raw)
    if not len(edges) or edges.min()<0 or edges.max()>=n or np.any(edges[:,0]==edges[:,1]):raise ValueError('Need nonempty graph without self loops')
    if any(int(b)!=b or b<0 for b in (budget_a,budget_b)):raise ValueError('Need nonnegative integer budgets')
    raw=raw.astype(np.int64);x=raw.copy() if start is None else np.asarray(start,dtype=np.int64).copy()
    if x.shape!=raw.shape or not np.isin(x,[0,1]).all() or np.any(x.sum(axis=1)>1):raise ValueError('Invalid start')
    budget=np.array([budget_a,budget_b],np.int64);used=np.sum(x!=raw,axis=0)
    if np.any(used>budget):raise ValueError('Start exceeds budgets')
    u,v=edges.T;m=len(edges);d=np.bincount(u,minlength=n).astype(np.int64);e=np.bincount(v,minlength=n).astype(np.int64)
    states=np.array([[0,0],[1,0],[0,1]],np.int64);sigma=1 if maximize else -1;moves=0
    for _ in range(max_moves):
        out=np.bincount(u,weights=x[v,1],minlength=n).astype(np.int64)
        inc=np.bincount(v,weights=x[u,0],minlength=n).astype(np.int64)
        left=int(d@x[:,0]);right=int(e@x[:,1]);best_gain=0;best=None
        for state in states:
            change=state-x;cost=used+(state!=raw)-(x!=raw)
            feasible=np.all(cost<=budget,axis=1)
            da,db=change.T;dl=d*da;dr=e*db
            # No self loops: contact delta contains no same-cell da*db term.
            gain=sigma*(m*(da*out+db*inc)-dl*right-left*dr-dl*dr)
            gain[~feasible]=0;index=int(np.argmax(gain))
            if int(gain[index])>best_gain:best_gain=int(gain[index]);best=(index,state.copy(),cost[index].copy())
        if best is None:break
        index,state,used=best;x[index]=state;moves+=1
    numerator,denominator=exact_numerator(edges,x)
    return {'numerator':numerator,'denominator':denominator,'moves':moves,
            'a_flips':np.flatnonzero(x[:,0]!=raw[:,0]).tolist(),'b_flips':np.flatnonzero(x[:,1]!=raw[:,1]).tolist()}


def repair_overlap(raw,labels):
    """Drop added overlaps, preserving raw identities and never increasing either budget."""
    raw=np.asarray(raw,int);x=np.asarray(labels,int).copy();overlap=x.sum(axis=1)>1
    x[overlap & (raw[:,1]==0),1]=0
    x[overlap & (raw[:,1]==1),0]=0
    if np.any(np.sum(x!=raw,axis=0)>np.sum(labels!=raw,axis=0)):raise AssertionError('Repair increased error budget')
    return x


def review_pair_bound(x,y,fields,error_cells):
    """Exact cheapest field set for one pair, a necessary condition for global sign certification."""
    x=np.asarray(x);y=np.asarray(y);fields=np.asarray(fields)
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=2 or fields.shape!=(len(x),):raise ValueError('Need aligned labels and fields')
    if not np.isin(x,[0,1]).all() or not np.isin(y,[0,1]).all() or np.any(x.sum(axis=1)>1) or np.any(y.sum(axis=1)>1):raise ValueError('Need categorical identities')
    if int(error_cells)!=error_cells or error_cells<0:raise ValueError('Need nonnegative integer review error budget')
    different=np.any(x!=y,axis=1);quota=2*int(error_cells)+1
    keys,counts=np.unique(fields[different],return_counts=True);order=np.lexsort((keys,-counts))
    keys=keys[order];counts=counts[order];cumulative=np.cumsum(counts)
    possible=int(different.sum())>=quota
    number=int(np.searchsorted(cumulative,quota)+1) if possible else None
    return {'differing_cells':int(different.sum()),'error_cells':int(error_cells),'required_distinct_disagreements':quota,
            'distinguishable_even_after_all_cells':possible,'minimum_unit_cost_fields_for_pair':number,
            'cheapest_pair_fields':keys[:number].tolist() if possible else [],
            'disagreements_by_field':dict(zip(keys.astype(str).tolist(),counts.tolist()))}
