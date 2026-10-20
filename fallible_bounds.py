"""Conditional bilinear bounds with two fallible annotation sources per marker."""
import numpy as np


def linear_optimum(weights,raw,reference,budget,errors,maximize=True):
    """Exact signed linear optimization over intersected partial Hamming balls."""
    raw=np.asarray(raw,int);reference=np.asarray(reference,int);weights=np.asarray(weights,float)
    if any(int(x)!=x or x<0 for x in (budget,errors)):raise ValueError('Nonnegative integer budgets required')
    if not np.isin(raw,[0,1]).all() or not np.isin(reference,[-1,0,1]).all():raise ValueError('Invalid calls')
    budget=min(int(budget),len(raw));errors=int(errors);known=reference>=0
    disagreement=known&(raw!=reference);agreement=known&~disagreement
    groups=[np.flatnonzero(mask) for mask in (disagreement,agreement,~known)]
    scores=weights*(1-2*raw)*(1 if maximize else -1)
    orders=[ids[np.argsort(-scores[ids],kind='stable')] for ids in groups]
    prefixes=[np.r_[0.,np.cumsum(scores[ids])] for ids in orders]
    best_u=np.zeros(len(prefixes[2]),int)
    for u in range(1,len(best_u)):
        best_u[u]=u if prefixes[2][u]>prefixes[2][best_u[u-1]] else best_u[u-1]
    best=-np.inf;counts=None;nd,ng,nu=map(len,groups)
    for d in range(min(nd,budget)+1):
        gmax=min(ng,budget-d,errors-nd+d)
        for g in range(gmax+1):
            u=int(best_u[min(nu,budget-d-g)])
            value=prefixes[0][d]+prefixes[1][g]+prefixes[2][u]
            if value>best:best=value;counts=(d,g,u)
    if counts is None:return None
    result=raw.copy()
    for ids,count in zip(orders,counts):result[ids[:count]]^=1
    return float(weights@result),result


def flippable(raw,reference,budget,errors,anchor):
    """Whether the label opposite the hybrid anchor is feasible at each cell."""
    known=reference>=0;disagreement=known&(raw!=reference)
    opposite=1-anchor;spent_raw=(opposite!=raw).astype(int)
    spent_ref=(known&(opposite!=reference)).astype(int)
    remaining_disagreements=int(disagreement.sum())-disagreement.astype(int)
    return ((budget-spent_raw>=0)&(errors-spent_ref>=0)&
            (remaining_disagreements<=budget-spent_raw+errors-spent_ref))


def review_bounds(edges,raw,reference,budgets,errors,witnesses=True):
    n=len(raw);u,v=np.asarray(edges).T;m=len(edges)
    raw=np.asarray(raw,int);reference=np.asarray(reference,int)
    anchor=np.where(reference>=0,reference,raw);oracles=[];feasible=[]
    for k in (0,1):
        def oracle(weights,maximize=True,k=k):
            return linear_optimum(weights,raw[:,k],reference[:,k],budgets[k],errors[k],maximize)
        initial=oracle(np.zeros(n))
        if initial is None:return {'empty':True}
        oracles.append(oracle);feasible.append(initial[1])
    wa=np.bincount(u,weights=anchor[v,1],minlength=n)
    wb=np.bincount(v,weights=anchor[u,0],minlength=n)
    base=int(np.sum(anchor[u,0]*anchor[v,1]))
    changes=[];possible=[];linear=[]
    for k,w in enumerate((wa,wb)):
        hybrid=anchor[:,k]
        changes.append(int(round(hybrid.sum()+oracles[k](1-2*hybrid)[0])))
        possible.append(flippable(raw[:,k],reference[:,k],budgets[k],errors[k],hybrid))
        linear.append([oracles[k](w,maximize=False)[0]-w@hybrid,oracles[k](w)[0]-w@hybrid])
    cross=[]
    for positive in (False,True):
        eligible=possible[0][u]&possible[1][v]&((anchor[u,0]==anchor[v,1])==positive)
        da=np.bincount(u[eligible],minlength=n);db=np.bincount(v[eligible],minlength=n)
        limits=[changes[0]*changes[1]]
        for k,degree in enumerate((da,db)):
            hybrid=anchor[:,k];limits.append(float(degree@hybrid+oracles[k](degree*(1-2*hybrid))[0]))
        cross.append(max(0.,min(limits)))
    outer=[max(0.,(base+linear[0][0]+linear[1][0]-cross[0])/m),
           min(1.,(base+linear[0][1]+linear[1][1]+cross[1])/m)]
    answer={'empty':False,'outer_bounds':outer,'cross_caps':cross,'max_changed_cells':changes}
    if not witnesses:return answer
    endpoints=[];details=[]
    for maximize in (False,True):
        candidates=[]
        for start in (0,1):
            a,b=[x.copy() for x in feasible]
            if start:a=oracles[0](np.ones(n),maximize=maximize)[1]
            for _ in range(12):
                old_a,old_b=a.copy(),b.copy()
                a=oracles[0](np.bincount(u,weights=b[v],minlength=n),maximize=maximize)[1]
                b=oracles[1](np.bincount(v,weights=a[u],minlength=n),maximize=maximize)[1]
                if np.array_equal(a,old_a) and np.array_equal(b,old_b):break
            value=float(np.mean(a[u]*b[v]));candidates.append((value,a,b))
        value,a,b=(max if maximize else min)(candidates,key=lambda x:x[0]);endpoints.append(value)
        details.append({'a_flips':np.flatnonzero(a!=raw[:,0]).tolist(),'b_flips':np.flatnonzero(b!=raw[:,1]).tolist()})
    if outer[0]>endpoints[0]+1e-12 or outer[1]<endpoints[1]-1e-12:raise AssertionError('Attained world outside bound')
    answer.update({'attained_bounds':endpoints,'endpoint_gaps':[endpoints[0]-outer[0],outer[1]-endpoints[1]],'details':details})
    return answer
