"""Classical count-conditioned MILP reference for an out-degree-regular graph."""
import time
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
from scipy.sparse import coo_array


def exact_numerator(edges,labels):
    """Integer arithmetic sign check, independent of floating interval calculations."""
    n=len(labels);m=len(edges);u,v=np.asarray(edges).T
    d=np.bincount(u,minlength=n);e=np.bincount(v,minlength=n)
    contacts=sum(int(labels[i,0])*int(labels[j,1]) for i,j in edges)
    left=sum(int(x)*int(y) for x,y in zip(d,labels[:,0]))
    right=sum(int(x)*int(y) for x,y in zip(e,labels[:,1]))
    return m*contacts-left*right,m*m


def maximize_by_counts(edges,labels,ba,bb,time_limit=1.):
    labels=np.asarray(labels,int);edges=np.asarray(edges,int);n=len(labels);m=len(edges);u,v=edges.T
    d=np.bincount(u,minlength=n);e=np.bincount(v,minlength=n)
    if not np.all(d==d[0]) or not m:raise ValueError('Need nonempty constant-outdegree graph')
    if any(int(b)!=b or b<0 for b in (ba,bb)):raise ValueError('Need nonnegative integer budgets')
    ba,bb=int(ba),int(bb)
    rows=[];cols=[];data=[];lo=[];hi=[]
    def add(ids,coeff,lower,upper):
        index=len(lo);rows.extend([index]*len(ids));cols.extend(ids);data.extend(coeff);lo.append(lower);hi.append(upper)
    for j,(i,k) in enumerate(edges):
        z=2*n+j
        add([z,int(i)],[1,-1],-np.inf,0);add([z,n+int(k)],[1,-1],-np.inf,0)
        add([z,int(i),n+int(k)],[1,-1,-1],-1,np.inf)
    for marker,budget in [(0,ba),(1,bb)]:
        a=labels[:,marker]
        add((marker*n+np.arange(n)).tolist(),(1-2*a).tolist(),-np.inf,int(budget-a.sum()))
    add(list(range(n)),[1]*n,0,0)
    matrix=coo_array((np.asarray(data,float),(np.asarray(rows,np.int32),np.asarray(cols,np.int32))),shape=(len(lo),2*n+m)).tocsc()
    raw_count=int(labels[:,0].sum());answers=[]
    for count in range(max(0,raw_count-ba),min(n,raw_count+ba)+1):
        lower=np.array(lo,float);upper=np.array(hi,float);lower[-1]=upper[-1]=count
        objective=np.r_[np.zeros(n),count*e/(n*m),np.full(m,-1/m)]
        start=time.perf_counter()
        result=milp(objective,integrality=np.r_[np.ones(2*n),np.zeros(m)],
              bounds=Bounds(np.zeros(2*n+m),np.ones(2*n+m)),
              constraints=LinearConstraint(matrix,lower,upper),
              options={'time_limit':time_limit,'mip_rel_gap':0.,'presolve':True})
        row={'count_a':count,'status':int(result.status),'seconds':time.perf_counter()-start,
             'upper_bound':None,'witness':None}
        if result.status==2:row['infeasible']=True
        else:
            dual=getattr(result,'mip_dual_bound',None)
            if dual is not None and np.isfinite(dual):row['upper_bound']=min(.25,float(-dual)+1e-7)
            if result.x is not None:
                x=np.stack([result.x[:n],result.x[n:2*n]],axis=1).round().astype(int)
                if np.sum(x[:,0]!=labels[:,0])<=ba and np.sum(x[:,1]!=labels[:,1])<=bb and x[:,0].sum()==count:
                    numerator,denominator=exact_numerator(edges,x)
                    row['witness']={'numerator':numerator,'denominator':denominator,
                        'a_flips':np.flatnonzero(x[:,0]!=labels[:,0]).tolist(),
                        'b_flips':np.flatnonzero(x[:,1]!=labels[:,1]).tolist()}
                    if row['upper_bound'] is not None and numerator/denominator>row['upper_bound']:
                        raise AssertionError('Independent integer witness exceeds solver bound')
        answers.append(row)
    feasible=[r for r in answers if not r.get('infeasible')]
    missing=any(r['upper_bound'] is None for r in feasible)
    return {'slices':answers,'upper_bound':None if missing else max(r['upper_bound'] for r in feasible),
            'all_slices_optimal':all(r['status'] in (0,2) for r in answers),
            'missing_dual_slices':sum(r['upper_bound'] is None for r in feasible)}
