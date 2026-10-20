"""General marker budgets: valid graph bounds, attainable witnesses, classical MILP."""
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
from scipy.sparse import coo_array
from measurement_error import top_budget_sum,degree_error_margin,pair_target


def _best_mask(scores,eligible,budget):
    ids=np.flatnonzero(eligible);count=min(budget,len(ids));result=np.zeros(len(scores),bool)
    if count:result[ids[np.argsort(-scores[ids],kind='stable')[:count]]]=True
    return result


def fast_bounds(edges,labels,budget_a,budget_b):
    """Simple directed unit-weight graph, census labels, externally bounded errors."""
    if any(int(b)!=b or b<0 for b in (budget_a,budget_b)):raise ValueError('Nonnegative integer budgets required')
    budget_a,budget_b=int(budget_a),int(budget_b)
    n=len(labels);u,v=edges.T;a=labels[:,0].astype(bool);b=labels[:,1].astype(bool)
    w=np.bincount(u,weights=b[v],minlength=n);z=np.bincount(v,weights=a[u],minlength=n)
    numerator=int(np.sum(a[u]&b[v]));details=[]
    for increase in (False,True):
        ea=~a if increase else a;eb=~b if increase else b
        cross=ea[u]&eb[v];cu=u[cross];cv=v[cross]
        da=np.bincount(cu,minlength=n);db=np.bincount(cv,minlength=n)
        separate=top_budget_sum(w[ea],budget_a)+top_budget_sum(z[eb],budget_b)
        if increase:
            upper=min(top_budget_sum((w+np.minimum(da,budget_b))[ea],budget_a)+top_budget_sum(z[eb],budget_b),
                      top_budget_sum(w[ea],budget_a)+top_budget_sum((z+np.minimum(db,budget_a))[eb],budget_b),
                      separate+min(top_budget_sum(da[ea],budget_a),top_budget_sum(db[eb],budget_b),budget_a*budget_b))
        else:upper=min(numerator,separate)
        best=-1.;witness=None;sign=1 if increase else -1
        for start in ('a','b'):
            x=_best_mask(w,ea,budget_a) if start=='a' else np.zeros(n,bool)
            y=_best_mask(z,eb,budget_b) if start=='b' else np.zeros(n,bool)
            for iteration in range(20):
                previous=(x.copy(),y.copy())
                if start=='a':
                    y=_best_mask(z+sign*np.bincount(cv,weights=x[cu],minlength=n),eb,budget_b)
                    x=_best_mask(w+sign*np.bincount(cu,weights=y[cv],minlength=n),ea,budget_a)
                else:
                    x=_best_mask(w+sign*np.bincount(cu,weights=y[cv],minlength=n),ea,budget_a)
                    y=_best_mask(z+sign*np.bincount(cv,weights=x[cu],minlength=n),eb,budget_b)
                if np.array_equal(x,previous[0]) and np.array_equal(y,previous[1]):break
            gain=float(w@x+z@y+sign*np.sum(x[cu]&y[cv]))
            if gain>best:best=gain;witness=(np.flatnonzero(x).tolist(),np.flatnonzero(y).tolist())
        # The degree bound supplies an independent valid cap.
        upper=min(upper,len(edges)*degree_error_margin(edges,n,budget_a,budget_b),len(edges)-numerator if increase else numerator)
        if best>upper+1e-8:raise AssertionError('Witness exceeds outer bound')
        details.append({'gain_lower':best,'gain_upper':float(upper),'a_flips':witness[0],'b_flips':witness[1]})
    outer=[(numerator-details[0]['gain_upper'])/len(edges),(numerator+details[1]['gain_upper'])/len(edges)]
    inner=[(numerator-details[0]['gain_lower'])/len(edges),(numerator+details[1]['gain_lower'])/len(edges)]
    return {'outer_bounds':outer,'attained_bounds':inner,'endpoint_gaps':[inner[0]-outer[0],outer[1]-inner[1]],'details':details}


def optimize_bounds(edges,labels,budget_a,budget_b,reference=None,reference_budgets=(0,0),relax=False,time_limit=5.):
    """Binary true labels and McCormick edge products; solver bounds, not rounded LP."""
    n=len(labels);m=len(edges);u,v=edges.T
    if m==0:raise ValueError('Nonempty graph required')
    rows=[];cols=[];values=[];lower=[];upper=[]
    def constraint(indices,coefficients,lo,hi):
        index=len(lower);rows.extend([index]*len(indices));cols.extend(indices);values.extend(coefficients);lower.append(lo);upper.append(hi)
    for j,(left,right) in enumerate(edges):
        edge=2*n+j
        constraint([edge,int(left)],[1.,-1.],-np.inf,0.)
        constraint([edge,n+int(right)],[1.,-1.],-np.inf,0.)
        constraint([edge,int(left),n+int(right)],[1.,-1.,-1.],-1.,np.inf)
    sources=[(np.asarray(labels), (budget_a,budget_b))]
    if reference is not None:sources.append((np.asarray(reference),reference_budgets))
    for calls,budgets in sources:
        for marker,budget in enumerate(budgets):
            mask=calls[:,marker]>=0;ids=np.flatnonzero(mask);observed=calls[ids,marker]
            constraint((marker*n+ids).tolist(),(1-2*observed).tolist(),-np.inf,float(budget-np.sum(observed)))
    matrix=coo_array((values,(np.asarray(rows,np.int32),np.asarray(cols,np.int32))),shape=(len(lower),2*n+m)).tocsc()
    objectives=[]
    for sign in (1.,-1.):
        objective=np.r_[np.zeros(2*n),np.full(m,sign)]
        result=milp(objective,integrality=np.r_[np.full(2*n,0 if relax else 1),np.zeros(m)],
                    bounds=Bounds(np.zeros(2*n+m),np.ones(2*n+m)),
                    constraints=LinearConstraint(matrix,lower,upper),
                    options={'time_limit':time_limit,'mip_rel_gap':0.})
        if result.status==2:return {'empty':True,'objectives':[]}
        primal=None if result.fun is None else float(result.fun)
        dual=primal if result.status==0 else getattr(result,'mip_dual_bound',None)
        if dual is not None and not np.isfinite(dual):dual=None
        feasible=None
        if not relax and result.x is not None:
            calls=np.stack([result.x[:n],result.x[n:2*n]],axis=1).round().astype(int)
            if all(np.sum((calls[:,k]!=observed[:,k])&(observed[:,k]>=0))<=budgets[k] for observed,budgets in sources for k in (0,1)):
                feasible=float(np.sum(calls[u,0]*calls[v,1]))*sign
        objectives.append({'status':int(result.status),'optimal':bool(result.success),'objective_primal':primal,
                           'objective_dual_bound':None if dual is None else float(dual),'feasible_integer_objective':feasible})
    # Expand solver tolerances; these are floating solver bounds, not formal rational certificates.
    lo=0. if objectives[0]['objective_dual_bound'] is None else max(0.,(objectives[0]['objective_dual_bound']-1e-7)/m)
    hi=1. if objectives[1]['objective_dual_bound'] is None else min(1.,(-objectives[1]['objective_dual_bound']+1e-7)/m)
    inner=None if relax or any(o['feasible_integer_objective'] is None for o in objectives) else [objectives[0]['feasible_integer_objective']/m,-objectives[1]['feasible_integer_objective']/m]
    return {'empty':False,'outer_bounds':[lo,hi],'attained_bounds':inner,'objectives':objectives,'relaxation':relax}
