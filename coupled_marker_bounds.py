"""Sharp census bounds for one A error and an arbitrary B error budget."""
import numpy as np


def top_budget(values,budget):
    if int(budget)!=budget or budget<0:raise ValueError('Budget must be a nonnegative integer')
    k=min(int(budget),len(values))
    if k==0:return 0.,0,0,0
    threshold=int(np.partition(values,len(values)-k)[len(values)-k])
    above=int(np.sum(values>threshold));at_least=int(np.sum(values>=threshold))
    total=float(values[values>threshold].sum()+threshold*(k-above))
    return total,threshold,above,at_least


def sharp_one_a_bounds(edges,labels,budget_b):
    """Simple directed graph, unit edge weights, at most one A flip and k B flips."""
    if int(budget_b)!=budget_b or budget_b<0:raise ValueError('Budget must be a nonnegative integer')
    n=len(labels);a=labels[:,0].astype(bool);b=labels[:,1].astype(bool);u,v=edges.T
    wa=np.bincount(u,weights=b[v],minlength=n).astype(int)
    zb=np.bincount(v,weights=a[u],minlength=n).astype(int)
    numerator=int(np.sum(a[u]&b[v]));extremes=[]
    for increase in (False,True):
        eligible_b=~b if increase else b;eligible_a=~a if increase else a
        k=min(budget_b,int(eligible_b.sum()))
        base,threshold,above,at_least=top_budget(zb[eligible_b],k)
        if k:
            mask=eligible_b[v];scores=zb[v[mask]];anchors=u[mask]
            high=np.bincount(anchors,weights=scores>threshold,minlength=n)
            equal=np.bincount(anchors,weights=scores==threshold,minlength=n)
            correction=high+np.minimum(equal,k-above) if increase else -high-np.maximum(0,equal-(at_least-k))
        else:correction=np.zeros(n)
        gains=wa+base+correction
        best=max(base,float(np.max(gains[eligible_a])) if eligible_a.any() else base)
        extremes.append(best)
    return np.array([(numerator-extremes[0])/len(edges),(numerator+extremes[1])/len(edges)])
