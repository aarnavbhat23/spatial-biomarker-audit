"""Sensitivity bounds for fixed graphs and externally bounded marker-call errors."""
import numpy as np


def pair_target(edges,labels):
    return float(np.mean(labels[edges[:,0],0]&labels[edges[:,1],1]))


def top_budget_sum(values,budget):
    values=np.asarray(values)
    if budget<0 or int(budget)!=budget:raise ValueError('Budget must be a nonnegative integer')
    count=min(int(budget),len(values))
    return float(np.sort(values)[-count:].sum()) if count else 0.


def degree_error_margin(edges,n,budget_a,budget_b):
    """Label-blind bound; available before observing either marker."""
    outdegree=np.bincount(edges[:,0],minlength=n)
    indegree=np.bincount(edges[:,1],minlength=n)
    return min(1.,(top_budget_sum(outdegree,budget_a)+top_budget_sum(indegree,budget_b))/len(edges))


def a_weights(edges,labels):
    return np.bincount(edges[:,0],weights=labels[edges[:,1],1],minlength=len(labels))


def sharp_a_bounds(edges,labels,budget_a):
    """Sharp full-census bounds for at most budget_a errors in A, B exact."""
    weights=a_weights(edges,labels);observed=labels[:,0];target=pair_target(edges,labels)
    lower=target-top_budget_sum(weights[observed],budget_a)/len(edges)
    upper=target+top_budget_sum(weights[~observed],budget_a)/len(edges)
    return np.array([max(0.,lower),min(1.,upper)])


def worst_a_flips(edges,labels,budget_a):
    """Maximize absolute target shift from this reference, with deterministic ties."""
    if budget_a<0 or int(budget_a)!=budget_a:raise ValueError('Budget must be a nonnegative integer')
    weights=a_weights(edges,labels);candidates=[]
    for current in (False,True):
        ids=np.flatnonzero(labels[:,0]==current)
        selected=ids[np.argsort(-weights[ids],kind='stable')[:budget_a]]
        candidates.append((float(weights[selected].sum()),selected,current))
    _,selected,current=max(candidates,key=lambda x:x[0])
    corrupted=labels.copy();corrupted[selected,0]=~corrupted[selected,0]
    return corrupted,selected,'decrease' if current else 'increase'


def enlarge(intervals,margins):
    intervals=np.asarray(intervals);margins=np.asarray(margins)
    return np.stack([np.maximum(0.,intervals[...,0]-margins),
                     np.minimum(1.,intervals[...,1]+margins)],axis=-1)
