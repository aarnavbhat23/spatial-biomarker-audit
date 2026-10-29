"""Conditional permutation sensitivity for a fixed measured-molecule population."""
import numpy as np
from scipy.spatial.distance import cdist


def scores(kernel,assignments):
    labels=np.atleast_2d(assignments).astype(float)
    na=labels.sum(axis=1);nb=labels.shape[1]-na
    if np.any(na==0) or np.any(nb==0):raise ValueError('Both channels required')
    return (labels@kernel.sum(axis=1)-np.sum((labels@kernel)*labels,axis=1))/(na*nb)


def kernel_bounds(positions,bandwidth,radius):
    d=cdist(positions,positions)
    nominal=np.exp(-d*d/(2*bandwidth*bandwidth))
    lower=np.exp(-(d+2*radius)**2/(2*bandwidth*bandwidth))
    upper=np.exp(-np.maximum(0,d-2*radius)**2/(2*bandwidth*bandwidth))
    return nominal,lower,upper


def permutations(labels,blocks,count,rng):
    output=np.tile(labels,(count,1))
    for block in np.unique(blocks):
        ids=np.flatnonzero(blocks==block)
        for row in output:row[ids]=rng.permutation(labels[ids])
    return output


def evaluate(positions,labels,permuted,bandwidth,radius,swaps):
    nominal,low,high=kernel_bounds(positions,bandwidth,radius)
    na=int(np.sum(labels));nb=len(labels)-na
    error=swaps/na+swaps/nb
    observed=float(scores(nominal,labels)[0])
    observed_lower=max(0.,float(scores(low,labels)[0])-error-1e-12)
    upper=np.minimum(1.,scores(high,permuted)+error+1e-12)
    null=scores(nominal,permuted)
    nominal_count=int(np.sum(null>=observed-1e-14))
    robust_count=int(np.sum(upper>=observed_lower))
    return {'score':observed,'observed_lower_bound':observed_lower,'swap_error':error,
            'nominal_exceedances':nominal_count,'robust_upper_exceedances':robust_count,
            'nominal_p':(1+nominal_count)/(len(permuted)+1),
            'robust_upper_p':(1+robust_count)/(len(permuted)+1)}
