"""Classical finite-population difference estimation under tile SRS."""
import numpy as np
from scipy.stats import t
from scipy.special import logsumexp


def betting_logwealth(mean, values, population, positive=True):
    """A fixed mixture of nonnegative mean-one bets under tile random ordering."""
    values=np.asarray(values,dtype=float)
    before=np.r_[0.,np.cumsum(values)[:-1]]
    remaining=(population*mean-before)/(population-np.arange(len(values)))
    stabilizers=np.array([0.01,0.05,0.1,0.5,1.0])[:,None]
    if positive:
        factors=(stabilizers+values)/(stabilizers+remaining)
    else:
        factors=(stabilizers+1-values)/(stabilizers+1-remaining)
    return float(logsumexp(np.log(factors).sum(axis=1))-np.log(len(stabilizers)))


def betting_interval(predicted, observed, indices, pair_counts, coefficient=1.0, alpha=0.05):
    """Invert two classical test martingales at a fixed sampling budget.

    The guarantee is conditional on a fixed finite population and uniform random
    tile ordering, with correct pair-total bounds. No spatial independence or
    correct prediction model is assumed. Bisection returns conservative outer
    endpoints. None denotes an empty confidence set on a rejection event.
    """
    predicted,observed,pair_counts=map(lambda x:np.asarray(x,dtype=float),(predicted,observed,pair_counts))
    indices=np.asarray(indices,dtype=int)
    population,size=len(predicted),len(indices)
    if not (0<alpha<1 and 0<=coefficient<=1 and 1<=size<=population) or len(np.unique(indices))!=size:
        raise ValueError('Invalid sampling design or coefficient')
    if len(observed)!=size or np.any(observed<0) or np.any(observed>pair_counts[indices]):
        raise ValueError('Observed totals violate their bounds')
    denominator=pair_counts.sum()
    if size==population:
        return [float(observed.sum()/denominator)]*2
    a=(-coefficient*predicted).min()
    b=(pair_counts-coefficient*predicted).max()
    width=b-a
    remaining=np.ones(population,dtype=bool)
    remaining[indices]=False
    target_low=observed.sum()/denominator
    target_high=(observed.sum()+pair_counts[remaining].sum())/denominator
    offset=coefficient*predicted.sum()+population*a
    scale=population*width
    lo,hi=(denominator*np.array([target_low,target_high])-offset)/scale
    normalized=(observed-coefficient*predicted[indices]-a)/width
    threshold=np.log(2/alpha)
    low,high=float(lo),float(hi)
    if betting_logwealth(high,normalized,population,True)>=threshold or betting_logwealth(low,normalized,population,False)>=threshold:
        return None
    if betting_logwealth(low,normalized,population,True)>=threshold:
        left,right=low,high
        for _ in range(24):
            middle=(left+right)/2
            if betting_logwealth(middle,normalized,population,True)>=threshold:
                left=middle
            else:
                right=middle
        low=left
    if betting_logwealth(high,normalized,population,False)>=threshold:
        left,right=low,high
        for _ in range(24):
            middle=(left+right)/2
            if betting_logwealth(middle,normalized,population,False)>=threshold:
                right=middle
            else:
                left=middle
        high=right
    if low>high:
        return None
    return np.clip((offset+scale*np.array([low,high]))/denominator,0,1).tolist()


def difference_estimate(predicted, observed, indices, pair_counts, coefficient=1.0, alpha=0.05):
    """Observed contains totals for sampled tiles only; predictions/counts are censused.

    Student intervals are approximate. The Hoeffding comparator is a conservative
    bounded-variable interval for uniform sampling without replacement. Neither
    interval addresses measurement error, graph detection, or new-patient shift.
    """
    predicted, observed, pair_counts = map(lambda a:np.asarray(a,dtype=float),
                                           (predicted,observed,pair_counts))
    indices = np.asarray(indices,dtype=int)
    population, size = len(predicted), len(indices)
    if not (2 <= size <= population) or len(np.unique(indices)) != size:
        raise ValueError('Need at least two distinct sampled tiles')
    if len(observed) != size or len(pair_counts) != population or pair_counts.sum() <= 0:
        raise ValueError('Incompatible observed totals or geometric pair counts')
    if not (0 <= coefficient <= 1 and 0 < alpha < 1):
        raise ValueError('Invalid coefficient or error probability')
    if np.any(observed<0) or np.any(observed>pair_counts[indices]):
        raise ValueError('Observed pair total is outside its known bounds')
    denominator = pair_counts.sum()
    residual = observed - coefficient*predicted[indices]
    estimate = (coefficient*predicted.sum() + population*residual.mean())/denominator
    variance = (population/denominator)**2 * (1-size/population) * residual.var(ddof=1)/size
    student_half = t.ppf(1-alpha/2,size-1)*np.sqrt(variance) if size<population else 0.0
    residual_range = (pair_counts-coefficient*predicted).max() - (-coefficient*predicted).min()
    bound_half = (population/denominator)*residual_range*np.sqrt(np.log(2/alpha)/(2*size)) if size<population else 0.0
    return {'estimate':float(estimate), 'variance':float(variance),
            'student_interval':np.clip([estimate-student_half,estimate+student_half],0,1).tolist(),
            'bounded_interval':np.clip([estimate-bound_half,estimate+bound_half],0,1).tolist()}


def union_area(rectangles):
    """Exact area of the union of axis-aligned measurement rectangles."""
    rectangles = np.asarray(rectangles)
    xs = np.unique(rectangles[:,[0,2]])
    area = 0.0
    for x0,x1 in zip(xs[:-1],xs[1:]):
        active = rectangles[(rectangles[:,0]<x1)&(rectangles[:,2]>x0)]
        spans = sorted((r[1],r[3]) for r in active)
        if not spans:
            continue
        bottom,top = spans[0]
        height = 0.0
        for lo,hi in spans[1:]:
            if lo>top:
                height+=top-bottom
                bottom,top=lo,hi
            else:
                top=max(top,hi)
        area+=(x1-x0)*(height+top-bottom)
    return float(area)
