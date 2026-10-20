"""Sharp partial-observation intervals: one uncertain marker, exact baseline B."""
import numpy as np
from measurement_error import top_budget_sum


def intervals(weights,totals,recorded,seen,reviewed,verified,budget):
    """Unknown recorded calls must be placeholders; only seen values are used."""
    errors=int(np.sum(reviewed&(recorded!=verified)));remaining=budget-errors
    if remaining<0:raise ValueError('Reviewed errors exceed the assumed error budget')
    known=weights[:,reviewed]@verified[reviewed].astype(float)
    pending=seen&~reviewed;unknown=~seen
    base=known+weights[:,pending]@recorded[pending].astype(float)
    low=[];high=[]
    for j,w in enumerate(weights):
        low.append((base[j]-top_budget_sum(w[pending&recorded],remaining))/totals[j])
        high.append((base[j]+top_budget_sum(w[pending&~recorded],remaining)+w[unknown].sum())/totals[j])
    return np.stack([low,high],axis=1)


def worst_added_ambiguity(positive,negative,new,budget):
    """Exact worst width contribution over all binary calls on a new field."""
    if budget==0:return 0.
    b=min(budget,len(positive)+len(negative)+len(new))
    pos=np.r_[0.,np.cumsum(np.sort(positive)[::-1])];neg=np.r_[0.,np.cumsum(np.sort(negative)[::-1])]
    fresh=np.r_[0.,np.cumsum(np.sort(new)[::-1])]
    a=np.arange(b+1)[:,None];c=np.arange(b+1)[None,:];count=a+c
    feasible=count<=len(new)
    values=pos[np.minimum(b-a,len(positive))]+neg[np.minimum(b-c,len(negative))]+fresh[np.minimum(count,len(new))]
    return float(np.max(values[feasible]))


def choose_action(weights,totals,owners,recorded,seen,reviewed,verified,budget,review_cost,policy,rng=None,field_cost=1.1):
    current=intervals(weights,totals,recorded,seen,reviewed,verified,budget)
    width=float(np.max(np.diff(current,axis=1)));remaining=budget-int(np.sum(reviewed&(recorded!=verified)))
    pending=seen&~reviewed;unseen=~seen;candidates=[]
    for field in np.unique(owners[unseen]):
        new=unseen&(owners==field)
        worst=[]
        for j,w in enumerate(weights):
            ambiguity=worst_added_ambiguity(w[pending&recorded],w[pending&~recorded],w[new],remaining)
            worst.append((w[unseen&~new].sum()+ambiguity)/totals[j])
        gain=width-max(worst)
        candidates.append((gain/field_cost, 'field',int(field)))
    # Only top-budget+1 members of each sign/claim can improve one-step ambiguity.
    review_candidates=set()
    for w in weights:
        for sign in (False,True):
            ids=np.flatnonzero(pending&(recorded==sign))
            review_candidates.update(ids[np.argsort(-w[ids],kind='stable')[:remaining+1]].tolist())
    if policy=='random_review':
        if candidates:return max(candidates,key=lambda x:(x[0],-x[2]))[1:]
        ids=np.flatnonzero(pending)
        return ('review',int(rng.choice(ids))) if len(ids) else None
    if policy=='census_first' and candidates:return max(candidates,key=lambda x:(x[0],-x[2]))[1:]
    for cell in sorted(review_candidates):
        updated=reviewed.copy();updated[cell]=True
        # Matching the recorded call is the worst review outcome: it spends no error budget.
        hypothetical=verified.copy();hypothetical[cell]=recorded[cell]
        future=intervals(weights,totals,recorded,seen,updated,hypothetical,budget)
        gain=width-float(np.max(np.diff(future,axis=1)))
        candidates.append((gain/review_cost,'review',int(cell)))
    if not candidates:return None
    return max(candidates,key=lambda x:(x[0],x[1]=='field',-x[2]))[1:]
