"""Residual PPS with nested observation reuse and predictable-scale betting."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from shapely import box, union_all

from study_region_scale import LAMBDAS, frame

ROOT=Path(__file__).resolve().parent


def residual_wealth(mean, values, lower, upper, positive):
    """Predictable bounded bets; each factor has conditional mean one at truth."""
    values=np.asarray(values);lower=np.asarray(lower);upper=np.asarray(upper)
    gap=np.maximum(mean-lower if positive else upper-mean,1e-15)
    bets=np.minimum(LAMBDAS[:,None],.5/gap[None,:])
    differences=values-mean if positive else mean-values
    factors=1+bets*differences[None,:]
    if np.any(factors<=0):
        raise ValueError('Candidate outside the retained predictable support')
    return float(logsumexp(np.log(factors).sum(axis=1))-np.log(len(LAMBDAS)))


def residual_interval(values, lower, upper, known_lower, known_upper, alpha=.05):
    """Invert monotone mixtures inside current deterministic bounds."""
    if known_lower==known_upper:return known_lower,known_upper
    threshold=np.log(2/alpha)
    a,b=known_lower,known_upper
    if residual_wealth(a,values,lower,upper,True)<threshold:
        lo=a
    else:
        for _ in range(24):
            mid=(a+b)/2
            if residual_wealth(mid,values,lower,upper,True)>=threshold:a=mid
            else:b=mid
        lo=a
    a,b=known_lower,known_upper
    if residual_wealth(b,values,lower,upper,False)<threshold:
        hi=b
    else:
        for _ in range(24):
            mid=(a+b)/2
            if residual_wealth(mid,values,lower,upper,False)>=threshold:b=mid
            else:a=mid
        hi=b
    return lo,hi


def nested_frames(data):
    fine=frame(data,128)
    frames={}
    for core in (128,512,2048):
        tiles,owner=np.unique(fine['tiles']//(core//128),axis=0,return_inverse=True)
        frames[core]={'owner':owner,'regions':len(tiles),
                      'rectangles':np.c_[tiles*core-80,(tiles+1)*core+80]}
    return fine['counts'],fine['positive'],frames


def choose_scale(remaining,frames,policy,overhead):
    candidates=(int(policy),) if policy!='adaptive' else tuple(frames)
    scores={};masses={}
    for core in candidates:
        masses[core]=np.bincount(frames[core]['owner'],weights=remaining,
                                minlength=frames[core]['regions'])
        expected_gain=np.dot(masses[core],masses[core])/remaining.sum()
        scores[core]=expected_gain/((core+160)**2+overhead)
    selected=max(scores,key=scores.get)
    return selected,masses[selected]


def residual_trial(counts,positive,frames,seed,policy,width,overhead,max_draws=2048):
    rng=np.random.default_rng(seed)
    total=int(counts.sum());remaining=counts.copy()
    known_positive=0.;values=[];lower=[];upper=[];rectangles=[];scales=[];checkpoints=[]
    truth=float(positive.sum()/total)
    covered_all=True;stopped=False
    lo,hi=0.,1.
    for step in range(1,max_draws+1):
        unknown=int(remaining.sum())
        if unknown==0:break
        core,masses=choose_scale(remaining,frames,policy,overhead)
        selected=int(rng.choice(len(masses),p=masses/unknown))
        members=(frames[core]['owner']==selected)&(remaining>0)
        acquired_count=int(remaining[members].sum())
        acquired_positive=float(positive[members].sum())
        lower.append(known_positive/total);upper.append((known_positive+unknown)/total)
        values.append(known_positive/total+unknown/total*acquired_positive/acquired_count)
        known_positive+=acquired_positive
        remaining[members]=0
        rectangles.append(frames[core]['rectangles'][selected]);scales.append(core)
        current_low=known_positive/total
        current_high=(known_positive+int(remaining.sum()))/total
        # Deterministic completion can safely be inspected after every acquisition.
        if step in (8,16,32,64,128,256,512,1024,2048) or step==max_draws or current_high-current_low<=width:
            lo,hi=residual_interval(values,lower,upper,current_low,current_high)
            covered=lo-1e-12<=truth<=hi+1e-12
            covered_all=covered_all and covered
            checkpoints.append({'draws':step,'interval':[lo,hi],'covered':bool(covered)})
            if 0<=hi-lo<=width:
                stopped=True;break
    rects=np.asarray(rectangles)
    area=float(union_all(box(*rects.T)).area)
    summed_area=float(np.sum((rects[:,2]-rects[:,0])*(rects[:,3]-rects[:,1])))
    return {'draws':len(values),'interval':[lo,hi],'stopped':stopped,
            'empty_confidence_set':lo>hi,'covered_all_checkpoints':bool(covered_all),
            'area_native_pixels':area,'summed_area_native_pixels':summed_area,
            'union_cost':area+overhead*len(values),'summed_cost':summed_area+overhead*len(values),
            'scale_counts':{str(core):scales.count(core) for core in frames},'checkpoints':checkpoints}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repetitions',type=int,default=10)
    args=parser.parse_args()
    result={'seed':20261027,'scope':'Exploratory development panels; stronger without-replacement baselines; no independent confirmation.',
            'repetitions':args.repetitions,'patients':{}}
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        counts,positive,frames=nested_frames(data)
        entries=[]
        for width in (.04,.005):
            for overhead in (0.,500000.):
                for policy in ('128','2048','adaptive'):
                    rows=[residual_trial(counts,positive,frames,20261027+1000*int(patient[3:])+repeat,
                                         policy,width,overhead) for repeat in range(args.repetitions)]
                    item={'policy':policy,'width_target':width,'overhead':overhead,'rows':rows,
                          'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                          'coverage_all_checkpoints':float(np.mean([r['covered_all_checkpoints'] for r in rows])),
                          'empty_sets':sum(r['empty_confidence_set'] for r in rows)}
                    for metric in ('draws','union_cost','summed_cost','area_native_pixels'):
                        item['mean_'+metric]=float(np.mean([r[metric] for r in rows]))
                    entries.append(item)
                    print(patient,width,overhead,policy,'cost',round(item['mean_union_cost']),
                          'draws',item['mean_draws'],'attain',item['stop_fraction'],flush=True)
        result['patients'][patient]=entries
        (ROOT/'results/residual_acquisition_study.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
