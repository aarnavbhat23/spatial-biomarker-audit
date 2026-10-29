"""Conventional molecular kernel proximity under conserved localization error.

This changes the observable; it is not a cell-type interaction estimator.
"""
import argparse
from itertools import product
import json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parent


def proximity(a,b,bandwidth):
    return float(np.exp(-cdist(a,b,'sqeuclidean')/(2*bandwidth**2)).mean())


def bounds(a,b,bandwidth,radius):
    distance=cdist(a,b)
    lower=float(np.exp(-(distance+2*radius)**2/(2*bandwidth**2)).mean())
    upper=float(np.exp(-np.maximum(0,distance-2*radius)**2/(2*bandwidth**2)).mean())
    score=proximity(a,b,bandwidth)
    error=2*radius*np.exp(-.5)/bandwidth
    # Outward guard for these numerical computations, not interval arithmetic.
    return {'score':score,'pairwise_bounds':[max(0.,lower-1e-12),min(1.,upper+1e-12)],
            'lipschitz_bounds':[max(0.,score-error),min(1.,score+error)]}


def oracle(a,b,bandwidth):
    values=[]
    for changes in product((-1,0,1),repeat=6):
        left=a.copy(); right=b.copy()
        left[:,0]+=changes[:3]; right[:,0]+=changes[3:]
        # Separate scalar pair calculation instead of cdist/proximity.
        total=0.
        for x in left:
            for y in right:
                total += np.exp(-sum(float(u-v)**2 for u,v in zip(x,y))/(2*bandwidth**2))
        values.append(total/9)
    return [float(min(values)),float(max(values))]


def run():
    small=[]; paired=[]
    for seed in (31,47,83):
        rng=np.random.default_rng(seed)
        a=np.rint(rng.uniform(0,20,(3,2))); b=np.rint(rng.uniform(0,20,(3,2)))
        for bandwidth in (4.,8.):
            outer=bounds(a,b,bandwidth,1.); exact=oracle(a,b,bandwidth)
            assert outer['pairwise_bounds'][0] <= exact[0] and outer['pairwise_bounds'][1] >= exact[1]
            small.append({'seed':seed,'bandwidth':bandwidth,'a':a.tolist(),'b':b.tolist(),
                          'bound':outer,'exact_discrete_extrema':exact,'worlds':729})
    for seed in (101,202,303):
        rng=np.random.default_rng(seed); noise_a=rng.normal(0,10,(256,2)); noise_b=rng.normal(0,10,(256,2))
        shared=(np.clip(noise_a+[50,50],0,100),np.clip(noise_b+[50,50],0,100))
        separate=(np.clip(noise_a+[30,50],0,100),np.clip(noise_b+[70,50],0,100))
        for bandwidth in (10.,20.):
            for radius in (1.,3.):
                first=bounds(*shared,bandwidth,radius); second=bounds(*separate,bandwidth,radius)
                # One additional unobserved molecule per channel changes a bounded
                # pair mean by at most the sum of the two missing-mass fractions.
                missing_error=2/257
                certified=first['pairwise_bounds'][0] > second['pairwise_bounds'][1]
                with_missing=first['pairwise_bounds'][0]-missing_error > second['pairwise_bounds'][1]+missing_error
                paired.append({'seed':seed,'bandwidth':bandwidth,'radius':radius,'per_channel_molecules':256,
                               'shared':first,'separated':second,'intervals_disjoint':certified,
                               'missing_mass_error_per_specimen':missing_error,
                               'disjoint_allowing_one_missing_per_channel':with_missing})
    # A deterministic budget violation: observed separated channels are actually
    # co-located, requiring a 20-unit shift per channel, outside rho=1.
    a=np.array([[30.,50.]]); b=np.array([[70.,50.]])
    declared=bounds(a,b,10.,1.)['pairwise_bounds']
    actual=proximity(np.array([[50.,50.]]),np.array([[50.,50.]]),10.)
    assert actual > declared[1]
    passed=sum(r['intervals_disjoint'] for r in paired)
    return {'protocol_commit':'1215d93','scope':'Known kernel observable; development feasibility only, no novelty or biological coverage claim.',
            'small':small,'paired_scaling':paired,
            'disjoint_pairs':passed,'continuation_gate_passed':passed >= 9,
            'out_of_budget_counterexample':{'declared_radius':1.,'actual_shift_each':20.,'declared_bounds':declared,'actual_proximity':actual}}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true'); args=parser.parse_args()
    result=run(); path=ROOT/'results/conserved_spatial_mass.json'
    if args.verify:
        saved=json.loads(path.read_text())
        assert saved['protocol_commit'] == result['protocol_commit']
        assert len(saved['small']) == len(result['small']) == 6
        assert len(saved['paired_scaling']) == len(result['paired_scaling']) == 12
        for a,b in zip(saved['small'],result['small']):
            assert (a['seed'],a['bandwidth'],a['worlds']) == (b['seed'],b['bandwidth'],b['worlds'])
            np.testing.assert_array_equal(a['a'],b['a']); np.testing.assert_array_equal(a['b'],b['b'])
            np.testing.assert_allclose(a['exact_discrete_extrema'],b['exact_discrete_extrema'],atol=1e-12,rtol=0)
            np.testing.assert_allclose(a['bound']['pairwise_bounds'],b['bound']['pairwise_bounds'],atol=1e-12,rtol=0)
        for a,b in zip(saved['paired_scaling'],result['paired_scaling']):
            assert (a['seed'],a['bandwidth'],a['radius']) == (b['seed'],b['bandwidth'],b['radius'])
            for key in ('shared','separated'):
                np.testing.assert_allclose(a[key]['pairwise_bounds'],b[key]['pairwise_bounds'],atol=1e-12,rtol=0)
                np.testing.assert_allclose(a[key]['lipschitz_bounds'],b[key]['lipschitz_bounds'],atol=1e-12,rtol=0)
                np.testing.assert_allclose(a[key]['score'],b[key]['score'],atol=1e-12,rtol=0)
            assert a['intervals_disjoint'] == b['intervals_disjoint']
            assert a['disjoint_allowing_one_missing_per_channel'] == b['disjoint_allowing_one_missing_per_channel']
            assert a['missing_mass_error_per_specimen'] == b['missing_mass_error_per_specimen']
        assert saved['out_of_budget_counterexample'] == result['out_of_budget_counterexample']
        assert saved['disjoint_pairs'] == result['disjoint_pairs']
        assert saved['continuation_gate_passed'] == result['continuation_gate_passed']
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Disjoint pairs:',result['disjoint_pairs'],'/12; gate:',result['continuation_gate_passed'])
    print('Allowing one missing molecule per channel:',sum(r['disjoint_allowing_one_missing_per_channel'] for r in result['paired_scaling']))
