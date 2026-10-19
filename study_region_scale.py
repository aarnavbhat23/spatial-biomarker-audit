"""Classical PPS region sampling with valid bounded-mean confidence sequences."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from shapely import box,union_all

ROOT=Path(__file__).resolve().parent
LAMBDAS=np.array([.125,.25,.5,1,2,4,8,16,32,64,128,256])


def logwealth(mean,values,positive):
    values=np.asarray(values)
    if positive:
        bet=np.minimum(LAMBDAS,.5/max(mean,1e-15))
        factors=1+bet[:,None]*(values[None,:]-mean)
    else:
        bet=np.minimum(LAMBDAS,.5/max(1-mean,1e-15))
        factors=1+bet[:,None]*(mean-values[None,:])
    return float(logsumexp(np.log(factors).sum(1))-np.log(len(LAMBDAS)))


def bounded_sequence(values,alpha=.05):
    """Invert two mixtures of bounded-mean test martingales.

    Each fixed bet is nonnegative and has conditional mean one at the true
    mean. The cap makes both wealths monotone in the candidate mean. Ville's
    inequality gives simultaneous coverage over time for IID index sampling.
    Bisection retains outer endpoints. These are established betting methods.
    """
    threshold=np.log(2/alpha)
    low,high=0.,1.
    if logwealth(0.,values,True)>=threshold:
        a,b=0.,1.
        for _ in range(24):
            mid=(a+b)/2
            if logwealth(mid,values,True)>=threshold:a=mid
            else:b=mid
        low=a
    if logwealth(1.,values,False)>=threshold:
        a,b=0.,1.
        for _ in range(24):
            mid=(a+b)/2
            if logwealth(mid,values,False)>=threshold:b=mid
            else:a=mid
        high=b
    return low,high


def frame(data,core):
    anchors=data['edges'][:,0]
    # Only anchor cells in the original fixed core panel define the target.
    cells=data['coords'][anchors]//core
    tiles,inverse=np.unique(cells,axis=0,return_inverse=True)
    counts=np.bincount(inverse)
    outcome=data['labels'][anchors,0]&data['labels'][data['edges'][:,1],1]
    positive=np.bincount(inverse,weights=outcome,minlength=len(tiles))
    rectangles=np.c_[tiles*core-80,(tiles+1)*core+80]
    return {'counts':counts,'positive':positive,'rectangles':rectangles,'tiles':tiles,'owners':inverse}


def acquisition_trial(f,draws,width=.04):
    counts=f['counts'];positive=f['positive'];total=counts.sum()
    seen=set();values=[];visited=[];stopped=False
    for k,selected in enumerate(draws):
        selected=int(selected);seen.add(selected);values.append(positive[selected]/counts[selected])
        # Inspect at fixed checkpoints; confidence sequence also permits stopping.
        if k+1 not in (8,16,32,64,128,256,512,1024,2048):continue
        ids=np.array(sorted(seen))
        known_positive=positive[ids].sum()
        lower=known_positive/total
        upper=(known_positive+total-counts[ids].sum())/total
        if len(ids)==len(counts):
            lo=hi=lower
        else:
            lo,hi=bounded_sequence(values)
            lo,hi=max(lo,lower),min(hi,upper)
        visited.append({'draws':k+1,'interval':[float(lo),float(hi)]})
        if 0<=hi-lo<=width:
            stopped=True;break
    ids=np.array(sorted(seen))
    rectangles=f['rectangles'][ids]
    area=float(union_all(box(*rectangles.T)).area)
    return {'draws':len(values),'unique_regions':len(ids),'area_native_pixels':area,
            'interval':[float(lo),float(hi)],'stopped':stopped,'empty_confidence_set':bool(lo>hi),
            'estimate':float(np.mean(values)), 'checkpoints':visited}


def evaluate(data,seed=20261022,repetitions=100,cores=(128,256,512,2048),width=.04):
    output=[];rng=np.random.default_rng(seed)
    truth=float(np.mean(data['labels'][data['edges'][:,0],0]&data['labels'][data['edges'][:,1],1]))
    for core in cores:
        f=frame(data,core)
        rows=[]
        for repeat in range(repetitions):
            draws=rng.choice(len(f['counts']),2048,p=f['counts']/f['counts'].sum())
            row=acquisition_trial(f,draws,width=width)
            row['covered']=bool(row['interval'][0]-1e-12<=truth<=row['interval'][1]+1e-12)
            row['repeat']=repeat
            rows.append(row)
        output.append({'core_native_pixels':core,'frame_regions':len(f['counts']),'truth':truth,
                       'target_width':width,
                       'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                       'coverage_at_stop':float(np.mean([r['covered'] for r in rows])),
                       'mean_draws':float(np.mean([r['draws'] for r in rows])),
                       'mean_unique_regions':float(np.mean([r['unique_regions'] for r in rows])),
                       'mean_area_native_pixels':float(np.mean([r['area_native_pixels'] for r in rows])),
                       'mean_width':float(np.mean([r['interval'][1]-r['interval'][0] for r in rows])),
                       'rows':rows})
        print('core',core,'area',output[-1]['mean_area_native_pixels'],'draws',output[-1]['mean_draws'],flush=True)
    return output


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repetitions',type=int,default=100)
    parser.add_argument('--precision-stress',action='store_true')
    args=parser.parse_args()
    if args.precision_stress:
        result={'scope':'Exploratory precision stress after confirmation; not a second confirmation.',
                'repetitions':args.repetitions,'patients':{}}
        for patient in ('CRC05','CRC06','CRC07','CRC09'):
            data=dict(np.load(ROOT/'data'/f'{patient}-acquisition.npz'))
            result['patients'][patient]={}
            for width in (.01,.005):
                result['patients'][patient][str(width)]=evaluate(data,seed=20261025+int(patient[3:]),
                    repetitions=args.repetitions,cores=(128,2048),width=width)
            (ROOT/'results/region_scale_precision_stress.json').write_text(json.dumps(result,indent=2)+'\n')
        return
    result={'seed':20261022,'target_width':.04,'alpha':.05,'repetitions':args.repetitions,
            'patients':{},'scope':'Development-only comparison; classical PPS and betting; geometry from measured segmentation.'}
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        result['patients'][patient]=evaluate(data,repetitions=args.repetitions)
        (ROOT/'results/region_scale_development.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
