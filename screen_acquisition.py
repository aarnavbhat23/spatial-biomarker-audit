"""Exploratory acquisition screens on already-used development panels."""
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta,binom
from shapely import box,union_all

ROOT=Path(__file__).resolve().parent


def exact_interval(positive,size,alpha=.05):
    positive=np.asarray(positive)
    low=beta.ppf(alpha/2,np.maximum(positive,1),size-positive+1)
    high=beta.ppf(1-alpha/2,positive+1,np.maximum(size-positive,1))
    return np.where(positive==0,0,low),np.where(positive==size,1,high)


def binomial_summary(probability,size):
    k=np.arange(size+1)
    low,high=exact_interval(k,size)
    mass=binom.pmf(k,size,probability)
    return {'sampled_pairs':size,'exact_coverage':float(mass@((low<=probability)&(probability<=high))),
            'expected_width':float(mass@(high-low)),
            'probability_width_below_004':float(mass@((high-low)<=.04)),
            'rmse':float(np.sqrt(probability*(1-probability)/size))}


def miss_probability(population,covered,size):
    return math.comb(population-covered,size)/math.comb(population,size) if population-covered>=size else 0.


def overlap_screen(data):
    coords=data['coords'];rectangles=data['rectangles'];edges=data['edges'];population=len(rectangles)
    masks=np.zeros(len(coords),dtype=np.uint64)
    for k,(x0,y0,x1,y1) in enumerate(rectangles):
        inside=(coords[:,0]>=x0)&(coords[:,0]<x1)&(coords[:,1]>=y0)&(coords[:,1]<y1)
        masks[inside]|=np.uint64(1)<<np.uint64(k)
    a,b=masks[edges[:,0]],masks[edges[:,1]]
    ac,bc,uc=map(np.bitwise_count,(a,b,a|b))
    signatures=np.unique(np.sort(np.stack([a,b],axis=1),axis=1),axis=0)
    counts,count_frequency=np.unique(np.stack([ac,bc,uc],axis=1),axis=0,return_counts=True)
    rows=[]
    for size in (4,8,16):
        if size>population:continue
        observed=0.
        for (x,y,u),frequency in zip(counts,count_frequency):
            inclusion=1-miss_probability(population,int(x),size)-miss_probability(population,int(y),size)+miss_probability(population,int(u),size)
            observed+=frequency*inclusion
        owned=len(edges)*size/population
        rows.append({'sampled_tiles':size,'ownership_expected_pairs':owned,'all_observed_expected_pairs':observed,
                     'relative_additional_pairs':float(observed/owned-1)})
    return {'observation_signatures':len(signatures),'geometric_pairs':len(edges),'rows':rows}


def physical_costs(data,size,rng,repetitions=5):
    edges=data['edges'];coords=data['coords']
    output=[]
    for side in (128,256,512):
        areas=[];counts=[]
        for _ in range(repetitions):
            selected=rng.integers(len(edges),size=size)
            center=coords[edges[selected]].mean(axis=1)
            # Identical physical squares need only be acquired once.
            center=np.unique(center,axis=0)
            lower=center-side/2;upper=center+side/2
            area=union_all(box(lower[:,0],lower[:,1],upper[:,0],upper[:,1])).area
            areas.append(area);counts.append(len(center))
        output.append({'side_native_pixels':side,'mean_union_area_native_pixels':float(np.mean(areas)),
                       'mean_regions':float(np.mean(counts)),'cost_repetitions':repetitions})
    return output


def main():
    rng=np.random.default_rng(20261020)
    report={'seed':20261020,'scope':'Exploratory existing-panel development; not independent confirmation.',
            'patients':{},'assumptions':['Known source-segmentation graph and error-free source positivity calls.',
              'Pair indices sampled independently with replacement; intervals condition on the fixed panel.',
              'Both endpoints lie within a square centered at their midpoint; radius is 80 native pixels.',
              'Area is only a cost proxy; per-region overhead is unmeasured.']}
    for patient in ('CRC02','CRC11','CRC19'):
        d=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        edges=d['edges'];owners=d['edge_tiles'];population=len(d['rectangles'])
        outcome=(d['labels'][edges[:,0],0]&d['labels'][edges[:,1],1])
        truth=float(outcome.mean())
        counts=np.bincount(owners,minlength=population)
        coarse=[]
        for size in (4,8,16):
            if size>population:continue
            areas=[]
            for _ in range(100):
                selected=rng.choice(population,size,replace=False)
                rectangles=d['rectangles'][selected]
                areas.append(union_all(box(*rectangles.T)).area)
            coarse.append({'sampled_tiles':size,'mean_union_area_native_pixels':float(np.mean(areas)),
                           'miss_single_hotspot_probability':1-size/population,
                           'largest_single_tile_pair_mass':float(counts.max()/counts.sum())})
        pairs=[]
        for size in (128,256,512,1024,2048,4096):
            row=binomial_summary(truth,size)
            row['costs']=physical_costs(d,size,rng)
            pairs.append(row)
        report['patients'][patient]={'truth':truth,'overlap':overlap_screen(d),'coarse_costs':coarse,'pair_sampling':pairs}
        print(patient,'truth',truth,'overlap',report['patients'][patient]['overlap']['rows'],flush=True)
    (ROOT/'results/acquisition_direction_screen.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
