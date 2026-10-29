"""Frozen VISp marker/scale family, with all missing markers retained."""
from itertools import combinations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from molecular_permutation import permutations,evaluate,kernel_bounds,scores

ROOT=Path(__file__).resolve().parent
GENES=('Gad2','Slc17a7','Pvalb','Sst','Vip','Cux2')
BANDWIDTHS=(5.,10.,20.)
ERRORS=((0.,0),(.5,0),(.5,1),(1.,1))


def main():
    source=ROOT/'data/visp/s3_spot_table.csv'
    checksum=hashlib.sha256(source.read_bytes()).hexdigest()
    manifest=json.loads((ROOT/'results/visp_source_manifest.json').read_text())
    assert checksum==manifest['sha256']
    df=pd.read_csv(source); xy=df[['rotated_x','rotated_y']].to_numpy()
    assert np.isfinite(xy).all()
    lower=xy.min(axis=0);upper=xy.max(axis=0);middle=(lower+upper)/2
    rng=np.random.default_rng(20261031);selection=[];rows=[]
    for half in ('left_development','right_confirmation'):
        region=(xy[:,0]<middle[0]) if half.startswith('left') else (xy[:,0]>=middle[0])
        chosen={}
        for gene in GENES:
            ids=np.flatnonzero(region & (df.gene.to_numpy()==gene))
            sampled=np.sort(rng.choice(ids,size=min(200,len(ids)),replace=False))
            chosen[gene]=sampled
            selection.append({'region':half,'gene':gene,'available':len(ids),'selected':len(sampled),'source_row_indices':sampled.tolist()})
        xmid=(lower[0]+middle[0])/2 if half.startswith('left') else (middle[0]+upper[0])/2
        for a,b in combinations(GENES,2):
            ia,ib=chosen[a],chosen[b]
            if min(len(ia),len(ib))<20:
                rows.append({'region':half,'genes':[a,b],'status':'unevaluable','reason':'Missing or fewer than 20 selected molecules','counts':[len(ia),len(ib)]})
                continue
            indices=np.r_[ia,ib];positions=xy[indices]
            labels=np.r_[np.ones(len(ia),int),np.zeros(len(ib),int)]
            blocks=2*(positions[:,0]>=xmid)+(positions[:,1]>=middle[1])
            blocked=permutations(labels,blocks,1999,rng)
            global_labels=permutations(labels,np.zeros(len(labels),int),1999,rng)
            same_cell=(df.cell.to_numpy()[ia,None]==df.cell.to_numpy()[ib][None,:])
            assigned=(df.cell.to_numpy()[ia,None]>0)&(df.cell.to_numpy()[ib][None,:]>0)
            for h in BANDWIDTHS:
                kernel=kernel_bounds(positions,h,0.)[0]
                observed=float(scores(kernel,labels)[0]);global_null=scores(kernel,global_labels)
                global_p=(1+int(np.sum(global_null>=observed-1e-14)))/2000
                results=[]
                for radius,swaps in ERRORS:
                    result=evaluate(positions,labels,blocked,h,radius,swaps)
                    result.update(radius=radius,swaps=swaps,
                                  bonferroni_upper_p=min(1.,45*result['robust_upper_p']),
                                  bonferroni_nominal_p=min(1.,45*result['nominal_p']))
                    results.append(result)
                cross=kernel[:len(ia),len(ia):]
                numerator=float(cross[same_cell & assigned].sum());denominator=float(cross.sum())
                rows.append({'region':half,'genes':[a,b],'status':'evaluated','bandwidth_um':h,'counts':[len(ia),len(ib)],
                             'stratum_counts':[{'block':int(c),'a':int(labels[blocks==c].sum()),'b':int(np.sum(blocks==c)-labels[blocks==c].sum())} for c in np.unique(blocks)],
                             'global_nominal_p':global_p,'global_bonferroni_p':min(1.,45*global_p),
                             'assigned_same_cell_kernel_fraction':numerator/denominator if denominator else None,
                             'same_cell_scope':'Source segmentation diagnostic only; not truth or a test filter.',
                             'sensitivity':results})
                print(half,a,b,h,'nominal',results[0]['bonferroni_nominal_p'],'robust',results[2]['bonferroni_upper_p'],flush=True)
    counts={}
    for half in ('left_development','right_confirmation'):
        selected=[r for r in rows if r['region']==half and r['status']=='evaluated']
        counts[half]={'evaluated_tests':len(selected),'missing_or_sparse_pairs':sum(r['region']==half and r['status']=='unevaluable' for r in rows),
                      'global_nominal_rejections':sum(r['global_bonferroni_p']<=.05 for r in selected),
                      'stratified_nominal_rejections':sum(r['sensitivity'][0]['bonferroni_nominal_p']<=.05 for r in selected),
                      'robust_rejections':[sum(r['sensitivity'][j]['bonferroni_upper_p']<=.05 for r in selected) for j in range(4)]}
    confirmation=counts['right_confirmation']
    gate=(confirmation['stratified_nominal_rejections']>0 and confirmation['robust_rejections'][2]>=1
          and confirmation['robust_rejections'][1]>=confirmation['stratified_nominal_rejections']/2)
    report={'protocol_commit':'a0c9526','source_sha256':checksum,'seed':20261031,'permutations':1999,
            'coordinate_bounds_um':[lower.tolist(),upper.tolist()],'split_x_um':float(middle[0]),
            'selection':selection,'rows':rows,'summary':counts,'continuation_gate_passed':bool(gate),
            'scope':'One public field, selected measured molecules; conditional permutation null and uncalibrated sensitivity budgets, not patient replication.'}
    (ROOT/'results/molecular_confirmation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(counts,indent=2));print('Gate:',gate)


if __name__=='__main__':main()
