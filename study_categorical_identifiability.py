"""Reuse the six enrichment development patients; never select fresh favorable cases."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from danenberg_confirmation import COLUMNS,retain,selected_labels,knn_edges
from categorical_identifiability import categorical_witness,repair_overlap,review_pair_bound
from contrast_count_slices import exact_numerator


def main():
    original=json.loads(Path('results/enrichment.json').read_text())
    primary=[r for r in original['rows'] if r['k']==12 and r['fraction']==.005]
    path=Path('data/danenberg/SingleCells.csv');expected=original['source']['files']['SingleCells.csv']['sha256']
    if hashlib.file_digest(path.open('rb'),'sha256').hexdigest()!=expected:raise ValueError('Source checksum changed')
    pairs={(r['patient'],r['image']) for r in primary};parts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);mask=np.array([(p,i) in pairs for p,i in zip(clean.metabric_id,clean.ImageNumber)])
        parts.append(clean.loc[mask])
    frame=pd.concat(parts);report={'protocol_commit':'217a54d','fresh_confirmation':False,
        'original_gate_changed':False,'source_sha256':expected,'rows':[]}
    for row in primary:
        clean=frame[(frame.metabric_id==row['patient'])&(frame.ImageNumber==row['image'])].sort_values('ObjectNumber')
        raw=selected_labels(clean).astype(int);xy=clean[['Location_Center_X','Location_Center_Y']].to_numpy();edges=knn_edges(xy,12)
        if np.any(raw.sum(axis=1)>1):raise ValueError('Source mapping is not mutually exclusive')
        if hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest()!=row['edges_sha256']:raise ValueError('Original graph changed')
        starts=[raw];old_overlap=[]
        for old in row['bounds']['witnesses']:
            x=raw.copy();x[old['a_flips'],0]^=1;x[old['b_flips'],1]^=1
            old_overlap.append(int(np.sum(x.sum(axis=1)>1)));starts.append(repair_overlap(raw,x))
        chosen=[]
        for maximize in (False,True):
            candidates=[categorical_witness(edges,raw,row['budget_each'],row['budget_each'],maximize,start=x) for x in starts]
            chosen.append((max if maximize else min)(candidates,key=lambda z:z['numerator']))
        worlds=[]
        for witness in chosen:
            x=raw.copy();x[witness['a_flips'],0]^=1;x[witness['b_flips'],1]^=1
            assert np.all(x.sum(axis=1)<=1) and np.all(np.sum(x!=raw,axis=0)<=row['budget_each'])
            assert exact_numerator(edges,x)==(witness['numerator'],witness['denominator']);worlds.append(x)
        opposite=chosen[0]['numerator']<0<chosen[1]['numerator']
        normalized=(xy-xy.min(axis=0))/np.maximum(xy.max(axis=0)-xy.min(axis=0),1e-12)
        grid=np.minimum((normalized*4).astype(int),3);fields=grid[:,0]*4+grid[:,1]
        result={'patient':row['patient'],'image':row['image'],'cells':len(raw),'budget_each':row['budget_each'],
            'original_witness_overlap_counts':old_overlap,'witnesses':chosen,'opposite_signs_proved':opposite,
            'pair_review_lower_bounds':[review_pair_bound(*worlds,fields,e) for e in (0,1)] if opposite else [],
            'cell_ids':clean.ObjectNumber.astype(int).tolist(),'fields':fields.tolist(),
            'raw_categories':(raw[:,0]+2*raw[:,1]).tolist(),'edges_sha256':row['edges_sha256']}
        report['rows'].append(result);print(row['patient'],opposite,[w['numerator'] for w in chosen],old_overlap,flush=True)
    report['summary']={'patients':len(primary),'categorical_opposite_sign_pairs':sum(r['opposite_signs_proved'] for r in report['rows']),
        'new_sign_certificates':0,'scope':'Feasible ambiguity witnesses only; local-search failure is inconclusive.'}
    Path('results/categorical_identifiability.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['summary'],flush=True)

if __name__=='__main__':main()
