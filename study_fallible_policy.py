"""Fresh patients and frozen corruption generators; policies never see hidden calls."""
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from danenberg_confirmation import COLUMNS,retain,selected_labels,knn_edges
from fallible_policy import fields_from_coordinates,run_policy
from measurement_error import pair_target

ROOT=Path(__file__).resolve().parent


def fresh_panels():
    path=ROOT/'data/danenberg/SingleCells.csv';manifest=json.loads((path.parent/'fetch_manifest.json').read_text())
    with path.open('rb') as source:digest=hashlib.file_digest(source,'sha256').hexdigest()
    if digest!=manifest['files']['SingleCells.csv']['sha256']:raise ValueError('Source changed')
    counts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        counts.append(retain(chunk).groupby(['metabric_id','ImageNumber']).size().rename('cells').reset_index())
    table=pd.concat(counts).groupby(['metabric_id','ImageNumber'],as_index=False)['cells'].sum()
    table=table[table['cells']>=500].sort_values(['metabric_id','ImageNumber']).groupby('metabric_id',sort=True).head(1).iloc[6:12]
    if len(table)!=6:raise ValueError('Insufficient fresh patients')
    selected=table.to_dict('records');pairs={(x['metabric_id'],x['ImageNumber']) for x in selected};parts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);mask=np.array([(p,i) in pairs for p,i in zip(clean['metabric_id'],clean['ImageNumber'])],bool)
        parts.append(clean.loc[mask])
    cells=pd.concat(parts)
    for entry in selected:
        frame=cells[(cells['metabric_id']==entry['metabric_id'])&(cells['ImageNumber']==entry['ImageNumber'])].sort_values('ObjectNumber')
        if frame.duplicated(['ImageNumber','ObjectNumber']).any():raise ValueError('Duplicate IDs')
        xy=frame[['Location_Center_X','Location_Center_Y']].to_numpy();labels=selected_labels(frame).astype(int)
        yield entry,knn_edges(xy,12),labels,fields_from_coordinates(xy),digest


def main():
    report={'protocol_commit':'1401618','scope':'Post-census fallible-review emulation, source annotations are proxy truth; hypothetical costs and budgets.',
            'rows':[],'selected':[]}
    for patient_index,(entry,edges,truth,fields,digest) in enumerate(fresh_panels()):
        report['selected'].append(entry);report['source_sha256']=digest;n=len(truth);budget=int(np.ceil(.005*n));u,v=edges.T
        influence=[np.bincount(u,weights=truth[v,1],minlength=n),np.bincount(v,weights=truth[u,0],minlength=n)]
        for repeat in range(2):
            seed=93000+1000*patient_index+repeat*10
            reference=truth.copy();ref_rng=np.random.default_rng(seed+100000)
            for k in (0,1):reference[int(ref_rng.integers(n)),k]^=1
            for profile in ('uniform','influential'):
                rng=np.random.default_rng(seed);raw=truth.copy()
                for k in (0,1):
                    ids=rng.choice(n,budget,replace=False) if profile=='uniform' else np.argsort(-influence[k],kind='stable')[:budget]
                    raw[ids,k]^=1
                for policy in ('uniform','degree','influence','lookahead','witness'):
                    start=time.perf_counter()
                    # Only this observation callback can access hidden reference responses.
                    result=run_policy(edges,raw,(budget,budget),(1,1),fields,lambda ids:reference[ids],policy,seed+200000)
                    target=pair_target(edges,truth)
                    valid=not result['empty'] and all(t['outer_bounds'][0]-1e-12<=target<=t['outer_bounds'][1]+1e-12 for t in result['trace'])
                    row={'patient':entry['metabric_id'],'image':int(entry['ImageNumber']),'cells':n,'a_positives':int(truth[:,0].sum()),
                         'b_positives':int(truth[:,1].sum()),'budget_each_marker':budget,'reference_error_budget_each':1,
                         'repeat':repeat,'seed':seed,'profile':profile,'policy':policy,'target':target,'valid_every_step':valid,
                         'result':result,'seconds':time.perf_counter()-start,
                         'graph_sha256':hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest(),
                         'raw_sha256':hashlib.sha256(raw.astype('i1').tobytes()).hexdigest(),
                         'reference_sha256':hashlib.sha256(reference.astype('i1').tobytes()).hexdigest()}
                    report['rows'].append(row)
                    (ROOT/'results/fallible_policy_study.json').write_text(json.dumps(report,indent=2)+'\n')
                    print(entry['metabric_id'],repeat,profile,policy,'valid',valid,'attained',result.get('attained'),
                          'cost',result.get('quoted_cost'),'seconds',row['seconds'],flush=True)


if __name__=='__main__':main()
