"""Frozen fresh-patient sensitivity screen for the changed enrichment estimand."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from danenberg_confirmation import COLUMNS,retain,selected_labels,knn_edges
from enrichment_bounds import contrast_bounds


def main():
    folder=Path('data/danenberg');path=folder/'SingleCells.csv'
    manifest=json.loads((folder/'fetch_manifest.json').read_text())
    if hashlib.file_digest(path.open('rb'),'sha256').hexdigest()!=manifest['files']['SingleCells.csv']['sha256']:
        raise ValueError('Source checksum changed')
    counts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);counts.append(clean.groupby(['metabric_id','ImageNumber']).size().rename('cells').reset_index())
    table=pd.concat(counts).groupby(['metabric_id','ImageNumber'],as_index=False)['cells'].sum()
    eligible=table[table['cells']>=500].sort_values(['metabric_id','ImageNumber']).groupby('metabric_id',sort=True).head(1)
    selected=eligible.iloc[12:18].to_dict('records')
    if len(selected)!=6:raise ValueError('Need six previously unused eligible patients')
    pairs={(x['metabric_id'],x['ImageNumber']) for x in selected};parts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);mask=np.array([(p,i) in pairs for p,i in zip(clean['metabric_id'],clean['ImageNumber'])])
        parts.append(clean.loc[mask])
    cells=pd.concat(parts);report={'protocol_commit':'d3fec2d','source':manifest,
        'scope':'Changed enrichment estimand, supplied label error budgets, fixed measured graph; no assay-error estimate.',
        'selected_images':selected,'rows':[]}
    for entry in selected:
        p,image=entry['metabric_id'],entry['ImageNumber']
        frame=cells[(cells['metabric_id']==p)&(cells['ImageNumber']==image)].sort_values('ObjectNumber')
        if frame.duplicated(['ImageNumber','ObjectNumber']).any():raise ValueError('Duplicate cell identifiers')
        labels=selected_labels(frame);n=len(frame);xy=frame[['Location_Center_X','Location_Center_Y']].to_numpy()
        for k in (6,12,24):
            edges=knn_edges(xy,k)
            for fraction in (.001,.005):
                budget=int(np.ceil(fraction*n));start=time.perf_counter()
                bounds=contrast_bounds(edges,labels,budget,budget)
                width=bounds['outer_bounds'][1]-bounds['outer_bounds'][0]
                control=bounds['separated_bounds'][1]-bounds['separated_bounds'][0]
                row={'patient':p,'image':int(image),'cells':n,'k':k,'edges':len(edges),
                     'marker_positives':labels.sum(axis=0).astype(int).tolist(),'fraction':fraction,'budget_each':budget,
                     'bounds':bounds,'width_reduction':1-width/control if control>0 else 0.,
                     'sign_certified':bounds['outer_bounds'][0]>0 or bounds['outer_bounds'][1]<0,
                     'seconds':time.perf_counter()-start,
                     'edges_sha256':hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest()}
                report['rows'].append(row)
                print(p,k,fraction,bounds['source_contrast'],bounds['outer_bounds'],row['width_reduction'],flush=True)
    primary=[r for r in report['rows'] if r['k']==12 and r['fraction']==.005]
    reduction=float(np.median([r['width_reduction'] for r in primary]));signs=sum(r['sign_certified'] for r in primary)
    report['summary']={'primary_patients':len(primary),'median_width_reduction':reduction,
                       'sign_certified':signs,'gate_passed':bool(reduction>=.25 and signs>=4)}
    Path('results/enrichment.json').write_text(json.dumps(report,indent=2)+'\n')
    print('SUMMARY',report['summary'],flush=True)

if __name__=='__main__':main()
