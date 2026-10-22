"""Post-screen diagnostic, with exact integer signs for all feasible witnesses."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from danenberg_confirmation import COLUMNS,retain,selected_labels,knn_edges
from contrast_count_slices import exact_numerator,maximize_by_counts
from enrichment_bounds import coupled_contrast_bounds

if __name__=='__main__':
    original=json.loads(Path('results/enrichment.json').read_text())
    primary=[r for r in original['rows'] if r['k']==12 and r['fraction']==.005]
    pairs={(r['patient'],r['image']) for r in primary};parts=[]
    for chunk in pd.read_csv('data/danenberg/SingleCells.csv',usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);mask=np.array([(p,i) in pairs for p,i in zip(clean['metabric_id'],clean['ImageNumber'])]);parts.append(clean.loc[mask])
    frame=pd.concat(parts);report={'refinement_protocol_commit':'2310d16','primary_gate_changed':False,'rows':[]}
    previous=json.loads(Path('results/enrichment_refinement.json').read_text()) if Path('results/enrichment_refinement.json').exists() else None
    for row in primary:
        clean=frame[(frame['metabric_id']==row['patient'])&(frame['ImageNumber']==row['image'])].sort_values('ObjectNumber')
        labels=selected_labels(clean).astype(int);edges=knn_edges(clean[['Location_Center_X','Location_Center_Y']].to_numpy(),12)
        witnesses=[]
        for witness in row['bounds']['witnesses']:
            x=labels.copy();x[witness['a_flips'],0]^=1;x[witness['b_flips'],1]^=1
            assert np.sum(x[:,0]!=labels[:,0])<=row['budget_each'] and np.sum(x[:,1]!=labels[:,1])<=row['budget_each']
            num,den=exact_numerator(edges,x)
            assert abs(num/den-witness['value'])<1e-12
            witnesses.append({'numerator':num,'denominator':den,'a_flips':witness['a_flips'],'b_flips':witness['b_flips']})
        result={'patient':row['patient'],'exact_witnesses':witnesses,
                'opposite_signs_proved':witnesses[0]['numerator']<0<witnesses[1]['numerator']}
        if row['patient']=='MB-0045':
            if previous is not None:
                result['count_slice_reference']=next(r['count_slice_reference'] for r in previous['rows'] if r['patient']=='MB-0045')
            else:
                result['count_slice_reference']=maximize_by_counts(edges,labels,row['budget_each'],row['budget_each'],time_limit=1.)
            print('REFERENCE',result['count_slice_reference']['upper_bound'],flush=True)
        report['rows'].append(result)
    report['structural_development']=[]
    for row in original['rows']:
        clean=frame[(frame['metabric_id']==row['patient'])&(frame['ImageNumber']==row['image'])].sort_values('ObjectNumber')
        labels=selected_labels(clean).astype(int);edges=knn_edges(clean[['Location_Center_X','Location_Center_Y']].to_numpy(),row['k'])
        result=coupled_contrast_bounds(edges,labels,row['budget_each'],row['budget_each'],base=row['bounds'])
        width=result['outer_bounds'][1]-result['outer_bounds'][0]
        control=row['bounds']['separated_bounds'][1]-row['bounds']['separated_bounds'][0]
        report['structural_development'].append({'patient':row['patient'],'k':row['k'],'fraction':row['fraction'],
            'bounds':result,'width_reduction':1-width/control,'sign_certified':result['outer_bounds'][0]>0 or result['outer_bounds'][1]<0})
    selected=[r for r in report['structural_development'] if r['k']==12 and r['fraction']==.005]
    report['development_summary']={'median_width_reduction':float(np.median([r['width_reduction'] for r in selected])),
        'sign_certified':sum(r['sign_certified'] for r in selected),'fresh_confirmation':False}
    Path('results/enrichment_refinement.json').write_text(json.dumps(report,indent=2)+'\n')
    print('exact opposite-sign witnesses',sum(r['opposite_signs_proved'] for r in report['rows']),flush=True)
