"""Independent all-pairs and rational check of the public expert-reference screen."""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd

base=Path(__file__).resolve().parent;source=base/'data/nucls/unbiased_anchors.csv'
result=json.loads((base/'results/nucls_pair_references.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest()==result['source_sha256']
x=pd.read_csv(source);classes={'tumor','mitotic_figure','fibroblast','lymphocyte','plasma_cell','macrophage','neutrophil',
'eosinophil','apoptotic_body','vascular_endothelium','myoepithelium','ductal_epithelium'}
experts=['SP.2','SP.3'];patient_values={};checks=[]
for row in result['fields']:
    f=x[x.fovname==row['fov']]
    keep=f[experts].isin(classes).all(axis=1)&np.isfinite(f[['xmin_relative','ymin_relative','xmax_relative','ymax_relative']]).all(axis=1)
    f=f[keep].sort_values('anchor_id');n=len(f);assert n==row['common_detected_labeled'];assert (n>=20)==row['eligible']
    if n<20:continue
    assert f.anchor_id.tolist()==row['anchor_ids']
    xy=f[['xmin_relative','ymin_relative']].to_numpy()/2+f[['xmax_relative','ymax_relative']].to_numpy()/2
    edges=[];ids=np.arange(n)
    for i in range(n):
        delta=xy-xy[i];distance=(delta*delta).sum(axis=1);order=np.lexsort((ids,distance));neighbors=order[order!=i][:6]
        edges.extend((i,int(j)) for j in neighbors)
    edges=np.array(edges,dtype='<i8');assert hashlib.sha256(edges.tobytes()).hexdigest()==row['edges_sha256']
    m=len(edges);u,v=edges.T;left=np.bincount(u,minlength=n);right=np.bincount(v,minlength=n);mapped=[]
    patient_values.setdefault(row['patient'],{e:[] for e in experts})
    for expert in experts:
        a=(f[expert]=='tumor').astype(int).to_numpy();b=(f[expert]=='lymphocyte').astype(int).to_numpy();mapped.append(np.column_stack([a,b]))
        contacts=sum(int(a[i])*int(b[j]) for i,j in edges)
        numerator=m*contacts-int(left@a)*int(right@b)
        assert numerator==row['experts'][expert]['numerator'] and m*m==row['experts'][expert]['denominator']
        patient_values[row['patient']][expert].append((m,Fraction(numerator,m*m)))
    pair=row['pairs'][0];distance=int(np.sum(np.any(mapped[0]!=mapped[1],axis=1)))
    assert distance==pair['mapped_category_disagreements']
    marker_distance=np.sum(mapped[0]!=mapped[1],axis=0).astype(int).tolist();assert marker_distance==pair['per_marker_disagreements']
    assert (distance<=2)==pair['one_wrong_mapped_cell_each_compatible']
    assert all(v<=2*int(np.ceil(.005*n)) for v in marker_distance)==pair['half_percent_each_marker_each_compatible']
    checks.append({'patient':row['patient'],'cells':n,'mapped_reference_disagreements':distance,'graph_and_integer_readout_verified':True})
for patient in result['patients']:
    for expert,value in patient['experts'].items():
        pairs=patient_values.get(patient['patient'],{}).get(expert,[])
        if not pairs:assert value is None;continue
        truth=sum((m*q for m,q in pairs),Fraction(0))/sum(m for m,q in pairs)
        assert (truth.numerator,truth.denominator)==(value['numerator'],value['denominator'])
old=json.loads((base/'results/nucls_references.json').read_text());assert old['summary']['eligible_fields']==0
assert (x['SP.1']=='DidNotAnnotateFOV').all()
verification={'fields':checks,'initial_three_expert_screen_unevaluable':True,'all_pair_readouts_verified':True,
              'scope':'Independent computational replay of external measured annotations, not biological truth or independent laboratory replication.'}
(base/'results/nucls_reference_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
print(json.dumps(verification,indent=2))
