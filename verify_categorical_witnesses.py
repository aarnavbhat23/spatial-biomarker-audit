"""Independent source labels, all-pairs graph and exact integer witness checks."""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from collections import defaultdict
import numpy as np
import pandas as pd

SP=Path(__file__).resolve().parent
report={}
artifact=json.loads((SP/'results/categorical_identifiability.json').read_text());patients={r['patient'] for r in artifact['rows']}
columns=['ImageNumber','ObjectNumber','metabric_id','cellPhenotype','is_epithelial','is_tumour','is_hotAggregate','Location_Center_X','Location_Center_Y']
source=SP/'data/danenberg/SingleCells.csv';assert hashlib.file_digest(source.open('rb'),'sha256').hexdigest()==artifact['source_sha256']
parts=[]
for chunk in pd.read_csv(source,usecols=columns,chunksize=100000):parts.append(chunk[chunk.metabric_id.isin(patients)])
allcells=pd.concat(parts);checks=[]
for row in artifact['rows']:
    f=allcells[(allcells.metabric_id==row['patient'])&(allcells.ImageNumber==row['image'])]
    valid=(f.is_tumour==1)&(f.is_hotAggregate==0)&np.isfinite(f[['Location_Center_X','Location_Center_Y']]).all(axis=1)
    valid&=f[['metabric_id','cellPhenotype','is_epithelial','ImageNumber','ObjectNumber']].notna().all(axis=1)
    valid&=f.cellPhenotype!='CD4^{+} T cells & APCs';f=f[valid].sort_values('ObjectNumber');n=len(f)
    raw=np.column_stack([f.is_epithelial.astype(int),f.cellPhenotype.isin(['CD4^{+} T cells','CD8^{+} T cells','T_{Reg} & T_{Ex}']).astype(int)])
    assert f.ObjectNumber.astype(int).tolist()==row['cell_ids'];assert (raw[:,0]+2*raw[:,1]).tolist()==row['raw_categories']
    xy=f[['Location_Center_X','Location_Center_Y']].to_numpy();ids=np.arange(n);edges=[]
    # All-pairs sorting is independent of the production KD-tree construction.
    for i in range(n):
        squared=np.sum((xy-xy[i])**2,axis=1);order=np.lexsort((ids,squared));neighbors=order[order!=i][:12]
        edges.extend((i,int(j)) for j in neighbors)
    edges=np.array(edges,dtype='<i8');assert hashlib.sha256(edges.tobytes()).hexdigest()==row['edges_sha256']
    m=len(edges);out=np.bincount(edges[:,0],minlength=n);inc=np.bincount(edges[:,1],minlength=n);nums=[];worlds=[]
    for w in row['witnesses']:
        x=raw.copy();x[w['a_flips'],0]^=1;x[w['b_flips'],1]^=1
        assert np.isin(x,[0,1]).all() and np.all(x.sum(axis=1)<=1)
        assert np.all(np.sum(x!=raw,axis=0)<=row['budget_each'])
        contacts=sum(int(x[i,0])*int(x[j,1]) for i,j in edges)
        numerator=m*contacts-sum(int(a)*int(b) for a,b in zip(out,x[:,0]))*sum(int(a)*int(b) for a,b in zip(inc,x[:,1]))
        assert numerator==w['numerator'] and m*m==w['denominator'];nums.append(numerator);worlds.append(x)
    assert row['opposite_signs_proved']==(nums[0]<0<nums[1])
    grid=np.minimum((((xy-xy.min(axis=0))/np.maximum(xy.max(axis=0)-xy.min(axis=0),1e-12))*4).astype(int),3)
    fields=grid[:,0]*4+grid[:,1];assert fields.tolist()==row['fields'];different=np.any(worlds[0]!=worlds[1],axis=1)
    for b in row['pair_review_lower_bounds']:
        assert b['differing_cells']==int(different.sum());quota=2*b['error_cells']+1
        counts=sorted([int(np.sum(different & (fields==k))) for k in np.unique(fields)],reverse=True)
        minimum=next((j for j in range(1,len(counts)+1) if sum(counts[:j])>=quota),None)
        assert minimum==b['minimum_unit_cost_fields_for_pair']
    checks.append({'patient':row['patient'],'categorical_opposite_signs':row['opposite_signs_proved'],'integer_numerators':nums,
                   'source_labels_graph_witnesses_verified':True})
report['spatial']=checks
report['scope']='Separate all-pairs implementation; no third-party biological replication.'
Path('results/categorical_identifiability_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
