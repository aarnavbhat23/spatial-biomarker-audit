"""Frozen external expert-reference sensitivity screen; no annotation is assumed true."""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import re
import ssl
import urllib.request
import numpy as np
import pandas as pd
from danenberg_confirmation import knn_edges
from contrast_count_slices import exact_numerator

FILE_ID='1RTyfmkZW8um3qmWqtu2gbu3Sgql-DQQj'
SHA='07edc0db091a51f90e58fcc18332c347f5f62a9356320c9aba7377e666469ff6'
EXPERTS=['SP.1','SP.2','SP.3']
BIOLOGICAL={'tumor','mitotic_figure','fibroblast','lymphocyte','plasma_cell','macrophage','neutrophil',
            'eosinophil','apoptotic_body','vascular_endothelium','myoepithelium','ductal_epithelium'}


def main(experts=EXPERTS,protocol='d666fa1',filename='nucls_references.json'):
    path=Path('data/nucls/unbiased_anchors.csv');path.parent.mkdir(parents=True,exist_ok=True)
    url='https://drive.google.com/uc?export=download&id='+FILE_ID
    if not path.exists():
        ca=Path('/etc/ssl/cert.pem');ctx=ssl.create_default_context(cafile=str(ca) if ca.exists() else None)
        blob=urllib.request.urlopen(url,context=ctx,timeout=90).read()
        if hashlib.sha256(blob).hexdigest()!=SHA:raise ValueError('Source checksum mismatch')
        path.write_bytes(blob)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHA:raise ValueError('Source checksum changed')
    data=pd.read_csv(path);encodings=sorted(set(data[experts].stack().dropna().astype(str)))
    if not set(encodings)<=BIOLOGICAL|{'undetected','unlabeled','DidNotAnnotateFOV'}:raise ValueError('Unexpected source encodings: '+repr(encodings))
    if data.anchor_id.duplicated().any():raise ValueError('Duplicate source anchor IDs')
    report={'protocol_commit':protocol,'experts':experts,'input_note':'NUCLS_INPUT_NOTE.md','source_url':url,'source_sha256':SHA,'license':'CC0 1.0',
            'source_doi':'10.1093/gigascience/giac037','source_rows':len(data),'class_encodings':encodings,
            'scope':'Independent external measured reference-choice study; tumor/lymphocyte identity, not protein truth.',
            'fields':[],'patients':[]};grouped={}
    for fov,f in data.groupby('fovname',sort=True):
        patient=re.search(r'TCGA-[^-]+-[^-]+',fov)[0]
        bounds=['xmin_relative','ymin_relative','xmax_relative','ymax_relative']
        keep=np.isfinite(f[bounds]).all(axis=1)&f[experts].isin(BIOLOGICAL).all(axis=1)
        frame=f[keep].sort_values('anchor_id');n=len(frame)
        row={'fov':fov,'patient':patient,'source_anchors':len(f),'common_detected_labeled':n,'excluded':len(f)-n,
             'eligible':n>=20};report['fields'].append(row)
        if n<20:continue
        xy=np.column_stack([(frame.xmin_relative+frame.xmax_relative)/2,(frame.ymin_relative+frame.ymax_relative)/2])
        edges=knn_edges(xy,6);m=len(edges);worlds=[];row['experts']={};row['pairs']=[]
        for expert in experts:
            labels=np.column_stack([frame[expert]=='tumor',frame[expert]=='lymphocyte']).astype(int);worlds.append(labels)
            numerator,denominator=exact_numerator(edges,labels)
            row['experts'][expert]={'numerator':numerator,'denominator':denominator,'contrast':numerator/denominator,
                'a_positives':int(labels[:,0].sum()),'b_positives':int(labels[:,1].sum())}
        for a,b in itertools.combinations(range(len(experts)),2):
            distance=int(np.sum(frame[experts[a]].to_numpy()!=frame[experts[b]].to_numpy()))
            mapped_distance=int(np.sum(np.any(worlds[a]!=worlds[b],axis=1)))
            per_marker=np.sum(worlds[a]!=worlds[b],axis=0).astype(int).tolist();budget=int(np.ceil(.005*n))
            row['pairs'].append({'experts':[experts[a],experts[b]],'full_class_disagreements':distance,
                'mapped_category_disagreements':mapped_distance,'per_marker_disagreements':per_marker,
                'one_wrong_mapped_cell_each_compatible':mapped_distance<=2,
                'half_percent_each_marker_each_compatible':all(d<=2*budget for d in per_marker)})
        row['edges']=m;row['anchor_ids']=frame.anchor_id.tolist();row['edges_sha256']=hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest()
        grouped.setdefault(patient,[]).append(row)
    for patient in sorted({re.search(r'TCGA-[^-]+-[^-]+',fov)[0] for fov in data.fovname.unique()}):
        fields=grouped.get(patient,[]);total=sum(r['edges'] for r in fields);values={}
        for expert in experts:
            value=sum((r['edges']*Fraction(r['experts'][expert]['numerator'],r['experts'][expert]['denominator']) for r in fields),Fraction(0))
            value=value/total if total else None
            values[expert]={'numerator':value.numerator,'denominator':value.denominator,'contrast':float(value)} if value is not None else None
        theta=[v['contrast'] for v in values.values() if v is not None]
        report['patients'].append({'patient':patient,'eligible_fields':len(fields),'experts':values,
            'observed_reference_sign_stable':bool(theta and (min(theta)>0 or max(theta)<0)),
            'opposite_reference_signs':bool(theta and min(theta)<0<max(theta))})
    fields=[r for r in report['fields'] if r['eligible']];pairs=[p for r in fields for p in r['pairs']]
    report['summary']={'evaluable':bool(fields),'patients_with_eligible_fields':len(grouped),'patients':len(report['patients']),'eligible_fields':len(fields),
        'common_cells_in_eligible_fields':sum(r['common_detected_labeled'] for r in fields),
        'reference_sign_stable_patients':sum(r['observed_reference_sign_stable'] for r in report['patients']),
        'opposite_reference_sign_patients':sum(r['opposite_reference_signs'] for r in report['patients']),
        'expert_field_pairs':len(pairs),'one_wrong_cell_pair_compatibility_failures':sum(not p['one_wrong_mapped_cell_each_compatible'] for p in pairs),
        'half_percent_marker_pair_compatibility_failures':sum(not p['half_percent_each_marker_each_compatible'] for p in pairs)}
    (Path('results')/filename).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['summary'],indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--available-pair',action='store_true');args=parser.parse_args()
    if args.available_pair:main(['SP.2','SP.3'],'5f30154','nucls_pair_references.json')
    else:main()
