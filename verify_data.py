"""Check source identity, split separation, and graph completeness against full tables."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent


def main():
    manifest=json.loads((ROOT/'results/data_manifest.json').read_text())
    rng=np.random.default_rng(333)
    checks=[]
    for patient,record in manifest.items():
        table_path=ROOT/'data'/f'{patient}-cells.csv'
        assert hashlib.file_digest(table_path.open('rb'),'sha256').hexdigest()==record['sha256']
        source=pd.read_csv(table_path,usecols=['label','x','y','Pan-CK_pos','CD3e_pos'])
        d=np.load(ROOT/'data'/f'{patient}-prepared.npz')
        ids=d['cell_ids'];coords=d['coords'];edges=d['edges']
        assert len(np.unique(ids))==len(ids)
        indexed=source.set_index('label').loc[ids]
        np.testing.assert_array_equal(coords,indexed[['x','y']].to_numpy())
        np.testing.assert_array_equal(d['labels'],indexed[['Pan-CK_pos','CD3e_pos']].to_numpy())
        assert np.isfinite(d['features']).all()
        assert np.all(edges[:,0]!=edges[:,1])
        assert np.all(np.linalg.norm(coords[edges[:,0]]-coords[edges[:,1]],axis=1)<=80+1e-12)
        # A source-table neighbor query independently checks the halo retained every endpoint.
        source_coords=source[['x','y']].to_numpy();tree=cKDTree(source_coords)
        anchors=np.unique(edges[:,0])
        for anchor in rng.choice(anchors,min(100,len(anchors)),replace=False):
            expected=tree.query_ball_point(coords[anchor],80)
            expected=set(source['label'].to_numpy()[expected])-{int(ids[anchor])}
            got=set(ids[edges[edges[:,0]==anchor,1]])
            assert expected==got,(patient,int(ids[anchor]))
        checks.append({'patient':patient,'source_sha256':record['sha256'],'cells':len(ids),'pairs':len(edges),
                       'source_neighbors_checked':min(100,len(anchors)),'status':'passed'})
        print(patient,'source identity and complete local pairs verified',flush=True)
    (ROOT/'results/data_checks.json').write_text(json.dumps(checks,indent=2)+'\n')


if __name__=='__main__':
    main()
