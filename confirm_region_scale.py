"""Execute the frozen acquisition comparison on four additional patients."""
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import fetch_data
from study_region_scale import evaluate

ROOT=Path(__file__).resolve().parent


def prepare_panel(patient):
    output=ROOT/'data'/f'{patient}-acquisition.npz'
    record=ROOT/'data'/f'{patient}-acquisition-preparation.json'
    if output.exists() and record.exists():
        return json.loads(record.read_text())
    table=pd.read_csv(ROOT/'data'/f'{patient}-cells.csv',usecols=['label','x','y','Pan-CK_pos','CD3e_pos'],
                      dtype={'label':'int32','x':'int32','y':'int32','Pan-CK_pos':'bool','CD3e_pos':'bool'})
    coords=table[['x','y']].to_numpy()
    candidates,counts=np.unique(coords//2048,axis=0,return_counts=True)
    lower=coords.min(0);upper=coords.max(0)
    allowed=(counts>=100)&((candidates*2048-80)>=lower).all(1)&(((candidates+1)*2048+80)<=upper).all(1)
    candidates=candidates[allowed]
    rng=np.random.default_rng(20261023+int(patient[3:]))
    tiles=candidates[rng.choice(len(candidates),32,replace=False)]
    tiles=tiles[np.lexsort((tiles[:,1],tiles[:,0]))]
    rectangles=np.c_[tiles*2048-80,(tiles+1)*2048+80]
    include=np.zeros(len(coords),dtype=bool)
    for x0,y0,x1,y1 in rectangles:
        include|=(coords[:,0]>=x0)&(coords[:,0]<x1)&(coords[:,1]>=y0)&(coords[:,1]<y1)
    selected=np.flatnonzero(include);xy=coords[selected]
    tile_map={tuple(t):i for i,t in enumerate(tiles)}
    owner=np.array([tile_map.get(tuple(t),-1) for t in xy//2048],dtype=np.int32)
    anchors=np.flatnonzero(owner>=0)
    neighbors=cKDTree(xy).query_ball_point(xy[anchors],80)
    left=np.repeat(anchors,[len(x) for x in neighbors]);right=np.concatenate(neighbors)
    keep=left!=right;edges=np.stack([left[keep],right[keep]],axis=1).astype(np.int32)
    # Independent checks against the complete source table, not just cached halos.
    source_tree=cKDTree(coords)
    source_ids=table['label'].to_numpy()
    chosen=rng.choice(anchors,min(100,len(anchors)),replace=False)
    for anchor in chosen:
        full=source_tree.query_ball_point(xy[anchor],80)
        local=cKDTree(xy).query_ball_point(xy[anchor],80)
        if set(source_ids[full])!=set(source_ids[selected][local]):
            raise AssertionError('Incomplete source neighborhood')
    if not np.all(np.sum((xy[edges[:,0]]-xy[edges[:,1]])**2,axis=1)<=80**2):
        raise AssertionError('Invalid pair distance')
    np.savez_compressed(output,coords=xy,labels=table[['Pan-CK_pos','CD3e_pos']].to_numpy()[selected],
                        cell_ids=source_ids[selected],edges=edges,edge_tiles=owner[edges[:,0]],
                        core_tiles=tiles,rectangles=rectangles)
    result={'patient':patient,'source_cells':len(table),'panel_cells':len(xy),'directed_pairs':len(edges),
            'candidate_cores':len(candidates),'selected_cores':tiles.tolist(),
            'complete_source_neighborhood_checks':len(chosen),
            'prepared_sha256':hashlib.file_digest(output.open('rb'),'sha256').hexdigest(),
            'scope':'Coordinate-only selected panel; source fluorescence segmentation; no image-model evaluation.'}
    record.write_text(json.dumps(result,indent=2)+'\n')
    print(patient,'prepared',len(xy),'cells',len(edges),'pairs',flush=True)
    return result


def main():
    fetch_data.SOURCES=json.loads((ROOT/'acquisition_sources.json').read_text())
    protocol_commit=subprocess.check_output(['git','log','-1','--format=%H','--','ACQUISITION_CONFIRMATION_PROTOCOL.md'],cwd=ROOT,text=True).strip()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        records=dict(pool.map(fetch_data.fetch_table,fetch_data.SOURCES['cell_tables'].items()))
    (ROOT/'data/acquisition_fresh_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    report={'protocol_commit':protocol_commit,'target_width':.04,'alpha':.05,'repetitions':100,
            'sources':records,'patients':{},'scope':'New-patient confirmation of a classical observation protocol, conditional on known source graph.'}
    for patient in fetch_data.SOURCES['cell_tables']:
        preparation=prepare_panel(patient)
        data=dict(np.load(ROOT/'data'/f'{patient}-acquisition.npz'))
        rows=evaluate(data,seed=20261024+int(patient[3:]),repetitions=100)
        baseline=next(row for row in rows if row['core_native_pixels']==2048)
        # Summed rectangle area charges all repeated overlap between distinct ROIs.
        for row in rows:
            side=row['core_native_pixels']+160
            row['mean_summed_area_native_pixels']=side**2*row['mean_unique_regions']
            row['union_area_ratio_to_coarse']=row['mean_area_native_pixels']/baseline['mean_area_native_pixels']
            extra=row['mean_unique_regions']-baseline['mean_unique_regions']
            row['overhead_break_even_native_pixel_area']=(baseline['mean_area_native_pixels']-row['mean_area_native_pixels'])/extra if extra>0 else None
        selected=next(row for row in rows if row['core_native_pixels']==128)
        passed=selected['stop_fraction']>=.95 and selected['mean_area_native_pixels']<baseline['mean_area_native_pixels']
        report['patients'][patient]={'preparation':preparation,'scales':rows,'primary_passed':bool(passed)}
        report['primary_passed']=all(v['primary_passed'] for v in report['patients'].values())
        (ROOT/'results/region_scale_confirmation.json').write_text(json.dumps(report,indent=2)+'\n')
        print(patient,'primary passed',passed,flush=True)


if __name__=='__main__':main()
