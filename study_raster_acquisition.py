"""Cost sensitivity to reusable fixed camera fields, not instrument validation."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from study_cell_reuse import cell_trial, grouped_indices, prepare_visibility

ROOT=Path(__file__).resolve().parent


def prepare_raster(data,core,field):
    prepared=prepare_visibility(data,core)
    rectangles=prepared['frame']['rectangles']
    expanded=np.c_[(rectangles[:,:2]//field)*field,
                   ((rectangles[:,2:]+field-1)//field)*field]
    tree=cKDTree(data['coords']);visible=[]
    for rectangle in expanded:
        center=(rectangle[:2]+rectangle[2:])/2
        radius=np.max(rectangle[2:]-rectangle[:2])/2
        ids=np.asarray(tree.query_ball_point(center,radius,p=np.inf),dtype=int)
        inside=np.all((data['coords'][ids]>=rectangle[:2])&(data['coords'][ids]<rectangle[2:]),axis=1)
        visible.append(ids[inside])
    prepared['frame']['rectangles']=expanded
    prepared['visible']=visible
    return prepared


def prepare_frame_pairs(data,field):
    """Partition edges by the one or two physical fields needed for their endpoints."""
    tiles,node_owner=np.unique(data['coords']//field,axis=0,return_inverse=True)
    endpoints=np.sort(node_owner[data['edges']],axis=1)
    pairs,owner=np.unique(endpoints,axis=0,return_inverse=True)
    rectangles=np.c_[tiles*field,(tiles+1)*field]
    node_groups=grouped_indices(node_owner,len(tiles))
    visible=[];physical=[]
    for first,second in pairs:
        ids=np.unique([first,second])
        visible.append(np.concatenate([node_groups[i] for i in ids]))
        physical.append(rectangles[ids])
    edges=data['edges'];edge_count=len(edges)
    incident=[ids%edge_count for ids in grouped_indices(np.r_[edges[:,0],edges[:,1]],len(data['coords']))]
    return {'frame':{'counts':np.bincount(owner),'owners':owner},
            'owned':grouped_indices(owner,len(pairs)),'incident':incident,'visible':visible,
            'physical_rectangles':physical}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--frame-pair',action='store_true')
    args=parser.parse_args()
    report={'scope':'Exploratory fixed-field caching sensitivity; hypothetical raster sizes, not validated instrument costs.',
            'seed':20261104,'target_width':.005,'repetitions':5,'patients':{}}
    if args.frame_pair:
        report=json.loads((ROOT/'results/raster_acquisition_study.json').read_text())
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        entries=[] if not args.frame_pair else [r for r in report['patients'][patient] if r['core_native_pixels'] is not None]
        for field in (512,1024,2048):
            for core in ((None,) if args.frame_pair else (128,2048)):
                prepared=prepare_frame_pairs(data,field) if core is None else prepare_raster(data,core,field)
                rows=[]
                for repeat in range(5):
                    row=cell_trial(data,prepared,20261104+1000*int(patient[3:])+repeat,.005)
                    count=row['area_native_pixels']/field**2
                    if abs(count-round(count))>1e-9:raise AssertionError('Raster union is not a whole number of fields')
                    row['physical_fields']=int(round(count));rows.append(row)
                entry={'field_native_pixels':field,'core_native_pixels':core,'query':'frame_pair' if core is None else 'rectangular_core','rows':rows,
                       'mean_physical_fields':float(np.mean([r['physical_fields'] for r in rows])),
                       'mean_logical_rounds':float(np.mean([r['draws'] for r in rows])),
                       'mean_area_native_pixels':float(np.mean([r['area_native_pixels'] for r in rows])),
                       'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                       'coverage_all_checkpoints':float(np.mean([r['covered_all_checkpoints'] for r in rows]))}
                entries.append(entry)
                print(patient,field,core,'fields',entry['mean_physical_fields'],'rounds',entry['mean_logical_rounds'],
                      'attain',entry['stop_fraction'],flush=True)
        report['patients'][patient]=entries
        (ROOT/'results/raster_acquisition_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
