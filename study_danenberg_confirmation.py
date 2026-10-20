"""Independent IMC tissue/platform confirmation of frozen conditional graph bounds."""
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from danenberg_confirmation import COLUMNS,retain,choose_images,selected_labels,knn_edges
from general_marker_bounds import fast_bounds
from measurement_error import degree_error_margin,pair_target

ROOT=Path(__file__).resolve().parent


def main():
    folder=ROOT/'data/danenberg';path=folder/'SingleCells.csv'
    manifest=json.loads((folder/'fetch_manifest.json').read_text())
    verified=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
    if verified!=manifest['files']['SingleCells.csv']['sha256']:raise ValueError('Source hash changed')
    counts=[];total=0;retained=0
    # Eligibility/selection precedes all target and bound evaluation.
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        total+=len(chunk);clean=retain(chunk);retained+=len(clean)
        counts.append(clean.groupby(['metabric_id','ImageNumber']).size().rename('cells').reset_index())
    table=pd.concat(counts).groupby(['metabric_id','ImageNumber'],as_index=False)['cells'].sum()
    selected=choose_images(table)
    if len(selected)!=6:raise ValueError('Insufficient eligible independent patients')
    pairs={(x['metabric_id'],x['ImageNumber']) for x in selected};parts=[]
    for chunk in pd.read_csv(path,usecols=COLUMNS,chunksize=100000):
        clean=retain(chunk);mask=np.array([(p,i) in pairs for p,i in zip(clean['metabric_id'],clean['ImageNumber'])],dtype=bool)
        parts.append(clean.loc[mask])
    cells=pd.concat(parts);report={'scope':'Independent cohort, source-annotation sensitivity; not biological label truth or assay validation.',
                                  'protocol_commit':'dd3041a','frozen_algorithm_commit':'9d8c6da','source':manifest,
                                  'rows_read':total,'retained_rows_all_images':retained,'excluded_rows_all_images':total-retained,
                                  'selected_images':selected,'rows':[]}
    for entry in selected:
        patient,image=entry['metabric_id'],entry['ImageNumber']
        frame=cells[(cells['metabric_id']==patient)&(cells['ImageNumber']==image)].sort_values('ObjectNumber')
        if frame.duplicated(['ImageNumber','ObjectNumber']).any():raise ValueError('Duplicate cell IDs')
        xy=frame[['Location_Center_X','Location_Center_Y']].to_numpy();labels=selected_labels(frame);n=len(frame)
        identity=np.c_[frame['ObjectNumber'].to_numpy(),xy,labels.astype(int)]
        for k in (6,12,24):
            edges=knn_edges(xy,k);target=pair_target(edges,labels)
            for fraction in (.001,.005):
                budget=int(np.ceil(n*fraction));start=time.perf_counter();bound=fast_bounds(edges,labels,budget,budget)
                elapsed=time.perf_counter()-start;margin=degree_error_margin(edges,n,budget,budget)
                outer=[max(0.,target-margin),min(1.,target+margin)]
                row={'patient':patient,'image':int(image),'cells':n,'a_positives':int(labels[:,0].sum()),'b_positives':int(labels[:,1].sum()),
                     'k':k,'edges':len(edges),'fraction_each_marker':fraction,'budget_each_marker':budget,'source_target':target,
                     'bounds':bound,'degree_outer_bounds':outer,'bound_seconds':elapsed,
                     'cells_sha256':hashlib.sha256(identity.tobytes()).hexdigest(),'edges_sha256':hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest()}
                report['rows'].append(row)
                print(patient,image,n,k,fraction,'gap pp',np.array(bound['endpoint_gaps'])*100,flush=True)
        (ROOT/'results/danenberg_confirmation.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
