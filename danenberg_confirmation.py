"""Frozen independent-cohort annotation mapping and per-image topological graphs."""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

COLUMNS=['ImageNumber','ObjectNumber','metabric_id','cellPhenotype','is_epithelial','is_tumour','is_hotAggregate',
         'Location_Center_X','Location_Center_Y']
TYPES={'CD4^{+} T cells','CD8^{+} T cells','T_{Reg} & T_{Ex}'}
AMBIGUOUS='CD4^{+} T cells & APCs'


def retain(frame):
    finite=np.isfinite(frame[['Location_Center_X','Location_Center_Y']]).all(axis=1)
    valid=(frame['is_tumour']==1)&(frame['is_hotAggregate']==0)&finite
    valid&=frame[['metabric_id','cellPhenotype','is_epithelial','ImageNumber','ObjectNumber']].notna().all(axis=1)
    valid&=frame['cellPhenotype']!=AMBIGUOUS
    result=frame.loc[valid].copy()
    if not result['is_epithelial'].isin([0,1]).all():raise ValueError('Unexpected epithelial annotation')
    return result


def choose_images(counts):
    eligible=counts[counts['cells']>=500].sort_values(['metabric_id','ImageNumber'])
    return eligible.groupby('metabric_id',sort=True).head(1).head(6).to_dict('records')


def selected_labels(frame):
    return np.stack([frame['is_epithelial'].to_numpy().astype(bool),frame['cellPhenotype'].isin(TYPES).to_numpy()],axis=1)


def knn_edges(coordinates,k):
    """Deterministic distance/index ties; no cross-image inputs are permitted."""
    xy=np.asarray(coordinates,float);n=len(xy)
    if k>=n or k<1:raise ValueError('Need 1 <= k < number of cells')
    if len(np.unique(xy,axis=0))!=n:raise ValueError('Duplicate coordinates require investigation')
    tree=cKDTree(xy);distances,_=tree.query(xy,k=k+1);result=[]
    for i in range(n):
        candidates=np.array(tree.query_ball_point(xy[i],np.nextafter(distances[i,-1],np.inf)),int)
        candidates=candidates[candidates!=i];squared=np.sum((xy[candidates]-xy[i])**2,axis=1)
        ordered=candidates[np.lexsort((candidates,squared))][:k]
        if len(ordered)!=k:raise AssertionError('Missing graph neighbours')
        result.extend((i,int(j)) for j in ordered)
    return np.array(result,int)
