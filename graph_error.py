"""Edge-budget propagation with uncertain graph size and image-level calibration."""
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree
from skimage.measure import label, regionprops_table


def pair_frequency_bounds(edges, count_low, count_high, deleted, added):
    """Enclose a [0,1]-marked edge mean after at most D deletions/A additions.

    Label bounds apply on the observed graph. Matched vertices must retain their
    marks; disagreements require a separate label uncertainty set. Empty target
    graphs have undefined means and conservatively return [0,1].
    """
    if any(not np.isfinite(v) or v < 0 for v in (edges, count_low, count_high, deleted, added)):
        raise ValueError('Finite nonnegative inputs required')
    if any(int(v) != v for v in (edges, deleted, added)) or not 0 <= count_low <= count_high <= edges:
        raise ValueError('Integer graph budgets and valid count bounds required')
    if edges == 0 or deleted >= edges:
        return 0., 1.
    low = max(0., count_low-deleted) / (edges-min(deleted,count_low)+added)
    high = (count_high+added) / (edges-min(deleted,edges-count_high)+added)
    return float(max(0., low)), float(min(1., high))


def conformal_upper(scores, alpha=.05):
    """Ordinary split-conformal order statistic, with the finite-sample +1."""
    values = np.asarray(scores, float)
    if values.ndim != 1 or not np.isfinite(values).all() or np.any(values<0) or not 0<alpha<1:
        raise ValueError('Nonnegative finite scores and alpha in (0,1) required')
    rank = math.ceil((len(values)+1)*(1-alpha))
    return float(np.sort(values)[rank-1]) if rank <= len(values) else math.inf


def decode_mask(image):
    """Authors' first-channel, 8-connected-component decoding for BBBC039."""
    if image.ndim == 3: image = image[:,:,0]
    return label(image, connectivity=2)


def match_instances(predicted, reference):
    """Maximum-cardinality one-to-one IoU>=.5 matching; labels must be contiguous."""
    if predicted.shape != reference.shape:
        raise ValueError('Mask shapes differ')
    npred, nref = int(predicted.max()), int(reference.max())
    if npred == 0 or nref == 0: return {}
    joint = np.bincount((predicted*(nref+1)+reference).ravel(), minlength=(npred+1)*(nref+1))
    joint = joint.reshape(npred+1,nref+1)
    union = np.bincount(predicted.ravel(),minlength=npred+1)[1:,None] + np.bincount(
        reference.ravel(),minlength=nref+1)[None,1:] - joint[1:,1:]
    iou = np.divide(joint[1:,1:],union,out=np.zeros_like(union,dtype=float),where=union>0)
    matching = maximum_bipartite_matching(csr_matrix(iou >= .5), perm_type='column')
    return {i+1:int(j+1) for i,j in enumerate(matching) if j>=0}


def centroids(mask):
    if not mask.max(): return np.empty((0,2))
    props = regionprops_table(mask, properties=('centroid',))
    return np.column_stack([props['centroid-0'],props['centroid-1']])


def radius_edges(points, radius=40):
    return set(cKDTree(points).query_pairs(radius)) if len(points)>1 else set()


def graph_errors(predicted, reference, radius=40):
    mapping = match_instances(predicted,reference)
    p, r = centroids(predicted), centroids(reference)
    pe, re = radius_edges(p,radius), radius_edges(r,radius)
    ids = [mapping.get(i+1, len(r)+i+1)-1 for i in range(len(p))]
    mapped = {tuple(sorted((ids[i],ids[j]))) for i,j in pe}
    deleted, added = len(mapped-re), len(re-mapped)
    denominator = max(1,len(pe))
    return {'predicted_nuclei':len(p),'reference_nuclei':len(r),'matched_nuclei':len(mapping),
            'instance_f1':2*len(mapping)/max(1,len(p)+len(r)),
            'predicted_edges':len(pe),'reference_edges':len(re),'deleted_edges':deleted,
            'added_edges':added,'score':max(deleted,added)/denominator}
