"""Safe position-conditioned outer bounds under a finite pixel-location contract."""
from fractions import Fraction
from itertools import product
import numpy as np
from scipy.spatial import cKDTree
from enrichment_bounds import coupled_contrast_bounds


def _sum_bounds(low, high, raw, budget):
    def top(x):
        return int(np.sort(x)[-budget:].sum()) if budget else 0
    return [int(low[raw].sum())-top(low[raw]), int(high[raw].sum())+top(high[~raw])]


def _tv_bound(m, deletions, additions):
    if m == 0 or deletions >= m:
        return 1.
    return max(deletions/m, additions/(m-deletions+additions))


def _enlarge(interval, amount):
    return [max(-.25, interval[0]-amount-1e-12), min(.25, interval[1]+amount+1e-12)]


def position_bounds(positions, labels, radius, budget_a, budget_b):
    """All cells choose any supplied position; no probability or coverage is assigned."""
    positions = np.asarray(positions); labels = np.asarray(labels)
    if positions.ndim != 3 or positions.shape[2] != 2 or not np.issubdtype(positions.dtype, np.integer):
        raise ValueError('Need integer positions shaped cells x choices x 2')
    n, choices, _ = positions.shape
    if n < 2 or choices == 0 or labels.shape != (n,) or not np.isin(labels, [0, 1, 2]).all():
        raise ValueError('Need aligned categorical identities and position choices')
    if not isinstance(radius, (int, np.integer)) or radius <= 0:
        raise ValueError('Need a positive integer radius')
    if any(int(b) != b or b < 0 for b in (budget_a, budget_b)):
        raise ValueError('Need nonnegative integer budgets')
    ba, bb = min(int(budget_a), n), min(int(budget_b), n)
    raw = positions[:, 0]
    displacement = np.linalg.norm(positions-raw[:, None], axis=2).max()
    pairs = cKDTree(raw).query_pairs(radius+2*displacement+1e-10, output_type='ndarray')
    delta = positions[pairs[:, 0], :, None]-positions[pairs[:, 1], None, :]
    possible = np.sum(delta*delta, axis=-1) <= radius*radius
    keep = possible.any(axis=(1, 2)); pairs, possible = pairs[keep], possible[keep]
    raw_present = possible[:, 0, 0]
    must = np.zeros((n, choices), np.int64); may = must.copy()
    bmust = must.copy(); bmay = must.copy(); removed = must.copy(); added = must.copy()
    for center, neighbor, tables in ((pairs[:, 0], pairs[:, 1], possible),
                                     (pairs[:, 1], pairs[:, 0], possible.transpose(0, 2, 1))):
        guaranteed = tables.all(axis=2); allowed = tables.any(axis=2)
        np.add.at(must, center, guaranteed); np.add.at(may, center, allowed)
        np.add.at(bmust, center, guaranteed*(labels[neighbor, None] == 2))
        np.add.at(bmay, center, allowed*(labels[neighbor, None] == 2))
        np.add.at(removed, center, raw_present[:, None] & ~guaranteed)
        np.add.at(added, center, ~raw_present[:, None] & allowed)
    dlow, dhigh = must.min(axis=1), may.max(axis=1)
    ml = int(2*((int(dlow.sum())+1)//2)); mh = int(2*(int(dhigh.sum())//2))
    m0 = int(2*raw_present.sum())
    independent_d = int(2*np.sum(raw_present & ~possible.all(axis=(1, 2))))
    independent_a = int(2*np.sum(~raw_present & possible.any(axis=(1, 2))))
    local_d = min(independent_d, int(2*(int(removed.max(axis=1).sum())//2)))
    local_a = min(independent_a, int(2*(int(added.max(axis=1).sum())//2)))
    raw_pairs = pairs[raw_present]
    edges = np.concatenate((raw_pairs, raw_pairs[:, ::-1]))
    binary = np.column_stack((labels == 1, labels == 2)).astype(int)
    fixed = coupled_contrast_bounds(edges, binary, ba, bb)['outer_bounds'] if m0 else [-.25, .25]
    independent = _enlarge(fixed, _tv_bound(m0, independent_d, independent_a))
    local = _enlarge(fixed, _tv_bound(m0, local_d, local_a))
    # A local B-neighbor bound can spend all B errors at each vertex; this relaxes
    # global consistency safely. A errors then use a global top-budget sum.
    wlow = np.maximum(0, bmust-bb).min(axis=1)
    whigh = np.minimum(may, bmay+bb).max(axis=1)
    contacts = _sum_bounds(wlow, whigh, labels == 1, ba)
    adegree = _sum_bounds(dlow, dhigh, labels == 1, ba)
    bdegree = _sum_bounds(dlow, dhigh, labels == 2, bb)
    moments = [-.25, .25]
    if ml > 0:
        al, ah = adegree[0]/mh, min(1., adegree[1]/ml)
        bl, bh = bdegree[0]/mh, min(1., bdegree[1]/ml)
        # Source identities are exclusive: pA+pB<=1 on the same undirected graph.
        if ah+bh <= 1:
            maximum_product = ah*bh
        else:
            feasible_low, feasible_high = max(al, 1-bh), min(ah, 1-bl)
            if feasible_low <= feasible_high:
                x = np.clip(.5, feasible_low, feasible_high)
                maximum_product = float(x*(1-x))
            else:
                maximum_product = .25  # Safe fallback for inconsistent relaxed marginals.
        moments = [max(-.25, contacts[0]/mh-maximum_product-1e-12),
                   min(.25, min(.5, contacts[1]/ml)-al*bl+1e-12)]
    outer = [max(local[0], moments[0]), min(local[1], moments[1])]
    return {'outer_bounds': outer, 'moment_bounds': moments, 'local_edit_bounds': local,
            'independent_edge_bounds': independent, 'fixed_graph_bounds': fixed,
            'directed_edge_count_bounds': [ml, mh], 'contact_count_bounds': contacts,
            'a_degree_bounds': adegree, 'b_degree_bounds': bdegree,
            'edit_caps': {'independent': [independent_d, independent_a], 'local': [local_d, local_a]},
            'may_include_empty_graph': ml == 0, 'candidate_pairs': len(pairs), 'raw_edges': m0}


def exact_position_oracle(positions, raw, radius, budget_a, budget_b):
    """All positions and categorical labels; rational extrema for nonempty graphs."""
    positions = np.asarray(positions); raw = np.asarray(raw); n, choices, _ = positions.shape
    labelings = np.array([x for x in product(range(3), repeat=n)
                         if sum((np.array(x) == 1) != (raw == 1)) <= budget_a
                         and sum((np.array(x) == 2) != (raw == 2)) <= budget_b], dtype=np.int64)
    a, b = (labelings == 1).astype(np.int64), (labelings == 2).astype(np.int64)
    low = high = None; low_world = high_world = None; graphs = {}; empty = 0
    for selected in product(range(choices), repeat=n):
        xy = positions[np.arange(n), selected]
        distance = np.sum((xy[:, None]-xy[None, :])**2, axis=2)
        graph = (distance <= radius*radius) & ~np.eye(n, dtype=bool)
        key = graph.tobytes()
        if key in graphs:
            continue
        graphs[key] = True; m = int(graph.sum())
        if m == 0:
            empty += 1; continue
        degree = graph.sum(axis=1)
        contact = ((a@graph)*b).sum(axis=1)
        numerator = m*contact-(a@degree)*(b@degree)
        for minimize, index in ((True, int(np.argmin(numerator))), (False, int(np.argmax(numerator)))):
            value = Fraction(int(numerator[index]), m*m)
            world = {'positions': list(selected), 'labels': labelings[index].tolist(),
                     'numerator': value.numerator, 'denominator': value.denominator}
            if minimize and (low is None or value < low):low, low_world = value, world
            if not minimize and (high is None or value > high):high, high_world = value, world
    return {'bounds': [float(low), float(high)] if low is not None else None,
            'witnesses': [low_world, high_world], 'unique_graphs': len(graphs),
            'empty_graphs': empty, 'admissible_labelings': len(labelings),
            'position_assignments': choices**n}
