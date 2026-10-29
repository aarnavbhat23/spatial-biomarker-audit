"""Complete small-world fixed measurement plans with persistent review errors."""
from itertools import combinations, product
import numpy as np


def enumerate_worlds(raw, edges, fragile, budgets=(1, 1)):
    """All categorical labelings and at most one fragile undirected edge toggle."""
    n = len(raw)
    if not n or any(x not in (0, 1, 2) for x in raw):
        raise ValueError('Need categorical raw calls')
    if len(budgets) != 2 or any(int(b) != b or b < 0 for b in budgets):
        raise ValueError('Need two nonnegative integer error budgets')
    for collection in (edges, fragile):
        if len(set(collection)) != len(collection) or any(not (0 <= u < v < n) for u, v in collection):
            raise ValueError('Need unique canonical undirected edges')
    graphs = [set(edges)] + [set(edges) ^ {edge} for edge in fragile]
    if any(not g for g in graphs):
        raise ValueError('Empty graph makes the covariance target undefined')
    worlds = []
    for labels in product(range(3), repeat=n):
        if any(sum((x == marker) != (y == marker) for x, y in zip(labels, raw)) > budget
               for marker, budget in zip((1, 2), budgets)):
            continue
        for graph in graphs:
            directed = [(u, v) for a, b in sorted(graph) for u, v in ((a, b), (b, a))]
            m = len(directed)
            contact = sum(labels[u] == 1 and labels[v] == 2 for u, v in directed)
            left = sum(labels[u] == 1 for u, _ in directed)
            right = sum(labels[v] == 2 for _, v in directed)
            numerator = m * contact - left * right
            worlds.append({'labels': list(labels), 'edges': [list(e) for e in sorted(graph)],
                           'fragile_present': [int(e in graph) for e in fragile],
                           'numerator': numerator, 'denominator': m*m,
                           'sign': (numerator > 0) - (numerator < 0)})
    return worlds


def measurement_actions(n, fragile):
    return ([{'kind': 'routine', 'index': i, 'cost': 1} for i in range(n)] +
            [{'kind': 'adjudication', 'index': i, 'cost': 4} for i in range(n)] +
            [{'kind': 'geometry', 'index': i, 'cost': 2} for i in range(len(fragile))])


def pair_requirements(worlds, actions, error_cells):
    """Quota 2E+1; one differing exact action supplies the whole quota."""
    if int(error_cells) != error_cells or error_cells < 0:
        raise ValueError('Need a nonnegative integer routine error budget')
    pairs = [(i, j) for i, j in combinations(range(len(worlds)), 2) if worlds[i]['sign'] != worlds[j]['sign']]
    quota = 2*int(error_cells)+1
    weights = np.zeros((len(pairs), len(actions)), dtype=np.int64)
    for row, (i, j) in enumerate(pairs):
        for col, action in enumerate(actions):
            key = 'fragile_present' if action['kind'] == 'geometry' else 'labels'
            different = worlds[i][key][action['index']] != worlds[j][key][action['index']]
            weights[row, col] = int(different) * (1 if action['kind'] == 'routine' else quota)
    return pairs, weights, quota


def optimal_fixed_plan(worlds, actions, error_cells, allowed_kinds=None):
    """Exhaust all subsets. No greedy optimality or population-scale guarantee."""
    pairs, weights, quota = pair_requirements(worlds, actions, error_cells)
    allowed = [i for i, a in enumerate(actions) if allowed_kinds is None or a['kind'] in allowed_kinds]
    plans = []
    for mask in range(1 << len(allowed)):
        selected = [a for bit, a in enumerate(allowed) if mask & (1 << bit)]
        plans.append((sum(actions[i]['cost'] for i in selected), len(selected), selected))
    plans.sort()
    checked = 0
    best = None
    for cost, _, selected in plans:
        checked += 1
        if np.all(weights[:, selected].sum(axis=1) >= quota):
            if best is None:
                best = {'cost': cost, 'actions': selected}
    if best is not None:
        return {**best, 'plans_checked': checked, 'opposite_decision_pairs': len(pairs)}
    return {'cost': None, 'actions': [], 'plans_checked': checked, 'opposite_decision_pairs': len(pairs)}


def greedy_fixed_plan(worlds, actions, error_cells):
    _, weights, quota = pair_requirements(worlds, actions, error_cells)
    remaining = np.full(len(weights), quota, dtype=np.int64)
    selected = []
    while np.any(remaining > 0):
        choices = []
        for i, action in enumerate(actions):
            if i not in selected:
                gain = int(np.minimum(remaining, weights[:, i]).sum())
                choices.append((gain/action['cost'], -action['cost'], -i, i))
        if not choices or max(choices)[0] == 0:
            return {'cost': None, 'actions': selected}
        index = max(choices)[-1]
        selected.append(index)
        remaining = np.maximum(0, remaining - weights[:, index])
    return {'cost': sum(actions[i]['cost'] for i in selected), 'actions': selected}


def prune_plan(worlds, actions, error_cells, plan):
    """Post-screen control: remove redundant actions in the frozen cost order."""
    _, weights, quota = pair_requirements(worlds, actions, error_cells)
    selected = list(plan['actions'])
    if plan['cost'] is None:
        return dict(plan)
    for index in sorted(selected, key=lambda i: (-actions[i]['cost'], i)):
        remaining = [i for i in selected if i != index]
        if np.all(weights[:, remaining].sum(axis=1) >= quota):
            selected = remaining
    return {'cost': sum(actions[i]['cost'] for i in selected), 'actions': selected}
