"""Independent complete-world, transcript and minimum-cost replay. No solver imports."""
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def answer(world, action):
    return world['fragile_present' if action['kind'] == 'geometry' else 'labels'][action['index']]


def separates(worlds, actions, selected, error):
    for x, y in combinations(worlds, 2):
        if x['sign'] == y['sign']:
            continue
        exact = any(a['kind'] != 'routine' and answer(x, a) != answer(y, a) for a in (actions[i] for i in selected))
        distance = sum(a['kind'] == 'routine' and answer(x, a) != answer(y, a) for a in (actions[i] for i in selected))
        if not exact and distance <= 2*error:
            return False
    return True


def transcripts(worlds, actions, selected, error):
    seen = {}
    for world in worlds:
        options = [range(3) if actions[i]['kind'] == 'routine' else (answer(world, actions[i]),) for i in selected]
        for transcript in product(*options):
            mistakes = sum(value != answer(world, actions[i]) for i, value in zip(selected, transcript)
                           if actions[i]['kind'] == 'routine')
            if mistakes > error:
                continue
            if transcript in seen and seen[transcript] != world['sign']:
                raise AssertionError('A reported certificate admits contradictory decisions')
            seen[transcript] = world['sign']
    return len(seen)


def main():
    report = json.loads((ROOT / 'results/joint_measurement.json').read_text())
    checked = []
    for row in report['rows']:
        worlds, actions = row['worlds'], row['actions']
        # Rebuild the complete admissible Cartesian product independently.
        raw = row['raw']; base = set(map(tuple, row['edges'])); fragile = list(map(tuple, row['fragile']))
        expected = set()
        for labels in product(range(3), repeat=len(raw)):
            if any(sum(int(x == marker) != int(y == marker) for x, y in zip(labels, raw)) > b
                   for marker, b in zip((1, 2), row['budgets'])):
                continue
            for toggle in (None, *fragile):
                graph = base if toggle is None else base.symmetric_difference((toggle,))
                expected.add((labels, tuple(sorted(graph))))
        actual = {(tuple(w['labels']), tuple(map(tuple, w['edges']))) for w in worlds}
        assert actual == expected and len(actual) == len(worlds)
        for w in worlds:
            graph = set(map(tuple, w['edges'])); labels = w['labels']
            assert w['fragile_present'] == [int(e in graph) for e in fragile]
            d = [sum(i in e for e in graph) for i in range(len(raw))]
            m = 2*len(graph)
            contact = sum((labels[u] == 1 and labels[v] == 2) or (labels[v] == 1 and labels[u] == 2) for u, v in graph)
            numerator = m*contact-sum(d[i] for i, c in enumerate(labels) if c == 1)*sum(d[i] for i, c in enumerate(labels) if c == 2)
            assert numerator == w['numerator'] and m*m == w['denominator']
            assert (numerator > 0)-(numerator < 0) == w['sign']
        for experiment in row['experiments']:
            error = experiment['error_cells']; optimum = experiment['optimal']; lower_count = 0
            for mask in range(1 << len(actions)):
                selected = [i for i in range(len(actions)) if mask & (1 << i)]
                if sum(actions[i]['cost'] for i in selected) < optimum['cost']:
                    assert not separates(worlds, actions, selected, error)
                    lower_count += 1
            transcript_counts = {}
            for key in ('optimal', 'greedy', 'pruned_greedy'):
                plan = experiment[key]
                assert plan['cost'] == sum(actions[i]['cost'] for i in plan['actions'])
                transcript_counts[key] = transcripts(worlds, actions, plan['actions'], error)
            restrictions = {'routine_only': {'routine'}, 'labels_only': {'routine', 'adjudication'},
                            'geometry_only': {'geometry'}}
            for name, kinds in restrictions.items():
                allowed = [i for i, a in enumerate(actions) if a['kind'] in kinds]
                costs = []
                for mask in range(1 << len(allowed)):
                    selected = [i for bit, i in enumerate(allowed) if mask & (1 << bit)]
                    if separates(worlds, actions, selected, error):
                        costs.append(sum(actions[i]['cost'] for i in selected))
                assert (min(costs) if costs else None) == experiment['restricted'][name]['cost']
            checked.append({'case': row['case'], 'error_cells': error, 'worlds': len(worlds),
                            'all_lower_cost_plans_rejected': lower_count, 'transcripts_verified': transcript_counts})
    (ROOT / 'results/joint_measurement_verification.json').write_text(json.dumps({'rows': checked}, indent=2) + '\n')
    print('Verified complete populations, exact targets, all selected-plan transcripts and all lower-cost alternatives.')


if __name__ == '__main__':
    main()
