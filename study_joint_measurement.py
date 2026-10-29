"""Frozen three-case, two-error-budget exhaustive measurement study."""
import json
from pathlib import Path
from joint_measurement import enumerate_worlds, measurement_actions, optimal_fixed_plan, greedy_fixed_plan, prune_plan

ROOT = Path(__file__).resolve().parent
CASES = (
    ((1, 2, 0, 0), ((0, 1), (1, 2), (2, 3)), ((0, 2), (1, 2))),
    ((1, 0, 2, 0), ((0, 1), (1, 2), (2, 3), (0, 3)), ((0, 2), (0, 3))),
    ((1, 2, 1, 2), ((0, 1), (0, 2), (0, 3)), ((1, 2), (0, 3))),
)


def main():
    rows = []
    for case, (raw, edges, fragile) in enumerate(CASES, 1):
        worlds = enumerate_worlds(raw, edges, fragile)
        actions = measurement_actions(len(raw), fragile)
        experiments = []
        for error in (0, 1):
            restrictions = {'routine_only': {'routine'}, 'labels_only': {'routine', 'adjudication'},
                            'geometry_only': {'geometry'}}
            result = {'error_cells': error, 'optimal': optimal_fixed_plan(worlds, actions, error),
                      'greedy': greedy_fixed_plan(worlds, actions, error),
                      'full_exact_cost': sum(a['cost'] for a in actions if a['kind'] != 'routine'),
                      'restricted': {name: optimal_fixed_plan(worlds, actions, error, kinds)
                                     for name, kinds in restrictions.items()}}
            result['pruned_greedy'] = prune_plan(worlds, actions, error, result['greedy'])
            experiments.append(result)
            print('case', case, 'E', error, 'worlds', len(worlds), 'optimal', result['optimal'],
                  'greedy', result['greedy'], flush=True)
        rows.append({'case': case, 'raw': raw, 'edges': edges, 'fragile': fragile, 'budgets': [1, 1],
                     'actions': actions, 'worlds': worlds, 'experiments': experiments})
    (ROOT / 'results/joint_measurement.json').write_text(json.dumps({
        'protocol_commit': '521524a', 'control_amendment_commit': '6d76bb3',
        'scope': 'Complete finite synthetic problems, hypothetical exact references and costs; fixed plans only.',
        'rows': rows}, indent=2) + '\n')


if __name__ == '__main__':
    main()
