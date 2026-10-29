"""Frozen small-instance exact checks and declared scaling cases."""
import json
from pathlib import Path
import time
import numpy as np
from position_uncertainty import position_bounds, exact_position_oracle

ROOT = Path(__file__).resolve().parent


def specimen(n, seed):
    rng = np.random.default_rng(seed); side = int(np.ceil(np.sqrt(n)))
    xy = np.array([(3*(i % side), 3*(i//side)) for i in range(n)])
    xy += rng.integers(-1, 2, size=(n, 2))
    labels = rng.choice(3, size=n, p=(.2, .4, .4))
    return xy[:, None]+np.array([[0, 0], [-1, 0], [1, 0]])[None], labels


def main():
    small, scale = [], []
    for n in (5, 6):
        for seed in (19, 43, 71):
            positions, labels = specimen(n, seed)
            for budget in (0, 1):
                bound = position_bounds(positions, labels, 5, budget, budget)
                exact = exact_position_oracle(positions, labels, 5, budget, budget)
                if exact['bounds'] is not None:
                    for key in ('outer_bounds', 'moment_bounds', 'local_edit_bounds', 'independent_edge_bounds'):
                        assert bound[key][0] <= exact['bounds'][0]+1e-10 and bound[key][1] >= exact['bounds'][1]-1e-10, (n, seed, budget, key)
                small.append({'n': n, 'seed': seed, 'budget': budget, 'positions': positions.tolist(),
                              'labels': labels.tolist(), 'bound': bound, 'exact': exact})
    for n in (128, 512, 2048):
        for seed in (101, 202):
            positions, labels = specimen(n, seed)
            for budget in (0, int(np.ceil(.005*n))):
                started = time.perf_counter(); bound = position_bounds(positions, labels, 5, budget, budget)
                elapsed = time.perf_counter()-started
                local_width = np.diff(bound['local_edit_bounds'])[0]; width = np.diff(bound['outer_bounds'])[0]
                scale.append({'n': n, 'seed': seed, 'budget': budget, 'bound': bound,
                              'seconds': elapsed, 'relative_width_reduction': float(1-width/local_width) if local_width else 0.})
                print(n, seed, budget, 'width', round(float(width), 6), 'improvement', scale[-1]['relative_width_reduction'], flush=True)
    median = float(np.median([r['relative_width_reduction'] for r in scale]))
    report = {'protocol_commit': 'b96dac4', 'scope': 'Complete discrete-position sensitivity model; no calibrated biological error or continuous geometry guarantee.',
              'small': small, 'scale': scale, 'median_relative_width_reduction': median, 'continuation_gate_passed': median >= .1}
    (ROOT/'results/position_uncertainty.json').write_text(json.dumps(report, indent=2)+'\n')
    print('All exhaustive cases covered; median improvement', median, 'gate', median >= .1)


if __name__ == '__main__':main()
