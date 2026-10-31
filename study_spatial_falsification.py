"""Frozen cross-channel coupling screen against explicit anatomy controls."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm

SIDE = 100.
RADIUS = 3.
SCENARIOS = ('uniform', 'shared_anatomy', 'within_cell', 'separated', 'coupled')


def sample_b(rng, scenario, count, n=128, centers=None):
    shape = (count, n, 2)
    if scenario == 'uniform':
        return rng.uniform(0, SIDE, shape)
    if scenario == 'separated':
        result = rng.uniform(0, SIDE, shape)
        result[:, :, 0] = rng.uniform(60, 90, (count, n))
        return result
    if scenario == 'within_cell':
        return (np.repeat(centers, n//len(centers), axis=0)[None, :, :]+rng.normal(0, 1.5, shape)) % SIDE
    centers = np.array([[25., 25.], [75., 70.], [35., 80.]])
    ids = rng.choice(3, size=(count, n), p=[.5, .3, .2])
    return (centers[ids]+rng.normal(0, 8., shape)) % SIDE


def generate(rng, scenario, n=128):
    centers = rng.uniform(5, 95, (16, 2)) if scenario == 'within_cell' else None
    a = sample_b(rng, scenario, 1, n, centers)[0]
    b = sample_b(rng, scenario, 1, n, centers)[0]
    if scenario == 'separated':
        a[:, 0] = rng.uniform(10, 40, n)
    if scenario == 'coupled':
        b[:96] = (a[:96]+rng.normal(0, .5, (96, 2))) % SIDE
    return a, b, centers


def nearest_distances(a, b):
    return cKDTree(a, boxsize=SIDE).query(np.asarray(b).reshape(-1, 2))[0].reshape(np.asarray(b).shape[:-1])


def score(a, b, radius=RADIUS):
    return (nearest_distances(a, b) <= radius).mean(axis=-1)


def mc_p(observed, null):
    return float((1+np.count_nonzero(np.asarray(null) >= observed-1e-12))/(len(null)+1))


def label_null(a, b, count, rng):
    points = np.concatenate([a, b])
    d = np.abs(points[:, None, :]-points[None, :, :])
    near = (np.sum(np.minimum(d, SIDE-d)**2, axis=-1) <= RADIUS**2).astype(float)
    assignments = np.stack([rng.permutation(np.arange(len(points)) < len(a)) for _ in range(count)])
    return (((assignments.astype(float) @ near) > 0) & ~assignments).sum(axis=1)/len(b)


def perturb(rng, positions, radius):
    angle = rng.uniform(0, 2*np.pi, len(positions))
    length = radius*np.sqrt(rng.uniform(size=len(positions)))
    return (positions+length[:, None]*np.column_stack([np.cos(angle), np.sin(angle)])) % SIDE


def robust_oracle(a_observed, b_observed, null_b, error_radius):
    # Round the two threshold comparisons outward: exact-boundary floating-point
    # distances must not inflate a lower bound or reduce an upper bound.
    lower = float((nearest_distances(a_observed, b_observed) <= RADIUS-2*error_radius-1e-10).mean())
    upper = (nearest_distances(a_observed, null_b) <= RADIUS+error_radius+1e-10).mean(axis=-1)
    return mc_p(lower, upper), lower


def wilson(successes, count):
    z = norm.ppf(.975); p = successes/count
    center = (p+z*z/(2*count))/(1+z*z/count)
    half = z*np.sqrt(p*(1-p)/count+z*z/(4*count*count))/(1+z*z/count)
    return [float(center-half), float(center+half)]


def run(output, repeats=200, permutations=199):
    rng = np.random.default_rng(2026102907)
    records = []; violations = 0
    for scenario in SCENARIOS:
        for replicate in range(repeats):
            a, b, centers = generate(rng, scenario)
            observed = float(score(a, b))
            null_b = sample_b(rng, scenario, permutations, centers=centers)
            oracle_p = mc_p(observed, score(a, null_b))
            tests = dict(oracle=oracle_p,
                         global_labels=mc_p(observed, label_null(a, b, permutations, rng)),
                         uniform_b=mc_p(observed, score(a, rng.uniform(0, SIDE, (permutations, len(b), 2)))))
            tile_b = np.floor(b/25)[None, :, :]*25+rng.uniform(0, 25, (permutations, len(b), 2))
            tests['conditional_tiles'] = mc_p(observed, score(a, tile_b))
            for rho in [.25, .75, 1.5]:
                ao, bo = perturb(rng, a, rho), perturb(rng, b, rho)
                robust_p, _ = robust_oracle(ao, bo, null_b, rho)
                violations += int(robust_p < oracle_p-1e-12)
                tests[f'robust_oracle_rho_{rho}'] = robust_p
            for method, p in tests.items():
                records.append(dict(scenario=scenario, replicate=replicate, method=method,
                                    p_value=p, rejected=p <= .05, latent_proximity=observed))
        print('Completed', scenario, flush=True)
    table = pd.DataFrame(records)
    summary = []
    for (scenario, method), group in table.groupby(['scenario', 'method'], sort=True):
        rejected = int(group.rejected.sum())
        summary.append(dict(scenario=scenario, method=method, repetitions=len(group), rejections=rejected,
                            rejection_rate=rejected/len(group), wilson_95=wilson(rejected, len(group))))
    gates = {}
    for method in ['oracle', 'robust_oracle_rho_0.25']:
        selected = [r for r in summary if r['method'] == method]
        null_pass = all(r['wilson_95'][1] <= .1 for r in selected if r['scenario'] != 'coupled')
        power_pass = all(r['wilson_95'][0] >= .7 for r in selected if r['scenario'] == 'coupled')
        gates[method] = dict(null_gate=bool(null_pass), power_gate=bool(power_pass), passed=bool(null_pass and power_pass))
    result = dict(seed=2026102907, repeats=repeats, monte_carlo_draws=permutations, n_per_channel=128,
                  radius=RADIUS, periodic_window_side=SIDE, summary=summary, gates=gates,
                  pointwise_domination_violations=violations,
                  measured_data=False, oracle_intensity_known=True,
                  novelty='Cross-nearest-neighbour statistic and conditional resampling are established controls; the bound uses conventional distance envelopes.')
    if violations:
        raise AssertionError('Robust p-value failed latent oracle domination')
    output.mkdir(exist_ok=True)
    table.to_csv(output/'spatial_falsification_trials.csv', index=False)
    (output/'spatial_falsification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(gates, indent=2))


if __name__ == '__main__':
    run(Path(__file__).resolve().parent/'results')
