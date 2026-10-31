"""Replay saved measurements without importing the study estimators or LP solver."""
import json
import gzip
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import betainc, gammainc, gammaincc
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent


def interval(k, n):
    z = 1.959963984540054
    p = k/n; den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [center-radius, center+radius]


def exact_bounds(count, exposure):
    low = np.zeros_like(count, dtype=float)
    positive = count > 0
    low[positive] = chi2.ppf(.0125, 2*count[positive])/(2*exposure)
    return low, chi2.ppf(.9875, 2*(count+1))/(2*exposure)


def replay_probe(study, counts_path):
    saved = np.load(counts_path)
    checked = 0
    for row in study['settings']:
        b, delta, eta, dose, truth = [row[k] for k in ['background', 'delta', 'efficiency_min', 'dose', 'truth']]
        x, y = saved[f'b{b}_delta{delta}_eta{eta}_E{dose}_z{truth}'].T
        method = row['method']
        if method == 'joint_rectangle':
            lx, ux = exact_bounds(x, dose); lb, ub = exact_bounds(y, 16*dose)
            present = lx > (1+delta)*ub
            absent = ux-(1-delta)*lb < eta
        elif method in ['paired_reference', 'conditional_reference']:
            n = x+y; p0 = (1+delta)/(17+delta)
            p1 = ((1-delta)*10+eta)/((17-delta)*10+eta)
            tail0 = np.ones(len(x)); take = x > 0
            tail0[take] = betainc(x[take], n[take]-x[take]+1, p0)
            tail1 = np.ones(len(x)); take = x < n
            tail1[take] = betainc(n[take]-x[take], x[take]+1, 1-p1)
            present, absent = tail0 <= .025, tail1 <= .025
        else:
            tail0 = np.ones(len(x)); take = x > 0
            tail0[take] = gammainc(x[take], dose*(1+delta)*10)
            tail1 = gammaincc(x+1, dose*((1-delta)*2+eta))
            present, absent = tail0 <= .025, tail1 <= .025
        decision = np.full(len(x), -1)
        decision[present&~absent] = 1; decision[absent&~present] = 0
        correct = decision == truth; wrong = (decision >= 0)&(decision != truth)
        np.testing.assert_allclose([correct.mean(), wrong.mean(), (decision < 0).mean()],
                                   [row['correct'], row['incorrect'], row['abstain']], atol=1e-14)
        np.testing.assert_allclose(interval(correct.sum(), len(x)), row['correct_wilson'], atol=1e-8)
        np.testing.assert_allclose(interval(wrong.sum(), len(x)), row['incorrect_wilson'], atol=1e-8)
        checked += len(x)
    return checked


def run():
    source = json.loads((ROOT/'results/fiducial_measurement.json').read_text())
    detail_bytes = (ROOT/'results'/source['geometry']['details_file']).read_bytes()
    assert hashlib.sha256(detail_bytes).hexdigest() == source['geometry']['details_sha256']
    details = json.loads(gzip.decompress(detail_bytes))
    codebook = pd.read_csv(ROOT/'datasets/spatial_v01/codebook.csv').set_index('name').iloc[:, 1:]
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    geometries = 0
    for w in details['witnesses']:
        bits = w['active_bits']
        first = w['first'].split(' + '); second = w['second'].split(' + ')
        np.testing.assert_array_equal(codebook.loc[first].sum(), codebook.loc[second].sum())
        h = ((codebook.loc[first[1]]-codebook.loc[first[0]]
              -codebook.loc[second[1]]+codebook.loc[second[0]])/2).iloc[bits].to_numpy()
        np.testing.assert_array_equal(h, w['required_center_difference'])
        selected = layout.iloc[bits]
        ar = np.eye(8)[selected.imagingRound.to_numpy()-1]
        ac = (selected.color.to_numpy()[:, None] == np.array([561, 650])).astype(float)
        matrices = dict(independent_bits=np.eye(len(bits)), shared_rounds=ar,
                        rounds_and_colors=np.column_stack([ar, ac]), known_rounds_colors_only=ac)
        for mode, matrix in matrices.items():
            rank_same = np.linalg.matrix_rank(matrix) == np.linalg.matrix_rank(np.column_stack([matrix, h]))
            assert rank_same == (w[mode+'_residual_linf'] < 1e-10)
        if w['lp_coefficients'] is not None:
            np.testing.assert_allclose(matrices['rounds_and_colors'] @ w['lp_coefficients'], h, atol=1e-8)
            assert max(abs(np.array(w['lp_coefficients'][:8]))) <= 2*w['minimum_per_world_round_radius_at_d1']+1e-8
        geometries += 1
    detections = pd.read_csv(ROOT/'results/fiducial_detections.csv').set_index(['round', 'detection_id'])
    matches = pd.read_csv(ROOT/'results/fiducial_matches.csv')
    for r in source['fiducials']['rounds']:
        rows = matches[matches['round'] == r['round']]
        shifts = []
        for _, match in rows.iterrows():
            reference = detections.loc[(1, match.reference_id)][['y', 'x']].to_numpy()
            moving = detections.loc[(r['round'], match.moving_id)][['y', 'x']].to_numpy()
            shifts.append(reference-moving)
            expected_train = int(np.floor(reference/128).sum())%2 == 0
            assert expected_train == match.training
        shifts = np.array(shifts)
        np.testing.assert_allclose(shifts, rows[['shift_y', 'shift_x']], atol=1e-10)
        train = rows.training.to_numpy(dtype=bool)
        fitted = np.median(shifts[train], axis=0)
        radii = np.sqrt(np.sum((shifts[~train]-fitted)**2, axis=1))
        np.testing.assert_allclose([np.median(radii), np.quantile(radii, .95), max(radii)],
                                   [r['heldout_median'], r['heldout_q95'], r['heldout_max']], atol=1e-10)
    sets = json.loads((ROOT/'results/fiducial_sets.json').read_text())
    certificates = 0
    for setting in sets['results']:
        ell = setting['allowance_pixels']; compatible = True
        for entry in setting['round_intervals']:
            rows = matches[matches['round'] == entry['round']][['shift_y', 'shift_x']].to_numpy()
            low = rows.max(axis=0)-ell; high = rows.min(axis=0)+ell
            np.testing.assert_allclose(low, entry['lower'], atol=1e-10)
            np.testing.assert_allclose(high, entry['upper'], atol=1e-10)
            compatible &= bool(np.all(low <= high))
        assert compatible == setting['contract_compatible']
        for row in setting['comparisons']:
            selected = layout.iloc[row['active_bits']]
            ar = np.eye(8)[selected.imagingRound.to_numpy()-1]
            ac = (selected.color.to_numpy()[:, None] == np.array([561, 650])).astype(float)
            a = np.column_stack([ar, ac, -np.array(row['difference'])])
            x = np.array(row['primal']); lam = np.array(row['equality_dual'])
            lower = np.array(row['lower_dual']); upper = np.array(row['upper_dual'])
            c = np.zeros(11); c[-1] = -1
            widths = np.array(row['round_difference_widths'])
            axis = ['y', 'x'].index(row['axis'])
            expected_widths = [0.]+[entry['upper'][axis]-entry['lower'][axis] for entry in setting['round_intervals']]
            np.testing.assert_allclose(widths, expected_widths, atol=1e-10)
            np.testing.assert_allclose(a @ x, 0, atol=1e-8)
            np.testing.assert_allclose(x[-1], row['maximum_alias_separation'], atol=1e-8)
            assert np.all(abs(x[:8]) <= widths+1e-8) and x[-1] >= -1e-8
            assert np.all(lower >= -1e-8) and np.all(upper <= 1e-8)
            np.testing.assert_allclose(a.T @ lam+lower+upper, c, atol=1e-8)
            np.testing.assert_allclose(lower[8:10], 0, atol=1e-8)
            np.testing.assert_allclose(upper[8:], 0, atol=1e-8)
            # Dual objective proves maximality without running an optimizer.
            dual_value = np.dot(-widths, lower[:8])+np.dot(widths, upper[:8])
            np.testing.assert_allclose(dual_value, -row['maximum_alias_separation'], atol=1e-8)
            certificates += 1
    decisions = replay_probe(source['probe'], ROOT/'results/paired_reference_counts.npz')
    joint = json.loads((ROOT/'results/paired_reference_joint.json').read_text())
    decisions += replay_probe(joint, ROOT/'results/paired_reference_joint_counts.npz')
    for method, gates in joint['gates'].items():
        rows = [r for r in joint['settings'] if r['method'] == method]
        strong = [r for r in rows if r['efficiency_min'] == 8 and r['dose'] == 64]
        assert gates['incorrect_decision_control'] == all(r['incorrect_wilson'][1] <= .075 for r in rows)
        assert gates['strong_signal_decisiveness'] == all(r['correct_wilson'][0] >= .8 for r in strong)
    result = dict(verified=True, oriented_geometry_cases=geometries, fiducial_matches=len(matches),
                  optimal_boundary_certificates=certificates, experiment_method_decisions=decisions,
                  biological_truth_verified=False, raw_peak_detection_replayed=False)
    (ROOT/'results/fiducial_measurement_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    run()
