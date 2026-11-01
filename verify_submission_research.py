"""Independent source-coordinate and full Gaussian-experiment replays."""
import gzip
import hashlib
import itertools
import json
from collections import defaultdict, deque
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import null_space, qr, solve_triangular
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent


def verify_reference():
    source = ROOT/'datasets/reference_calibration_v01/cross_channel_matches.csv'
    frame = pd.read_csv(source)
    study = json.loads((ROOT/'results/cross_channel_calibration.json').read_text())
    payload = (ROOT/'results'/study['details_file']).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == study['details_sha256']
    study.update(json.loads(gzip.decompress(payload)))
    assert hashlib.sha256(source.read_bytes()).hexdigest() == study['source_sha256']
    graph = defaultdict(set); edges = []
    frozen = json.loads((ROOT/'datasets/reference_calibration_v01/frozen_feature_groups.json').read_text())
    anchors = {tuple(sorted((int(c), round(x, 8), round(y, 8)) for c, x, y in r['nodes'])):r for r in frozen['groups']}
    for row in frame.to_numpy():
        first = (int(row[16]), round(float(row[0]), 8), round(float(row[1]), 8))
        second = (int(row[17]), round(float(row[8]), 8), round(float(row[9]), 8))
        graph[first].add(second); graph[second].add(first); edges.append((first, second))
    assigned, nodes = {}, {}
    for first in graph:
        if first in assigned:
            continue
        seen, queue = {first}, deque([first])
        while queue:
            for other in graph[queue.popleft()]:
                if other not in seen:
                    seen.add(other); queue.append(other)
        anchor = anchors[tuple(sorted(seen))]
        encoded = ';'.join(f'{c}:{float(x).hex()}:{float(y).hex()}' for c, x, y in sorted(anchor['nodes']))
        key = hashlib.sha256(encoded.encode()).hexdigest()
        assert key == anchor['group_id']
        nodes[key] = anchor['nodes']
        for node in seen:
            assigned[node] = key
    assert len(nodes) == study['physical_feature_components']
    predicted = {(r['mode'], r['method'], r['source_row']):r for r in study['predictions']}
    checked = 0
    for mode in ['hashed', 'spatial_extrapolation']:
        ordered = sorted(nodes) if mode == 'hashed' else sorted(nodes, key=lambda key:(np.mean([n[1] for n in nodes[key]]), key))
        split = {key:('train' if i < len(nodes)//2 else 'calibration' if i < len(nodes)//2+len(nodes)//4 else 'test') for i, key in enumerate(ordered)}
        for saved in [s for s in study['splits'] if s['mode'] == mode]:
            assert split[saved['group_id']] == saved['split']
        for method in ['translation', 'affine']:
            score = defaultdict(list)
            for (ref, comp), subset in frame.groupby(['ref_ch', 'comp_ch']):
                indices = subset.index.to_numpy()
                groups = [assigned[edges[i][0]] for i in indices]
                train = np.array([split[g] == 'train' for g in groups])
                x = subset[['comp_x', 'comp_y']].to_numpy(); y = subset[['ref_x', 'ref_y']].to_numpy()
                if method == 'translation':
                    estimate = x+np.median(y[train]-x[train], axis=0)
                else:
                    matrix = np.column_stack([x[train], np.ones(train.sum())])
                    q, r = qr(matrix, mode='economic')
                    coefficients = solve_triangular(r, q.T @ y[train])
                    estimate = np.column_stack([x, np.ones(len(x))]) @ coefficients
                for k, index in enumerate(indices):
                    saved = predicted[mode, method, int(index)]
                    np.testing.assert_allclose(estimate[k], [saved['predicted_x'], saved['predicted_y']], atol=1e-8, rtol=0)
                    residual = float(np.linalg.norm(estimate[k]-y[k]))
                    assert abs(residual-saved['residual']) < 1e-8
                    assert saved['group_id'] == groups[k] and saved['split'] == split[groups[k]]
                    score[groups[k]].append(residual); checked += 1
            maximum = {key:max(value) for key, value in score.items()}
            calibration = sorted(maximum[key] for key in maximum if split[key] == 'calibration')
            radius = calibration[int(np.ceil(.95*(len(calibration)+1)))-1]
            test = [maximum[key] for key in maximum if split[key] == 'test']
            summary = next(s for s in study['settings'] if s['mode'] == mode and s['method'] == method)
            np.testing.assert_allclose([radius, np.median(test), max(test)],
                                       [summary['radius'], summary['median_test_group_residual'], summary['maximum_test_group_residual']], atol=1e-8)
            assert sum(value <= radius for value in test) == summary['test_covered']
    return dict(physical_feature_components=len(nodes), source_predictions=checked, independent_identity_truth=False)


def verify_gaussian():
    study = json.loads((ROOT/'results/photon_ceiling.json').read_text())
    base = json.loads((ROOT/'results/fiducial_measurement.json').read_text())
    details = json.loads(gzip.decompress((ROOT/'results'/base['geometry']['details_file']).read_bytes()))
    cases = [w for w in details['witnesses'] if w['first'] == 'Adamdec1 + Adgrg6']
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    saved = np.load(ROOT/'results/photon_ceiling_statistics.npz')
    rng = np.random.default_rng(2026110102)
    scored = 0
    for row in study['cases']:
        case = cases[row['case']]
        selected = layout.iloc[case['active_bits']]
        ar = np.eye(8)[selected.imagingRound.to_numpy()-1]
        ac = (selected.color.to_numpy()[:, None] == np.array([561, 650])).astype(float)
        h = row['separation']*np.array(case['required_center_difference'])
        vg = .09/row['exposure']; vf = row['fiducial_sd']**2/row['fiducial_references']
        std = np.r_[np.full(len(h), np.sqrt(vg)), np.full(8, np.sqrt(vf))]
        matrix = np.block([[ar, ac], [np.eye(8), np.zeros((8, 2))]])/std[:, None]
        signal = np.r_[h, np.zeros(8)]/std
        orthogonal = null_space(matrix.T)
        projected = orthogonal @ (orthogonal.T @ signal)
        information = float(projected @ projected)
        np.testing.assert_allclose(information, row['squared_information'], atol=1e-8)
        color_null = null_space(ac.T)
        drift = np.linalg.pinv(color_null.T @ ar) @ (color_null.T @ h)
        ceiling = float(drift @ drift/vf)
        np.testing.assert_allclose(ceiling, row['ceiling_information'], atol=1e-8)
        np.testing.assert_allclose(norm.cdf(-np.sqrt(information)/2), row['optimal_error'], atol=1e-12)
        if row['case'] != 0:
            continue
        errors, extra_errors = [], []
        for truth in [0, 1]:
            noise = rng.normal(size=(2000, len(h)))*np.sqrt(vg)
            fiducials = rng.normal(size=(2000, 8))*np.sqrt(vf)
            joint_observation = np.column_stack([truth*h+noise, fiducials])
            statistic = joint_observation @ (projected/std)-information/2
            probe = rng.normal(4*truth, 1., size=2000)
            key = f"d{row['separation']}_tau{row['fiducial_sd']}_n{row['fiducial_references']}_E{row['exposure']}_truth{truth}"
            np.testing.assert_allclose(np.column_stack([statistic, probe]), saved[key], atol=1e-8, rtol=0)
            errors.append(float(((statistic > 0) != truth).mean()))
            extra_errors.append(float(((statistic+4*probe-8 > 0) != truth).mean()))
            scored += len(statistic)
        np.testing.assert_allclose(errors, row['empirical_errors'])
        np.testing.assert_allclose(extra_errors, row['empirical_probe_errors'])
    return dict(theory_settings=len(study['cases']), saved_full_experiment_statistics=scored, real_covariance_validated=False)


if __name__ == '__main__':
    result = dict(verified=True, reference=verify_reference(), gaussian=verify_gaussian())
    (ROOT/'results/submission_research_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
