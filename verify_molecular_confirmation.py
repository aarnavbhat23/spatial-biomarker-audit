"""Source-to-result replay using separate distance/statistic calculations."""
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parent


def cross_means(kernel, labels):
    a = np.atleast_2d(labels).astype(float)
    b = 1-a
    return np.einsum('ij,ij->i', a @ kernel, b)/(a.sum(axis=1)*b.sum(axis=1))


def main():
    report = json.loads((ROOT/'results/molecular_confirmation.json').read_text())
    path = ROOT/'data/visp/s3_spot_table.csv'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == report['source_sha256']
    df = pd.read_csv(path)
    xy = df[['rotated_x', 'rotated_y']].to_numpy()
    low, high = xy.min(axis=0), xy.max(axis=0)
    middle = (low+high)/2
    np.testing.assert_allclose(report['coordinate_bounds_um'], [low, high], rtol=0, atol=0)
    genes = ('Gad2', 'Slc17a7', 'Pvalb', 'Sst', 'Vip', 'Cux2')
    rng = np.random.default_rng(report['seed'])
    checked = 0
    max_discrepancy = 0.
    summaries = {}
    for half in ('left_development', 'right_confirmation'):
        region = xy[:, 0] < middle[0] if half.startswith('left') else xy[:, 0] >= middle[0]
        chosen = {}
        for gene in genes:
            ids = np.flatnonzero(region & (df.gene.to_numpy() == gene))
            chosen[gene] = np.sort(rng.choice(ids, min(200, len(ids)), replace=False))
            saved = next(x for x in report['selection'] if (x['region'], x['gene']) == (half, gene))
            assert saved['available'] == len(ids) and saved['selected'] == len(chosen[gene])
            np.testing.assert_array_equal(chosen[gene], saved['source_row_indices'])
        summary = dict(evaluated_tests=0, missing_or_sparse_pairs=0, global_nominal_rejections=0,
                       stratified_nominal_rejections=0, robust_rejections=[0]*4)
        xmid = (low[0]+middle[0])/2 if half.startswith('left') else (middle[0]+high[0])/2
        for a, b in itertools.combinations(genes, 2):
            matching = [x for x in report['rows'] if x['region'] == half and x['genes'] == [a, b]]
            ia, ib = chosen[a], chosen[b]
            if min(len(ia), len(ib)) < 20:
                assert len(matching) == 1 and matching[0]['status'] == 'unevaluable'
                summary['missing_or_sparse_pairs'] += 1
                continue
            positions = xy[np.r_[ia, ib]]
            labels = np.r_[np.ones(len(ia), int), np.zeros(len(ib), int)]
            blocks = 2*(positions[:, 0] >= xmid) + (positions[:, 1] >= middle[1])
            shuffles = []
            for grouping in (blocks, np.zeros(len(labels), int)):
                shuffled = np.tile(labels, (1999, 1))
                for block in np.unique(grouping):
                    ids = np.flatnonzero(grouping == block)
                    for row in shuffled:
                        row[ids] = rng.permutation(labels[ids])
                shuffles.append(shuffled)
            distance = squareform(pdist(positions))
            for h in (5., 10., 20.):
                saved = next(x for x in matching if x['bandwidth_um'] == h)
                kernel = np.exp(-distance**2/(2*h*h))
                observed = float(cross_means(kernel, labels)[0])
                global_p = (1 + np.sum(cross_means(kernel, shuffles[1]) >= observed-1e-14))/2000
                assert global_p == saved['global_nominal_p']
                summary['global_nominal_rejections'] += int(45*global_p <= .05)
                summary['evaluated_tests'] += 1
                for j, (rho, k) in enumerate(((0., 0), (.5, 0), (.5, 1), (1., 1))):
                    sensitivity = saved['sensitivity'][j]
                    error = k/len(ia)+k/len(ib)
                    lower = max(0., cross_means(np.exp(-(distance+2*rho)**2/(2*h*h)), labels)[0]-error-1e-12)
                    upper = np.minimum(1., cross_means(np.exp(-np.maximum(0., distance-2*rho)**2/(2*h*h)), shuffles[0])+error+1e-12)
                    nominal_count = int(np.sum(cross_means(kernel, shuffles[0]) >= observed-1e-14))
                    upper_count = int(np.sum(upper >= lower))
                    assert nominal_count == sensitivity['nominal_exceedances']
                    assert upper_count == sensitivity['robust_upper_exceedances']
                    assert (1+nominal_count)/2000 == sensitivity['nominal_p']
                    assert (1+upper_count)/2000 == sensitivity['robust_upper_p']
                    assert min(1., 45*sensitivity['nominal_p']) == sensitivity['bonferroni_nominal_p']
                    assert min(1., 45*sensitivity['robust_upper_p']) == sensitivity['bonferroni_upper_p']
                    max_discrepancy = max(max_discrepancy, abs(observed-sensitivity['score']), abs(lower-sensitivity['observed_lower_bound']))
                    summary['robust_rejections'][j] += int(45*sensitivity['robust_upper_p'] <= .05)
                    if j == 0:
                        summary['stratified_nominal_rejections'] += int(45*sensitivity['nominal_p'] <= .05)
                    checked += 1
        assert summary == report['summary'][half]
        summaries[half] = summary
    confirmation = summaries['right_confirmation']
    gate = (confirmation['stratified_nominal_rejections'] > 0 and confirmation['robust_rejections'][2] >= 1
            and confirmation['robust_rejections'][1] >= confirmation['stratified_nominal_rejections']/2)
    assert gate == report['continuation_gate_passed']
    assert max_discrepancy < 1e-12
    output = dict(verified_sensitivity_results=checked, maximum_statistic_discrepancy=max_discrepancy,
                  source_sha256=report['source_sha256'], summary=summaries, continuation_gate_passed=gate,
                  scope='Separate source replay and arithmetic; not independent biological replication.')
    (ROOT/'results/molecular_confirmation_verification.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
