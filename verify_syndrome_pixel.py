"""Check fitted-point readouts against exhaustive small-network optimization."""
from collections import Counter
import itertools
import gzip
import hashlib
import json
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def verify():
    study = json.loads((ROOT/'results/syndrome_pixel.json').read_text())
    book = pd.read_csv(ROOT/'datasets/spatial_v01/codebook.csv')
    codes = book.iloc[:, 2:].to_numpy(dtype=int)
    names = book.name.tolist()
    comparison = ROOT/'comparators/white_pixel'
    assert hashlib.sha256(gzip.decompress((comparison/'inputs.json.gz').read_bytes())).hexdigest() == study['inputs_sha256']
    manifest = tomllib.loads((comparison/'Manifest.toml').read_text())['deps']
    for name, commit in study['package_commits'].items():
        assert manifest[name][0]['repo-rev'] == commit
        saved = next(p for p in study['packages'] if p['name'] == name)
        assert saved['tree_hash'] == manifest[name][0]['git-tree-sha1']
    check = np.loadtxt(comparison/'parity_check.csv', delimiter=',', dtype=int)
    assert check.shape == (5, 16) and np.all(check @ codes.T % 2 == 0)
    checked, candidate_counts = 0, []
    for row in study['outputs']:
        counts = Counter(p['block'] for p in row['fitted_points'])
        cap = 40 if row['mode'] == 'measured' else 20
        saturated = [bit for bit, n in sorted(counts.items()) if n >= cap]
        assert row['saturated_bits'] == saturated
        if row['mode'] == 'measured':
            assert saturated == list(range(1, 17))
            continue
        assert len(row['fitted_points']) == 8 and max(counts.values()) == 1
        points = {p['block']-1:p for p in row['fitted_points']}
        candidates = []
        for name, code in zip(names, codes):
            active = frozenset(np.flatnonzero(code))
            if not active.issubset(points):
                continue
            xy = np.array([[points[b]['x'], points[b]['y']] for b in active])
            if np.linalg.norm(xy[:, None]-xy[None, :], axis=2).max() > 4+1e-7:
                continue
            brightness = np.log2([points[b]['w'] for b in active])
            cost = 7*np.var(xy, axis=0, ddof=1).sum()+3.2*np.var(brightness, ddof=1)
            candidates.append((name, active, cost))
        # Each source word uses four of the eight detected dots, so at most
        # two disjoint words can be chosen. Enumerate the complete feasible set.
        choices = [((), 8.)]
        for size in [1, 2]:
            for selected in itertools.combinations(candidates, size):
                used = set().union(*(c[1] for c in selected))
                if len(used) == 4*size:
                    choices.append((tuple(sorted(c[0] for c in selected)), 8-4*size+sum(c[2] for c in selected)))
        chosen = tuple(sorted(row['decoded_genes']))
        cost = next(value for genes, value in choices if genes == chosen)
        minimum = min(value for _, value in choices)
        np.testing.assert_allclose(cost, minimum, atol=1e-5)
        if row['mode'] == 'alias':
            assert sum(row['paired_errors']) == 1
        else:
            assert row['exact_pair_recovered']
        checked += 1; candidate_counts.append(len(candidates))
    assert checked == 12
    convergence = json.loads((ROOT/'results/syndrome_pixel_convergence.json').read_text())
    assert convergence['completed'] is False and convergence['convergence_gate'] is False
    raw = gzip.decompress((comparison/'convergence_inputs.json.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest() == convergence['input_sha256']
    original_input = json.loads(gzip.decompress((comparison/'inputs.json.gz').read_bytes()))
    changed_input = json.loads(raw)
    assert changed_input['protocol_commit'] == '4c10382'
    for before, after in zip(original_input['stacks'], changed_input['stacks']):
        expected = dict(before)
        if before['mode'] == 'measured':
            expected['max_iters'] = 200
        assert expected == after
    log = (ROOT/convergence['log_file']).read_bytes()
    assert hashlib.sha256(log).hexdigest() == convergence['log_sha256']
    logged = gzip.decompress(log).decode()
    assert 'There are currently 0 solution(s)' in logged
    assert logged.count('Hit max iters in frank-wolfe!') == convergence['reported_fit_cap_messages']
    result = dict(exhaustively_verified_synthetic_networks=checked, candidate_counts=candidate_counts,
                  measured_bits_hitting_cap=16, measured_accuracy_evaluated=False,
                  extended_run_completed=False, outside_expert_feedback=False)
    (ROOT/'results/syndrome_pixel_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    verify()
