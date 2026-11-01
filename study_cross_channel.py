"""Frozen whole-feature registration and split-conformal reference benchmark."""
import hashlib
import gzip
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from study_barcode_collisions import wilson

ROOT = Path(__file__).resolve().parent


def source_groups(frame, frozen=None):
    parent = {}
    def find(node):
        parent.setdefault(node, node)
        if parent[node] != node:
            parent[node] = find(parent[node])
        return parent[node]
    edges = []
    for row in frame.to_numpy():
        first = (int(row[16]), round(float(row[0]), 8), round(float(row[1]), 8))
        second = (int(row[17]), round(float(row[8]), 8), round(float(row[9]), 8))
        parent[find(second)] = find(first)
        edges.append((first, second))
    groups = defaultdict(list)
    for node in parent:
        groups[find(node)].append(node)
    ids, group_nodes = {}, {}
    anchors = None if frozen is None else {tuple(sorted((int(c), round(x, 8), round(y, 8)) for c, x, y in r['nodes'])):r for r in frozen}
    for nodes in groups.values():
        if anchors is None:
            encoded = ';'.join(f'{c}:{x.hex()}:{y.hex()}' for c, x, y in sorted(nodes))
            ident = hashlib.sha256(encoded.encode()).hexdigest()
            group_nodes[ident] = sorted(nodes)
        else:
            anchor = anchors[tuple(sorted(nodes))]
            ident = anchor['group_id']
            group_nodes[ident] = anchor['nodes']
        for node in nodes:
            ids[node] = ident
    return [ids[a] for a, _ in edges], group_nodes


def partition(nodes, mode):
    if mode == 'hashed':
        ordered = sorted(nodes)
    else:
        ordered = sorted(nodes, key=lambda key:(np.mean([n[1] for n in nodes[key]]), key))
    n = len(ordered)
    return {key:('train' if i < n//2 else 'calibration' if i < n//2+n//4 else 'test')
            for i, key in enumerate(ordered)}


def conformal_radius(scores, coverage=.95):
    k = int(np.ceil((len(scores)+1)*coverage))
    return float(np.sort(scores)[k-1]) if k <= len(scores) else float('inf')


def fit_transform(x, y, method):
    if method == 'translation':
        return np.median(y-x, axis=0)
    matrix = np.column_stack([x, np.ones(len(x))])
    if np.linalg.matrix_rank(matrix) < 3:
        raise ValueError('Affine training geometry is rank deficient')
    return np.linalg.lstsq(matrix, y, rcond=None)[0]


def predict_transform(x, parameters, method):
    return x+parameters if method == 'translation' else np.column_stack([x, np.ones(len(x))]) @ parameters


def run():
    path = ROOT/'datasets/reference_calibration_v01/cross_channel_matches.csv'
    if hashlib.sha256(path.read_bytes()).hexdigest() != '7ebb373a323bdc3254745c3083dd62ccaf559876596c704dacfdfc2fceaf32ab':
        raise ValueError('Source mismatch')
    frame = pd.read_csv(path)
    np.testing.assert_array_equal(frame.iloc[:, :3], frame.iloc[:, 4:7])
    frozen = json.loads((ROOT/'datasets/reference_calibration_v01/frozen_feature_groups.json').read_text())
    if frozen['source_sha256'] != hashlib.sha256(path.read_bytes()).hexdigest():
        raise ValueError('Frozen feature manifest belongs to a different source')
    group_ids, nodes = source_groups(frame, frozen['groups'])
    invalid = [key for key, group in nodes.items() if len(group) != len(set(n[0] for n in group))]
    if invalid or len(nodes) < 60:
        raise ValueError('Frozen identity/group eligibility failed')
    frame['group_id'] = group_ids
    summaries, predictions, splits = [], [], []
    for mode in ['hashed', 'spatial_extrapolation']:
        split = partition(nodes, mode)
        for key in sorted(nodes):
            splits.append(dict(mode=mode, group_id=key, split=split[key], nodes=nodes[key]))
        for method in ['translation', 'affine']:
            scored = []
            for (ref, comp), subset in frame.groupby(['ref_ch', 'comp_ch'], sort=True):
                training = subset.group_id.map(split).to_numpy() == 'train'
                if training.sum() < 10:
                    raise ValueError('Insufficient training matches')
                x = subset[['comp_x', 'comp_y']].to_numpy()
                y = subset[['ref_x', 'ref_y']].to_numpy()
                parameters = fit_transform(x[training], y[training], method)
                estimate = predict_transform(x, parameters, method)
                residual = np.linalg.norm(estimate-y, axis=1)
                for k, (index, row) in enumerate(subset.iterrows()):
                    record = dict(mode=mode, method=method, source_row=int(index), group_id=row.group_id,
                                  split=split[row.group_id], ref_channel=int(ref), comp_channel=int(comp),
                                  predicted_x=float(estimate[k, 0]), predicted_y=float(estimate[k, 1]),
                                  residual=float(residual[k]))
                    scored.append(record); predictions.append(record)
            group_scores = defaultdict(list)
            for row in scored:
                group_scores[row['group_id']].append(row['residual'])
            scores = {key:max(value) for key, value in group_scores.items()}
            calibration = [scores[key] for key in scores if split[key] == 'calibration']
            test = np.array([scores[key] for key in scores if split[key] == 'test'])
            radius = conformal_radius(calibration)
            covered = int((test <= radius).sum())
            summaries.append(dict(mode=mode, method=method, train_groups=sum(v == 'train' for v in split.values()),
                                  calibration_groups=len(calibration), test_groups=len(test), radius=radius,
                                  test_covered=covered, test_coverage=covered/len(test),
                                  coverage_wilson=wilson(covered, len(test)),
                                  median_test_group_residual=float(np.median(test)),
                                  maximum_test_group_residual=float(test.max())))
    primary = {r['method']:r for r in summaries if r['mode'] == 'hashed'}
    affine, translation = primary['affine'], primary['translation']
    gates = dict(coverage=affine['coverage_wilson'][0] >= .80, precision=affine['radius'] <= 1.,
                 improvement=affine['median_test_group_residual'] <= .8*translation['median_test_group_residual'])
    output = dict(protocol_commit='f4fd37b', source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  source_rows=len(frame), physical_feature_components=len(nodes), invalid_components=invalid,
                  source_position='0', source_fields=1, independently_identified_specimens=None, selected_matches_only=True,
                  source_units='source xy coordinate units; no physical scale inferred',
                  gates=gates, usefulness_pass=all(gates.values()), settings=summaries,
                  splits=splits, predictions=predictions, biological_truth_verified=False,
                  probe_efficiency_calibrated=False)
    details = {key:output.pop(key) for key in ['splits', 'predictions']}
    detail_path = ROOT/'results/cross_channel_calibration_details.json.gz'
    detail_path.write_bytes(gzip.compress(json.dumps(details, separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    output['details_file'] = detail_path.name
    output['details_sha256'] = hashlib.sha256(detail_path.read_bytes()).hexdigest()
    (ROOT/'results/cross_channel_calibration.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(groups=len(nodes), gates=gates, settings=summaries), indent=2))


if __name__ == '__main__':
    run()
