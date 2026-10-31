"""Structured registration ambiguity and one-image paired-reference feasibility."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.optimize import linprog
from scipy.spatial import cKDTree
from scipy.stats import binom, poisson
from skimage.feature import peak_local_max

from study_barcode_collisions import wilson

ROOT = Path(__file__).resolve().parent
FIDUCIAL_SHA = '3219c1f7dc1b25ae4471ce6453f7119b22c17603053cca4ff6c766915827a247'
ORGANIZATION_SHA = '3292fb375f6a192fd992dbf7d85bc07586b9bf818995710ceb687dfeea25999c'


def write_csv(path, rows):
    if not rows:
        raise ValueError('No records to write')
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def design(rounds, colors, mode):
    n = len(rounds)
    round_matrix = np.array([[float(r == value) for value in range(1, 9)] for r in rounds])
    color_matrix = np.array([[float(c == value) for value in [561, 650]] for c in colors])
    if mode == 'independent_bits':
        return np.eye(n)
    if mode == 'shared_rounds':
        return round_matrix
    if mode == 'rounds_and_colors':
        return np.column_stack([round_matrix, color_matrix])
    if mode == 'known_rounds_colors_only':
        return color_matrix
    raise ValueError(mode)


def projection_witness(matrix, difference):
    coefficients = np.linalg.lstsq(matrix, difference, rcond=None)[0]
    residual = difference-matrix @ coefficients
    return coefficients, residual


def minimum_round_radius(matrix, difference):
    """Per-world l-infinity drift radius; nuisance differences have radius 2r."""
    # Eight round columns, two free color columns, and scalar per-world radius.
    constraints = np.zeros((16, 11))
    for index in range(8):
        constraints[2*index, index] = 1
        constraints[2*index+1, index] = -1
        constraints[2*index:2*index+2, -1] = -2
    objective = np.zeros(11); objective[-1] = 1
    solution = linprog(objective, A_ub=constraints, b_ub=np.zeros(16),
                       A_eq=np.column_stack([matrix, np.zeros(len(matrix))]), b_eq=difference,
                       bounds=[(None, None)]*10+[(0, None)], method='highs')
    if solution.status == 2:
        return None, None
    if not solution.success:
        raise RuntimeError(solution.message)
    assert np.max(abs(matrix @ solution.x[:-1]-difference)) < 1e-8
    assert np.max(abs(solution.x[:8])) <= 2*solution.x[-1]+1e-8
    return float(solution.x[-1]), solution.x[:-1].tolist()


def geometry_audit():
    metadata_path = ROOT/'datasets/spatial_v01/dataorganization.csv'
    if hashlib.sha256(metadata_path.read_bytes()).hexdigest() != ORGANIZATION_SHA:
        raise ValueError('Source acquisition layout checksum mismatch')
    metadata = pd.read_csv(metadata_path).iloc[:16]
    codebook = pd.read_csv(ROOT/'datasets/spatial_v01/codebook.csv')
    if metadata.readoutName.tolist() != codebook.columns.tolist()[2:]:
        raise ValueError('Codebook/readout order mismatch')
    codes = dict(zip(codebook.name, codebook.iloc[:, 2:].to_numpy(dtype=int)))
    previous = json.loads((ROOT/'results/barcode_collisions.json').read_text())['codebook']['genes_only']['collisions']
    rows, witnesses = [], []
    selected_classes = 0
    for group_index, group in enumerate(previous):
        signature = np.array(group['signature'])
        if signature.max() > 1:
            continue
        selected_classes += 1
        active = signature == 1
        rounds = metadata.imagingRound.to_numpy()[active]
        colors = metadata.color.to_numpy()[active]
        for first, second in itertools.combinations(group['worlds'], 2):
            center = (-.5*codes[first[0]]+.5*codes[first[1]])[active]
            for reverse in [False, True]:
                alternative = second[::-1] if reverse else second
                other_center = (-.5*codes[alternative[0]]+.5*codes[alternative[1]])[active]
                difference = center-other_center
                record = dict(class_index=group_index, first=' + '.join(first),
                              second=' + '.join(alternative), reversed=reverse)
                for mode in ['independent_bits', 'shared_rounds', 'rounds_and_colors', 'known_rounds_colors_only']:
                    coefficients, residual = projection_witness(design(rounds, colors, mode), difference)
                    record[mode+'_residual_linf'] = float(abs(residual).max())
                radius, coefficients = minimum_round_radius(design(rounds, colors, 'rounds_and_colors'), difference)
                record['minimum_per_world_round_radius_at_d1'] = radius
                rows.append(record)
                witnesses.append(dict(**record, active_bits=np.flatnonzero(active).tolist(),
                                      required_center_difference=difference.tolist(), lp_coefficients=coefficients))
    modes = ['independent_bits', 'shared_rounds', 'rounds_and_colors', 'known_rounds_colors_only']
    summary = dict(binary_sum_classes=selected_classes, excluded_nonbinary_classes=len(previous)-selected_classes,
                   oriented_world_comparisons=len(rows),
                   exact_aliases_by_model={mode:sum(r[mode+'_residual_linf'] < 1e-10 for r in rows) for mode in modes},
                   comparisons=rows, witnesses=witnesses,
                   target_comparisons=[r for r in rows if r['first'] == 'Adamdec1 + Adgrg6'],
                   fiducial_color=metadata.fiducialColor.unique().tolist(), readout_colors=sorted(metadata.color.unique().tolist()),
                   calibrated_target_registration=False)
    return summary


def detect_centroids(image):
    image = image.astype(float)
    band = np.maximum(gaussian_filter(image, 1)-gaussian_filter(image, 4), 0)
    peaks = peak_local_max(band, min_distance=8, threshold_abs=float(np.quantile(band, .999)), exclude_border=12)
    ordering = np.lexsort((peaks[:, 1], peaks[:, 0], -band[peaks[:, 0], peaks[:, 1]]))
    retained = peaks[ordering[:300]]
    result = []
    for y, x in retained:
        patch = band[y-3:y+4, x-3:x+4]
        yy, xx = np.indices(patch.shape)
        result.append([y-3+float((yy*patch).sum()/patch.sum()), x-3+float((xx*patch).sum()/patch.sum())])
    return np.array(result).reshape(-1, 2), len(peaks)


def fiducial_audit(dataset_root, output):
    path = dataset_root/'spatial/fiducials_oriented_unaligned.npy'
    with path.open('rb') as stream:
        if hashlib.file_digest(stream, 'sha256').hexdigest() != FIDUCIAL_SHA:
            raise ValueError('Frozen fiducial source checksum mismatch')
    images = np.load(path, mmap_mode='r')
    original = json.loads((ROOT/'datasets/spatial_v01/registration.json').read_text())
    detections = [detect_centroids(image) for image in images]
    reference = detections[0][0]
    points, matches, summaries = [], [], []
    for index, (moving, count_before_cap) in enumerate(detections):
        for ident, (y, x) in enumerate(moving):
            points.append(dict(round=index+1, detection_id=ident, y=float(y), x=float(x)))
        shift = np.array([original[index]['shift_y_pixels'], original[index]['shift_x_pixels']])
        distance, ref_ids = cKDTree(reference).query(moving+shift)
        _, reverse = cKDTree(moving+shift).query(reference)
        kept = np.flatnonzero((distance <= 2)&(reverse[ref_ids] == np.arange(len(moving))))
        reference_ids = ref_ids[kept]
        training = ((np.floor(reference[reference_ids, 0]/128)+np.floor(reference[reference_ids, 1]/128)).astype(int)%2) == 0
        observed = reference[reference_ids]-moving[kept]
        evaluable = int(training.sum()) >= 10 and int((~training).sum()) >= 10
        fitted = np.median(observed[training], axis=0) if evaluable else None
        residuals = observed-fitted if evaluable else np.full_like(observed, np.nan)
        for k, moving_id in enumerate(kept):
            matches.append(dict(round=index+1, reference_id=int(reference_ids[k]), moving_id=int(moving_id),
                                training=bool(training[k]), shift_y=float(observed[k, 0]), shift_x=float(observed[k, 1]),
                                residual_y=float(residuals[k, 0]) if evaluable else None,
                                residual_x=float(residuals[k, 1]) if evaluable else None,
                                coarse_match_distance=float(distance[moving_id])))
        validation = np.linalg.norm(residuals[~training], axis=1) if evaluable else None
        summaries.append(dict(round=index+1, detected_before_cap=count_before_cap, retained=len(moving),
                              matched=len(kept), unmatched=len(moving)-len(kept),
                              train_matches=int(training.sum()), heldout_matches=int((~training).sum()),
                              evaluable=evaluable, identity_reference=index == 0,
                              fitted_shift=None if fitted is None else fitted.tolist(),
                              fitted_minus_original=None if fitted is None else (fitted-shift).tolist(),
                              heldout_median=None if not evaluable else float(np.median(validation)),
                              heldout_q95=None if not evaluable else float(np.quantile(validation, .95)),
                              heldout_max=None if not evaluable else float(validation.max())))
    write_csv(output/'fiducial_detections.csv', points)
    write_csv(output/'fiducial_matches.csv', matches)
    return dict(source_sha256=FIDUCIAL_SHA, source_doi='10.5281/zenodo.15678603', field='000',
                rounds=summaries, residuals_are_calibrated_error_bounds=False,
                chromatic_error_measured=False, independent_validation=False,
                matching_uses_prior_full_image_shift=True, residuals_truncated_by_matching=True)


def paired_decision(target, background, delta, efficiency_min, area=16):
    total = target+background
    p0 = (1+delta)/(area+1+delta)
    p1 = ((1-delta)*10+efficiency_min)/((area+1-delta)*10+efficiency_min)
    present = binom.sf(target-1, total, p0) <= .025
    absent = binom.cdf(target, total, p1) <= .025
    return np.where(present&~absent, 1, np.where(absent&~present, 0, -1))


def target_only_decision(target, dose, delta, efficiency_min):
    max_absent = (1+delta)*10
    min_present = (1-delta)*2+efficiency_min
    present = poisson.sf(target-1, dose*max_absent) <= .025
    absent = poisson.cdf(target, dose*min_present) <= .025
    return np.where(present&~absent, 1, np.where(absent&~present, 0, -1))


def probe_screen(output):
    rng = np.random.default_rng(2026103003)
    rows, saved, repeats = [], {}, 2000
    for b, delta, eta, dose, truth in itertools.product([2., 10.], [0., .1, .25], [.5, 2., 8.], [1, 4, 16, 64], [0, 1]):
        mismatch = delta if truth == 0 else -delta
        x = rng.poisson(dose*(b*(1+mismatch)+eta*truth), repeats)
        y = rng.poisson(dose*16*b, repeats)
        key = f'b{b}_delta{delta}_eta{eta}_E{dose}_z{truth}'
        saved[key] = np.column_stack([x, y]).astype(np.int32)
        for name, decision in [('paired_reference', paired_decision(x, y, delta, eta)),
                               ('target_only', target_only_decision(x, dose, delta, eta))]:
            correct = decision == truth
            incorrect = (decision >= 0)&(decision != truth)
            rows.append(dict(background=b, delta=delta, efficiency_min=eta, dose=dose, truth=truth, method=name,
                             correct=float(correct.mean()), incorrect=float(incorrect.mean()), abstain=float((decision < 0).mean()),
                             correct_wilson=wilson(int(correct.sum()), repeats),
                             incorrect_wilson=wilson(int(incorrect.sum()), repeats)))
    x = rng.poisson(64*15., repeats); y = rng.poisson(64*16*10., repeats)
    saved['uncalibrated_mismatch_null'] = np.column_stack([x, y]).astype(np.int32)
    stress = paired_decision(x, y, 0., 8.)
    paired = [r for r in rows if r['method'] == 'paired_reference']
    strong = [r for r in paired if r['efficiency_min'] == 8 and r['dose'] == 64]
    gates = dict(incorrect_decision_control=all(r['incorrect_wilson'][1] <= .075 for r in paired),
                 strong_signal_decisiveness=all(r['correct_wilson'][0] >= .80 for r in strong))
    np.savez_compressed(output/'paired_reference_counts.npz', **saved)
    write_csv(output/'paired_reference_summary.csv', rows)
    return dict(seed=2026103003, repetitions=repeats, settings=rows, gates=gates,
                conditional_feasibility_pass=all(gates.values()),
                out_of_contract_null_false_presence=float((stress == 1).mean()),
                real_probe_acquired=False, calibrated_photon_units=False, distinct_new_test=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset-root', type=Path, required=True)
    args = parser.parse_args()
    output = ROOT/'results'; output.mkdir(exist_ok=True)
    geometry = geometry_audit()
    fiducials = fiducial_audit(args.dataset_root, output)
    probe = probe_screen(output)
    detail_path = output/'fiducial_geometry_witnesses.json.gz'
    details = {k:geometry[k] for k in ['comparisons', 'witnesses']}
    detail_path.write_bytes(gzip.compress(json.dumps(details, sort_keys=True, allow_nan=False).encode(), mtime=0))
    geometry_summary = {k:v for k,v in geometry.items() if k not in details}
    geometry_summary['details_file'] = detail_path.name
    geometry_summary['details_sha256'] = hashlib.sha256(detail_path.read_bytes()).hexdigest()
    result = dict(protocol_commit='ee47e68', geometry=geometry_summary, fiducials=fiducials, probe=probe,
                  real_one_image_claim_settled=False)
    (output/'fiducial_measurement.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(geometry={k:v for k,v in geometry.items() if k not in ['comparisons', 'witnesses']},
                          fiducials=fiducials, probe={k:v for k,v in probe.items() if k != 'settings'}), indent=2))
