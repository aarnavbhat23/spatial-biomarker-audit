"""Conditional Gaussian localization information and a calibration ceiling."""
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from study_fiducial_measurement import design

ROOT = Path(__file__).resolve().parent


def calibrated_information(h, ar, ac, gene_variance, fiducial_variance):
    covariance = gene_variance*np.eye(len(h))+fiducial_variance*ar @ ar.T
    inv_h = np.linalg.solve(covariance, h)
    inv_a = np.linalg.solve(covariance, ac)
    nuisance = np.linalg.lstsq(ac.T @ inv_a, ac.T @ inv_h, rcond=None)[0]
    residual = h-ac @ nuisance
    weight = np.linalg.solve(covariance, residual)
    return float(residual @ weight), weight


def joint_information(h, ar, ac, gene_variance, fiducial_variance):
    # Independent implementation in the full joint observation space.
    signal = np.r_[h, np.zeros(ar.shape[1])]
    nuisance = np.block([[ar, ac], [np.eye(ar.shape[1]), np.zeros((ar.shape[1], ac.shape[1]))]])
    std = np.r_[np.full(len(h), np.sqrt(gene_variance)), np.full(ar.shape[1], np.sqrt(fiducial_variance))]
    whitened = signal/std
    matrix = nuisance/std[:, None]
    fitted = matrix @ np.linalg.lstsq(matrix, whitened, rcond=None)[0]
    return float(np.sum((whitened-fitted)**2))


def information_ceiling(h, ar, ac, fiducial_variance):
    projection = np.eye(len(h))-ac @ np.linalg.pinv(ac)
    matrix = projection @ ar
    u = np.linalg.lstsq(matrix, projection @ h, rcond=None)[0]
    if np.linalg.norm(matrix @ u-projection @ h) > 1e-9:
        return np.inf
    return float(u @ u/fiducial_variance)


def cases():
    source = json.loads((ROOT/'results/fiducial_measurement.json').read_text())
    details = json.loads(gzip.decompress((ROOT/'results'/source['geometry']['details_file']).read_bytes()))
    return [w for w in details['witnesses'] if w['first'] == 'Adamdec1 + Adgrg6']


def run():
    rng = np.random.default_rng(2026110102)
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    rows, scores = [], {}
    for case_index, case in enumerate(cases()):
        selected = layout.iloc[case['active_bits']]
        ar = design(selected.imagingRound, selected.color, 'shared_rounds')
        ac = design(selected.imagingRound, selected.color, 'known_rounds_colors_only')
        for d, tau, n, exposure in itertools.product([.25, .5, 1.], [.1, .3], [1, 16], [1, 16, 256]):
            h = d*np.array(case['required_center_difference'])
            vg, vf = .09/exposure, tau*tau/n
            squared, weight = calibrated_information(h, ar, ac, vg, vf)
            joint = joint_information(h, ar, ac, vg, vf)
            ceiling = information_ceiling(h, ar, ac, vf)
            row = dict(case=case_index, first=case['first'], alternative=case['second'], separation=d,
                       fiducial_sd=tau, fiducial_references=n, exposure=exposure,
                       squared_information=squared, joint_squared_information=joint, ceiling_information=ceiling,
                       optimal_error=float(norm.cdf(-np.sqrt(squared)/2)),
                       ceiling_error=float(norm.cdf(-np.sqrt(ceiling)/2)),
                       probe_error=float(norm.cdf(-np.sqrt(squared+16)/2)),
                       practical_error_threshold=.05)
            if case_index == 0:
                errors, probe_errors = [], []
                for truth in [0, 1]:
                    noise = rng.normal(size=(2000, len(h)))*np.sqrt(vg)
                    fiducials = rng.normal(size=(2000, 8))*np.sqrt(vf)
                    observation = truth*h+noise-fiducials @ ar.T
                    statistic = observation @ weight-squared/2
                    probe = rng.normal(4*truth, 1., size=2000)
                    decision = statistic > 0
                    extra_decision = statistic+4*probe-8 > 0
                    errors.append(float((decision != truth).mean()))
                    probe_errors.append(float((extra_decision != truth).mean()))
                    key = f'd{d}_tau{tau}_n{n}_E{exposure}_truth{truth}'
                    scores[key] = np.column_stack([statistic, probe])
                row['empirical_errors'] = errors
                row['empirical_probe_errors'] = probe_errors
            rows.append(row)
    agreement = max(abs(r['squared_information']-r['joint_squared_information']) for r in rows)
    deviations = [abs(error-r[field]) for r in rows for empirical, field in [('empirical_errors', 'optimal_error'), ('empirical_probe_errors', 'probe_error')]
                  for error in r.get(empirical, [])]
    result = dict(protocol_commit='fdece4d', seed=2026110102, cases=rows,
                  gates=dict(joint_projection_agreement=agreement <= 1e-8, monte_carlo_agreement=max(deviations) <= .05),
                  maximum_joint_discrepancy=agreement, maximum_empirical_error_deviation=max(deviations),
                  practical_counts=dict(total=len(rows), current_photons=sum(r['optimal_error'] <= .05 for r in rows),
                                        infinite_photons=sum(r['ceiling_error'] <= .05 for r in rows),
                                        added_calibrated_probe=sum(r['probe_error'] <= .05 for r in rows)),
                  calibrated_real_variances=False, new_decoder=False, independent_biological_validation=False)
    (ROOT/'results/photon_ceiling.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    np.savez_compressed(ROOT/'results/photon_ceiling_statistics.npz', **scores)
    print(json.dumps({k:v for k,v in result.items() if k != 'cases'}, indent=2))


if __name__ == '__main__':
    run()
