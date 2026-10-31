"""Frozen conventional joint-count control after the conditional-test failure."""
import itertools
import json
import numpy as np
from scipy.stats import gamma

from study_fiducial_measurement import ROOT, paired_decision, target_only_decision, write_csv
from study_barcode_collisions import wilson


def poisson_bounds(count, exposure):
    count = np.asarray(count)
    lower = gamma.ppf(.0125, np.maximum(count, 1))/exposure
    lower = np.where(count == 0, 0., lower)
    upper = gamma.ppf(.9875, count+1)/exposure
    return lower, upper


def joint_flags(x, y, dose, delta, efficiency_min):
    lx, ux = poisson_bounds(x, dose)
    lb, ub = poisson_bounds(y, 16*dose)
    return lx-(1+delta)*ub > 0, ux-(1-delta)*lb < efficiency_min


def joint_decision(x, y, dose, delta, efficiency_min):
    present, absent = joint_flags(x, y, dose, delta, efficiency_min)
    return np.where(present&~absent, 1, np.where(absent&~present, 0, -1))


def decisions(x, y, dose, delta, eta):
    return {'joint_rectangle': joint_decision(x, y, dose, delta, eta),
            'conditional_reference': paired_decision(x, y, delta, eta),
            'target_only': target_only_decision(x, dose, delta, eta)}


def screen():
    rng = np.random.default_rng(2026103004)
    repeats, rows, saved = 2000, [], {}
    for b, delta, eta, dose, truth in itertools.product([2., 10.], [0., .1, .25], [.5, 2., 8.], [1, 4, 16, 64], [0, 1]):
        mismatch = delta if truth == 0 else -delta
        x = rng.poisson(dose*(b*(1+mismatch)+eta*truth), repeats)
        y = rng.poisson(dose*16*b, repeats)
        saved[f'b{b}_delta{delta}_eta{eta}_E{dose}_z{truth}'] = np.column_stack([x, y]).astype(np.int32)
        for method, decision in decisions(x, y, dose, delta, eta).items():
            correct = decision == truth
            wrong = (decision >= 0)&(decision != truth)
            rows.append(dict(background=b, delta=delta, efficiency_min=eta, dose=dose, truth=truth, method=method,
                             correct=float(correct.mean()), incorrect=float(wrong.mean()), abstain=float((decision < 0).mean()),
                             correct_wilson=wilson(int(correct.sum()), repeats), incorrect_wilson=wilson(int(wrong.sum()), repeats)))
    x = rng.poisson(64*15., repeats); y = rng.poisson(64*16*10., repeats)
    saved['uncalibrated_mismatch_null'] = np.column_stack([x, y]).astype(np.int32)
    stress = {}
    for method, decision in decisions(x, y, 64, 0., 8.).items():
        stress[method] = dict(false_presence=float((decision == 1).mean()),
                              correct_absence=float((decision == 0).mean()), abstain=float((decision < 0).mean()))
    present, absent = joint_flags(x, y, 64, 0., 8.)
    stress['joint_rectangle']['both_flags'] = float((present&absent).mean())
    gates = {}
    for method in decisions(x[:1], y[:1], 64, 0., 8.):
        selected = [r for r in rows if r['method'] == method]
        strong = [r for r in selected if r['efficiency_min'] == 8 and r['dose'] == 64]
        gates[method] = dict(incorrect_decision_control=all(r['incorrect_wilson'][1] <= .075 for r in selected),
                             strong_signal_decisiveness=all(r['correct_wilson'][0] >= .80 for r in strong),
                             minimum_strong_correct=min(r['correct'] for r in strong),
                             minimum_strong_correct_wilson_lower=min(r['correct_wilson'][0] for r in strong),
                             maximum_incorrect=max(r['incorrect'] for r in selected))
    output = ROOT/'results'
    np.savez_compressed(output/'paired_reference_joint_counts.npz', **saved)
    write_csv(output/'paired_reference_joint_summary.csv', rows)
    result = dict(protocol_commit='0a55e8b', seed=2026103004, repetitions=repeats, settings=rows,
                  gates=gates, out_of_contract_stress=stress, real_probe_acquired=False,
                  independent_biological_validation=False, distinct_new_test=False)
    (output/'paired_reference_joint.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(gates=gates, stress=stress), indent=2))


if __name__ == '__main__':
    screen()
