"""Conditional mean-center ambiguity given actual selected bead displacements."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from study_fiducial_measurement import ROOT, design


def translation_intervals(observed, allowance):
    return np.max(observed, axis=0)-allowance, np.min(observed, axis=0)+allowance


def maximum_separation(matrix, difference, widths):
    # The last variable is the source separation. Two feasible translations
    # from the same observed interval can differ by at most its full width.
    constraint = np.column_stack([matrix, -np.asarray(difference)])
    bounds = [(-float(w), float(w)) for w in widths]+[(None, None)]*2+[(0, None)]
    objective = np.zeros(11); objective[-1] = -1
    answer = linprog(objective, A_eq=constraint, b_eq=np.zeros(len(difference)), bounds=bounds, method='highs')
    if not answer.success:
        raise RuntimeError(answer.message)
    return dict(maximum_alias_separation=float(answer.x[-1]), primal=answer.x.tolist(),
                equality_dual=answer.eqlin.marginals.tolist(), lower_dual=answer.lower.marginals.tolist(),
                upper_dual=answer.upper.marginals.tolist())


def run():
    previous = json.loads((ROOT/'results/fiducial_measurement.json').read_text())
    matched = pd.read_csv(ROOT/'results/fiducial_matches.csv')
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    witnesses = [w for w in previous['geometry']['witnesses'] if w['first'] == 'Adamdec1 + Adgrg6']
    observations = [matched[matched['round'] == r][['shift_y', 'shift_x']].to_numpy() for r in range(2, 9)]
    consistency_minimum = max(float(np.ptp(values, axis=0).max()/2) for values in observations)
    rows = []
    for allowance in [.1, .25, .5, 1.]:
        intervals = [translation_intervals(values, allowance) for values in observations]
        compatible = all(np.all(lo <= hi) for lo, hi in intervals)
        record = dict(allowance_pixels=allowance, contract_compatible=bool(compatible),
                      round_intervals=[dict(round=r, lower=lo.tolist(), upper=hi.tolist())
                                       for r, (lo, hi) in enumerate(intervals, start=2)], comparisons=[])
        if compatible:
            for axis_index, axis in enumerate(['y', 'x']):
                widths = [0.]+[float(hi[axis_index]-lo[axis_index]) for lo, hi in intervals]
                for witness in witnesses:
                    bits = witness['active_bits']
                    selected = layout.iloc[bits]
                    matrix = design(selected.imagingRound, selected.color, 'rounds_and_colors')
                    row = dict(axis=axis, first=witness['first'], second=witness['second'],
                               active_bits=bits, difference=witness['required_center_difference'],
                               round_difference_widths=widths,
                               **maximum_separation(matrix, witness['required_center_difference'], widths))
                    record['comparisons'].append(row)
        rows.append(record)
    result = dict(protocol_commit='2e91333', minimum_allowance_for_selected_match_consistency=consistency_minimum,
                  calibrated_allowance=False, independent_validation=False, results=rows)
    (ROOT/'results/fiducial_sets.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(consistency_minimum=consistency_minimum,
                         settings=[dict(allowance=r['allowance_pixels'], compatible=r['contract_compatible'],
                                        thresholds=[c['maximum_alias_separation'] for c in r['comparisons']]) for r in rows]), indent=2))


if __name__ == '__main__':
    run()
