"""Post-outcome algebraic diagnostic; no replacement hypothesis tests."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parent


def main():
    report = json.loads((ROOT/'results/molecular_confirmation.json').read_text())
    xy = pd.read_csv(ROOT/'data/visp/s3_spot_table.csv')[['rotated_x', 'rotated_y']].to_numpy()
    rows = []
    for result in report['rows']:
        if result['status'] != 'evaluated':
            continue
        selected = [next(x['source_row_indices'] for x in report['selection']
                         if x['region'] == result['region'] and x['gene'] == g) for g in result['genes']]
        assert len(selected[0]) == len(selected[1])
        m = len(selected[0]); h = result['bandwidth_um']
        distance = squareform(pdist(xy[np.concatenate(selected)]))
        kernel = np.exp(-distance**2/(2*h*h))
        sign = np.r_[np.ones(m), -np.ones(m)]
        mmd2 = float(sign @ kernel @ sign/m**2)
        score = float(kernel[:m, m:].mean())
        total = float(kernel.sum())
        assert abs(score - (total/(4*m*m)-mmd2/4)) < 1e-12
        assert mmd2 >= -1e-12
        mean = (total-2*m)/(2*m*(2*m-1))
        assert score-mean <= 1/(2*m-1)+1e-12
        rows.append(dict(region=result['region'], genes=result['genes'], bandwidth_um=h,
                         biased_mmd_squared=mmd2, score=score, global_permutation_mean=mean,
                         excess_over_global_mean=score-mean, positive_excess_ceiling=1/(2*m-1),
                         one_swap_allowance=2/m,
                         zero_lower_at_half_um_one_swap=result['sensitivity'][2]['observed_lower_bound'] == 0.))
    output = dict(scope='Post-outcome diagnostic based on established MMD algebra; no new p-values or confirmation claim.',
                  rows=rows, below_global_mean=sum(r['excess_over_global_mean'] < 0 for r in rows),
                  zero_lower_at_half_um_one_swap=sum(r['zero_lower_at_half_um_one_swap'] for r in rows))
    (ROOT/'results/molecular_null_diagnostic.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
