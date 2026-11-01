"""Render fixed calibration and Gaussian-information results."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
reference = json.loads((root/'results/cross_channel_calibration.json').read_text())
information = json.loads((root/'results/photon_ceiling.json').read_text())
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
primary = [r for r in reference['settings'] if r['mode'] == 'hashed']
x = np.arange(2)
axes[0].bar(x-.17, [r['median_test_group_residual'] for r in primary], .34, label='Median held-out feature error', color='#397c8c')
axes[0].bar(x+.17, [r['radius'] for r in primary], .34, label='95% calibration radius', color='#be9153')
axes[0].axhline(1, ls='--', color='black', lw=1, label='Frozen radius limit')
axes[0].set_xticks(x, ['Translation', 'Affine'])
axes[0].set_ylabel('Source coordinate units')
axes[0].set_ylim(0, 3.6)
axes[0].set_title('A  Lower typical error; primary radius gate fails', loc='left', fontsize=10)
axes[0].legend(fontsize=8, loc='upper left')
for n, color in [(1, '#be9153'), (16, '#397c8c')]:
    chosen = [r for r in information['cases'] if r['case'] == 0 and r['separation'] == .25 and r['fiducial_sd'] == .3 and r['fiducial_references'] == n]
    axes[1].plot([r['exposure'] for r in chosen], [100*r['optimal_error'] for r in chosen], 'o-', color=color, label=f'{n} independent reference(s)/round')
    axes[1].axhline(100*chosen[0]['ceiling_error'], color=color, ls=':', lw=1)
axes[1].axhline(5, color='black', ls='--', lw=1, label='5% ideal decision error')
axes[1].set_xscale('log'); axes[1].set_xticks([1, 16, 256], ['1', '16', '256'])
axes[1].set_xlabel('Hypothetical gene-photon multiplier')
axes[1].set_ylabel('Minimax classification error (%)')
axes[1].set_title('B  Fixed calibration can limit photon gains', loc='left', fontsize=10)
axes[1].legend(fontsize=8)
fig.suptitle('One selected reference field; Gaussian localization limit uses assumed known variances', fontsize=10)
fig.savefig(root/'results/submission_research.png', dpi=180)
