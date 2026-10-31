"""Plot retained outputs; no model selection or re-estimation."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
study = json.loads((root/'results/fiducial_measurement.json').read_text())
joint = json.loads((root/'results/paired_reference_joint.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), constrained_layout=True)
geometry = study['geometry']
values = list(geometry['exact_aliases_by_model'].values())
labels = ['Independent\nbit shifts', 'Shared\nround drift', 'Round drift\n+ color offsets', 'Known rounds\n+ color offsets']
axes[0].bar(np.arange(4), np.array(values)/geometry['oriented_world_comparisons']*100, color='#527f9f')
axes[0].set_xticks(np.arange(4), labels, fontsize=8)
axes[0].set_ylabel('Exact aliases among oriented comparisons (%)')
axes[0].set_ylim(0, 111)
for i, value in enumerate(values):
    axes[0].text(i, value/5274*100+2, f'{value:,}', ha='center', fontsize=9)
axes[0].set_title('A  Shared acquisition constraints', loc='left', fontsize=11)
methods = ['joint_rectangle', 'conditional_reference', 'target_only']
correct, errors = [], [[], []]
for method in methods:
    rows = [r for r in joint['settings'] if r['method'] == method and r['efficiency_min'] == 8 and r['dose'] == 64]
    worst = min(rows, key=lambda r:r['correct'])
    correct.append(worst['correct']*100)
    errors[0].append(max(0, worst['correct']-worst['correct_wilson'][0])*100)
    errors[1].append(max(0, worst['correct_wilson'][1]-worst['correct'])*100)
axes[1].bar(np.arange(3), correct, yerr=errors, capsize=4, color=['#287c74', '#a9824b', '#a9824b'])
axes[1].axhline(80, ls='--', color='black', lw=1, label='Required Wilson lower bound')
axes[1].set_xticks(np.arange(3), ['Joint\nrectangle', 'Conditional\nreference', 'Target\nonly'])
axes[1].set_ylim(0, 111)
axes[1].set_ylabel('Worst strong-setting correct decision (%)')
for i, value in enumerate(correct):
    axes[1].text(i, value+3, f'{value:.2f}%', ha='center', fontsize=9)
axes[1].legend(fontsize=7, loc='center right')
axes[1].set_title('B  One hypothetical added image', loc='left', fontsize=11)
rounds = study['fiducials']['rounds'][1:]
axes[2].plot([r['round'] for r in rounds], [r['heldout_q95'] for r in rounds], 'o-', color='#527f9f', label='95th percentile')
axes[2].plot([r['round'] for r in rounds], [r['heldout_median'] for r in rounds], 's-', color='#287c74', label='Median')
axes[2].set_xticks(range(2, 10)); axes[2].set_ylim(0, .25)
axes[2].set_xlabel('Fiducial round (9 = morphology)')
axes[2].set_ylabel('Selected-match held-out residual (pixels)')
axes[2].legend(fontsize=8)
axes[2].set_title('C  Measured 473-nm consistency', loc='left', fontsize=11)
fig.suptitle('Conditional feasibility: error allowance, color transfer and probe calibration remain unvalidated', fontsize=11)
fig.savefig(root/'results/fiducial_measurement.png', dpi=180)
