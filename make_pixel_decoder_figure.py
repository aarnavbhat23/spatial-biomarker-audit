"""Scientific figure for the fixed full-pixel comparison."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from study_photon_ceiling import cases
from study_pixel_decoder import inputs

ROOT = Path(__file__).resolve().parent
study = json.loads((ROOT/'results/pixel_decoder.json').read_text())
names, codes, _, mapping = inputs()
case = cases()[0]
first, second = [case[k].split(' + ') for k in ['first', 'second']]
bits = case['active_bits']
coef = case['lp_coefficients']
center0 = np.array([-.5 if codes[names.index(first[0]), b] else .5 for b in bits])
center1 = np.array([-.5 if codes[names.index(second[0]), b] else .5 for b in bits])
shift = np.array([coef[mapping[b][0]]+coef[8+mapping[b][1]] for b in bits])
observed = center0-shift/2
assert np.allclose(observed, center1+shift/2)
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout='constrained')
for values, color, label, marker in [(center0, '#2364aa', 'World A, before shifts', 'o'),
                                     (center1, '#b85b26', 'World B, before shifts', 's'),
                                     (observed, '#222222', 'Both worlds, observed', 'x')]:
    axes[0].scatter(values, np.arange(8), color=color, marker=marker, label=label, s=45)
axes[0].set(yticks=np.arange(8), yticklabels=[f'bit {b+1}' for b in bits], xlabel='Emitter x position (d=1 pixel)',
            title='A  Different explanations, identical gene images')
axes[0].legend(fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.2), frameon=False)
errors = [.5, 1-study['registered_control_accuracy'], 1-study['decoder_probe_accuracy']]
axes[1].bar(range(3), errors, color=['#777777', '#2364aa', '#49906e'], width=.6)
axes[1].scatter(range(3), errors, color='black', s=15)
for x, error in enumerate(errors):
    axes[1].text(x, error+.02, f'{100*error:.0f}%', ha='center')
axes[1].set(xticks=range(3), xticklabels=['Ambiguous\ngene pixels', 'Known registration\n(high signal)', 'Extra specific\nprobe'],
            ylim=(0, .64), ylabel='Forced claim error', title='B  Published ISTDECO on controlled images')
fig.suptitle('Restricted equal-brightness construction; these are not tissue accuracy estimates', fontsize=11)
fig.savefig(ROOT/'results/pixel_decoder.png', dpi=190)
fig.savefig(ROOT/'results/pixel_decoder.svg')
