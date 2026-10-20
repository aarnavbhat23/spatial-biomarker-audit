"""Expose the physical-field versus feedback-round tradeoff."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
report=json.loads((ROOT/'results/batch_acquisition_study.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
for ax,(patient,entries) in zip(axes,report['patients'].items()):
    for policy,color,marker in (('batch_pps','#146a87','o'),('greedy_certificate','#ba692d','s')):
        rows=[e for e in entries if e['field_native_pixels']==1024 and e['policy']==policy]
        ax.plot([e['mean_rounds'] for e in rows],[e['mean_fields'] for e in rows],marker=marker,
                color=color,label='Sampling' if policy=='batch_pps' else 'Deterministic')
        for e in rows:
            batch=e['batch_size'] or e['greedy_batch']
            ax.annotate(str(batch),(e['mean_rounds'],e['mean_fields']),xytext=(4,5),textcoords='offset points',fontsize=8)
    ax.set(title=patient,xscale='log',xlabel='Feedback rounds',ylabel='Physical fields')
    ax.grid(alpha=.2)
axes[0].legend(fontsize=8)
fig.suptitle('1024-pixel fields; width 0.005 — labels are batch sizes, lower-left is cheaper',fontsize=12)
fig.savefig(ROOT/'results/batch_acquisition_study.png',dpi=170)
