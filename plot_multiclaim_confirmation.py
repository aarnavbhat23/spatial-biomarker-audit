"""Display frozen field-count gains and the latency boundary."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/multiclaim_confirmation.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
patients=list(r['patients']);policies=['envelope','mixture','focused_reuse','greedy']
colors=['#146a87','#ba692d','#6a8598','#737373']
for k,(policy,color) in enumerate(zip(policies,colors)):
    heights=[next(e for e in r['patients'][p]['methods'] if e['policy']==policy)['mean_fields'] for p in patients]
    axes[0].bar(np.arange(3)+(k-1.5)*.2,heights,width=.2,color=color,label=policy.replace('_',' '))
axes[0].set(xticks=range(3),xticklabels=patients,ylabel='Mean physical fields',title='All three claims reach width .005',ylim=(0,220))
axes[0].legend(fontsize=8,ncol=2,loc='upper center');axes[0].grid(axis='y',alpha=.2)
h=np.linspace(0,40,200)
for patient in patients:
    es=r['patients'][patient]['methods'];a=next(e for e in es if e['policy']=='envelope');g=next(e for e in es if e['policy']=='greedy')
    difference=a['mean_fields']+h*a['mean_rounds']-g['mean_fields']-h*g['mean_rounds']
    axes[1].plot(h,difference,label=patient)
axes[1].axhline(0,color='black',lw=1);axes[1].set(xlabel='Feedback-round cost (field exposures)',ylabel='Envelope cost minus greedy cost',title='Below zero favors envelope; above favors greedy')
axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
fig.suptitle('Frozen confirmation on three new ORION patients; hypothetical 1024-pixel fields',fontsize=12)
fig.savefig(ROOT/'results/multiclaim_confirmation.png',dpi=170)
