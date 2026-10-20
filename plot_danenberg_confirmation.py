"""Display every frozen independent-cohort graph/budget setting without refitting."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/danenberg_confirmation.json').read_text())['rows']
patients=list(dict.fromkeys(x['patient'] for x in rows))
fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharey=True)
for axis,fraction in zip(axes,(.001,.005)):
    for index,patient in enumerate(patients):
        selected=[x for x in rows if x['patient']==patient and x['fraction_each_marker']==fraction]
        for offset,row in zip((-.2,0,.2),selected):
            x=index+offset
            degree=np.diff(row['degree_outer_bounds'])[0]*100
            outer=np.diff(row['bounds']['outer_bounds'])[0]*100
            attained=np.diff(row['bounds']['attained_bounds'])[0]*100
            first=index==0 and offset==-.2
            axis.plot([x,x],[attained,outer],color='#126c8a',lw=3)
            axis.scatter(x,degree,color='#b17635',marker='x',s=27,label='Degree outer width' if first else None)
            axis.scatter(x,outer,color='#126c8a',s=18,label='Fast outer width' if first else None)
            axis.scatter(x,attained,facecolors='white',edgecolors='#126c8a',s=18,label='Attainable width' if first else None)
    axis.axhline(.5,color='#933c45',ls='--',lw=1,label='Target width: 0.5 pp')
    axis.set_title(f'{fraction*100:g}% error allowance per marker')
    axis.set_xticks(range(len(patients)),[x.replace('MB-','') for x in patients])
    axis.set_xlabel('Patient (within each: k=6, 12, 24)')
    axis.grid(axis='y',alpha=.2)
axes[0].set_ylim(0,max(np.diff(x['degree_outer_bounds'])[0]*100 for x in rows)*1.15)
axes[0].set_ylabel('Pair-fraction interval width (percentage points)')
axes[1].legend(loc='upper left',fontsize=8)
fig.suptitle('Independent IMC cohort: fixed graph bounds with hypothetical label errors',fontsize=12)
fig.text(.5,.012,'Source phenotypes are proxy annotations; no assay accuracy or acquisition saving is validated.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,.96))
fig.savefig(ROOT/'results/danenberg_confirmation.png',dpi=180)
