"""Distinguish attainable ambiguity from looseness of the certified outer bound."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/general_marker_study.json').read_text())['panels']
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,fraction in zip(axes,(.001,.005)):
    selected=[r for r in rows if r['fraction']==fraction];x=np.arange(len(selected))
    outer=[np.diff(r['fast']['outer_bounds'])[0]*100 for r in selected]
    attained=[np.diff(r['fast']['attained_bounds'])[0]*100 for r in selected]
    degree=[np.diff(r['degree_outer_bounds'])[0]*100 for r in selected]
    ax.bar(x,outer,color='lightsteelblue',label='Valid width upper bound')
    ax.bar(x,attained,color='tab:blue',label='Attainable width lower bound')
    ax.plot(x,degree,marker='.',color='0.35',label='Degree outer bound')
    ax.axhline(.5,color='tab:red',linestyle='--',linewidth=1,label='Requested width: 0.5 pp')
    ax.set_xticks(x,[r['patient']+' / '+str(r['radius']) for r in selected],rotation=45,ha='right',fontsize=8)
    ax.set_ylabel('Full-census ambiguity width (percentage points)');ax.set_title(f'Both marker budgets: {fraction*100:g}% each')
    ax.grid(axis='y',alpha=.2)
axes[0].legend(fontsize=7,loc='upper left')
fig.suptitle('Development source graphs: error budgets are sensitivity assumptions',fontsize=12)
fig.savefig(ROOT/'results/general_marker_bounds.png',dpi=160)
