import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

r=json.loads(Path('results/graph_error.json').read_text());test=r['splits']['test']
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
scores=np.array([x['score'] for x in test]);f1=np.array([x['instance_f1'] for x in test])
ax[0].scatter(f1,scores+.01,color='#21618c',s=24)
ax[0].axhline(r['calibration_quantile']+.01,color='#b03a2e',ls='--',label=f"Calibration q={r['calibration_quantile']:.2f}")
ax[0].set(yscale='log',xlabel='Instance F1 on held-out image',ylabel='Graph-error score + 0.01',title='Good segmentation scores do not certify the graph')
ax[0].legend(fontsize=8)
width=np.array([x['quarter_mark_diagnostic'][1]-x['quarter_mark_diagnostic'][0] for x in test])
ax[1].plot(np.arange(1,51),np.sort(width),color='#b03a2e',lw=2,label='Calibrated diagnostic width')
ax[1].axhline(.1,color='#777777',ls='--',label='Frozen useful-precision gate')
ax[1].set(xlabel='Held-out image (sorted)',ylabel='Interval width',ylim=(0,1.08),title='48/50 covered; precision gate fails')
ax[1].legend(fontsize=8,loc='center right')
fig.savefig('results/graph_error.png',dpi=180)
