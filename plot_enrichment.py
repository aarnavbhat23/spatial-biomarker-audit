import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=json.loads(Path('results/enrichment.json').read_text());d=json.loads(Path('results/enrichment_refinement.json').read_text())
primary=[x for x in r['rows'] if x['k']==12 and x['fraction']==.005]
refined=[x for x in d['structural_development'] if x['k']==12 and x['fraction']==.005]
fig,ax=plt.subplots(1,2,figsize=(10.5,4),layout='constrained');x=np.arange(6)
for key,offset,label,color in [('separated_bounds',-.18,'Separated control','#999999'),('outer_bounds',0,'Frozen joint bound','#21618c')]:
 widths=[p['bounds'][key][1]-p['bounds'][key][0] for p in primary]
 ax[0].bar(x+offset,widths,width=.18,label=label,color=color)
widths=[p['bounds']['outer_bounds'][1]-p['bounds']['outer_bounds'][0] for p in refined]
ax[0].bar(x+.18,widths,width=.18,label='Post-screen development',color='#b03a2e')
ax[0].set(xticks=x,xticklabels=[p['patient'].replace('MB-','') for p in primary],xlabel='Patient',ylabel='Covariance interval width',title='34% frozen / 53% development median reduction')
ax[0].legend(fontsize=8)
for i,p in enumerate(primary):
 lo,hi=p['bounds']['outer_bounds'];t=p['bounds']['source_contrast']
 ax[1].plot([lo,hi],[i,i],color='#21618c',lw=3);ax[1].scatter([t],[i],color='#21618c',s=20)
ax[1].axvline(0,color='black',ls='--')
ax[1].set(yticks=x,yticklabels=[p['patient'] for p in primary],xlabel='Covariance (fixed measured graph)',title='Frozen sign certificate: 1 of 6 patients')
fig.savefig('results/enrichment.png',dpi=180)
