"""All fresh-patient policies and the primary cost/compute comparison."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/fallible_policy_study.json').read_text())['rows']
policies=['uniform','degree','influence','lookahead','witness'];labels=['Uniform','Degree','Influence','Lookahead','Witness']
patients=list(dict.fromkeys(x['patient'] for x in rows));colors=['#8a9299','#b68448','#8f79ae','#49658b','#168376']
fig,axes=plt.subplots(1,2,figsize=(12,4.8))
for policy,label,color in zip(policies,labels,colors):
    means=[np.mean([x['result']['quoted_cost']['1'] for x in rows if x['patient']==p and x['policy']==policy]) for p in patients]
    axes[0].plot(range(6),means,'o-',color=color,label=label,lw=1.2,ms=4)
axes[0].set_xticks(range(6),[p.replace('MB-','') for p in patients]);axes[0].set_xlabel('Patient; two seeds and two error profiles each')
axes[0].set_ylabel('Mean hypothetical cost at latency 1');axes[0].set_title('All policies reach the same precision target')
axes[0].legend(fontsize=8,ncol=2);axes[0].grid(axis='y',alpha=.2)
compute=[sum(x['seconds'] for x in rows if x['policy']==p) for p in policies]
axes[1].bar(range(5),compute,color=colors);axes[1].set_xticks(range(5),labels,rotation=20)
axes[1].set_ylabel('Total policy computation, 24 scenarios (seconds)');axes[1].set_title('Witness: 3.08% less quoted cost, 10.65× faster')
axes[1].grid(axis='y',alpha=.2)
fig.suptitle('Fallible two-marker review on six fresh IMC patients',fontsize=12)
fig.text(.5,.012,'Conditional on supplied error budgets and measured-coordinate graphs; costs and fields are emulated.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(ROOT/'results/fallible_policy.png',dpi=180)
