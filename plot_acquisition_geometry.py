"""Render confirmation, overhead sensitivity and labeled precision stress."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
report=json.loads((ROOT/'results/region_scale_confirmation.json').read_text())
stress=json.loads((ROOT/'results/region_scale_precision_stress.json').read_text())
patients=list(report['patients'])
fig,axes=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
for core,label in [(128,'128-pixel core'),(256,'256'),(512,'512'),(2048,'2048')]:
    ratios=[]
    for patient in patients:
        row=next(x for x in report['patients'][patient]['scales'] if x['core_native_pixels']==core)
        ratios.append(row['union_area_ratio_to_coarse'])
    axes[0].plot(patients,ratios,marker='o',label=label)
axes[0].set(ylabel='Area / coarse comparator area',title='Frozen 4 pp precision target');axes[0].legend(fontsize=8)
for patient in patients:
    rows=report['patients'][patient]['scales']
    small=next(x for x in rows if x['core_native_pixels']==128)
    coarse=next(x for x in rows if x['core_native_pixels']==2048)
    overhead=np.linspace(0,1_000_000,100)
    ratio=(small['mean_area_native_pixels']+overhead*small['mean_unique_regions'])/(coarse['mean_area_native_pixels']+overhead*coarse['mean_unique_regions'])
    axes[1].plot(overhead/1e6,ratio,label=patient)
axes[1].axhline(1,color='black',ls='--',lw=1)
axes[1].set(xlabel='Per-ROI overhead (million pixel-area equivalents)',ylabel='Small / coarse total cost',title='Area benefit can disappear with overhead')
axes[1].legend(fontsize=8)
x=np.arange(len(patients));small=[];coarse=[]
for patient in patients:
    rows=stress['patients'][patient]['0.005']
    small.append(rows[0]['stop_fraction']);coarse.append(rows[1]['stop_fraction'])
axes[2].bar(x-.18,small,width=.36,label='128-pixel core');axes[2].bar(x+.18,coarse,width=.36,label='2048-pixel core')
axes[2].set(xticks=x,xticklabels=patients,ylim=(0,1.15),ylabel='Fraction reaching precision within budget',title='Exploratory 0.5 pp stress test')
axes[2].legend(fontsize=8)
fig.savefig(ROOT/'results/acquisition_geometry.png',dpi=180)
