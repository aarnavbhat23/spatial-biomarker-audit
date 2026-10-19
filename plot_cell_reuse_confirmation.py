"""Show both the confirmed precision gain and the cost-model limitation."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    data=json.loads((ROOT/'results/cell_reuse_confirmation.json').read_text())
    patients=list(data['patients']);union=[];summed=[];cell=[];core=[]
    for patient in patients:
        fine,coarse,ablation=data['patients'][patient]['methods']
        union.append(fine['mean_area_native_pixels']/coarse['mean_area_native_pixels'])
        summed.append(fine['mean_summed_area_native_pixels']/coarse['mean_summed_area_native_pixels'])
        cell.append(fine['stop_fraction']);core.append(ablation['stop_fraction'])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.8));x=np.arange(4);width=.35
    axes[0].bar(x-width/2,union,width,label='Union area',color='#247c8d')
    axes[0].bar(x+width/2,summed,width,label='Summed rectangle area',color='#b35437')
    axes[0].axhline(1,color='#444444',ls='--',lw=1)
    axes[0].set(xticks=x,xticklabels=patients,ylabel='Fine / coarse acquisition area',ylim=(0,1.4),
                title='Savings depend on physical cost')
    axes[0].legend(frameon=False,fontsize=9)
    axes[1].bar(x-width/2,np.array(cell)*100,width,label='Cell-level reuse',color='#247c8d')
    axes[1].bar(x+width/2,np.array(core)*100,width,label='Core-total reuse',color='#b35437')
    axes[1].axhline(95,color='#444444',ls='--',lw=1)
    axes[1].set(xticks=x,xticklabels=patients,ylabel='Trials attaining target precision (%)',ylim=(0,119),
                title='Cell reuse rescues tight precision')
    axes[1].legend(frameon=False,fontsize=9,loc='upper center',ncol=2)
    fig.text(.025,.025,'Four new patients; 20 trials per method/patient; width target 0.5 percentage points. Instrument overhead is excluded here.',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,1))
    fig.savefig(ROOT/'results/cell_reuse_confirmation.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
