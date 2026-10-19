"""Render the complete confirmation results and a fixed-region alignment check."""
import json
from pathlib import Path

import fsspec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import zarr

ROOT=Path(__file__).resolve().parent


def main():
    report=json.loads((ROOT/'results/pair_study.json').read_text())
    figure,axes=plt.subplots(2,3,figsize=(13,7),sharex=True)
    methods={'measured_only':('Measured only','#30353c'),
             'marker_product_tuned':('Marker product','#c65c31'),
             'pair_constant_tuned':('Pair prevalence','#7d6899'),
             'pair_distance_tuned':('Distance only','#327d9a'),
             'pair_image_tuned':('Direct H&E pair','#27816e')}
    for row,patient in enumerate(('CRC04','CRC12')):
        for method,(label,color) in methods.items():
            results=[r for r in report['rows'] if r['patient']==patient and r['method']==method]
            x=[r['sampled_tiles'] for r in results]
            axes[row,0].plot(x,[100*r['rmse'] for r in results],'o-',label=label,color=color,markersize=4)
            axes[row,1].plot(x,[100*r['student_coverage'] for r in results],'o-',color=color,markersize=4)
            axes[row,2].plot(x,[100*r['betting_width'] for r in results],'o-',color=color,markersize=4)
        axes[row,0].set_ylabel(patient+'\nPair-fraction RMSE (pp)')
        axes[row,1].axhline(95,color='#333333',ls='--',lw=1)
        axes[row,1].set_ylim(0,103)
        for ax in axes[row]:
            ax.set_xticks([4,8,16,32]);ax.grid(alpha=0.18)
            ax.set_xlabel('Acquired core tiles out of 32')
    axes[0,0].set_title('Point estimation error')
    axes[0,1].set_title('Approximate 95% interval coverage (%)')
    axes[0,2].set_title('Rigorous interval width (pp)')
    figure.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=5,fontsize=9)
    figure.suptitle('Independent patient confirmation: direct pair prediction fails the fixed success criterion',fontsize=13)
    figure.tight_layout(rect=(0,.06,1,.94))
    figure.savefig(ROOT/'results/confirmation.png',dpi=180)
    plt.close(figure)

    # First predetermined CRC02 core tile, not selected for appearance or marker content.
    manifest=json.loads((ROOT/'results/data_manifest.json').read_text())['CRC02']
    d=np.load(ROOT/'data/CRC02-prepared.npz')
    tx,ty=d['core_tiles'][0]
    x0,y0=int(tx)*2048,int(ty)*2048
    with fsspec.open(manifest['image_url'],'rb',block_size=262144,cache_type='blockcache') as file:
        with tifffile.TiffFile(file) as image:
            store=image.series[0].levels[2].aszarr()
            array=zarr.open(store,mode='r')
            crop=np.array(array[y0//4:(y0+2048)//4,x0//4:(x0+2048)//4,:])
            store.close()
    coords=d['coords'];inside=((coords[:,0]>=x0)&(coords[:,0]<x0+2048)&(coords[:,1]>=y0)&(coords[:,1]<y0+2048))
    figure,axes=plt.subplots(1,2,figsize=(10,5))
    for ax in axes:
        ax.imshow(crop);ax.axis('off')
    for marker,label,color in [(0,'Pan-CK gate','#e9523f'),(1,'CD3e gate','#12b7c8')]:
        use=inside&d['labels'][:,marker]
        axes[1].scatter((coords[use,0]-x0)/4,(coords[use,1]-y0)/4,s=4,c=color,label=label,alpha=.6)
    axes[0].set_title('Registered H&E; fixed first core tile')
    axes[1].set_title('Source measured-marker gates at cell centroids')
    axes[1].legend(loc='lower right',fontsize=8)
    figure.suptitle('CRC02 alignment check; graph depends on source fluorescence segmentation',fontsize=11)
    figure.tight_layout()
    figure.savefig(ROOT/'results/region_check.png',dpi=160)
    plt.close(figure)


if __name__=='__main__':
    main()
