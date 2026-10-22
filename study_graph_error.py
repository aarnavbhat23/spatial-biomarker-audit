"""Frozen manual-annotation calibration study; official image partitions."""
import io
import json
import math
from pathlib import Path
import time
import zipfile

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.stats import beta
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
import tifffile

from fetch_geometry import fetch
from graph_error import decode_mask, graph_errors, conformal_upper, pair_frequency_bounds


def segment(image,multiplier,distance):
    smooth=ndimage.gaussian_filter(image.astype(float),1)
    foreground=smooth>threshold_otsu(smooth)*multiplier
    cc,n=ndimage.label(foreground);sizes=np.bincount(cc.ravel())
    foreground&=sizes[cc]>=30
    foreground=ndimage.binary_fill_holes(foreground)
    edt=ndimage.distance_transform_edt(foreground)
    peaks=peak_local_max(edt,min_distance=distance,labels=foreground,exclude_border=False)
    markers=np.zeros(image.shape,int)
    if len(peaks):markers[tuple(peaks.T)]=np.arange(1,len(peaks)+1)
    return watershed(-edt,markers,mask=foreground).astype(int)


def run():
    root,manifest=fetch();start=time.perf_counter()
    with zipfile.ZipFile(root/'metadata.zip') as meta:
        splits={k:meta.read('metadata/'+k+'.txt').decode().splitlines() for k in ['training','validation','test']}
        plates=dict(line.split(',') for line in meta.read('metadata/filenames_and_plates.csv').decode().splitlines())
    assert [len(splits[k]) for k in splits]==[100,50,50]
    assert len(set(sum(splits.values(),[])))==200
    with zipfile.ZipFile(root/'images.zip') as images,zipfile.ZipFile(root/'masks.zip') as masks:
        def load(name):
            image=tifffile.imread(io.BytesIO(images.read('images/'+str(Path(name).with_suffix('.tif')))))
            mask=decode_mask(np.asarray(Image.open(io.BytesIO(masks.read('masks/'+name)))))
            return image,mask
        training=[load(name) for name in splits['training']]
        settings=[]
        for multiplier in [.8,1.,1.2]:
            for distance in [5,8,12]:
                f1=[graph_errors(segment(im,multiplier,distance),gt)['instance_f1'] for im,gt in training]
                settings.append({'multiplier':multiplier,'distance':distance,'mean_f1':float(np.mean(f1))})
                print('training',settings[-1],flush=True)
        chosen=max(settings,key=lambda x:x['mean_f1'])
        # Record selection before accessing validation/test outcomes.
        (root/'selection.json').write_text(json.dumps({'settings':settings,'chosen':chosen},indent=2)+'\n')
        result={'protocol_commit':'80f5d8f','source':manifest,'selection':chosen,'training_settings':settings,
                'radius_pixels':40,'alpha':.05,'splits':{},
                'selection_annotation_cost':{'fields':len(training),'reference_nuclei':sum(int(gt.max()) for _,gt in training),
                    'pixels':sum(im.size for im,_ in training)}}
        for split in ['validation','test']:
            rows=[]
            for name in splits[split]:
                im,gt=load(name);pred=segment(im,chosen['multiplier'],chosen['distance'])
                rows.append({'image':name,'plate':plates[name],**graph_errors(pred,gt)})
            result['splits'][split]=rows
            if split=='validation':
                q=conformal_upper([r['score'] for r in rows]);result['calibration_quantile']=q
                result['calibration_cost']={'fields':len(rows),'reference_nuclei':sum(r['reference_nuclei'] for r in rows),
                     'pixels':sum(load(r['image'])[0].size for r in rows),'timings':'Annotation and instrument timings unavailable'}
        rows=result['splits']['test'];covered=0;widths=[];naive=0
        for row in rows:
            m=row['predicted_edges'];budget=math.ceil(q*max(1,m))
            row['budget_each']=budget
            row['covered']=row['deleted_edges']<=budget and row['added_edges']<=budget
            row['naive_covered']=row['deleted_edges']==0 and row['added_edges']==0
            row['quarter_mark_diagnostic']=pair_frequency_bounds(m,.25*m,.25*m,budget,budget)
            covered+=row['covered'];naive+=row['naive_covered']
            widths.append(row['quarter_mark_diagnostic'][1]-row['quarter_mark_diagnostic'][0])
        n=len(rows)
        result['summary']={'test_images':n,'covered':covered,'naive_zero_error_covered':naive,
              'descriptive_binomial_95':[float(beta.ppf(.025,covered,n-covered+1)) if covered else 0.,
                    float(beta.ppf(.975,covered+1,n-covered)) if covered<n else 1.],
              'median_diagnostic_width':float(np.median(widths)),
              'mean_instance_f1':float(np.mean([r['instance_f1'] for r in rows])),
              'gate_passed':bool(covered>=47 and np.median(widths)<=.10),
              'compute_seconds':time.perf_counter()-start}
        result['plate_counts']={s:len(set(plates[n] for n in names)) for s,names in splits.items()}
        Path('results/graph_error.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in result.items() if k not in ['splits','training_settings']},indent=2),flush=True)

if __name__=='__main__':run()
