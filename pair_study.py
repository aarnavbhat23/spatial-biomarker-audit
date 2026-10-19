"""Direct H&E pair prediction, with untouched patient confirmation."""
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import brier_score_loss,roc_auc_score

from estimators import betting_interval,difference_estimate
from run_pilot import totals

ROOT=Path(__file__).resolve().parent


def pair_features(data,indices,image=True):
    edges=data['edges'][indices]
    left,right=edges[:,0],edges[:,1]
    distance=np.linalg.norm(data['coords'][left]-data['coords'][right],axis=1)[:,None]/80
    if not image:
        return distance
    a,b=data['features'][left],data['features'][right]
    return np.concatenate([a,b,np.abs(a-b),distance],axis=1)


def true_pairs(data,indices):
    edges=data['edges'][indices]
    return data['labels'][edges[:,0],0]*data['labels'][edges[:,1],1]


def predict_edges(data,models,marker_models,prevalence):
    probabilities=np.column_stack([model.predict_proba(data['features'])[:,1] for model in marker_models])
    edges=data['edges']
    output={'marker_product':probabilities[edges[:,0],0]*probabilities[edges[:,1],1],
            'pair_constant':np.full(len(edges),prevalence,dtype=np.float32)}
    for name,model in models.items():
        probability=np.empty(len(edges),dtype=np.float32)
        for start in range(0,len(edges),20000):
            indices=slice(start,start+20000)
            probability[indices]=model.predict_proba(pair_features(data,indices,image=name=='pair_image'))[:,1]
        output[name]=probability
    return output


def main():
    manifest=json.loads((ROOT/'results/data_manifest.json').read_text())
    data={p:dict(np.load(ROOT/'data'/f'{p}-prepared.npz')) for p in manifest}
    train=['CRC01','CRC03'];validation='CRC19';confirmation=['CRC04','CRC12']
    lengths=[len(data[p]['edges']) for p in train]
    selected=np.random.default_rng(20261019).choice(sum(lengths),100000,replace=False)
    boundary=lengths[0]
    indices={train[0]:selected[selected<boundary],train[1]:selected[selected>=boundary]-boundary}
    features=np.concatenate([pair_features(data[p],indices[p]) for p in train])
    labels=np.concatenate([true_pairs(data[p],indices[p]) for p in train])
    prevalence=float(labels.mean())
    models={name:HistGradientBoostingClassifier(max_iter=100,max_leaf_nodes=15,l2_regularization=1.,random_state=20261019)
            for name in ('pair_image','pair_distance')}
    models['pair_image'].fit(features,labels)
    models['pair_distance'].fit(features[:,-1:],labels)
    marker_models=joblib.load(ROOT/'results/baseline.joblib')
    val_predictions=predict_edges(data[validation],models,marker_models,prevalence)
    val_data=data[validation]
    actual=true_pairs(val_data,slice(None))
    coefficients={};validation_scores={}
    val_y=np.bincount(val_data['edge_tiles'],weights=actual,minlength=len(val_data['rectangles']))
    for name,probability in val_predictions.items():
        val_f=np.bincount(val_data['edge_tiles'],weights=probability,minlength=len(val_y))
        coefficients[name]=float(np.clip(np.cov(val_y,val_f,ddof=1)[0,1]/np.var(val_f,ddof=1),0,1)) if np.var(val_f)>0 else 0.
        validation_scores[name]={'brier':float(brier_score_loss(actual,probability)),
                                 'auc':float(roc_auc_score(actual,probability))}
    choice={'training_patients':train,'validation_patient':validation,'confirmation_patients':confirmation,
            'training_edge_count':len(labels),'training_edge_prevalence':prevalence,
            'validation_coefficients':coefficients,'validation_scores':validation_scores,
            'primary':'Validation-tuned coefficients for every predictor; no choice on confirmation outcomes.'}
    joblib.dump(models,ROOT/'results/pair_models.joblib')
    choice['model_sha256']=hashlib.file_digest((ROOT/'results/pair_models.joblib').open('rb'),'sha256').hexdigest()
    (ROOT/'results/pair_frozen_choices.json').write_text(json.dumps(choice,indent=2)+'\n')
    np.savez_compressed(ROOT/'results/pair_training_indices.npz',**indices)
    print('Pair models and choices frozen before confirmation.',flush=True)
    output=[];patients={}
    rng=np.random.default_rng(20261019)
    for patient in confirmation:
        d=data[patient]
        predictions=predict_edges(d,models,marker_models,prevalence)
        actual=true_pairs(d,slice(None))
        owners=d['edge_tiles'];population=len(d['rectangles'])
        counts=np.bincount(owners,minlength=population)
        y=np.bincount(owners,weights=actual,minlength=population)
        target=float(y.sum()/counts.sum())
        totals_by_model={name:np.bincount(owners,weights=p,minlength=population) for name,p in predictions.items()}
        scores={name:{'brier':float(brier_score_loss(actual,p)), 'auc':float(roc_auc_score(actual,p)) if np.std(p)>0 else 0.5,
                      'virtual_pair_fraction':float(p.mean())} for name,p in predictions.items()}
        np.savez_compressed(ROOT/'results'/f'{patient}-pair-predictions.npz',**predictions)
        patients[patient]={'true_pair_fraction':target,'edge_metrics':scores,'tile_pair_counts':counts.tolist(),
                           'tile_measured_totals':y.tolist(),'tile_predicted_totals':{name:v.tolist() for name,v in totals_by_model.items()}}
        for size in (4,8,16,population):
            draws=[rng.choice(population,size,replace=False) for _ in range(1000)]
            methods=[('measured_only',np.zeros(population),0.)]
            for name,f in totals_by_model.items():
                methods.extend([(name+'_raw',f,1.),(name+'_tuned',f,coefficients[name])])
            for method,f,coefficient in methods:
                estimates=[];coverage=[];widths=[];rigorous=[];rigorous_widths=[];empty=[]
                for selected in draws:
                    result=difference_estimate(f,y[selected],selected,counts,coefficient)
                    estimates.append(result['estimate'])
                    lo,hi=result['student_interval']
                    coverage.append(lo-1e-12<=target<=hi+1e-12);widths.append(hi-lo)
                    band=betting_interval(f,y[selected],selected,counts,coefficient)
                    empty.append(band is None)
                    rigorous.append(band is not None and band[0]-1e-12<=target<=band[1]+1e-12)
                    rigorous_widths.append(band[1]-band[0] if band is not None else 0.)
                row={'patient':patient,'method':method,'sampled_tiles':size,
                     'rmse':float(np.sqrt(np.mean((np.array(estimates)-target)**2))),
                     'student_coverage':float(np.mean(coverage)),'student_width':float(np.mean(widths)),
                     'betting_coverage':float(np.mean(rigorous)),'betting_width':float(np.mean(rigorous_widths)),
                     'empty_betting_intervals':int(sum(empty))}
                output.append(row)
            print(patient,'completed confirmation budget',size,flush=True)
    primary={p:{r['method']:r['rmse'] for r in output if r['patient']==p and r['sampled_tiles']==8} for p in confirmation}
    direct=np.mean([primary[p]['pair_image_tuned'] for p in confirmation])
    passes=(all(primary[p]['pair_image_tuned']<primary[p]['marker_product_tuned'] for p in confirmation) and
            all(direct<np.mean([primary[p][control+'_tuned'] for p in confirmation]) for control in ('marker_product','pair_constant','pair_distance')))
    report={'frozen_choices':choice,'patients':patients,'rows':output,'primary_success':bool(passes),
            'scope':'Two-patient internal confirmation. Prediction and sampling methods are established; no architectural novelty claim.'}
    (ROOT/'results/pair_study.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Primary success criterion:',passes,flush=True)


if __name__=='__main__':
    main()
