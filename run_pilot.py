"""Fit image-only baselines, freeze validation choices, audit held-out tile sampling."""
import hashlib
import json
import platform
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import binomtest
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from estimators import difference_estimate, union_area

ROOT = Path(__file__).resolve().parent
MARKERS = ['Pan-CK_pos','CD3e_pos']


def score(labels, probabilities):
    return {marker:{'brier':float(brier_score_loss(labels[:,i],probabilities[:,i])),
                    'auc':float(roc_auc_score(labels[:,i],probabilities[:,i])) if len(np.unique(labels[:,i]))==2 else None,
                    'prevalence':float(labels[:,i].mean())} for i,marker in enumerate(MARKERS)}


def totals(data, probabilities):
    edges,owners=data['edges'],data['edge_tiles']
    n=len(data['rectangles'])
    measured=data['labels'][edges[:,0],0]*data['labels'][edges[:,1],1]
    predicted=probabilities[edges[:,0],0]*probabilities[edges[:,1],1]
    return (np.bincount(owners,weights=measured,minlength=n),
            np.bincount(owners,weights=predicted,minlength=n),
            np.bincount(owners,minlength=n))


def main():
    manifest=json.loads((ROOT/'results/data_manifest.json').read_text())
    data={patient:dict(np.load(ROOT/'data'/f'{patient}-prepared.npz')) for patient in manifest}
    train=[patient for patient,record in manifest.items() if record['split']=='train']
    validation=[patient for patient,record in manifest.items() if record['split']=='validation']
    tests=[patient for patient,record in manifest.items() if record['split']=='test']
    x=np.concatenate([data[p]['features'] for p in train])
    y=np.concatenate([data[p]['labels'] for p in train])
    candidates={
        'logistic':[make_pipeline(StandardScaler(),LogisticRegression(C=1.0,max_iter=500)) for _ in MARKERS],
        'histogram_boosting':[HistGradientBoostingClassifier(max_iter=100,max_leaf_nodes=15,l2_regularization=1.,random_state=1729) for _ in MARKERS],
    }
    records={}
    validation_predictions={}
    for name,models in candidates.items():
        for i,model in enumerate(models):
            model.fit(x,y[:,i])
        validation_predictions[name]={p:np.column_stack([model.predict_proba(data[p]['features'])[:,1] for model in models]) for p in validation}
        records[name]={p:score(data[p]['labels'],validation_predictions[name][p]) for p in validation}
        print('Fitted',name,'validation',records[name],flush=True)
    def validation_loss(name):
        return np.mean([metric['brier'] for p in validation for metric in records[name][p].values()])
    chosen=min(candidates,key=validation_loss)
    models=candidates[chosen]
    val_totals=[totals(data[p],validation_predictions[chosen][p]) for p in validation]
    actual=np.concatenate([row[0] for row in val_totals])
    predicted=np.concatenate([row[1] for row in val_totals])
    coefficient=float(np.clip(np.cov(actual,predicted,ddof=1)[0,1]/np.var(predicted,ddof=1),0,1)) if np.var(predicted)>0 else 0.0
    frozen={'chosen_model':chosen,'validation_prediction_coefficient':coefficient,
            'validation_metrics':records,'training_patients':train,'validation_patients':validation,
            'n_training_cells':len(x),'test_patients':tests}
    # Persist model/choices before any test metrics are computed.
    joblib.dump(models,ROOT/'results/baseline.joblib')
    (ROOT/'results/frozen_model_choices.json').write_text(json.dumps(frozen,indent=2)+'\n')
    predictions={p:np.column_stack([model.predict_proba(data[p]['features'])[:,1] for model in models]) for p in tests}
    results=[]
    sampling={}
    rng=np.random.default_rng(20261018)
    for patient in tests:
        d=data[patient]
        probability=predictions[patient]
        np.savez_compressed(ROOT/'results'/f'{patient}-predictions.npz',cell_ids=d['cell_ids'],probabilities=probability)
        truth,guess,counts=totals(d,probability)
        population=len(counts)
        target=float(truth.sum()/counts.sum())
        patient_record={'marker_metrics':score(d['labels'],probability),'true_pair_fraction':target,
                        'virtual_pair_fraction':float(guess.sum()/counts.sum()),
                        'tile_measured_totals':truth.tolist(),'tile_predicted_totals':guess.tolist(),
                        'tile_pair_counts':counts.tolist(),'sampling':[]}
        complete_area=union_area(d['rectangles'])
        for size in (2,4,8,16,population):
            sampled=np.array([rng.choice(population,size,replace=False) for _ in range(1000)])
            area=np.array([union_area(d['rectangles'][idx]) for idx in sampled])
            for method,weight in [('measured_only',0.0),('difference',1.0),('validation_tuned',coefficient)]:
                runs=[difference_estimate(guess,truth[idx],idx,counts,weight) for idx in sampled]
                estimates=np.array([r['estimate'] for r in runs])
                row={'patient':patient,'method':method,'sampled_tiles':size,
                     'mean_acquired_panel_area_fraction':float(area.mean()/complete_area),
                     'bias':float(estimates.mean()-target),'rmse':float(np.sqrt(np.mean((estimates-target)**2)))}
                for interval in ('student_interval','bounded_interval'):
                    bands=np.array([r[interval] for r in runs])
                    covered=(bands[:,0]<=target+1e-12)&(bands[:,1]>=target-1e-12)
                    hits=int(covered.sum())
                    monte_carlo=binomtest(hits,len(covered)).proportion_ci()
                    row[interval]={'coverage':float(covered.mean()),'mean_width':float(np.diff(bands,axis=1).mean()),
                                   'monte_carlo_95_interval':[float(monte_carlo.low),float(monte_carlo.high)]}
                results.append(row)
                patient_record['sampling'].append(row)
            print(patient,'completed budget',size,flush=True)
        sampling[patient]=patient_record
    report={'scope':'Five-patient feasibility study; two held-out patients. Classical tile estimators, not a novel method.',
            'python':platform.python_version(),'frozen_choices':frozen,'patients':sampling,'rows':results}
    (ROOT/'results/pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    for patient in tests:
        print(patient,'truth / virtual',sampling[patient]['true_pair_fraction'],sampling[patient]['virtual_pair_fraction'],flush=True)


if __name__=='__main__':
    main()
