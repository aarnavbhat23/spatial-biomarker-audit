"""Exploratory stronger controls after the initial pilot; no model retraining."""
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.metrics import brier_score_loss

from estimators import betting_interval,difference_estimate,union_area
from run_pilot import totals

ROOT=Path(__file__).resolve().parent


def main():
    manifest=json.loads((ROOT/'results/data_manifest.json').read_text())
    frozen=json.loads((ROOT/'results/frozen_model_choices.json').read_text())
    data={p:dict(np.load(ROOT/'data'/f'{p}-prepared.npz')) for p in manifest}
    train=np.concatenate([data[p]['labels'] for p in frozen['training_patients']])
    prevalence=train.mean(axis=0)
    models=joblib.load(ROOT/'results/baseline.joblib')
    validation=frozen['validation_patients'][0]
    d=data[validation]
    probability=np.tile(prevalence,(len(d['labels']),1))
    actual,guess,counts=totals(d,probability)
    density_coefficient=float(np.clip(np.cov(actual,guess,ddof=1)[0,1]/np.var(guess,ddof=1),0,1)) if np.var(guess)>0 else 0.
    protocol={'status':'Exploratory extension after inspecting initial pilot; not independent confirmation.',
              'added_control':'Training marker prevalence times known geometric pair count; no H&E image features.',
              'density_coefficient_from_validation':density_coefficient,
              'betting_stabilizers':[0.01,0.05,0.1,0.5,1.0],
              'confidence':'Fixed-budget two-sided test-martingale inversion with alpha=0.05.',
              'sampling_seed':20261018,'repetitions':1000}
    (ROOT/'results/confidence_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rng=np.random.default_rng(20261018)
    output=[]
    for patient in frozen['test_patients']:
        d=data[patient]
        predictions=np.load(ROOT/'results'/f'{patient}-predictions.npz')['probabilities']
        measured,predicted,counts=totals(d,predictions)
        density=counts*np.prod(prevalence)
        population=len(counts)
        target=measured.sum()/counts.sum()
        for size in (2,4,8,16,population):
            indices=[rng.choice(population,size,replace=False) for _ in range(1000)]
            for method,guess,coefficient in [('measured_only',predicted,0.),('difference',predicted,1.),
                                             ('validation_tuned',predicted,frozen['validation_prediction_coefficient']),
                                             ('density_only',density,density_coefficient)]:
                estimates=[];intervals=[]
                for selected in indices:
                    observed=measured[selected]
                    estimates.append(difference_estimate(guess,observed,selected,counts,coefficient)['estimate'])
                    intervals.append(betting_interval(guess,observed,selected,counts,coefficient))
                covered=[band is not None and band[0]-1e-12<=target<=band[1]+1e-12 for band in intervals]
                widths=[band[1]-band[0] if band is not None else 0. for band in intervals]
                output.append({'patient':patient,'sampled_tiles':size,'method':method,
                               'rmse':float(np.sqrt(np.mean((np.array(estimates)-target)**2))),
                               'betting_coverage':float(np.mean(covered)),'betting_mean_width':float(np.mean(widths)),
                               'empty_intervals':sum(band is None for band in intervals)})
                print(patient,size,method,'coverage',np.mean(covered),'width',np.mean(widths),flush=True)
    (ROOT/'results/confidence_audit.json').write_text(json.dumps({'protocol':protocol,'rows':output},indent=2)+'\n')


if __name__=='__main__':
    main()
