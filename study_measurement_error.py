"""Development-only injected marker-call errors; no biological truth claims."""
import json
from pathlib import Path
import numpy as np
from measurement_error import pair_target,degree_error_margin,sharp_a_bounds,worst_a_flips,enlarge
from study_multiclaim import claim_masks,multiclaim_trial
from study_raster_acquisition import prepare_frame_pairs

ROOT=Path(__file__).resolve().parent


def covers(intervals,truth):
    return bool(np.all((intervals[...,0]-1e-12<=truth)&(truth<=intervals[...,1]+1e-12)))


def main():
    report={'scope':'Development-only injected errors relative to source calls, not biological validation.',
            'seed':20261114,'radii_native_pixels':[40,60,80],'target_width':.005,
            'joint_alpha':.05,'patients':{}}
    rng=np.random.default_rng(20261114)
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));labels=data['labels'];n=len(labels)
        masks=claim_masks(data);edges=[data['edges'][mask] for mask in masks]
        truth=np.array([pair_target(e,labels) for e in edges]);prepared=prepare_frame_pairs(data,1024)
        census=[]
        for fraction in (.0001,.001,.005,.01):
            budget=int(np.ceil(n*fraction));bounds=np.array([sharp_a_bounds(e,labels,budget) for e in edges])
            census.append({'fraction':fraction,'budget_a':budget,'bounds':bounds.tolist(),
                           'widths':np.diff(bounds,axis=1).ravel().tolist(),
                           'geometry_margins':[degree_error_margin(e,n,budget,0) for e in edges]})
        rows=[]
        for fraction in (.001,.005):
            budget=int(np.ceil(n*fraction));margins=np.array([degree_error_margin(e,n,budget,0) for e in edges])
            for corruption in ('random','adversarial'):
                if corruption=='random':
                    changed=labels.copy();selected=rng.choice(n,budget,replace=False)
                    changed[selected,0]=~changed[selected,0];direction='mixed'
                else:changed,selected,direction=worst_a_flips(edges[-1],labels,budget)
                corrupted=dict(data);corrupted['labels']=changed
                recorded=np.array([pair_target(e,changed) for e in edges]);trials=[]
                for repeat in range(5):
                    row=multiclaim_trial(corrupted,prepared,masks,1024,20261114+1000*int(patient[3:])+repeat,'envelope')
                    checkpoints=np.array([c['intervals'] for c in row['checkpoints']])
                    row['source_covered_all_checkpoints']=covers(checkpoints,truth)
                    row['enlarged_source_covered_all_checkpoints']=covers(enlarge(checkpoints,margins),truth)
                    row['enlarged_final_widths']=np.diff(enlarge(row['intervals'],margins),axis=1).ravel().tolist()
                    trials.append(row)
                entry={'fraction':fraction,'budget_a':budget,'actual_flips':len(selected),'corruption':corruption,
                       'direction':direction,'source_targets':truth.tolist(),'recorded_targets':recorded.tolist(),
                       'bias':(recorded-truth).tolist(),'geometry_margins':margins.tolist(),'rows':trials}
                rows.append(entry)
                print(patient,fraction,corruption,'bias',entry['bias'],
                      'recorded coverage',sum(r['joint_covered_all_checkpoints'] for r in trials),
                      'source coverage',sum(r['source_covered_all_checkpoints'] for r in trials),
                      'robust coverage',sum(r['enlarged_source_covered_all_checkpoints'] for r in trials),flush=True)
        report['patients'][patient]={'cells':n,'source_targets':truth.tolist(),'census':census,'experiments':rows}
        (ROOT/'results/measurement_error_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
