"""Frozen general-error screen; retain solver gaps and precision failures."""
import hashlib,json,time
from pathlib import Path
import numpy as np
from general_marker_bounds import fast_bounds,optimize_bounds
from measurement_error import pair_target,degree_error_margin
from study_multiclaim import claim_masks

ROOT=Path(__file__).resolve().parent


def geometry(n,seed):
    rng=np.random.default_rng(11000+100*n+seed);xy=rng.uniform(size=(n,2))
    distance=np.linalg.norm(xy[:,None]-xy[None,:],axis=-1)
    edges=np.argwhere((distance>0)&(distance<=.35))
    labels=(rng.uniform(size=(n,2))<np.array([.3,.4])).astype(int)
    return edges,labels


def main():
    report={'scope':'Conditional sensitivity, assumed budgets, hypothetical reference observations.',
            'protocol_commit':'c10c31e','synthetic':[],'panels':[],'references':[],'input_sha256':{}}
    for n in (12,24,48):
        for seed in range(3):
            edges,labels=geometry(n,seed)
            for fraction in (.1,.25):
                budget=int(np.ceil(n*fraction));start=time.perf_counter();fast=fast_bounds(edges,labels,budget,budget)
                fast_seconds=time.perf_counter()-start;start=time.perf_counter();integer=optimize_bounds(edges,labels,budget,budget)
                integer_seconds=time.perf_counter()-start;start=time.perf_counter();lp=optimize_bounds(edges,labels,budget,budget,relax=True)
                row={'n':n,'seed':seed,'edges':len(edges),'budget_each_marker':budget,'fraction':fraction,
                     'fast':fast,'integer':integer,'lp':lp,'fast_seconds':fast_seconds,'integer_seconds':integer_seconds,
                     'lp_seconds':time.perf_counter()-start}
                report['synthetic'].append(row)
                print('synthetic',n,seed,budget,'integer optimal',[x['optimal'] for x in integer['objectives']],flush=True)
    (ROOT/'results/general_marker_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for patient in ('CRC02','CRC11','CRC19'):
        path=ROOT/'data'/f'{patient}-prepared.npz';report['input_sha256'][patient]=hashlib.sha256(path.read_bytes()).hexdigest()
        data=dict(np.load(path));labels=data['labels'];n=len(labels)
        for radius,mask in zip((40,60,80),claim_masks(data)):
            edges=data['edges'][mask]
            if len(np.unique(edges,axis=0))!=len(edges):raise ValueError('Simple graph required')
            target=pair_target(edges,labels)
            for fraction in (.001,.005):
                budget=int(np.ceil(n*fraction));start=time.perf_counter();fast=fast_bounds(edges,labels,budget,budget)
                seconds=time.perf_counter()-start;margin=degree_error_margin(edges,n,budget,budget)
                degree=[max(0.,target-margin),min(1.,target+margin)]
                row={'patient':patient,'n':n,'edges':len(edges),'radius':radius,'fraction':fraction,'budget_each_marker':budget,
                     'fast':fast,'degree_outer_bounds':degree,'fast_seconds':seconds}
                report['panels'].append(row)
                print('panel',patient,radius,fraction,'outer width',np.diff(fast['outer_bounds'])[0],
                      'gaps',fast['endpoint_gaps'],flush=True)
        (ROOT/'results/general_marker_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for seed in range(3):
        edges,truth=geometry(8,seed);raw=truth.copy();raw[:2]=1-raw[:2]
        for reviewed,errors in ((0,0),(4,1),(8,1),(8,0)):
            reference=np.full_like(raw,-1);reference[:reviewed]=truth[:reviewed]
            answer=optimize_bounds(edges,raw,2,2,reference=reference,reference_budgets=(errors,errors))
            report['references'].append({'seed':seed,'n':8,'edges':len(edges),'reviewed_cells_both_markers':reviewed,
                                         'reference_error_budget_each_marker':errors,'truth_target':pair_target(edges,truth),'answer':answer})
    (ROOT/'results/general_marker_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
