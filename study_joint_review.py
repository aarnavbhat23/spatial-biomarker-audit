"""Online observation prototype on fixed development geometry and injected errors."""
import json
from pathlib import Path
import numpy as np
from acquisition_emulator import ObservationService,HttpClient,http_service,replay
from joint_review import intervals
from measurement_error import worst_a_flips
from study_multiclaim import claim_masks

ROOT=Path(__file__).resolve().parent


def fixture(patient,corruption):
    data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));labels=data['labels'];coords=data['coords'];edges=data['edges']
    fields,owners=np.unique(coords//1024,axis=0,return_inverse=True);selected=owners<min(8,len(fields))
    anchors=np.flatnonzero(selected);mapping=np.full(len(coords),-1);mapping[anchors]=np.arange(len(anchors))
    masks=claim_masks(data);weights=[];totals=[]
    for mask in masks:
        chosen=mask&selected[edges[:,0]];e=edges[chosen];totals.append(len(e))
        weights.append(np.bincount(mapping[e[:,0]],weights=labels[e[:,1],1],minlength=len(anchors)))
    weights=np.array(weights);totals=np.array(totals);reference=labels[anchors,0];budget=int(np.ceil(.005*len(anchors)))
    recorded=reference.copy();rng=np.random.default_rng(20261115+int(patient[3:]))
    if corruption=='random':ids=rng.choice(len(anchors),budget,replace=False)
    else:
        choices=[]
        for sign in (False,True):
            candidates=np.flatnonzero(reference==sign);ids=candidates[np.argsort(-weights[-1,candidates],kind='stable')[:budget]]
            choices.append((float(weights[-1,ids].sum()),ids))
        _,ids=max(choices,key=lambda x:x[0])
    recorded[ids]=~recorded[ids]
    return weights,totals,owners[anchors],recorded,reference,budget


def main():
    report={'scope':'Development-only exact-B prerequisite and ideal review; hypothetical incremental costs.','rows':[]}
    for patient in ('CRC02','CRC11','CRC19'):
        for corruption in ('random','adversarial'):
            args=fixture(patient,corruption);truth=args[0]@args[4]/args[1]
            for cost in (.05,.5):
                for policy in ('joint','census_first','random_review'):
                    service=ObservationService(*args,review_cost=cost)
                    row=replay(service,policy,seed=20261115)
                    covered=all(np.all((np.array(c['intervals'])[:,0]-1e-12<=truth)&(truth<=np.array(c['intervals'])[:,1]+1e-12)) for c in row['checkpoints'])
                    row.update({'patient':patient,'corruption':corruption,'review_cost':cost,'source_cells':len(args[2]),'budget_a':args[-1],
                                'truth':truth.tolist(),'joint_reference_covered_all_checkpoints':bool(covered),'transport':'in_process_observation_interface'})
                    report['rows'].append(row)
                    print(patient,corruption,cost,policy,'fields',row['fields'],'reviews',row['reviews'],'cost',round(row['simulated_time'],3),'covered',covered,flush=True)
                    (ROOT/'results/joint_review_study.json').write_text(json.dumps(report,indent=2)+'\n')
    args=fixture('CRC19','adversarial');service=ObservationService(*args,review_cost=.05)
    with http_service(service) as url:row=replay(HttpClient(url),'joint',seed=20261115)
    row.update({'patient':'CRC19','corruption':'adversarial','review_cost':.05,'transport':'localhost_HTTP'})
    report['http_replay']=row
    (ROOT/'results/joint_review_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
