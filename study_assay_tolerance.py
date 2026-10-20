"""Bracket error budgets using sufficient and necessary precision certificates."""
import argparse,json
from pathlib import Path
import numpy as np
from general_marker_bounds import fast_bounds
from study_multiclaim import claim_masks

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refine',action='store_true');args=parser.parse_args()
    if args.refine:
        refine();return
    rows=[]
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));n=len(data['labels'])
        for radius,mask in zip((40,60,80),claim_masks(data)):
            edges=data['edges'][mask];lo=0;hi=int(np.ceil(n*.001));checks=[]
            def check(budget):
                result=fast_bounds(edges,data['labels'],budget,budget)
                outer=result['outer_bounds'][1]-result['outer_bounds'][0]
                attained=result['attained_bounds'][1]-result['attained_bounds'][0]
                row={'budget_each_marker':budget,'outer_width':outer,'attained_width':attained,
                     'guarantees_precision':outer<=.005,'proves_precision_impossible':attained>.005,
                     'details':result['details']}
                checks.append(row);return row
            assert check(lo)['guarantees_precision'];assert check(hi)['proves_precision_impossible']
            while hi-lo>1:
                mid=(lo+hi)//2;result=check(mid)
                if result['guarantees_precision']:lo=mid
                elif result['proves_precision_impossible']:hi=mid
                else:break
            row={'patient':patient,'radius':radius,'cells':n,'target_width':.005,
                 'largest_admissible_budget_bracket':[lo,hi-1],
                 'guaranteed_budget_fraction':lo/n,'first_certified_impossible_budget':hi,
                 'threshold_determined_exactly':hi-lo==1,'checks':checks}
            rows.append(row);print(patient,radius,'budget bracket',[lo,hi-1],flush=True)
    (ROOT/'results/assay_tolerance_study.json').write_text(json.dumps({'scope':'Exploratory conditional assay requirements; source annotations are the reference.','protocol_commit':'b9062fb','rows':rows},indent=2)+'\n')


def refine():
    path=ROOT/'results/assay_tolerance_study.json';report=json.loads(path.read_text());rows=[]
    for patient in ('CRC02','CRC11','CRC19'):
        data=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'));n=len(data['labels'])
        for radius,mask in zip((40,60,80),claim_masks(data)):
            edges=data['edges'][mask];checks=[]
            for budget in range(int(np.ceil(n*.001))+1):
                answer=fast_bounds(edges,data['labels'],budget,budget)
                checks.append({'budget_each_marker':budget,'outer_width':answer['outer_bounds'][1]-answer['outer_bounds'][0],
                               'attained_width':answer['attained_bounds'][1]-answer['attained_bounds'][0],'details':answer['details']})
            good=max(x['budget_each_marker'] for x in checks if x['outer_width']<=.005)
            bad=min(x['budget_each_marker'] for x in checks if x['attained_width']>.005)
            assert good<bad
            rows.append({'patient':patient,'radius':radius,'cells':n,'largest_admissible_budget_bracket':[good,bad-1],
                         'threshold_determined_exactly':bad-good==1,'checks':checks})
            print('refined',patient,radius,[good,bad-1],flush=True)
    report['refinement_protocol_commit']='7cdcad9';report['refined_rows']=rows
    path.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
