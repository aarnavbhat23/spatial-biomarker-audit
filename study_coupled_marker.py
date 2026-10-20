"""Development-only exact census bounds, one A error and arbitrary B budget."""
import hashlib,json,time
from pathlib import Path
import numpy as np
from coupled_marker_bounds import sharp_one_a_bounds
from measurement_error import pair_target,degree_error_margin
from study_multiclaim import claim_masks

ROOT=Path(__file__).resolve().parent


def main():
    report={'scope':'Conditional fixed-graph sensitivity; source calls are not biological truth.',
            'budget_a':1,'b_fractions':[0,.001,.005],'rows':[],'input_sha256':{}}
    for patient in ('CRC02','CRC11','CRC19'):
        path=ROOT/'data'/f'{patient}-prepared.npz'
        report['input_sha256'][patient]=hashlib.sha256(path.read_bytes()).hexdigest()
        data=dict(np.load(path));n=len(data['labels'])
        for radius,mask in zip((40,60,80),claim_masks(data)):
            edges=data['edges'][mask]
            if len(np.unique(edges,axis=0))!=len(edges):raise ValueError('Simple graph required')
            reference=pair_target(edges,data['labels'])
            for fraction in report['b_fractions']:
                budget=int(np.ceil(n*fraction));start=time.perf_counter()
                bounds=sharp_one_a_bounds(edges,data['labels'],budget)
                seconds=time.perf_counter()-start
                margin=degree_error_margin(edges,n,1,budget)
                outer=[max(0.,reference-margin),min(1.,reference+margin)]
                row={'patient':patient,'cells':n,'edges':len(edges),'radius':radius,
                     'budget_a':1,'budget_b':budget,'b_fraction':fraction,'source_target':reference,
                     'bounds':bounds.tolist(),'width':float(bounds[1]-bounds[0]),
                     'degree_outer_bounds':outer,'degree_outer_width':outer[1]-outer[0],
                     'bound_seconds':seconds}
                report['rows'].append(row);print(row,flush=True)
        (ROOT/'results/coupled_marker_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
