"""Execute the frozen joint-claim study on three additional ORION patients."""
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
import fetch_data
from confirm_region_scale import prepare_panel
from study_raster_acquisition import prepare_frame_pairs
from study_multiclaim import claim_masks,multiclaim_trial

ROOT=Path(__file__).resolve().parent


def main():
    sources=json.loads((ROOT/'multiclaim_sources.json').read_text());fetch_data.SOURCES=sources
    protocol_commit=subprocess.check_output(['git','log','-1','--format=%H','--','MULTICLAIM_CONFIRMATION_PROTOCOL.md'],cwd=ROOT,text=True).strip()
    files=['study_multiclaim.py','study_raster_acquisition.py','study_cell_reuse.py',
           'study_residual_acquisition.py','confirm_region_scale.py','confirm_multiclaim.py']
    report={'protocol_commit':protocol_commit,'joint_alpha':.05,'width':.005,'field':1024,'batch':16,
            'radii':[40,60,80],'repetitions':20,'seed':20261110,'patients':{},
            'scope':'Patient-disjoint same-cohort confirmation; fixed source graph and exact source calls.',
            'implementation_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        records=dict(pool.map(fetch_data.fetch_table,sources['cell_tables'].items()))
    report['sources']=records
    (ROOT/'data/multiclaim_fresh_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    for patient in sources['cell_tables']:
        preparation=prepare_panel(patient);data=dict(np.load(ROOT/'data'/f'{patient}-acquisition.npz'))
        prepared=prepare_frame_pairs(data,1024);masks=claim_masks(data);entries=[]
        for policy in ('envelope','mixture','focused_reuse','greedy'):
            rows=[]
            for repeat in range(1 if policy=='greedy' else 20):
                row=multiclaim_trial(data,prepared,masks,1024,20261110+1000*int(patient[3:])+repeat,policy)
                row['repeat']=repeat;rows.append(row)
            entry={'policy':policy,'rows':rows,'attainment':float(np.mean([r['stopped'] for r in rows])),
                   'joint_coverage_all_checkpoints':float(np.mean([r['joint_covered_all_checkpoints'] for r in rows])),
                   'empty_sets':sum(r['empty_set_any_checkpoint'] for r in rows)}
            for key in ('fields','rounds','queries'):entry['mean_'+key]=float(np.mean([r[key] for r in rows]))
            entry['cost_scenarios']={str(h):entry['mean_fields']+h*entry['mean_rounds'] for h in (0,1,10,100)}
            entries.append(entry);print(patient,policy,'fields',entry['mean_fields'],'rounds',entry['mean_rounds'],
                                        'attain',entry['attainment'],'coverage',entry['joint_coverage_all_checkpoints'],flush=True)
        envelope,mixture,focused,greedy=entries
        passed=(envelope['attainment']>=.95 and all(c['attainment']>=.95 and envelope['mean_fields']<c['mean_fields'] for c in (mixture,focused)))
        report['patients'][patient]={'preparation':preparation,'claim_edge_counts':masks.sum(axis=1).tolist(),
                                     'methods':entries,'primary_passed':bool(passed)}
        report['primary_passed']=all(x['primary_passed'] for x in report['patients'].values())
        (ROOT/'results/multiclaim_confirmation.json').write_text(json.dumps(report,indent=2)+'\n')
        print(patient,'primary passed',passed,flush=True)


if __name__=='__main__':main()
