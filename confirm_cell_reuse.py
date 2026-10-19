"""Frozen tight-precision comparison on four additional ORION patients."""
import concurrent.futures
import json
import subprocess
from pathlib import Path

import numpy as np

import fetch_data
from confirm_region_scale import prepare_panel
from study_cell_reuse import cell_trial, prepare_visibility
from study_residual_acquisition import nested_frames, residual_trial

ROOT=Path(__file__).resolve().parent


def main():
    sources=json.loads((ROOT/'cell_reuse_sources.json').read_text())
    fetch_data.SOURCES=sources
    protocol_commit=subprocess.check_output(['git','log','-1','--format=%H','--','CELL_REUSE_CONFIRMATION_PROTOCOL.md'],cwd=ROOT,text=True).strip()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        records=dict(pool.map(fetch_data.fetch_table,sources['cell_tables'].items()))
    (ROOT/'data/cell_reuse_fresh_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    report={'protocol_commit':protocol_commit,'target_width':.005,'alpha':.05,'repetitions':20,
            'seed':20261102,'sources':records,'patients':{},
            'scope':'Patient-disjoint confirmation within the same ORION cohort; known source graph and exact thresholded labels.'}
    for patient in sources['cell_tables']:
        preparation=prepare_panel(patient)
        data=dict(np.load(ROOT/'data'/f'{patient}-acquisition.npz'))
        counts,positive,frames=nested_frames(data)
        entries=[]
        for core,reuse in ((128,'cell_visibility'),(2048,'cell_visibility'),(128,'nested_core')):
            prepared=prepare_visibility(data,core) if reuse=='cell_visibility' else None
            rows=[]
            for repeat in range(20):
                seed=20261102+1000*int(patient[3:])+repeat
                row=(cell_trial(data,prepared,seed,.005) if reuse=='cell_visibility' else
                     residual_trial(counts,positive,frames,seed,str(core),.005,0))
                row['repeat']=repeat;rows.append(row)
            entry={'core':core,'reuse':reuse,'rows':rows,
                   'stop_fraction':float(np.mean([r['stopped'] for r in rows])),
                   'coverage_all_checkpoints':float(np.mean([r['covered_all_checkpoints'] for r in rows])),
                   'empty_sets':sum(r['empty_confidence_set'] for r in rows)}
            for metric in ('draws','area_native_pixels','summed_area_native_pixels'):
                entry['mean_'+metric]=float(np.mean([r[metric] for r in rows]))
            entries.append(entry)
            print(patient,core,reuse,'area',round(entry['mean_area_native_pixels']),
                  'draws',entry['mean_draws'],'attain',entry['stop_fraction'],flush=True)
        fine,coarse,ablation=entries
        extra=fine['mean_draws']-coarse['mean_draws']
        fine['overhead_break_even_native_pixel_area']=(coarse['mean_area_native_pixels']-fine['mean_area_native_pixels'])/extra if extra>0 else None
        passed=(fine['stop_fraction']>=.95 and coarse['stop_fraction']>=.95 and
                fine['mean_area_native_pixels']<coarse['mean_area_native_pixels'])
        report['patients'][patient]={'preparation':preparation,'methods':entries,'primary_passed':bool(passed)}
        report['primary_passed']=all(r['primary_passed'] for r in report['patients'].values())
        (ROOT/'results/cell_reuse_confirmation.json').write_text(json.dumps(report,indent=2)+'\n')
        print(patient,'primary passed',passed,flush=True)


if __name__=='__main__':main()
