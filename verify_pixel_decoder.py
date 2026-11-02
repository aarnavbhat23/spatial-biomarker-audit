"""Independent image identities, mass readouts and comparator provenance checks."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def verify():
    study = json.loads((ROOT/'results/pixel_decoder.json').read_text())
    payload = (ROOT/'results'/study['details_file']).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == study['details_sha256']
    detail = json.loads(gzip.decompress(payload))
    expected = {'istdeco.py':'85bca2f1d5d0e9d6e0b5a8f6c2c5fd65524403a5ec8ef122aef10664a611e65a',
                'codebook.py':'cd0e6ce89e78d1d5b07e0cd47e013543c6dfd023aaa33fbd90a3cfcab2caff39',
                'psf.py':'b84b92256440cb761ce7ecaab774642ad51fce4e25da146573057df2810588da'}
    for name, checksum in expected.items():
        assert hashlib.sha256((ROOT/'third_party/istdeco'/name).read_bytes()).hexdigest() == checksum
    assert study['decoder']['commit'] == '264c9383a2d73d5669561bff179fd4c359d2e103'
    previous = json.loads((ROOT/'results/fiducial_measurement.json').read_text())
    witnesses = json.loads(gzip.decompress((ROOT/'results'/previous['geometry']['details_file']).read_bytes()))['witnesses']
    witnesses = [w for w in witnesses if w['first'] == 'Adamdec1 + Adgrg6']
    book = pd.read_csv(ROOT/'datasets/spatial_v01/codebook.csv')
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    names = book.name.tolist()
    code = dict(zip(names, book.iloc[:, 2:].to_numpy(dtype=int)))
    identity_count = 0
    for case in witnesses:
        signature = code['Adamdec1']+code['Adgrg6']
        first, second = [case[k].split(' + ') for k in ['first', 'second']]
        np.testing.assert_array_equal(signature, code[second[0]]+code[second[1]])
        assert max(signature) == 1
        coefficients = case['lp_coefficients']
        for d in [1., 3.]:
            for bit in np.flatnonzero(signature):
                row = layout.iloc[bit]
                displacement = d*(coefficients[int(row.imagingRound)-1]+coefficients[8+[561, 650].index(int(row.color))])
                center0 = (-d/2 if code[first[0]][bit] else d/2)-displacement/2
                center1 = (-d/2 if code[second[0]][bit] else d/2)+displacement/2
                np.testing.assert_allclose(center0, center1, atol=1e-12)
                # Equality of positions implies equality for ANY common translated
                # PSF; test another non-Gaussian profile as an independent check.
                grid = np.arange(-12, 13)
                np.testing.assert_allclose(1/(1+(grid-center0)**2), 1/(1+(grid-center1)**2), atol=1e-12)
                identity_count += 1
    masses = np.load(ROOT/'results/pixel_decoder_masses.npz')
    correct, probe_correct, registered = [], [], []
    for row in detail['synthetic_rows']:
        case = witnesses[row['case']]
        first, second = [case[k].split(' + ') for k in ['first', 'second']]
        suffix = '_alias' if row['mode'] == 'alias_paired' else ('_registered' if row['mode'] == 'known_registration' else '_probe')+f'_t{row["truth"]}'
        mass = masses[row['key']+suffix]
        score = sum(mass[names.index(g)] for g in first)-sum(mass[names.index(g)] for g in second)
        np.testing.assert_allclose(row['score'], score, atol=1e-8)
        assert row['decision'] == int(score < 0)
        if row['mode'] == 'alias_paired':
            assert sum(row['errors_across_two_worlds']) == 1
        else:
            assert row['correct'] == (row['decision'] == row['truth'])
            if row['mode'] == 'added_probe':
                correct.append(row['correct']); probe_correct.append(row['conventional_correct'])
            elif row['distance'] == 3 and row['exposure'] == 16:
                registered.append(row['correct'])
    np.testing.assert_allclose([study['registered_control_accuracy'], study['decoder_probe_accuracy'], study['conventional_probe_accuracy']],
                               [np.mean(registered), np.mean(correct), np.mean(probe_correct)], atol=1e-12)
    assert len(detail['measured']['rows']) == 20
    assert len(masses.files) == 960
    output = dict(checked_active_image_identities=identity_count, checked_decoder_readouts=len(detail['synthetic_rows']),
                  decoder_runs=study['decoder_calls'], upstream_source_files=3,
                  independent_gene_truth=False, external_expert_review=False)
    (ROOT/'results/pixel_decoder_verification.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    verify()
