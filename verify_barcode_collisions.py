"""Independent sorted-row codebook audit and Poisson-tail count replay."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import poisson


def verify(root):
    result = json.loads((root/'results/barcode_collisions.json').read_text())
    path = root/'datasets/spatial_v01/codebook.csv'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == result['source_sha256']
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    names = [r['name'] for r in rows]
    codes = np.array([[int(v) for v in list(r.values())[2:]] for r in rows])
    verified_worlds = 0
    for subset, summary in result['codebook'].items():
        indices = [i for i, name in enumerate(names) if subset == 'all_codewords' or not name.startswith('blank_')]
        entries = sorted((tuple(codes[i]+codes[j]), (names[i], names[j]))
                         for ii, i in enumerate(indices) for j in indices[ii:])
        distinct, ambiguous_classes, ambiguous_worlds = 0, 0, 0
        group, previous = [], None
        independently_found = set()
        for signature, pair in entries+[(None, None)]:
            if signature != previous and group:
                distinct += 1
                if len(group) > 1:
                    ambiguous_classes += 1; ambiguous_worlds += len(group)
                    independently_found.add((previous, tuple(sorted(group))))
                group = []
            group.append(pair); previous = signature
        recorded = {(tuple(g['signature']), tuple(sorted(tuple(w) for w in g['worlds'])))
                    for g in summary['collisions']}
        assert independently_found == recorded
        assert (distinct, ambiguous_classes, ambiguous_worlds) == (
            summary['distinct_signatures'], summary['ambiguous_classes'], summary['ambiguous_worlds'])
        verified_worlds += len(entries)
    checked = 0
    with np.load(root/'results/barcode_probe_counts.npz', allow_pickle=False) as saved:
        assert 2*len(saved.files) == len(result['photon_settings'])
        for row in result['photon_settings']:
            b, eta, n, z = row['background'], row['efficiency'], row['exposures'], row['target']
            counts = saved[f'b{b}_eta{eta}_n{n}_z{z}']
            threshold = 1.2*b if row['method'] == 'calibrated_background' else b
            # Invert by direct Poisson upper-tail probability, not Gamma quantiles.
            reject = poisson.sf(counts-1, n*threshold) < .025
            assert abs(reject.mean()-row['rejection_rate']) < 1e-12
            checked += len(counts)
    verification = dict(status='passed', codebook_worlds_checked=verified_worlds,
                        method_replicate_decisions_replayed=checked,
                        unique_photon_replicates=checked//2,
                        reference='Sorted integer sums and direct Poisson survival probabilities',
                        independent_molecular_ground_truth=False)
    (root/'results/barcode_collision_verification.json').write_text(json.dumps(verification, indent=2)+'\n')
    print(json.dumps(verification, indent=2))


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
