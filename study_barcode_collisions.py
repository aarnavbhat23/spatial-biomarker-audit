"""Exact additive barcode ambiguity and a calibrated independent-probe screen."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import gamma

CODEBOOK_SHA = '416db756170566d6e46ed21a8272d41c75114d8c95e384d4f0289df92bcd8c34'


def collision_classes(codes):
    codes = np.asarray(codes)
    if codes.ndim != 2 or not np.isin(codes, [0, 1]).all():
        raise ValueError('Expected binary codewords')
    if len(np.unique(codes, axis=0)) != len(codes):
        raise ValueError('Codewords must be unique')
    groups = defaultdict(list)
    for i, j in itertools.combinations_with_replacement(range(len(codes)), 2):
        groups[tuple((codes[i]+codes[j]).tolist())].append((i, j))
    return groups


def poisson_interval(count, exposures, alpha=.05):
    counts = np.asarray(count)
    lower = np.zeros_like(counts, dtype=float)
    positive = counts > 0
    lower[positive] = gamma.ppf(alpha/2, counts[positive])/exposures
    upper = gamma.ppf(1-alpha/2, counts+1)/exposures
    return lower, upper


def wilson(successes, total):
    z = 1.959963984540054
    p = successes/total
    middle = (p+z*z/(2*total))/(1+z*z/total)
    radius = z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return [float(middle-radius), float(middle+radius)]


def audit_codebook(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != CODEBOOK_SHA:
        raise ValueError('Frozen source codebook checksum mismatch')
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    names = [r['name'] for r in rows]
    columns = list(rows[0])[2:]
    codes = np.array([[int(r[c]) for c in columns] for r in rows])
    output = {}
    for subset in ['all_codewords', 'genes_only']:
        indices = [i for i, name in enumerate(names) if subset == 'all_codewords' or not name.startswith('blank_')]
        groups = collision_classes(codes[indices])
        ambiguous = []
        verified_probes = 0
        for signature, pairs in sorted(groups.items()):
            if len(pairs) <= 1:
                continue
            original_pairs = [(indices[i], indices[j]) for i, j in pairs]
            # Independent integer-coordinate summation, rather than array-key lookup.
            for i, j in original_pairs:
                assert all(int(codes[i, bit])+int(codes[j, bit]) == signature[bit]
                           for bit in range(codes.shape[1]))
            for world in original_pairs:
                for gene in set(world):
                    assert all(gene not in other for other in original_pairs if other != world)
                    verified_probes += 1
            ambiguous.append(dict(signature=list(signature),
                                  worlds=[[names[i], names[j]] for i, j in original_pairs]))
        ambiguous.sort(key=lambda g:g['worlds'])
        output[subset] = dict(codewords=len(indices), pair_worlds=len(indices)*(len(indices)+1)//2,
                              distinct_signatures=len(groups), ambiguous_classes=len(ambiguous),
                              ambiguous_worlds=sum(len(g['worlds']) for g in ambiguous),
                              class_size_distribution=dict(sorted(Counter(len(v) for v in groups.values()).items())),
                              exact_probe_separations_verified=verified_probes,
                              first_collision=ambiguous[0] if ambiguous else None,
                              collisions=ambiguous)
    return output


def run(codebook, output):
    audit = audit_codebook(codebook)
    rng = np.random.default_rng(2026103002)
    repetitions, rows, saved = 2000, [], {}
    for background, efficiency, exposures in itertools.product([2., 10.], [.5, 2., 8.], [1, 4, 16, 64]):
        minimum, maximum = .8*background, 1.2*background
        for target in [0, 1]:
            rate = maximum if target == 0 else minimum+efficiency
            counts = rng.poisson(exposures*rate, repetitions)
            lower, upper = poisson_interval(counts, exposures)
            key = f'b{background}_eta{efficiency}_n{exposures}_z{target}'
            saved[key] = counts.astype(np.int32)
            for name, threshold in [('calibrated_background', maximum), ('nominal_background', background)]:
                reject = lower > threshold
                rows.append(dict(background=background, efficiency=efficiency, exposures=exposures,
                                 target=target, method=name, true_rate=rate,
                                 rejection_rate=float(reject.mean()),
                                 rejection_wilson=wilson(int(reject.sum()), repetitions),
                                 mean_interval_coverage=float(((lower <= rate)&(rate <= upper)).mean()),
                                 uniform_identification_possible=efficiency > maximum-minimum))
    robust = [r for r in rows if r['method'] == 'calibrated_background']
    nulls = [r for r in robust if r['target'] == 0]
    strong = [r for r in robust if r['target'] == 1 and r['efficiency'] == 8 and r['exposures'] == 64]
    overlap_witnesses = []
    for background, efficiency in itertools.product([2., 10.], [.5, 2., 8.]):
        if efficiency <= .4*background:
            overlap_witnesses.append(dict(background=background, efficiency=efficiency,
                                          positive_background=.8*background,
                                          negative_background=.8*background+efficiency,
                                          common_rate=.8*background+efficiency))
    gates = dict(source_gene_collision_exists=bool(audit['genes_only']['ambiguous_classes']),
                 robust_null_control=all(r['rejection_wilson'][1] <= .075 for r in nulls),
                 strong_probe_power=all(r['rejection_wilson'][0] >= .80 for r in strong))
    result = dict(seed=2026103002, protocol='BARCODE_COLLISION_PROTOCOL.md', protocol_commit='0315fe7',
                  source_doi='10.5281/zenodo.15678603', source_sha256=CODEBOOK_SHA,
                  replicates_per_setting=repetitions, codebook=audit, photon_settings=rows,
                  exact_background_overlap_witnesses=overlap_witnesses, gates=gates,
                  conditional_feasibility_pass=all(gates.values()),
                  measured_molecular_collisions_established=False, biological_validation=False,
                  novel_decoder_established=False)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'barcode_probe_counts.npz', **saved)
    (output/'barcode_collisions.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with (output/'barcode_probe_summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    compact = {k:v for k,v in result.items() if k not in ['codebook', 'photon_settings']}
    compact['codebook'] = {k:{x:y for x,y in v.items() if x != 'collisions'} for k,v in audit.items()}
    print(json.dumps(compact, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--codebook', type=Path, default=Path(__file__).resolve().parent/'datasets/spatial_v01/codebook.csv')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent/'results')
    args = parser.parse_args()
    run(args.codebook, args.output)
