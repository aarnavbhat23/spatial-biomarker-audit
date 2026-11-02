"""Full-pixel counterexamples and an unchanged published ISTDECO comparator."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from study_photon_ceiling import cases

ROOT = Path(__file__).resolve().parent
VENDOR = ROOT/'third_party/istdeco'


def source_check():
    metadata = json.loads((VENDOR/'upstream.json').read_text())
    for name, checksum in metadata['files'].items():
        if hashlib.sha256((VENDOR/name).read_bytes()).hexdigest() != checksum:
            raise ValueError(f'Modified external comparator: {name}')
    return metadata


def inputs():
    book = pd.read_csv(ROOT/'datasets/spatial_v01/codebook.csv')
    layout = pd.read_csv(ROOT/'datasets/spatial_v01/dataorganization.csv').iloc[:16]
    if layout.readoutName.tolist() != list(book.columns[2:]):
        raise ValueError('Acquisition/codebook bit mismatch')
    mapping = [(int(r)-1, [561, 650].index(int(c))) for r, c in zip(layout.imagingRound, layout.color)]
    codes = book.iloc[:, 2:].to_numpy(dtype=float)
    dictionary = np.zeros((len(book), 8, 2), dtype=np.float32)
    for bit, (r, c) in enumerate(mapping):
        dictionary[:, r, c] = codes[:, bit]
    return book.name.tolist(), codes, dictionary, mapping


def psf(center, size=25, sigma=1.2):
    yy, xx = np.indices((size, size), dtype=float)
    shape = np.exp(-((yy-center[0])**2+(xx-center[1])**2)/(2*sigma*sigma))
    return shape/shape.sum()


def render(pair, names, codes, mapping, distance, offsets=None, exposure=1., extra_probe=False):
    rounds = 9 if extra_probe else 8
    image = np.full((rounds, 2, 25, 25), 2.*exposure)
    offsets = np.zeros(16) if offsets is None else np.asarray(offsets)
    for side, gene in zip([-.5, .5], pair):
        code = codes[names.index(gene)]
        for bit in np.flatnonzero(code):
            r, c = mapping[bit]
            image[r, c] += 200.*exposure*psf((12., 12.+side*distance+offsets[bit]))
    if extra_probe and 'Adamdec1' in pair:
        image[8, 0] += 200.*exposure*psf((12., 12.))
    return image


def decode(image, dictionary, background):
    source_check()
    import torch
    torch.set_num_threads(1)
    sys.path.insert(0, str(VENDOR)) if str(VENDOR) not in sys.path else None
    from istdeco import ISTDeco
    model = ISTDeco(np.asarray(image, dtype=np.float32), dictionary, (1.2, 1.2), b=background)
    x, q, loss = model.run(niter=100, acc=1., suppress_radius=None)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(loss.numpy())):
        raise ValueError('Nonfinite published decoder output')
    return x.sum(axis=(1, 2)).astype(float), float(loss[-1])


def claim(mass, names, first, alternative):
    score = sum(mass[names.index(g)] for g in first)-sum(mass[names.index(g)] for g in alternative)
    return int(score < 0), float(score)


def synthetic():
    rng = np.random.default_rng(2026110104)
    names, codes, dictionary, mapping = inputs()
    probe_dictionary = np.pad(dictionary, ((0, 0), (0, 1), (0, 0)))
    probe_dictionary[names.index('Adamdec1'), 8, 0] = 1
    rows, means, masses = [], {}, {}
    for case_index, case in enumerate(cases()):
        first, alternative = [case[s].split(' + ') for s in ['first', 'second']]
        coefficients = np.array(case['lp_coefficients'])
        shifts = np.array([coefficients[r]+coefficients[8+c] for r, c in mapping])
        for d, exposure in itertools.product([1., 3.], [1., 16.]):
            left = render(first, names, codes, mapping, d, -d*shifts/2, exposure)
            right = render(alternative, names, codes, mapping, d, d*shifts/2, exposure)
            key = f'c{case_index}_d{d}_E{exposure}'
            means[key] = dict(max_pixel_mean_difference=float(abs(left-right).max()),
                              paired_worlds=[first, alternative])
            for repeat in range(12):
                common = rng.poisson(left)
                mass, loss = decode(common, dictionary, 2*exposure)
                choice, score = claim(mass, names, first, alternative)
                sample_key = f'{key}_r{repeat}'
                masses[sample_key+'_alias'] = mass
                rows.append(dict(key=sample_key, case=case_index, distance=d, exposure=exposure,
                                 mode='alias_paired', decision=choice, score=score,
                                 errors_across_two_worlds=[choice != 0, choice != 1], loss=loss))
                for truth, pair in enumerate([first, alternative]):
                    registered = rng.poisson(render(pair, names, codes, mapping, d, exposure=exposure))
                    control_mass, control_loss = decode(registered, dictionary, 2*exposure)
                    decision, control_score = claim(control_mass, names, first, alternative)
                    masses[f'{sample_key}_registered_t{truth}'] = control_mass
                    rows.append(dict(key=sample_key, case=case_index, distance=d, exposure=exposure,
                                     mode='known_registration', truth=truth, decision=decision,
                                     score=control_score, correct=decision == truth, loss=control_loss))
                    # Reuse gene pixels: only the independent probe round is new.
                    probe_mean = render(pair, names, codes, mapping, d, exposure=exposure, extra_probe=True)[8:]
                    probe = rng.poisson(probe_mean)
                    augmented = np.concatenate([common, probe])
                    probe_mass, probe_loss = decode(augmented, probe_dictionary, 2*exposure)
                    probe_decision, probe_score = claim(probe_mass, names, first, alternative)
                    null = np.full((25, 25), 2.*exposure)
                    positive = null+200.*exposure*psf((12., 12.))
                    loglr = float((probe[0, 0]*np.log(positive/null)-(positive-null)).sum())
                    conventional_decision = int(loglr < 0)
                    masses[f'{sample_key}_probe_t{truth}'] = probe_mass
                    rows.append(dict(key=sample_key, case=case_index, distance=d, exposure=exposure,
                                     mode='added_probe', truth=truth, decision=probe_decision,
                                     score=probe_score, correct=probe_decision == truth, loss=probe_loss,
                                     conventional_decision=conventional_decision,
                                     conventional_correct=conventional_decision == truth,
                                     probe_log_likelihood_ratio=loglr))
    settings = []
    for case, d, exposure, mode in itertools.product(range(4), [1., 3.], [1., 16.],
                                                    ['alias_paired', 'known_registration', 'added_probe']):
        selected = [r for r in rows if (r['case'], r['distance'], r['exposure'], r['mode']) == (case, d, exposure, mode)]
        error = np.mean([e for r in selected for e in r['errors_across_two_worlds']]) if mode == 'alias_paired' else 1-np.mean([r['correct'] for r in selected])
        settings.append(dict(case=case, distance=d, exposure=exposure, mode=mode, error=float(error)))
    hard_controls = [r for r in rows if r['mode'] == 'known_registration' and r['distance'] == 3 and r['exposure'] == 16]
    probes = [r for r in rows if r['mode'] == 'added_probe']
    result = dict(seed=2026110104, settings=settings, pixel_mean_checks=means,
                  gates=dict(pixel_equal=all(r['max_pixel_mean_difference'] < 1e-9 for r in means.values()),
                             paired_error=all(r['error'] == .5 for r in settings if r['mode'] == 'alias_paired'),
                             registered_control=bool(np.mean([r['correct'] for r in hard_controls]) >= .8),
                             published_decoder_with_probe=bool(np.mean([r['correct'] for r in probes]) >= .95),
                             conventional_probe=bool(np.mean([r['conventional_correct'] for r in probes]) >= .95)),
                  registered_control_accuracy=float(np.mean([r['correct'] for r in hard_controls])),
                  decoder_probe_accuracy=float(np.mean([r['correct'] for r in probes])),
                  conventional_probe_accuracy=float(np.mean([r['conventional_correct'] for r in probes])),
                  rows=rows)
    return result, masses


def measured(dataset_root):
    image_path = dataset_root/'spatial/middle_z14_oriented_unaligned.npy'
    with image_path.open('rb') as stream:
        checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
    if checksum != 'fbc4d3d0d034373ee896863f540cf11f7ad089f65128a8bb1a7bfa5490afb223':
        raise ValueError('Source image checksum mismatch')
    images = np.load(image_path, mmap_mode='r')[:16]
    transforms = json.loads((ROOT/'datasets/spatial_v01/registration.json').read_text())
    scales = json.loads((ROOT/'results/measured_barcodes.json').read_text())['bit_scales']
    signals = np.array([np.maximum(gaussian_filter(im.astype(np.float32), 1)-gaussian_filter(im.astype(np.float32), 3), 0)/scale
                        for im, scale in zip(images, scales)])
    names, codes, dictionary, mapping = inputs()
    variants = {'nominal': [[r['shift_y_pixels'], r['shift_x_pixels']] for r in transforms]}
    variants.update({f'tile_{t}':[r['tile_shifts_yx'][t] for r in transforms] for t in range(4)})
    rows = []
    for y, x in [(600, 600), (600, 1100), (1100, 600), (1100, 1100)]:
        grid = np.indices((48, 48)).reshape(2, -1).T+[y, x]
        for variant, shifts in variants.items():
            stack = np.zeros((8, 2, 48, 48), dtype=np.float32)
            for bit, (r, c) in enumerate(mapping):
                stack[r, c] = map_coordinates(signals[bit], (grid-np.asarray(shifts[r])).T,
                                              order=1, mode='constant', cval=0., prefilter=False).reshape(48, 48)
            mass, loss = decode(stack, dictionary, 1e-5)
            rows.append(dict(crop_y=y, crop_x=x, variant=variant, gene_masses=mass.tolist(),
                             dominant_gene=names[int(mass.argmax())], loss=loss))
    return dict(source_sha256=checksum, field='000', z=14, crops=4, variants=5, rows=rows,
                known_gene_truth=False, independent_biological_validation=False)


def run(dataset_root):
    source = source_check()
    import torch
    result, masses = synthetic()
    real = measured(dataset_root)
    details = dict(synthetic_rows=result.pop('rows'), measured=real)
    path = ROOT/'results/pixel_decoder_details.json.gz'
    path.write_bytes(gzip.compress(json.dumps(details, separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    np.savez_compressed(ROOT/'results/pixel_decoder_masses.npz', **masses)
    result.update(protocol_commit='29a7eb4', decoder=source, torch_version=torch.__version__,
                  decoder_calls=len(masses)+len(real['rows']), gene_names=inputs()[0],
                  measured_summary=dict(crops=4, variants=5, known_gene_truth=False,
                                        dominant_gene_changes=sum(r['dominant_gene'] != next(n['dominant_gene'] for n in real['rows']
                                            if n['variant'] == 'nominal' and (n['crop_y'], n['crop_x']) == (r['crop_y'], r['crop_x']))
                                            for r in real['rows'] if r['variant'] != 'nominal')),
                  details_file=path.name, details_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  independent_biological_validation=False, new_decoder=False)
    (ROOT/'results/pixel_decoder.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['gates', 'decoder_calls', 'registered_control_accuracy', 'decoder_probe_accuracy',
                                           'conventional_probe_accuracy', 'measured_summary']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset-root', type=Path, required=True)
    run(parser.parse_args().dataset_root)
