"""Conventional fixed-candidate barcode sensitivity on one development plane."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from skimage.feature import peak_local_max


def assign(vectors, codes, threshold=.5):
    norm = np.linalg.norm(vectors, axis=1)
    unit = vectors/np.maximum(norm[:, None], 1e-30)
    cosine = np.clip(unit @ codes.T, -1., 1.)
    ordering = np.argsort(-cosine, axis=1, kind='stable')[:, :2]
    distance = np.sqrt(np.maximum(0., 2-2*np.take_along_axis(cosine, ordering, axis=1)))
    return ordering[:, 0], (distance[:, 0] <= threshold)&(norm > 0), distance[:, 0], distance[:, 1]-distance[:, 0]


def sample_vectors(signals, positions, shifts, rounds):
    result = np.empty((len(positions), len(signals)), dtype=np.float64)
    for bit, signal in enumerate(signals):
        displaced = positions-np.asarray(shifts[rounds[bit]-1])
        result[:, bit] = map_coordinates(signal, displaced.T, order=1, mode='constant', cval=0., prefilter=False)
    return result


def run(root, output):
    metadata = pd.read_csv(root/'spatial/middle_plane_map.csv').iloc[:16]
    codebook = pd.read_csv(root/'spatial/metadata/codebook.csv')
    if metadata.readout.tolist() != codebook.columns.tolist()[2:]:
        raise ValueError('Bit order is inconsistent with the source codebook')
    images_path = root/'spatial/middle_z14_oriented_unaligned.npy'
    inventory = json.loads((root/'release_checksums.json').read_text())['files']
    expected = next(r['sha256'] for r in inventory if r['path'] == 'spatial/middle_z14_oriented_unaligned.npy')
    with images_path.open('rb') as handle:
        if hashlib.file_digest(handle, 'sha256').hexdigest() != expected:
            raise ValueError('Source compact image checksum mismatch')
    images = np.load(images_path, mmap_mode='r')[:16]
    transforms = json.loads((root/'spatial/registration.json').read_text())
    nominal = [[r['shift_y_pixels'], r['shift_x_pixels']] for r in transforms]
    rounds = metadata['round'].to_numpy()
    signals = np.empty(images.shape, dtype=np.float32); scales = []
    for index, image in enumerate(images):
        value = np.maximum(gaussian_filter(image.astype(np.float32), 1)-gaussian_filter(image.astype(np.float32), 3), 0.)
        scale = float(np.quantile(value, .999))
        if scale <= 0:
            raise ValueError('A source bit has no positive scale; no silent bit exclusion')
        signals[index] = value/scale; scales.append(scale)
    # Candidate discovery occurs once, under nominal transformations.
    grid = np.indices(images.shape[1:]).reshape(2, -1).T
    strength_squared = np.zeros(len(grid), dtype=np.float32)
    for bit, signal in enumerate(signals):
        shifted = grid-np.asarray(nominal[rounds[bit]-1])
        values = map_coordinates(signal, shifted.T, order=1, mode='constant', cval=0., prefilter=False)
        strength_squared += values*values
    strength = np.sqrt(strength_squared).reshape(images.shape[1:])
    threshold = float(np.quantile(strength, .995))
    positions = peak_local_max(strength, min_distance=3, threshold_abs=threshold, exclude_border=12)
    before_cap = len(positions)
    ordering = np.lexsort((positions[:, 1], positions[:, 0], -strength[positions[:, 0], positions[:, 1]]))
    positions = positions[ordering[:5000]]
    codes = codebook.iloc[:, 2:].to_numpy(dtype=float)
    codes /= np.linalg.norm(codes, axis=1)[:, None]
    variants = {'nominal': nominal}
    for tile in range(4):
        variants[f'fiducial_tile_{tile}'] = [r['tile_shifts_yx'][tile] for r in transforms]
    for sign in [-1, 1]:
        changed = np.array(nominal, dtype=float)
        for index in range(8):
            changed[index, 0] += sign*((-1)**index)
        variants[f'imposed_one_pixel_{sign:+d}'] = changed.tolist()
    tables = []; decoded = {}; summaries = []
    for name, shifts in variants.items():
        vectors = sample_vectors(signals, positions, shifts, rounds)
        identities, accepted, distances, margins = assign(vectors, codes)
        labels = codebook.name.to_numpy()[identities]
        blanks = np.char.startswith(labels.astype(str), 'blank_')
        decoded[name] = (identities, accepted)
        tables.append(pd.DataFrame(dict(candidate_id=np.arange(len(positions)), y=positions[:, 0], x=positions[:, 1],
                                        variant=name, nearest_codeword=labels, accepted=accepted,
                                        blank=blanks, distance=distances, distance_margin=margins)))
        summaries.append(dict(variant=name, candidates=len(positions), accepted=int(accepted.sum()),
                              accepted_gene_calls=int((accepted&~blanks).sum()), accepted_blank_calls=int((accepted&blanks).sum())))
    nominal_ids, nominal_accepted = decoded['nominal']
    stable = nominal_accepted.copy(); tile_stable = nominal_accepted.copy()
    for summary in summaries:
        name = summary['variant']; identities, accepted = decoded[name]
        identity_changes = nominal_accepted&accepted&(identities != nominal_ids)
        summary['acceptance_changes_from_nominal'] = int((accepted != nominal_accepted).sum())
        summary['identity_changes_among_jointly_accepted'] = int(identity_changes.sum())
        same = accepted&(identities == nominal_ids)
        stable &= same
        if name.startswith('fiducial_tile'):
            tile_stable &= same
    result = dict(source_doi='10.5281/zenodo.15678603', source_images_sha256=expected,
                  field_id='000', z_index=14, split='development', codewords=140, blank_codewords=6,
                  bit_scales=scales, candidate_strength_threshold=threshold, acceptance_distance=.5,
                  candidates_before_cap=before_cap, retained_candidates=len(positions), truncated=before_cap > 5000,
                  variants=summaries, nominal_accepted=int(nominal_accepted.sum()),
                  stable_across_fiducial_tiles=int(tile_stable.sum()), stable_across_all_variants=int(stable.sum()),
                  gene_call_truth_available=False, detection_recall_evaluated=False,
                  calibrated_error_bound=False, independent_biological_validation=False,
                  novelty='Transparent conventional nearest-codeword sensitivity screen, not a new decoder.')
    output.mkdir(exist_ok=True)
    pd.concat(tables).to_csv(output/'measured_barcode_candidates.csv', index=False)
    (output/'measured_barcodes.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'bit_scales'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset-root', type=Path, required=True)
    args = parser.parse_args()
    run(args.dataset_root, Path(__file__).resolve().parent/'results')
