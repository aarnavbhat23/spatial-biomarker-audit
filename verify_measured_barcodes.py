"""Replay fixed candidate assignments via manual interpolation and direct distances."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist


def bilinear(image, y, x):
    iy, ix = np.floor(y).astype(int), np.floor(x).astype(int)
    fy, fx = y-iy, x-ix
    return ((1-fy)*(1-fx)*image[iy, ix]+fy*(1-fx)*image[iy+1, ix]
            +(1-fy)*fx*image[iy, ix+1]+fy*fx*image[iy+1, ix+1])


def verify(root, output):
    result = json.loads((output/'measured_barcodes.json').read_text())
    table = pd.read_csv(output/'measured_barcode_candidates.csv')
    image_path = root/'spatial/middle_z14_oriented_unaligned.npy'
    with image_path.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == result['source_images_sha256']
    images = np.load(image_path, mmap_mode='r')[:16]
    metadata = pd.read_csv(root/'spatial/middle_plane_map.csv').iloc[:16]
    codebook = pd.read_csv(root/'spatial/metadata/codebook.csv')
    codes = codebook.iloc[:, 2:].to_numpy(dtype=float)/2
    assert np.allclose(np.linalg.norm(codes, axis=1), 1.)
    transforms = json.loads((root/'spatial/registration.json').read_text())
    nominal = np.array([[r['shift_y_pixels'], r['shift_x_pixels']] for r in transforms])
    reference = table[table.variant == 'nominal'].sort_values('candidate_id')
    assert len(reference) == result['retained_candidates'] and reference.candidate_id.is_unique
    signals = []
    for index, image in enumerate(images):
        value = np.clip(gaussian_filter(image.astype(np.float32), 1)-gaussian_filter(image.astype(np.float32), 3), 0, None)
        assert abs(float(np.quantile(value, .999))-result['bit_scales'][index]) < 1e-6
        signals.append(value/result['bit_scales'][index])
    maximum_distance_error = 0.; replays = {}; checked = 0
    for variant, group in table.groupby('variant'):
        group = group.sort_values('candidate_id')
        np.testing.assert_array_equal(group[['candidate_id', 'x', 'y']], reference[['candidate_id', 'x', 'y']])
        shifts = nominal.copy()
        if variant.startswith('fiducial_tile_'):
            shifts = np.array([r['tile_shifts_yx'][int(variant[-1])] for r in transforms])
        elif variant.startswith('imposed_one_pixel_'):
            sign = int(variant.split('_')[-1])
            shifts[:8, 0] += sign*np.array([1, -1]*4)
        values = []
        for index, signal in enumerate(signals):
            dy, dx = shifts[int(metadata['round'].iloc[index])-1]
            values.append(bilinear(signal, group.y.to_numpy()-dy, group.x.to_numpy()-dx))
        vectors = np.stack(values, axis=1)
        norms = np.linalg.norm(vectors, axis=1)
        distances = cdist(vectors/np.maximum(norms[:, None], 1e-30), codes)
        choices = np.argmin(distances, axis=1)
        best = distances[np.arange(len(choices)), choices]
        accepted = (best <= .5)&(norms > 0)
        assert np.array_equal(codebook.name.to_numpy()[choices], group.nearest_codeword.to_numpy())
        assert np.array_equal(accepted, group.accepted.to_numpy())
        maximum_distance_error = max(maximum_distance_error, float(np.max(abs(best-group.distance))))
        replays[variant] = (choices, accepted)
        checked += len(group)
    assert maximum_distance_error < 1e-6
    ids, accepted = replays['nominal']; all_stable = accepted.copy(); tile_stable = accepted.copy()
    for name, (other_ids, other_accepted) in replays.items():
        all_stable &= other_accepted&(ids == other_ids)
        if name.startswith('fiducial_tile'):
            tile_stable &= other_accepted&(ids == other_ids)
    assert int(tile_stable.sum()) == result['stable_across_fiducial_tiles']
    assert int(all_stable.sum()) == result['stable_across_all_variants']
    report = dict(status='passed', candidate_variant_rows_checked=checked,
                  maximum_distance_replay_error=maximum_distance_error,
                  reference_computation='manual bilinear interpolation and direct Euclidean distances',
                  verified_detection_recall=False, verified_biological_accuracy=False)
    (output/'measured_barcode_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--dataset-root', type=Path, required=True)
    verify(parser.parse_args().dataset_root, Path(__file__).resolve().parent/'results')
