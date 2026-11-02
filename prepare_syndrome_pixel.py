"""Export fixed pixel inputs and an exact binary parity check for Julia."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from study_photon_ceiling import cases
from study_pixel_decoder import ROOT, inputs, render


def parity_check(codes):
    matrix = np.array(codes, dtype=np.uint8).copy()
    pivots, row = [], 0
    for column in range(matrix.shape[1]):
        choices = np.flatnonzero(matrix[row:, column])
        if len(choices) == 0:
            continue
        pivot = row+choices[0]
        matrix[[row, pivot]] = matrix[[pivot, row]]
        for other in range(len(matrix)):
            if other != row and matrix[other, column]:
                matrix[other] ^= matrix[row]
        pivots.append(column); row += 1
    vectors = []
    for column in range(matrix.shape[1]):
        if column not in pivots:
            vector = np.zeros(matrix.shape[1], dtype=np.uint8)
            vector[column] = 1
            for index, pivot in enumerate(pivots):
                vector[pivot] = matrix[index, column]
            vectors.append(vector)
    check = np.array(vectors)
    assert np.all(check @ codes.T % 2 == 0)
    return check


def prepare(dataset_root, output):
    names, codes, dictionary, mapping = inputs()
    check = parity_check(codes.astype(int))
    stacks = []
    for index, case in enumerate(cases()):
        first, second = [case[k].split(' + ') for k in ['first', 'second']]
        coefficients = case['lp_coefficients']
        shifts = np.array([coefficients[r]+coefficients[8+c] for r, c in mapping])
        alias = render(first, names, codes, mapping, 3., -3*shifts/2, 16.)
        controls = [render(pair, names, codes, mapping, 3., exposure=16.) for pair in [first, second]]
        for mode, image in zip(['alias', 'registered_first', 'registered_second'], [alias]+controls):
            frames = [image[r, c].tolist() for r, c in mapping]
            stacks.append(dict(key=f'case{index}_{mode}', case=index, mode=mode, first=first, alternative=second,
                               images=frames, background=32., min_weight=1., max_iters=20))
    image_path = dataset_root/'spatial/middle_z14_oriented_unaligned.npy'
    with image_path.open('rb') as stream:
        checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
    assert checksum == 'fbc4d3d0d034373ee896863f540cf11f7ad089f65128a8bb1a7bfa5490afb223'
    images = np.load(image_path, mmap_mode='r')[:16]
    scales = json.loads((ROOT/'results/measured_barcodes.json').read_text())['bit_scales']
    transforms = json.loads((ROOT/'datasets/spatial_v01/registration.json').read_text())
    grid = np.indices((48, 48)).reshape(2, -1).T+[600, 600]
    frames = []
    for bit, (r, c) in enumerate(mapping):
        signal = np.maximum(gaussian_filter(images[bit].astype(np.float32), 1)-gaussian_filter(images[bit].astype(np.float32), 3), 0)/scales[bit]
        shift = [transforms[r]['shift_y_pixels'], transforms[r]['shift_x_pixels']]
        frame = map_coordinates(signal, (grid-np.array(shift)).T, order=1, mode='constant', cval=0., prefilter=False).reshape(48, 48)*1000
        frames.append(frame.tolist())
    stacks.append(dict(key='measured_600_600', case=None, mode='measured', images=frames,
                       background=0., min_weight=20., max_iters=40))
    output.mkdir(parents=True, exist_ok=True)
    book = pd.DataFrame(codes.astype(int), columns=[f'block{i+1}' for i in range(16)])
    book.insert(0, 'gene', names)
    book.to_csv(output/'codebook.csv', index=False)
    np.savetxt(output/'parity_check.csv', check, delimiter=',', fmt='%d')
    (output/'inputs.json').write_text(json.dumps(dict(protocol_commit='a380e64', stacks=stacks), separators=(',', ':'), allow_nan=False)+'\n')
    print(json.dumps(dict(stacks=len(stacks), parity_shape=list(check.shape), source_sha256=checksum)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.dataset_root, args.output)
