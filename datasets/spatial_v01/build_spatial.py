"""Build a measured, single-field development pack with reversible provenance."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, shift as image_shift
from skimage.registration import phase_cross_correlation
import tifffile

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'spatial'
META = OUT/'metadata'


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def parse_array(value, dtype=float):
    tokens = str(value).strip().strip('[]').split()
    return np.array([dtype(x) for x in tokens])


def orient(image, microscope):
    """TIFF row/column -> source-configured row/column, before registration."""
    if microscope['transpose']:
        image = image.T
    if microscope['flip_vertical']:
        image = image[::-1, :]
    if microscope['flip_horizontal']:
        image = image[:, ::-1]
    return image


def plane_map(organization):
    rows = []
    for row in organization.itertuples(index=False):
        frames, z = parse_array(row.frame, int), parse_array(row.zPos)
        if len(frames) != 28 or len(z) != 28 or len(set(frames)) != 28:
            raise ValueError('Unexpected z-to-frame mapping')
        for z_index, (frame, z_value) in enumerate(zip(frames, z)):
            rows.append(dict(round=int(row.imagingRound), frame=int(frame),
                             z_index=z_index, z_position_source=float(z_value),
                             z_position_unit='unspecified_in_source_configuration',
                             readout=row.readoutName, channel=row.channelName,
                             organization_channel_number=int(row.bitNumber),
                             bit_number=int(row.bitNumber) if row.bitNumber <= 16 else None,
                             color=int(row.color), role='readout' if row.bitNumber <= 16 else 'morphology'))
    frame_table = pd.DataFrame(rows)
    if frame_table.duplicated(['round', 'frame']).any():
        raise ValueError('Two readouts assigned to one source plane')
    if not frame_table.frame.between(0, 83).all():
        raise ValueError('Out-of-range TIFF frame')
    # Keep all remaining pages. Their channel/z mapping is not explicitly listed
    # as an imaging readout, so no z/color identity is manufactured for them.
    for imaging_round in range(1, 10):
        group = organization[organization.imagingRound == imaging_round]
        fiducial_frames = set(group.fiducialFrame.astype(int))
        if len(fiducial_frames) != 1 or set(group.fiducialImagingRound) != {imaging_round}:
            raise ValueError('Ambiguous fiducial metadata')
        fiducial_frame = next(iter(fiducial_frames))
        for frame in sorted(set(range(84)) - set(frame_table.loc[frame_table['round'] == imaging_round, 'frame'])):
            rows.append(dict(round=imaging_round, frame=frame, z_index=None,
                             z_position_source=None, z_position_unit=None,
                             readout=None, channel=None, bit_number=None,
                             color=int(group.fiducialColor.iloc[0]) if frame == fiducial_frame else None,
                             role='designated_fiducial' if frame == fiducial_frame else 'unassigned_source_plane'))
    result = pd.DataFrame(rows).sort_values(['round', 'frame']).reset_index(drop=True)
    if len(result) != 756 or result.duplicated(['round', 'frame']).any():
        raise ValueError('Incomplete plane provenance')
    return result


def highpass(image):
    image = image.astype(np.float32)
    return gaussian_filter(image, 1) - gaussian_filter(image, 8)


def correlation(first, second, valid):
    a, b = first[valid].astype(float), second[valid].astype(float)
    a -= a.mean(); b -= b.mean()
    denominator = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.dot(a, b)/denominator) if denominator else None


def estimate_registration(reference, moving):
    ref, mov = highpass(reference), highpass(moving)
    displacement, error, _ = phase_cross_correlation(ref, mov, upsample_factor=10, normalization=None)
    aligned = image_shift(mov, displacement, order=1, mode='constant', cval=np.nan, prefilter=False)
    valid = np.isfinite(aligned)
    h, w = ref.shape
    tile_shifts = []
    for sy, sx in [(slice(0, h//2), slice(0, w//2)), (slice(0, h//2), slice(w//2, w)),
                   (slice(h//2, h), slice(0, w//2)), (slice(h//2, h), slice(w//2, w))]:
        local, _, _ = phase_cross_correlation(ref[sy, sx], mov[sy, sx], upsample_factor=10, normalization=None)
        tile_shifts.append(local.tolist())
    disagreement = np.linalg.norm(np.asarray(tile_shifts)-displacement, axis=1)
    return dict(shift_y_pixels=float(displacement[0]), shift_x_pixels=float(displacement[1]),
                registration_error=float(error), correlation_before=correlation(ref, mov, np.ones(ref.shape, bool)),
                correlation_after=correlation(ref, aligned, valid),
                valid_fraction=float(valid.mean()), tile_shifts_yx=tile_shifts,
                max_tile_disagreement_pixels=float(disagreement.max()),
                interpretation='estimated translation; tile disagreement is a diagnostic, not a calibrated error bound')


def main():
    source_metadata = json.loads((META/'zenodo_record.json').read_text())
    for item in source_metadata['files']:
        if item['key'].endswith('.zip'):
            continue
        with (META/item['key']).open('rb') as handle:
            digest = hashlib.file_digest(handle, 'md5').hexdigest()
        if f'md5:{digest}' != item['checksum']:
            raise ValueError('Publisher metadata checksum mismatch')
    organization = pd.read_csv(META/'dataorganization.csv')
    microscope = json.loads((META/'microscope.json').read_text())
    codebook = pd.read_csv(META/'codebook.csv')
    readouts = organization[organization.bitNumber <= 16].sort_values('bitNumber').readoutName.tolist()
    if codebook.columns.tolist()[2:] != readouts or len(codebook) != 140:
        raise ValueError('Codebook and imaging readouts do not agree')
    if not np.isin(codebook.iloc[:, 2:].to_numpy(), [0, 1]).all():
        raise ValueError('Nonbinary codebook')
    codes = codebook.iloc[:, 2:].to_numpy()
    if codebook.name.nunique() != len(codebook) or len(np.unique(codes, axis=0)) != len(codebook):
        raise ValueError('Duplicate codebook labels or binary words')
    distances = (codes[:, None, :] != codes[None, :, :]).sum(axis=2)
    np.fill_diagonal(distances, 999)
    planes = plane_map(organization)
    stats, fiducials, middle, middle_metadata = [], [], [], []
    receipts = []
    for imaging_round in range(1, 10):
        path = OUT/'raw'/f'merFISH_merged_{imaging_round:02d}_000.tiff'
        receipt = json.loads(path.with_suffix('.receipt.json').read_text())
        if path.stat().st_size != receipt['bytes'] or sha256(path) != receipt['sha256']:
            raise ValueError('Raw TIFF does not match its verified acquisition receipt')
        receipts.append(receipt)
        with tifffile.TiffFile(path) as tiff:
            if len(tiff.pages) != 84 or tiff.series[0].shape != (84, 1608, 1608):
                raise ValueError('Unexpected raw image geometry')
            for frame, page in enumerate(tiff.pages):
                image = page.asarray()
                if image.dtype != np.uint16 or image.shape != (1608, 1608):
                    raise ValueError('Unexpected raw plane dtype/shape')
                stats.append(dict(round=imaging_round, frame=frame, minimum=int(image.min()),
                                  maximum=int(image.max()), mean=float(image.mean()),
                                  std=float(image.std()), zero_fraction=float((image == 0).mean()),
                                  uint16_saturation_fraction=float((image == 65535).mean())))
            fiducial_frame = int(organization.loc[organization.imagingRound == imaging_round, 'fiducialFrame'].iloc[0])
            fiducials.append(orient(tiff.pages[fiducial_frame].asarray(), microscope))
            selection = planes[(planes['round'] == imaging_round) & (planes.z_index == 14)].sort_values('organization_channel_number')
            for row in selection.itertuples(index=False):
                middle.append(orient(tiff.pages[row.frame].asarray(), microscope))
                middle_metadata.append(dict(array_index=len(middle)-1, round=imaging_round,
                                            frame=row.frame, z_index=14, readout=row.readout,
                                            channel=row.channel, role=row.role,
                                            raw_file=str(path.relative_to(ROOT)),
                                            raw_sha256=receipt['sha256'], field_id='000', split='development'))
        print('Read and checked all planes:', path.name, flush=True)
    planes['field_id'] = '000'
    planes['split'] = 'development'
    planes['raw_file'] = planes['round'].map(lambda r: f'spatial/raw/merFISH_merged_{r:02d}_000.tiff')
    planes['sha256'] = planes['round'].map({i+1: x['sha256'] for i, x in enumerate(receipts)})
    planes = planes.merge(pd.DataFrame(stats), on=['round', 'frame'], validate='one_to_one')
    planes.to_parquet(OUT/'planes.parquet', index=False)
    middle, fiducials = np.stack(middle), np.stack(fiducials)
    np.save(OUT/'middle_z14_oriented_unaligned.npy', middle)
    np.save(OUT/'fiducials_oriented_unaligned.npy', fiducials)
    pd.DataFrame(middle_metadata).to_csv(OUT/'middle_plane_map.csv', index=False)
    transforms = []
    for index, moving in enumerate(fiducials):
        estimate = estimate_registration(fiducials[0], moving)
        estimate.update(round=index+1, reference_round=1, source_frame=3,
                        transform_space='source-configured oriented row/column',
                        method='Gaussian bandpass sigma 1 minus sigma 8; phase cross-correlation; 0.1-pixel grid',
                        split='development', field_id='000')
        transforms.append(estimate)
    (OUT/'registration.json').write_text(json.dumps(transforms, indent=2, allow_nan=False)+'\n')
    # Nan-filled boundaries explicitly distinguish unavailable support from zero intensity.
    aligned = np.lib.format.open_memmap(OUT/'middle_z14_registered.npy', mode='w+', dtype=np.float32, shape=middle.shape)
    for index, metadata in enumerate(middle_metadata):
        transform = transforms[metadata['round']-1]
        displacement = (transform['shift_y_pixels'], transform['shift_x_pixels'])
        aligned[index] = image_shift(middle[index].astype(np.float32), displacement,
                                     order=1, mode='constant', cval=np.nan, prefilter=False)
    aligned.flush()
    del aligned
    audit = dict(raw_rounds=9, raw_planes=len(planes), height=1608, width=1608,
                 z_planes_per_round=28, raw_bytes=sum(x['bytes'] for x in receipts),
                 compact_z_index=14, compact_readout_and_morphology_planes=len(middle),
                 codewords=len(codebook), blank_codewords=int(codebook.name.str.startswith('blank_').sum()),
                 readout_bits=16, verified_file_crc32_count=9, local_sha256_verified_count=9,
                 minimum_codebook_hamming_distance=int(distances.min()),
                 codebook_bit_weight_values=np.unique(codes.sum(axis=1)).astype(int).tolist(),
                 archive_whole_file_md5_verified=False, all_raw_planes_retained=True,
                 spatial_resolution_um_per_pixel=microscope['microns_per_pixel'],
                 fields=1, split='development', independent_biological_units=None,
                 molecular_truth_available=False, molecular_coordinates_decoded=False,
                 registration_is_estimated=True, calibrated_localization_error_available=False,
                 model_evaluations_run=0)
    (OUT/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    (OUT/'source_manifest.json').write_text(json.dumps(receipts, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
