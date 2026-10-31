# Spatial measurement dataset — development release 0.1

## What this is

A bounded raw-image curation of field `000` from
[MERFISH 4T1 Tissue Raw Unaligned Images](https://doi.org/10.5281/zenodo.15678603),
by Jenkin Tsui and colleagues, 2025. Full source creators are retained in
`metadata/zenodo_record.json`. The source describes mouse tissue bearing 4T1
tumors, with 54 fields, nine imaging rounds and 28 z slices per field.

This release selects the lowest field identifier, before inspecting its pixel
values or testing molecular associations. It retains that field's nine complete
raw TIFF stacks: 756 uint16 planes, each 1608 × 1608 pixels, totaling
3,909,650,274 uncompressed bytes. It includes 16 readout bits, two morphology
channels, the source codebook's 134 nonblank entries and six blank entries.

**One field is one development field.** Nine rounds, 28 z slices, image crops,
gene pairs and derived perturbations are not nine or hundreds of independent
specimens. Animal/specimen identities are unavailable in the selected source
metadata. Nothing in this release is designated an independent test cohort.

## What is new in this curation

- Exact source-plane provenance, including channel/readout/z mappings where
  explicitly specified; remaining pages are preserved without invented labels.
- Bounded range acquisition of the selected ZIP members, ZIP CRC32 validation
  and local SHA-256 receipts. Truncated downloads cannot be accepted as images.
- A compact middle-z view with explicit orientation and registered/unaligned
  alternatives, plus the designated fiducial plane from every round.
- Full raw-plane intensity QC and per-round registration diagnostics, with
  transformations kept separate from the source measurements.

The underlying images are previously published data. Curation, schema and audit
are our contributions here; they are not a new experimental collection.
Source and derived data retain CC BY 4.0 attribution. No novelty claim is made
for phase correlation, image alignment, decoding, segmentation-free statistics
or conventional robustness bounds.

## Provenance and geometry

`raw/merFISH_merged_0N_000.tiff` contains the original ZIP member bytes, unchanged.
Each accompanying receipt records source DOI, archive member, byte range, CRC32,
size and SHA-256. The 171 GB parent archive was not downloaded in full, so its
whole-archive publisher MD5 has **not** been verified. Member CRC32 verifies the
retrieved decompressed bytes against the archive index; local SHA-256 supports
subsequent integrity checks and is not a publisher-supplied per-member hash.

`planes.parquet` inventories all 756 source pages with plane QC. Source `frame`
indices are zero-based. The data-organization CSV explicitly maps 28 frames
per readout and uses frame 3 as the fiducial for each round. Other unmatched
pages are `unassigned_source_plane`; their identity is not guessed from array
stride. Source z coordinates 0–27 are retained with unspecified physical unit;
the configuration alone does not establish z spacing in micrometres.

The microscope metadata gives x/y resolution **0.174129 micrometres per pixel**,
horizontal and vertical flips, and no transpose. Raw TIFFs remain unchanged.
The compact arrays apply those source-configured flips before registration.
Coordinate arrays and transformations use `(row, column)` / `(y, x)` order.
Source stage positions are preserved; no unverified global stage mapping or
specimen identity is constructed from their row numbers.

## Compact arrays and registration

- `middle_z14_oriented_unaligned.npy`: `(18, 1608, 1608)`, uint16. Z index 14 is
  the fixed middle-plane convention chosen before viewing data. The first 16
  planes follow codebook bit order; the final two follow source Nucleus and Cell
  channels. Always use `middle_plane_map.csv` for identity.
- `fiducials_oriented_unaligned.npy`: `(9, 1608, 1608)`, uint16, source frame 3.
- `middle_z14_registered.npy`: same 18-plane geometry, float32, registered into
  round-1 coordinates. Missing interpolation support is NaN, never fabricated
  zero signal. Linear interpolation changes individual intensities; use the
  unaligned array or raw stacks when exact integer values matter.
- `registration.json`: translations applied to the oriented images; Gaussian
  bandpass followed by phase cross-correlation on a 0.1-pixel grid. Positive
  shifts move image content toward increasing row/column indices. Every round
  is retained, including any with poor diagnostics.

Fiducial alignment diagnostics include before/after correlation and agreement
among four image tiles. These are useful error indicators, **not** empirical
confidence intervals, molecular localization truth, or a calibrated bound on
chromatic/nonrigid error. A successful fiducial translation does not validate
the positions of every gene signal.

## Absent labels and prohibited interpretations

There are no curated molecular coordinates, independently validated gene calls,
true cell boundaries, repeated independent captures, detection probabilities,
real acquisition-cost measurements, clinical labels or independent-cohort results
in this release. The source blank codewords are useful controls, not automatically
known truth for all detection errors. No decoder or downstream association model
was evaluated during the build. A single z plane is a convenience view, not a
replacement for the retained 3D data.

## Valid next experiments

1. Run a pinned conventional decoder on all z planes, preserving spots, barcode
   evidence, blank assignments, localization scores and transformation provenance.
   Its calls are estimates; never label them ground truth.
2. Freeze the spatial estimand and decoder/registration sensitivity protocol
   using this development field. Report failures and abstentions as outcomes.
3. Separate measured uncertainty (fiducial residuals, backgrounds and blank
   responses) from explicitly simulated stress tests. Artificial molecule deletion
   does not establish the real missed-detection distribution.
4. Obtain independently annotated spots or a repeat-capture reference to calibrate
   uncertainty coverage. Validate across specimen/laboratory groups after those
   identities are verified. Additional fields from the same specimen are not
   independent biological validation.
5. For acquisition efficiency, collect actual scanner jobs and feedback latency
   in a lab pilot. Byte sizes and CPU timings cannot substitute for scanner costs.

The metadata folder also retains the earlier metadata-only ingestion audits.
Those describe that historical stage; `spatial/audit.json` is the completed
image-build inventory. See the root validation report for release status.
