# Spatial falsification and measured barcode sensitivity

The [protocol](MEASURED_SPATIAL_PROTOCOL.md) was committed as `d0f2bc8` before
outcomes. The benchmark tests the intended distinction between shared anatomy
and additional spot coupling. The raw-image study uses only the previously
selected single development field. Neither establishes independent biology.

## Falsification benchmark

Each row below contains 200 independent synthetic datasets, with 199 Monte Carlo
draws per method. Entries are rejection percentages at nominal alpha=.05.
The separated-domain case tests absence of positive coupling; negative association
is not the target. Same-cell coexpression is also not the coupling alternative.

| Scenario | Uniform B | Coarse conditional tiles | Global labels | Known-anatomy oracle | Robust oracle, rho=.25 |
|---|---:|---:|---:|---:|---:|
| Independent uniform | 7.0 | 6.5 | 6.0 | 6.5 | 0 |
| Shared heterogeneous anatomy | 100 | 100 | 3.0 | 2.5 | 0 |
| Independent within shared cells | 100 | 100 | 4.0 | .5 | 0 |
| Separated domains | 0 | 0 | 0 | 0 | 0 |
| Coupled spots | 100 | 100 | 100 | 100 | 100 |

Uniform and coarse-tile resampling give complete false-positive failure in the
two anatomically structured null cases. Global label permutation already works
well in these particular synthetic cases; it must not be hidden to make the
oracle construction appear uniquely capable. It targets its own exchangeability
null and is not universally validated by this small benchmark.

The **robust oracle at rho=.25 passes** the frozen null/power screen: zero
rejections in each null case (Wilson upper bound 1.885%) and 200/200 coupling
rejections (lower bound 98.115%). It assumes the true anatomy distribution is
known and the imposed position errors are bounded correctly. These inputs are
not supplied by the real field.

The **unadjusted oracle fails the conservative screen**: 13/200 uniform-null
rejections have a Wilson upper bound of 10.802%, above the fixed 10% gate.
Its interval includes 5%; this is not evidence disproving conditional Monte Carlo
validity. No extra repetitions were added to change the declared decision.

Robust-oracle power falls to **2/200** at rho=.75 and **0/200** at rho=1.5,
with proximity radius fixed at 3. These are simulation coordinate units, not
calibrated tissue micrometres. Small-error success does not establish useful
precision at realistic error levels.

The [conditional proof](CONDITIONAL_PROXIMITY_THEORY.md) explains validity and
the missing assumptions. All 3,000 simulated error configurations satisfy
pointwise domination of the latent oracle p value. Exhaustive small-case checks
cover 512 admissible perturbation configurations and caught a floating-point
boundary error; inward/outward comparison guards fix it. The complete fixed
study was rerun after that correction. No scientific gate or scenario changed.

## Raw-image barcode screen

The source is field 000, middle z index 14, from
[Tsui et al., Zenodo 15678603](https://doi.org/10.5281/zenodo.15678603), CC BY 4.0.
The analysis uses a transparent conventional nearest-codeword implementation,
not a full MERlin/Serval replication. All 1,295 selected candidates and all seven
variants are retained; the 5,000-candidate cap was not reached.

| Alignment variant | Accepted candidates | Acceptance changes from nominal |
|---|---:|---:|
| Nominal | 112 | 0 |
| Fiducial tile 0 | 112 | 0 |
| Fiducial tile 1 | 109 | 5 |
| Fiducial tile 2 | 111 | 3 |
| Fiducial tile 3 | 108 | 6 |
| Imposed alternating -1-pixel shifts | 71 | 49 |
| Imposed alternating +1-pixel shifts | 76 | 58 |

**105/112 (93.75%)** nominally accepted candidates retain acceptance and identity
under all four fiducial estimates. **41/112 (36.61%)** do so under all seven
variants. A pixel is .174129 micrometres in this source; the one-pixel changes
are deliberate stress tests, not an estimated localization-error distribution.

There are **no identity changes among jointly accepted candidates** in these
variants. The observed instability concerns threshold acceptance, not demonstrated
wrong-gene switching. All accepted assignments are to nonblank codewords; zero
blank assignments does not establish zero FDR or correct gene identities.
Candidate discovery is fixed to nominal alignment, so this does not test end-to-end
detection recall. Single-plane analysis does not replace a full 3D decoder.

## Verification and artifacts

- [All 7,000 simulation method/trial results](results/spatial_falsification_trials.csv) and [aggregate results/intervals](results/spatial_falsification.json).
- [All 9,065 candidate/variant records](results/measured_barcode_candidates.csv) and [barcode summary](results/measured_barcodes.json).
- [Independent interpolation/distance replay](results/measured_barcode_verification.json): every candidate/variant label and acceptance decision agrees; maximum distance discrepancy 2.47e-8.
- Six new unit tests cover the geometric error envelope, periodic geometry,
  resolution failure, registration direction, ambiguity rejection and barcode
  normalization. The full local suite has 90 passing tests, including the
  existing localhost acquisition emulator.

Reproduce the simulation with `python study_spatial_falsification.py`.
The [bounded source-acquisition package](datasets/spatial_v01/README.md) downloads
only the required field and rebuilds the image inputs. Then run
`study_measured_barcodes.py` and `verify_measured_barcodes.py` with the documented
`--dataset-root`. Large source images remain untracked.

## Paper decision

Defer submission. The benchmark and raw data sharpen the requirements but do not
establish a new spatial statistic or decoder. The most defensible direction is
a calibrated decision/abstention method that handles uncertainty in anatomy AND
measurement, with a nonvacuous guarantee. That requires a full conventional decoder
as a control, independent spot/reference measurements or repeat captures, and
identified specimens. The current oracle proof does not cover unknown intensity,
missed detections, false spots, gene-call errors or real scanner costs.

Source acquisition and sensitivity are reproducible. Useful real-data calibration,
independent validation and originality beyond established spatial-process methods
remain unmet submission gates.
