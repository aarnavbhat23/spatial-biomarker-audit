# Acquisition units change the preferred sampling strategy

The tighter confirmation established an information benefit from reusing measured
cells, but its physical-cost interpretation depended on variable-size rectangles.
This exploratory study instead charges a reusable raster of fixed camera fields.
That is a plausible acquisition model, not a verified feature of a particular
instrument.

## Two acquisition rules and a stronger baseline

For logical cores of width 128 or 2048, acquire every aligned camera field touched
by the core and its halo. Reuse all labels in those fields and never charge a field
twice. Logical query rounds remain separate from physical camera exposures.

The stronger baseline partitions target edges by the unordered pair of camera
fields containing their endpoints. Residual PPS chooses such a group, whose
measurement needs only one or two fields. It retains every visible cell and every
logically determined zero edge. A diagonal cross-field pair needs two fields,
not all four fields in its bounding rectangle. Ownership, visibility and this
cost distinction are independently tested.

All three methods use the same residual conditional-mean estimator and the same
valid sequential uncertainty law. Five trials per condition were run on three
previously used development patients, at width .005. These are exploratory cost
sensitivity tests, not new-patient confirmation.

## Results

Each entry is mean physical fields / mean logical rounds. Every trial reached
the target precision and all inspected intervals covered truth in this small run.

| Patient | Field width | Fine logical core | Coarse logical core | Frame-pair sampling |
|---|---:|---:|---:|---:|
| CRC02 | 512 | 478.4 / 256.0 | 954.4 / 27.2 | 508.8 / 474.6 |
| CRC02 | 1024 | 179.2 / 128.0 | 414.4 / 26.8 | 133.2 / 128.0 |
| CRC02 | 2048 | 67.8 / 56.4 | 201.6 / 23.4 | 60.6 / 59.4 |
| CRC11 | 512 | 471.6 / 256.0 | 934.4 / 26.8 | 504.6 / 471.4 |
| CRC11 | 1024 | 176.8 / 128.0 | 397.6 / 26.2 | 134.2 / 128.0 |
| CRC11 | 2048 | 76.2 / 59.6 | 187.8 / 22.4 | 65.8 / 64.0 |
| CRC19 | 512 | 222.4 / 124.4 | 395.2 / 11.2 | 186.6 / 173.0 |
| CRC19 | 1024 | 68.8 / 47.0 | 168.0 / 11.0 | 57.6 / 55.6 |
| CRC19 | 2048 | 23.2 / 17.8 | 79.0 / 9.0 | 18.6 / 18.4 |

Frame-pair sampling used the fewest fields in seven of nine conditions. It was
strictly worse than fine rectangular queries in both fields and rounds on CRC02
and CRC11 at field width 512. It is not a universal winner or a proven optimal
acquisition policy. The sparse selected panels can also exaggerate the penalty
for coarse fields that extend beyond the target panel.

A consistent cost model is

    C = (field_area + field_overhead) * number_of_fields
        + round_overhead * number_of_rounds.

At field width 1024, frame-pair sampling reduced fields versus the fine logical
core on CRC02 and CRC11 without increasing rounds. On CRC19 it saved 11.2 fields
but needed 8.6 more rounds, so it is cheaper only when the per-round cost is below
about 1.30 times the per-field cost. Coarse logical queries use far fewer rounds;
high adaptive latency can therefore reverse their exposure disadvantage.

## Practical meaning and novelty limits

The [manufacturer's public specifications](https://www.rarecyte.com/orion-instruments.html)
list whole-slide scan speed and resolution. They do not establish the hypothetical
field widths used here, an arbitrary ROI API, field reuse, or feedback latency.
Those quantities cannot be inferred from aggregate whole-slide throughput alone.
No measured monetary or acquisition-time savings are claimed.

The next useful experiment needs an actual field size and a timing model, then a
frozen comparison of frame-pair, fine-rectangle and coarse-rectangle policies on
untouched panels or another cohort. Fixed staining costs, a nuclear prescan and
source-graph errors must also enter a deployment calculation. Acquiring fewer
fields can reduce imaging work under this model; it does not demonstrate that
fewer antibodies or staining rounds are needed.

This extends established survey/graph sampling and observation caching to a
specific spatial-readout acquisition problem. It is a serious stronger baseline,
not evidence that probability sampling itself is novel. The remaining methodological
question is how to choose acquisition groups and batch sizes to minimize validated
physical cost while preserving the inference guarantee.

Reproduce the rectangle comparison with `python study_raster_acquisition.py`, then
append the stronger baseline with `python study_raster_acquisition.py --frame-pair`.
Both sets of trial records are retained in `results/raster_acquisition_study.json`.
The optional physical-field path in `cell_trial` leaves the earlier rectangle
confirmation's default path unchanged.
