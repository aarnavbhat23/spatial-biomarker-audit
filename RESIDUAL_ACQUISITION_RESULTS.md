# Reusing observations before optimizing acquisition policy

## What changed

The original geometry study used PPS sampling with replacement. The new baseline
samples only unresolved contributions and shares known fine-core totals across
nested scales. A second version stores cell-level visibility, reuses halo labels,
and resolves a zero pair whenever either relevant measured marker is negative.
Every method uses the same fixed directed-pair target and a valid betting sequence.
The original confirmation and its tighter-precision failures remain in the record.

## Conditional-mean argument

Before acquisition let K be the known positive count, U the unresolved edge count,
and M the fixed total. Choose a scale using only the current history. At that scale,
regions partition the unresolved owned edges into counts M_t and unknown totals Y_t.
Select t with probability M_t/U and form

    Z = K/M + (U/M) Y_t/M_t.

Conditional on the history, E[Z]=T, even when earlier outcomes influenced scale
selection. Its known support is [K/M,(K+U)/M]. After constructing this observation,
additional newly identifiable edges can be included in the next step's known set.
Repeated sampling from a depleted region is unnecessary.

For candidate mean mu, each positive bet uses
lambda_t(mu)=min(lambda_grid,0.5/max(mu-lower_t,epsilon)), and the negative bet
uses upper_t-mu. At the true mean, the factors are nonnegative with conditional
mean one. Mixing fixed component strategies preserves the test-martingale law.
Two one-sided thresholds 2/alpha and Ville's inequality give time-uniform coverage.
Candidate means outside current deterministic bounds can be discarded exactly.
The implementation inverts the monotone wealths by bisection. This is an
application of established predictable betting, not a new inference theorem.

For full cell reuse, known zeros require only one measured negative marker;
known positive edges require both relevant positive markers. With no additional
constraints on unmeasured binary marker bits, assigning all those bits zero or
one attains the lower and upper total bounds. Exhaustive completion tests verify
this statement on small graphs. It does not assert that every such assignment is
biologically plausible.

## Development results

Ten repetitions per condition were run on previously used CRC02, CRC11 and CRC19.
These are three patients, not hundreds of independent biological replicates.
The geometry-only adaptive rule maximizes expected newly learned pair count per
rectangle-area-plus-overhead cost. It failed to beat the best fixed scale
consistently and is not selected for confirmation.

At zero per-ROI overhead and width 0.04, its union area was approximately 75.9M,
87.1M and 39.2M native pixels, versus 17.2M, 18.2M and 13.9M for fixed 128-pixel
cores. At overhead 500,000 pixel-area equivalents, coarse cores were generally
competitive or cheaper. Maximizing acquired count is not the same as minimizing
cost to a specified precision.

Full cell reuse was compared separately, using a fresh simulation seed and the
same fixed scale. At width 0.005:

| Patient | Fine core-only area | Fine cell-reuse area | Fine draws before / after | Attainment before / after | Coarse cell-reuse area |
|---|---:|---:|---:|---:|---:|
| CRC02 | 69.58M | 60.80M | 2048 / 1024 | 100% / 100% | 133.56M |
| CRC11 | 71.82M | 95.64M | 2048 / 1944.8 | 90% / 100% | 128.43M |
| CRC19 | 34.15M | 30.22M | 1024 / 512 | 100% / 100% | 55.04M |

CRC11's area increase is retained. The lower-cost predecessor did not reliably
attain the same precision. At width 0.04, full reuse increased fine-core area on
all three patients. Sharing more observations changes the residual sampling
probabilities and does not guarantee lower physical cost. All inspected intervals
covered their target in these trials; the theoretical guarantee remains 95%.

The selected candidate is therefore **fine-core cell-level reuse for the tighter
0.5-percentage-point target when per-acquisition overhead is low**. Its next
confirmation must use untouched patients and the stronger coarse cell-reuse
comparator. It is not a universal winner or an optimal adaptive policy.

## Feasibility and prior work

The public ORION processing source at commit
`43387e3eebdf3d40103c4b54502f357561798fc9` reads channel 0 (Hoechst) for nuclear
segmentation and duplicates it into both inputs of its segmentation model.
This supports investigating a nuclear-counterstain-first acquisition workflow;
it does not demonstrate an H&E-only graph or prospectively measured assay savings.
Source: [processing script](https://github.com/labsyspharm/orion-crc/blob/43387e3eebdf3d40103c4b54502f357561798fc9/scripts-processing/processing/unmicst-scale-parallel.py).
The full graph still comes from the authors' processed fluorescence tables.

[Adaptive Monte Carlo via Bandit Allocation](https://proceedings.mlr.press/v32/neufeld14.html)
already studies sequential choice among unbiased estimators with different costs.
[Adaptive drawing probabilities](https://www150.statcan.gc.ca/n1/pub/12-001-x/2023001/article/00009-eng.htm)
and [without-replacement confidence sequences](https://arxiv.org/abs/2006.04347)
are also established. The scientific question is the joint effect of physical
visibility, shared cell labels, precision and actual measurement cost in this
spatial assay. Novelty priority remains unproven; hardware costs and an independent
cohort are still missing.

Reproduce with `python study_residual_acquisition.py` and
`python study_cell_reuse.py`. Full trial records are retained in the corresponding
JSON files under `results/`. Neither script trains a marker predictor.

The [frozen confirmation](CELL_REUSE_CONFIRMATION_RESULTS.md) is now complete.
All four patients passed the union-area criterion; summed-area and overhead
results substantially limit any physical-cost claim.
