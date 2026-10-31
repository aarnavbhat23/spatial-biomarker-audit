# Spatial falsification and raw-image screen — frozen 2026-10-29

Two development experiments, with all outcomes retained. Neither is an independent
biological confirmation. No new gene-pair discovery or clinical/acquisition claim.

## A. Distinguish anatomy from cross-channel coupling

Fixed seed 2026102907; 200 independent simulated datasets per scenario; 128 A and
128 B positions; periodic 100 x 100 coordinate window; proximity radius 3.
Scenarios: independent uniform channels; independent channels sharing a fixed
three-component heterogeneous intensity; independent within-cell positions with
eight molecules/channel in each of 16 shared cells; separated channel domains;
and 96/128 B spots paired to A with Gaussian displacement SD .5 on top of shared
anatomy. Same-cell coexpression is not the spot-coupling alternative here.

Statistic: fraction of B positions within radius 3 of an A position (an established
cross-nearest-neighbour statistic, not a new statistic). Compare 199-draw Monte
Carlo tests at alpha .05: global label permutation; uniform-B resampling;
conditional within-25-unit-tile B resampling preserving tile counts; and oracle
B resampling from the known generating intensity/known cell membership. Oracle
information is unavailable in real tissue. Uniform and coarse-tile controls may
fail under heterogeneous anatomy; failures must not be hidden.

Error sensitivity: independently perturb observed points within disks of radii
.25, .75 and 1.5. For the oracle test, lower-bound observed proximity using
r-2*rho, and upper-bound each null score using r+rho (null B draws are exact,
observed A is uncertain). p_upper=(1+#upper_null>=lower_observed)/(200).
This dominates the same oracle Monte Carlo test on latent positions when every
recorded error is within rho. It assumes the true null intensity is known.
It is not a solution for unknown anatomy or a newly calibrated real-image bound.

Numerical implementation correction: exhaustive boundary checks exposed floating
point exclusion at a null upper-bound radius. Use a 1e-10 inward guard for the
observed lower comparison and outward guard for the null upper comparison in
this 100-unit window. This correction preserves the declared conservative bound;
it does not change the scenarios, scientific thresholds or selection gates.

Report all 200-replicate rejection counts, Wilson intervals, power and failures.
Screen: oracle and rho=.25 robust oracle have 95% upper rejection-rate bounds
<=.10 on each no-coupling case and 95% lower power bound >=.70 on coupling.
At 2*rho>=r, generic observed lower bounds are zero; retain this resolution limit.
Exact finite-configuration checks and direct-distance replay verify the bound
and implementation. No positive outcome is called algorithmic novelty.

## B. Measured barcode sensitivity

Use only the previously selected development field 000, fixed middle z index 14,
from Zenodo 15678603; preserve its source checksums. No additional field is opened.
This is a small transparent conventional nearest-codeword screen, NOT a full
MERlin/Serval replication, validated molecule detector or end-to-end decoder benchmark.

Use the 16 source-ordered readout planes. Background filter Gaussian sigma1 minus
sigma3, clip the filtered signal at zero; scale each bit by its 99.9th percentile.
Candidate strength is the norm of the nominally registered scaled 16-bit vector.
Find local maxima with minimum separation 3 pixels, above its 99.5th percentile,
exclude a 12-pixel border; retain strongest at most 5000, with row/column tie order.
Record any truncation. Normalize intensity vectors and match the 140 source
codewords (including all six blanks) by Euclidean distance; fixed acceptance .5.

Keep the SAME candidates and bit scales when replacing the nominal registration
by each of four fiducial-tile translation estimates. Also apply two explicitly
synthetic alternating-round +/-1-pixel y shifts. Report accepted gene/blank calls,
acceptance changes, identity changes, distance margins and all-variant stability.
This tests conditional assignment stability, not detection recall or gene-call
accuracy. Blank counts are NOT converted to an FDR without an appropriate model.
Tile agreement is not a calibrated coordinate-error bound. Preserve all candidate
records, all variants and failure states. No downstream molecular association is
declared significant from this single-field sensitivity exercise.

Next gate: a full pinned conventional decoder, reference spots or repeat captures,
and independently identified specimens before claiming calibrated biological
robustness. A stable candidate is not automatically a true molecule.
