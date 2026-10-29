# Position-conditioned bounds: development protocol

Freeze before evaluation on 2026-10-28. This tests usefulness of a bound, not
biological truth, scanner savings or a new generic optimization algorithm.

Use a radius graph and mutually exclusive other/A/B identities. Each cell can
occupy one of three integer positions: its recorded coordinate, or one pixel
left/right. The uncertainty set includes every combination of those positions
and every categorical identity vector within the supplied per-marker Hamming
budgets. This is a complete discrete position contract, not coverage of arbitrary
continuous displacement, missing/merged cells or calibrated assay errors.

For each cell position, determine neighbors that must/may be connected across
all allowed positions of the other cell. Propagate local degree, contact and
marker-degree bounds into the degree-adjusted covariance. Intersect these with
existing fixed-graph label bounds enlarged by a valid total-variation bound on
graph edits. Compare independent-edge, local edit-budget and final joint-moment
bounds. Every comparison uses the same target and admissible errors. Do not
claim an advantage over the weak independent-edge control when the local edit
control already explains it.

Small cases: n=5,6, seeds 19,43,71, radius 5; coordinates are the first n sites of
a square grid with spacing 3, plus independent integer jitter in [-1,1]. Labels
are seeded categorical draws with probabilities (.2,.4,.4). Budgets (0,0),(1,1).
Enumerate all 3^n position assignments and all admissible categorical labels.
Compute extrema with integer numerators and denominators. Retain undefined
empty-graph cases explicitly. Any incorrect outer bound is a hard failure.

Scaling cases: n=128,512,2048, seeds 101,202, radius 5; spacing 3 and the same
jitter. Labels are drawn with the same probabilities. Budgets per marker are
0 and ceil(.005*n). Measure all bound variants on one CPU and save every case.
These random labels deliberately do not guarantee a nonzero biological signal.
Report widths, vacuity, sign certification and runtime, without selecting only
positive cases. No fresh patient outcome is consumed by this development screen.

Continue the claimed new-bound direction only if exhaustive checks pass and the
final bound improves median width by at least 10% over the stronger local-edit
control on the scaling cases. This is a development usefulness gate, not a
significance test or a biomedical tolerance. If it fails, retain the oracle and
the failure; do not promote the method as a paper's algorithmic contribution.
