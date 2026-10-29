# Position-conditioned robustness: valid but insufficiently useful

The [protocol](POSITION_UNCERTAINTY_PROTOCOL.md) was committed as `b96dac4`
before evaluation. This round implements a complete finite set of allowed
positions and categorical identities, then tests safe bounds for the spatial
contrast. It does not add arbitrary latent states or assume annotation truth.

## Results and decision

Every outer interval contains the exact attainable contrast range in all 12
small cases. A separate verifier checks **3,669 distinct graph/label worlds**,
component count bounds and rational endpoint witnesses without importing the
production bounds or oracle.

On the 12 declared scaling cases, the median width reduction against the
stronger local-edit control is **9.3775%**, below the frozen **10%** gate. The
gate fails. Widths are .4323–.5000 out of the entire possible .5 range. None of
the 12 intervals certifies a positive or negative contrast.

| Cells | Width, no label errors (two seeds) | Width, positive label budget (two seeds) |
|---:|---|---|
| 128 | .43582, .43553 | .46910, .46934 |
| 512 | .43575, .43713 | .49832, .49795 |
| 2,048 | .43442, .43225 | .50000, .50000 |

The positive budget is ceil(.005*n) per marker; the two marker error sets need
not coincide. Position choices are the recorded integer location or one pixel
left/right, with radius five. All combinations are allowed. Labels are random
with probabilities (.2,.4,.4) for other/A/B, so zero sign certification alone is
not a test of power against a specified biological signal. The broad intervals
and the failed comparative gate are the reasons to stop this bound variant.

The modest no-label-error gains must not be presented as the overall result.
Do not tune the gate after these outcomes. This is development evidence and not
an independent cohort evaluation.

## Why the bounds are outward

For an undirected graph, let m be the number of directed edges, C the number of
A-to-B edges, and DA, DB the sums of degrees over A and B cells. The target is
F=C/m-DA*DB/m^2 for m>0. Categories A and B are mutually exclusive. Empty graphs
have no defined edge-distribution contrast and are flagged, never assigned zero.

For each allowed central-cell position, a neighbor is mandatory if it is within
radius for all of its positions, and possible if it is within radius for any.
Taking minima/maxima over the center's positions bounds each degree. Summing
these bounds gives [mL,mU]; parity is safely rounded inward to even integers
because directed edges occur in reverse pairs. Analogous per-position counts
bound deletions/additions relative to the recorded graph. Summing per-vertex
maxima relaxes the shared-position constraint and therefore remains outward.

For fixed labels, two uniform distributions on directed edge sets of sizes m0
and m0-D+A have total variation

    max(D/m0, A/(m0-D+A)).

Replacing D and A by their upper bounds is conservative; use one when the
denominator may vanish. The covariance of two binary variables is 1-Lipschitz
in total variation. To see this, its four probability-coordinate derivatives,
after subtracting a common constant, are 0, -b, -a, 1-a-b. Their range is at most
one for a,b in [0,1]. Integrating the derivative along the segment between two
probability laws proves the bound. Apply this to the same admissible labeling
on the recorded and displaced graphs, then union over labelings. Enlarging the
fixed-graph label interval by the TV bound is therefore safe.

The additional moment bound separately encloses C, DA and DB. At each vertex,
at most B_B B-label flips alter its neighbor count. Spending that whole budget
at every vertex is a conservative relaxation. At most B_A A-label flips alter
the sum over central vertices; removing/adding the largest admissible terms
bounds that sum. The same construction bounds degree-weighted marker totals.
For mL>0, divide lower numerators by mU and upper numerators by mL. Bound the
marginal product using pA+pB<=1, and use C/m<=1/2 for exclusive labels on an
undirected graph. Intersect this with the local-edit interval. All count
arithmetic in the reported cases is integral; final floating endpoints receive
an outward 1e-12 numerical guard, not a general arbitrary-precision certificate.

## Scope and reproduction

The geometry contract is discrete and exhaustive, but not calibrated to an
imaging assay. It excludes arbitrary continuous displacements, missing/merged
cells and cell-type-dependent detection. It cannot be sold as a biological
coverage guarantee or causal cell-interaction estimator.

```sh
python -m unittest discover -s tests -v
python verify_position_uncertainty.py
python study_position_uncertainty.py
```

Full outcomes: [position_uncertainty.json](results/position_uncertainty.json).
Independent checks: [position_uncertainty_verification.json](results/position_uncertainty_verification.json).

The useful retained output is an exact small-case oracle and a tested failure
boundary. This candidate should not be a paper's claimed algorithmic advance.
The broader project still supports a measurement-audit paper direction, but
requires target-matched adjudicated identities/geometry and useful held-out
certificates. Existing expert disagreement is not ground truth; actual scanner
feedback costs remain unmeasured. No fresh patient holdout was consumed here.
