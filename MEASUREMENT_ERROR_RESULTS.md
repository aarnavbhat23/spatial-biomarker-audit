# Sampling precision does not remove marker-call uncertainty

**In all three development panels, allowing 0.5% of PanCK calls to be wrong
leaves a sharp full-census target range wider than the requested 0.5 percentage
points.** The corresponding range is 1.00-2.97 percentage points across the
three patients and three distance claims. These are adversarial sensitivity
budgets, not estimated assay error rates or observations of biological truth.

[Protocol](MEASUREMENT_ERROR_PROTOCOL.md), [bounds](measurement_error.py),
[study](study_measurement_error.py), [all results](results/measurement_error_study.json).
The protocol was written locally before execution, not separately committed
as a prospective confirmation. Only development CRC02/CRC11/CRC19 were used.
The existing acquisition policy and three radii were unchanged.

## A bound available before acquisition

Fix a directed graph E, M=|E| and a pair target T=M^-1 sum_(u,v) A_u B_v.
Suppose recorded binary calls differ from true calls at at most b_A cells for
A and b_B for B. Coordinates and graph membership remain exact in this model.
Let d_out and d_in be graph degrees. Then

\[
|T-\widetilde T|\le\Delta=
\min\left\{1,\frac{\operatorname{top}_{b_A}(d_{out})+
 \operatorname{top}_{b_B}(d_{in})}{M}\right\}.
\]

Here top_b denotes the sum of the b largest entries, capped at the number of
cells. For each edge, |AB-A'B'| is at most the indicator that its A endpoint
changed plus the indicator that its B endpoint changed. Summing affected edges
and bounding by the largest degrees proves the inequality. Double counting
only makes this an outer bound. It need not be sharp.

Enlarge every valid interval [L,U] for the fixed recorded-call target to
[max(0,L-Delta), min(1,U+Delta)]. The original simultaneous coverage guarantee
then applies to T under the deterministic error budgets. There is no new
independence assumption: any fixed corrupted call table is a finite population
for the existing sampler. Margins use geometry alone, with no hidden labels.
If budgets are estimated with failure probability delta, that uncertainty
must also be accounted for, for example by an alpha+delta union bound. The
study does not estimate these budgets.

## Sharp full-census ambiguity when only A may be wrong

With all B calls exact and observed, let w_u=sum_(u,v) B_v. From a fully observed
recorded table, the attainable target interval with at most b_A wrong A calls is

\[
\left[\widetilde T-
 \frac{\operatorname{top}_{b_A}\{w_u:\widetilde A_u=1\}}{M},
 \widetilde T+
 \frac{\operatorname{top}_{b_A}\{w_u:\widetilde A_u=0\}}{M}\right].
\]

The endpoints are attained by flipping the largest-weight positive or negative
calls. This is sharp because each flip changes exactly w_u edge outcomes,
and opposite-direction flips cannot improve the respective extremum. This
label-weighted calculation is a full-census diagnostic, not an online proposal.
Exhaustive small-graph tests check both endpoint attainability and two-marker
geometry bounds over all allowed binary completions.

If this sharp range is wider than epsilon, the same observed table is
compatible with two true populations separated by more than epsilon. An
always-width-at-most-epsilon interval based on that table cannot cover both.
Even if it is randomized, its worst-case coverage over those two populations
is at most 1/2. Thus uniform 95% coverage at that width is impossible without
new assumptions or independent label information. This is an identifiability
limit under the bounded-error model, not a claim about the actual assay.

| Patient | Cells | Census width at 0.1% errors, radii 40/60/80 (pp) | Width at 0.5% errors (pp) |
|---|---:|---|---|
| CRC02 | 56,220 | .357 / .303 / .279 | 1.284 / 1.092 / .998 |
| CRC11 | 72,922 | .380 / .319 / .291 | 1.543 / 1.303 / 1.184 |
| CRC19 | 11,654 | .788 / .703 / .657 | 2.968 / 2.695 / 2.605 |

Budgets are ceil(n*f); in the smallest panel this is 12 cells at f=.001.
Error-location influence matters, so a percentage alone is insufficient.

## Injected-error acquisition results

At fractions .001 and .005, inject either uniformly random A flips or an exact
maximal absolute radius-80 shift with B fixed. Apply the same cell set to all
three radii. Five repeated acquisitions share each fixed corrupted table;
there are 12 corrupted tables, not 60 independent corruption replicates.

| Outcome over 60 acquisition paths | Paths jointly covering all three targets at every checkpoint |
|---|---:|
| Original intervals, recorded corrupted-call targets | 60/60 |
| Original intervals, original source-call reference targets | 10/60 |
| Geometry-enlarged intervals, original source-call reference targets | 60/60 |

All 30 adversarial paths failed source-reference coverage somewhere. At .001,
radius-80 adversarial bias was .199/.201/.527 pp for CRC02/CRC11/CRC19;
at .005 it was .732/.832/2.161 pp. Even a small bias can defeat a narrow
interval near its boundary, so these results do not imply a minimum bias of
.5 pp is necessary for failure. Enlarged intervals were wider than the .5 pp
stopping target; the original sampler's stop flag certifies recorded calls
only. Empirical 60/60 coverage is a sanity check, not the proof of coverage.

## Decision

Continue shared spatial certification, with a sharper question: **which
measurements and independent label reviews are sufficient for a valid spatial
claim at a real cost?** More imaging alone cannot solve the demonstrated
bounded-error ambiguity. A next method should allocate between imaging and
reviewing influential marker calls, retaining a valid joint certificate and
comparing with fixed QC plus greedy/full acquisition. It must count review
cost and specify an independent review standard; relabeling the same image
with another threshold is not independent ground truth.

Segmentation-aware phenotyping already has substantial precedent, including
[STARLING](https://pmc.ncbi.nlm.nih.gov/articles/PMC11700195/).
The simple perturbation and identification bounds above are not themselves
claimed as a new general statistical method. Cell-coordinate, missing-cell
and graph-membership errors remain outside this study.

Reproduce: `python -m unittest discover -s tests -v`, then
`python study_measurement_error.py` with the existing development arrays.
