# Sharp joint-marker bounds in a restricted scalable case

**Both markers can be uncertain without enumerating every possible population.**
For a fixed simple directed graph with unit edge weights, at most one wrong A
call and any specified B error budget, the full-census pair fraction has an
attainable lower and upper bound computable in O(n+M) operations using selection
and edge counts. This is a restricted algorithmic result, not a general solution
for large error budgets in both markers or an estimated assay error rate.

The [development protocol](COUPLED_MARKER_PROTOCOL.md) fixes three previously
used patients, three radii and three B budgets before reading these outcomes.
These are sensitivity calculations around source calls. They do not validate
biological truth, an acquisition policy or an independent cohort. The algorithm
uses elementary monotonicity and top-k sensitivity; no priority claim over the
combinatorial optimization literature is made.

## Derivation

Write N=sum_(i,j in E) A_i B_j, w_i=sum_(i,j in E) B_j and
z_j=sum_(i,j in E) A_i. A and B are the recorded binary calls. Since all edge
weights are nonnegative, a maximizer needs only 0-to-1 flips; a minimizer needs
only 1-to-0 flips. Opposite-direction flips can be removed without worsening the
objective or exceeding either budget. Thus consider no A flip or one eligible
A flip, followed by the optimal B flips.

For the maximum, eligible B cells have B=0. Let k be the B budget capped at their
count, T the sum of their top k scores z_j, t the kth score, c_> the number of
eligible scores above t, and s=k-c_>. With no A flip the gain is T. Flipping
A_i=0 to 1 adds w_i and increases each outgoing eligible neighbor score by one.
Let h_i count those neighbors with z_j>t, and e_i those with z_j=t. Because the
graph is simple and scores are integers, neighbors below t cannot improve the
top-k sum: a score t-1 can only tie the existing threshold. Every increment above
t improves the sum, and at most s threshold increments fit. The exact gain is

    w_i + T + h_i + min(e_i, s).

For the minimum, eligible B cells have B=1; T and t now summarize their top-k
removal scores. Let c_>= count eligible scores at least t and s=c_>=-k. Flipping
A_i=1 to 0 subtracts w_i and decreases each outgoing eligible score by one. With
h_i/e_i defined relative to this threshold, the largest decrease is

    w_i + T - h_i - max(0, e_i-s).

Above-threshold decrements each reduce the top-k sum by one. Threshold
reductions matter only after exhausting the s tied replacements. Take the
maximum of these decreases and the no-A-flip T. The lower endpoint is
(N-largest_decrease)/M; the upper endpoint is (N+largest_increase)/M. When k=0,
use T=h=e=0. Each endpoint is attained by the selected A flip (or none) and the
corresponding top-k B flips.

Top-k selection takes linear work, and h/e/w/z for every possible A cell are
computed together by edge counts. Additional working memory is O(n+M) in this
vectorized implementation. The simple, unit-weight graph requirement matters:
parallel edges or arbitrary weights can change a neighbor's score by more than
one and invalidate this threshold formula. Prepared graphs are checked for
uniqueness. This function does not solve two large error budgets, partial
acquisition or review allocation.

## Validation and development results

The unit test exhausts all 64 recorded three-cell binary populations, all 64
possible true populations, and B budgets 0..3. Direct extrema equal both endpoints.
A separate counterexample shows why adding observed-label one-marker sensitivity
bounds can fail: on the single edge 0->1 with all recorded calls zero, changing
A0 alone or B1 alone changes nothing, but changing both gives a target of one.
The existing label-blind degree outer bound remains valid.

At radius 80 native pixels, full-census interval widths in percentage points:

| Patient | B budget at 0.1% | Sharp width | Degree outer width | B budget at 0.5% | Sharp width | Degree outer width |
|---|---:|---:|---:|---:|---:|---:|
| CRC02 | 57 | 0.6234 | 0.9248 | 282 | 2.3684 | 3.2495 |
| CRC11 | 73 | 0.4025 | 0.4905 | 365 | 1.7502 | 2.2680 |
| CRC19 | 12 | 0.8047 | 1.0982 | 59 | 3.0332 | 3.9532 |

A budget is **one cell** in every row. Across all 27 panel/radius/budget
calculations, the exact bound takes 0.35-36.21 milliseconds in this local run,
excluding input loading and graph uniqueness checks. These timings are not
instrument latencies. The bound is 12%-66% narrower than the degree outer bound
across the entire screen, including B=0 cases. At 0.1% B uncertainty, sharp widths
across the nine graphs are 0.4025-0.8766 percentage points: six exceed the target
0.5 percentage points even though A uncertainty is only one cell. Better
computation cannot remove this information limitation. At 0.5%, all exceed it.
Full bounds, times and input hashes are in
[`coupled_marker_study.json`](results/coupled_marker_study.json).

## Next gate

Extend exact or certified optimization to meaningful budgets in both markers,
with lower/upper solver bounds and explicit gaps. Compare a direct integer
optimization formulation on small graphs before using any scalable relaxation.
Only then combine partial acquisition and review, including correlated review
errors. The [finite review study](FINITE_REVIEW_RESULTS.md) explains why repeated
responses alone cannot substitute for orthogonal evidence. Budget calibration,
a new cohort and the [lab pilot](LAB_PILOT_PROTOCOL.md) remain necessary for an
empirical instrument-saving claim.

Reproduce after the existing data preparation:

```bash
python -m unittest discover -s tests -v
python study_coupled_marker.py
```
