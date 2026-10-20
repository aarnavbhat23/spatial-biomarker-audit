# Conditional spatial bounds with fallible partial references

The graph is fixed, directed, simple and unit weighted. Raw binary marker
vectors rA,rB are known at every cell. A partial reference z has values -1
(unreviewed), 0 or 1. For each marker separately the true vector x satisfies

    sum_i [x_i != r_i] <= B
    sum_{i:z_i>=0} [x_i != z_i] <= E.

Budgets are supplied assumptions. They are not estimated from the data. Marker
errors may be dependent; no stochastic independence of errors is required.
Repeated responses from one cell do not create new constraints when errors are
counted by systematically wrong cells. This is review after raw-call census,
conditioned on measured-image coordinates, not an H&E-only acquisition model.

## Exact signed linear oracle

Let f=x XOR r. Partition cells into D (reviewed raw/reference disagreement),
G (reviewed agreement), and U (unreviewed). For d,g,u flips in those groups,

    d+g+u <= B,        |D|-d+g <= E.

For any signed objective c'x, write c'r + sum_i c_i*(1-2*r_i)*f_i. For fixed
counts, selecting the largest d/g/u signed weights within the respective groups
is exact by exchange. Sort the three lists and form prefix sums. Scan d,g;
the best admissible u is a prefix maximum up to B-d-g, including u=0. Forced
negative-weight flips are retained when needed to meet the disagreement bound.
Negating c gives the minimum. No feasible count triple means an empty set.
Time is O(n log n + B^2), memory O(n+B), with B capped at n.

An individual proposed value x_i=v is feasible exactly when the remaining raw
and reference budgets are nonnegative and the remaining disagreement count is
at most their sum. This supplies a coordinate possibility mask without n oracle
calls. These are classical cardinality/Hamming-constraint calculations; no new
general optimization primitive is claimed.

## Bilinear outer bound

The target numerator is T=a'Wb, divided by m graph edges. Let h equal the
reference where reviewed and raw calls elsewhere. h need not itself be feasible.
Write a=hA+deltaA and b=hB+deltaB. Then

    T = hA'W hB + (W hB)'deltaA + (W' hA)'deltaB
                    + deltaA'W deltaB.

The exact linear oracle bounds each linear term over its own feasible marker
set. Each delta_i=(1-2*h_i)*t_i, where t_i indicates a change from h. An edge's
interaction sign is positive when its hybrid endpoint labels agree and negative
when they differ. Ignore edges whose endpoint cannot change in any feasible
world. For each sign, count incident degrees dA,dB in this remaining edge set.

Let kA,kB be exact maxima of the number of changes from h. The unsigned
interaction contribution of either sign is at most

    min(kA*kB, max_feasible_a sum_i dA_i*tA_i,
                 max_feasible_b sum_j dB_j*tB_j).

The first bound uses graph simplicity; the others count each interaction edge
at one of its changed endpoints. The latter maxima are linear-oracle calls
since t_i=h_i+(1-2*h_i)*x_i. Subtract the negative cap from the lower linear
bound and add the positive cap to the upper linear bound, then intersect [0,m].
This proves the outer interval for every feasible pair of marker populations.

Alternating exact linear optimization over a and b gives feasible lower/upper
witness populations. It is not a global optimizer. We retain two fixed starts
and at most 12 alternating iterations; validity of witnesses and outer bounds
does not depend on convergence. A witness and its outer endpoint bracket the
exact endpoint, exposing any remaining optimization gap. A fixed number of
oracle calls and graph passes gives O(m+n log n+B_A^2+B_B^2) work per iteration.

## Adaptive review and validity

Every review adds constraints to the feasible set under a global systematic
reference budget E. Thus sets are nested. Each bound covers every member of
its current set, regardless of how the next field was selected. Intersecting
the interval with all earlier valid intervals preserves validity and enforces
monotone bounds. No optional-stopping correction is needed for these deterministic
conditional bounds; their validity still depends entirely on the stated error
budgets and fixed graph being appropriate.

The policy receives only raw calls, graph/field identifiers and already queried
responses. The hidden reference array is accessible solely through a field
observation callback. Endpoint witnesses may guide which field to review, but
they are not a complete minimax set of possible reference responses. The
lookahead control uses three feasible-population proxy responses for ranking;
the final bound does not assume one of those proxies is correct.

## Verification and precedent

Tests compare signed linear optima and coordinate feasibility against every
binary vector on five cells across 180 generated constraint configurations.
Graph bounds, witness feasibility and empty sets are checked against all 256
two-marker worlds on four cells across 120 configurations. Additional tests
retain the systematic-reference precision floor and check all five policies
for valid monotone transcripts on a small graph.

Active review and collection already have substantial precedent:
[active label cleaning](https://doi.org/10.1038/s41467-022-28818-3),
[active statistical inference](https://proceedings.mlr.press/v235/zrnic24a.html),
and [robust active sampling](https://arxiv.org/abs/2511.08991). Our target here is
a conditional bilinear finite-graph claim under two fallible annotation sources,
not classifier retraining or a new general active-learning framework. A claim
of priority for this exact bound would require a broader literature review.
