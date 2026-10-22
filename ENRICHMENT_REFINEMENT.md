# Post-screen computational reference — frozen before solver outcomes

The primary screen is complete and failed its combined gate: median width
reduction .340148 but only 1/6 signs certified. Four primary patients already
have feasible witnesses on both sides of zero. MB-0045 remains unresolved by
the fast bound despite both alternating witnesses being negative.

Run a computational reference on that unresolved case only; this is a post-
outcome refinement, not a new independent confirmation or changed primary gate.
For constant outgoing degree k, m=nk and the A abundance is sum(a)/n. Condition
on each possible A count K in [sum(raw_a)-B_A,sum(raw_a)+B_A] intersect [0,n].
Then theta=sum_edges(a_u*b_v)/m - K*sum_j(e_j*b_j)/(n*m).

For each K, maximize theta using ordinary binary a,b variables and continuous
edge products with exact binary McCormick constraints, both Hamming budgets,
and sum(a)=K. No dense n-by-n matrix is required. Use SciPy/HiGHS MILP with
1-second per-slice limit, zero relative gap target, presolve enabled. Combine
maximization upper bounds by their maximum over all count slices; any missing
dual bound gives an unresolved global bound, never a fabricated certificate.
Retain statuses, incumbents, dual bounds and runtime. Pad bounds by 1e-7 for
solver tolerance; these are numerical solver bounds, not rational certificates.
Verify all returned integer witnesses independently with integer numerator
m*sum_edges(a_u*b_v)-(sum_i d_i*a_i)*(sum_j e_j*b_j).

Also independently verify every primary fast-bound witness and its sign using
that exact integer numerator. Opposite-sign feasible witnesses prove ambiguity
within the supplied fixed-graph Hamming model. A loose outer interval alone
cannot prove ambiguity. Count slicing, MILP and McCormick linearization are
classical, not new optimization primitives. This diagnoses where further
optimization could help versus where additional measurements are necessary.

## Structural bound after the limited MILP reference

All 21 requested one-second MB-0045 slices stopped without dual bounds; measured
wall time was about 111 seconds. This does not certify either sign. Before
its outcomes, freeze a separate development tightening: for flip signs
sA=1-2a, sB=1-2b and objective sign sigma, the positive sparse interaction has
sigma*sA_i*sB_j>0; the positive rank-one interaction has the opposite sign.
For each A vertex, cap the first by min(B_B,eligible outgoing degree)/m and
the second by d_i/m^2 times the top-B_B incoming-degree sum in the appropriate
B-label class. Add that row cap to the signed A linear gain BEFORE top-B_A
selection; bound B linear gains separately. Repeat by columns and for both
objective signs. Intersect with the original direct and separated bounds.

This exploits row-specific coupling lost by adding a global interaction cap
after separate linear maximization. It is a classical cardinality relaxation
specialized to sparse-plus-rank-one structure, not a new optimizer. Re-evaluate
all 36 already-read cases as DEVELOPMENT ONLY. Preserve original outcomes and
the original failed gate. No additional patients are screened to hunt for a
passing confirmation. Exhaustive binary-world checks must validate the new
bound before interpreting its outcome.
