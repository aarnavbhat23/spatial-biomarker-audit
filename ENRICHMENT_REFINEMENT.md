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
