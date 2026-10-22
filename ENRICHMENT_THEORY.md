# Sparse conditional bounds for endpoint-marker covariance

For a simple directed graph W with m>0 edges, d=W1 and e=W'1, define
C=W/m-de'/m^2. Its rows and columns sum to zero. For binary marker vectors a,b,

    theta = a'C b = P(A_source=1,B_target=1) - P(A_source=1)P(B_target=1),

where the random directed edge is uniform. The marginals are degree weighted.
This is an established covariance/modularity-style contrast; no new biological
interaction definition is claimed. It measures association above this particular
abundance baseline, not causality, molecular communication or clinical benefit.
The binary covariance lies in [-1/4,1/4] by Cauchy–Schwarz.

Let raw a,b be binary. Compatible true vectors differ in at most B_A and B_B
coordinates. No true prevalence is fixed. Let delta_a=s_A*x, delta_b=s_B*y,
where s_A=1-2a, s_B=1-2b and binary x,y indicate flips. Then

    theta_true-theta_raw = g_A'delta_a + g_B'delta_b + delta_a'C delta_b,
    g_A=C b, g_B=C'a.

## Initial direct bound (frozen screen)

Each linear gain is optimized exactly by taking up to B positive signed flip
weights; negative gains give the minimum. Zero flips are always allowed.
For the interaction, the sparse term has magnitude at most
min(B_A*B_B, top_B_A(d), top_B_B(e))/m. The rank-one term has magnitude at most
top_B_A(d)*top_B_B(e)/m^2. An independent cap is
sqrt(B_A*B_B*||C||_1*||C||_infinity), from the spectral norm bound and the
Euclidean norms of the binary flip vectors. Take the smaller cap.

Absolute row sums can be computed without materializing C:

    sum_j |C_ij| = d_i/m + sum_(i,j in E) (|1/m-d_i*e_j/m^2|-d_i*e_j/m^2).

Column sums follow symmetrically. Compute Cb and C'a by sparse graph sums and
rank-one corrections. After graph validation the bound uses O(m+n log n) work
and O(m+n) memory; the implementation's duplicate-edge validation additionally
sorts edges. A fixed number of alternating exact linear optimizations supplies
feasible witnesses, without claiming a global solution.

The separated control propagates unchanged positive-contact bounds and exact
degree-weighted A/B marginal bounds independently through contact-p_A*p_B.
Both methods bound THE SAME covariance. Their intersection remains valid.
The compared widths must never be interpreted as widths for the original raw
contact-frequency target.

## Row/column coupling (post-screen development)

For an objective sign sigma in {-1,+1}, write signed linear gains as
p_i=sigma*s_A_i*g_A_i and q_j=sigma*s_B_j*g_B_j. A positive sparse interaction
can occur where sigma*s_A_i*s_B_j>0. A positive contribution from the negative
rank-one term requires the opposite sign. Thus for every fixed i, its possible
interaction gain is at most

    R_i = min(B_B, number of eligible outgoing neighbors of i)/m
          + d_i/m^2 * top_B_B(e among the eligible B-label class).

This follows by dropping negative contributions and separately optimizing the
two nonnegative pieces. It relaxes the requirement that both maxima use the
same B flips, so it is an upper bound, not an exact optimizer. Consequently

    sigma*(theta_true-theta_raw)
        <= top_B_A(p+R) + top_B_B(q).

Compute the analogous column bound and take the smaller one. Select positive
weights only, since budgets are upper limits. Crucially, add R_i BEFORE selecting
the A coordinates; adding a global interaction cap afterward loses this coupling.
Apply both objective signs and intersect with the previous bounds. Sparse
neighbors plus two degree-sorted label classes implement each sign in
O(m+n log n) work. Exhaustive small-world tests check the resulting enclosure.
This is a specialization of classical cardinality relaxations, not a claim to
have invented a new general optimization method.

## Independent sign witnesses and count reference

For every feasible binary population, the exact rational numerator is

    m*sum_(i,j in E)(a_i*b_j) - (sum_i d_i*a_i)*(sum_j e_j*b_j),

with denominator m^2. Python integer arithmetic checks witness signs without
rounding a small covariance. Opposite-sign witnesses prove sign ambiguity ONLY
in the supplied Hamming model on this fixed graph. That model allows overlap
between the two positive marker classes; imposing biological class exclusivity,
additional references or correlated error restrictions can change the result.
The witnesses are not claims that either population is the true tissue state.

If all outgoing degrees equal k, conditioning on the integer true A count K
makes the rank-one correction linear in b. Standard binary McCormick constraints
then give a sparse MILP for each K. Taking the maximum of valid maximization
bounds over every admissible K bounds the global maximum. A missing bound in
even one feasible slice prevents a global numerical certificate. Time-limited
solver incumbents are witnesses only; they are never substituted for dual bounds.
Solver tolerance padding does not make floating solver bounds rational proofs.

## Limits and precedents

Marker error budgets remain supplied assumptions. Geometry, segmentation and
cell correspondence are fixed; the earlier measured graph failure is unresolved.
An enrichment sign says nothing about human outcome prediction or protein
function. Covariance, Newman degree correction, binary metadata-mixing bounds,
KAMP spatial adjustment, graph robustness, McCormick linearization and signed
cardinality optimization are precedents linked in the protocol. Any priority
claim for this exact sparse certificate needs a broader independent review.
