# Complete joint measurement audit, frozen 2026-10-28

Question: what is the cheapest fixed measurement set sufficient to decide the
sign of graph covariance when both cell identities and graph edges are unknown?
This is an exact small-population reference problem, not patient validation,
scanner-cost measurement, or a scalable acquisition method.

Use four cells. Categories 0,1,2 mean other, epithelial, T cell. Enumerate every
category assignment within one raw-call error per marker. Enumerate every graph
obtained by toggling at most one specified fragile undirected edge of the raw
graph. Represent each undirected edge by both directed orientations. The target
is the three-way sign (-1,0,+1) of

    N = m * sum_edges a[u]*b[v] - sum_edges a[u] * sum_edges b[v].

Use all three cases below, without replacing failures:

1. raw labels (1,2,0,0), path edges 01,12,23, fragile edges 02,12;
2. raw labels (1,0,2,0), cycle edges 01,12,23,03, fragile edges 02,03;
3. raw labels (1,2,1,2), star edges 01,02,03, fragile edges 12,03.

Available actions: one persistent categorical review per cell (cost 1), exact
independent adjudication per cell (cost 4), exact geometry verification per
fragile edge (cost 2). Use E=0 and E=1 wrong routine-reviewed cells separately.
Exact measurements and costs are hypothetical assumptions, not empirical facts.
Do not interpret two agreeing human readers as an exact adjudication.

For each pair of admissible worlds with different target signs, a fixed plan
must either contain an exact measurement on which they disagree, or more than
2E routine measurements on which their categorical identities disagree. This
is the disjoint-Hamming-balls criterion applied to a decision, not a new coding
theorem. It is necessary and sufficient: a shared transcript exists exactly
when neither condition holds. Check every plan (2^10 subsets), so minimum cost
is global only for this complete finite problem and a fixed, nonadaptive plan.

Compare optimal cost with full exact measurement and with a cost-normalized
greedy pair-cover rule. Report failure under routine-only, label-only, and
geometry-only action restrictions. Save complete worlds, plans and pair counts.
Independently verify each selected plan by enumerating every feasible noisy
transcript and confirming that no transcript has two signs. Independently
enumerate lower-cost plans using the pair criterion to certify optimality.

Existing work already treats persistent-noise testing and decision trees:
[Jia, Navidi, Nagarajan and Ravi](https://jmlr.org/papers/v25/23-1484.html).
The purpose here is to determine whether the spatial target provides useful
structure beyond that general framework. No positive novelty conclusion from a
four-cell exhaustive solver is permitted. A next contribution would require a
scalable certified optimization method, measured error/cost assumptions, and
independent data; all three remain open.
