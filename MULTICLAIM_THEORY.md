# Shared spatial certificates: range, overlap and complementarity

**One acquisition can update several spatial claims with a joint error guarantee.**
The best distribution for the largest normalized unresolved-group weight has a
closed form. This is a range-optimal proposal, not a proof of minimum physical
measurement cost. Natural pair-resolution utility has complementary effects,
so a standard submodular-coverage greedy guarantee does not apply directly.

## Fixed targets and observations

There are Q geometry-fixed directed edge sets E_j, of sizes M_j>0. An edge's
outcome is A_u B_v, where A and B are binary marker calls on cells. Acquisition
reveals both marker calls in one or two physical camera fields. A shared cache
therefore resolves edges across every claim. An edge is resolved if both endpoint
calls needed by the product are observed, or one observed factor is zero.

Before acquisition, let K_j count known positive edges, U_j unresolved edges,
and m_jg unresolved edges owned by physical field-pair group g. Then

`L_j = K_j/M_j`, `Ubound_j = (K_j+U_j)/M_j`.

These bounds are jointly valid for every acquisition history. Completing all
unknown marker bits with zero attains every lower endpoint simultaneously;
completing them with one attains every upper endpoint simultaneously. Thus each
coordinate bound is sharp under unrestricted binary completions. The Cartesian
box need not equal the full joint feasible set. Biology may impose extra constraints
that this model does not use.

## A minimax range proposal

For active claims define `a_jg=m_jg/M_j`, `e_g=max_j a_jg`, and `H=sum_g e_g`.
Assume at least one unresolved active edge. Draw a group using

`p_g = e_g/H`.

The selected physical fields reveal Y_jg positive PRE-ACQUISITION unresolved
edges of claim j in that group. The residual importance value is

`Z_j = K_j/M_j + Y_jg/(M_j p_g)`.

Its conditional expectation is T_j, the target positive fraction. Before drawing
g, its support lies in `[K_j/M_j, K_j/M_j+H]`. The support can exceed [0,1],
although the target is a fraction. The betting calculation must use this support.

**Proposition.** Among proposals supporting all unresolved active claims, this
proposal minimizes `B(p)=max_(j,g) a_jg/p_g`, with optimal value H.

**Proof.** If B is a feasible maximum weight, `p_g>=e_g/B`. Summing gives
`1>=H/B`, so B>=H. For `p_g=e_g/H`, every ratio is at most H, and at least one
claim attains it at every nonempty group. The lower bound is attained.

This controls the worst possible positive-group residual without assuming any
probabilistic relationship between cell labels. It need not minimize variance,
the number of camera fields, or feedback latency. Multiple importance sampling
and minimax proposal design have substantial precedent; this specialization is
not presented as a new general importance-sampling principle.

## What overlap buys, and what it cannot buy

`max_j(U_j/M_j) <= H <= sum_j(U_j/M_j)`.

The lower bound follows by comparing a maximum with each row; the upper bound
follows by replacing a column maximum with its sum. If every normalized unresolved
group profile is identical, sharing costs no additional range as Q increases.
If profiles have disjoint group supports, H equals the sum. Thus overlap has a
precise range benefit; it does not remove the multiple-claim error budget.

For comparison, at a fixed number n of all-claim updates, adaptive Hoeffding
concentration with the deterministic initial H_0 yields simultaneous intervals
of radius `H_0 sqrt(log(2Q/alpha)/(2n))`. The known positive offset may change,
but each predictable range length is at most H_0 because unresolved counts only
decrease. Union bounding the Q fixed-time tails gives this formula. Full width
epsilon is sufficient once `n >= 2 H_0^2 log(2Q/alpha)/epsilon^2`.
This deliberately conservative fixed-time bound is not the stopping rule used
in the experiments and does not establish an optimal sample complexity.

There is also a useful limitation of that surrogate. Let c_g>=0 be a known
additive query cost, ignoring batching/caching interactions. The same envelope
proposal minimizes `B(p)^2 sum_g p_g c_g`: since `B>=H` and `Bp_g>=e_g`,

`B^2 sum_g p_g c_g >= B sum_g e_g c_g >= H sum_g e_g c_g`,

with equality at p=e/H. Consequently this range-based cost bound cannot adapt
the proposal to actual camera cost. A better physical-cost policy needs a
sharper variance criterion or an explicit model of field reuse and batching.

## Joint sequential validity

Freeze the proposal, remaining groups and importance bounds for each batch.
Sample indices independently with replacement. For each active claim the
existing nonnegative betting factors have conditional mean one at its target;
the product across a completed batch also has mean one. Correlation between
claims and between biological cells does not change this argument: the fixed
population outcomes are conditioned on, and random acquisition indices provide
the sampling randomness.

Give claim j an error budget alpha_j, with sum alpha_j<=alpha. Its confidence
sequence is valid relative to the full acquisition history, including the other
claims' observed labels. A union bound gives simultaneous coverage for all claims
at every inspected boundary. No independence between confidence sequences is
needed. Stopped claims retain their intervals, intersected with later valid
deterministic bounds. Their betting processes receive no unsupported updates.

The bound used for a different proposal is `max_g a_jg/p_g`, computed BEFORE
the group draw. Using the selected group's smaller bound to choose the bet is
not justified by the pre-draw conditional-mean argument. A focused proposal may
miss some claims' unresolved groups; those claims cannot receive an unbiased
importance update. The focused-reuse baseline updates every claim with full
support, and caches all observations regardless.

## Why ordinary greedy coverage theory is insufficient

Take one edge whose endpoints lie in different camera fields, with both true
marker factors positive. Observing either field alone resolves nothing;
observing the second after the first resolves the edge. Marginal utility rises
from zero to one, violating submodularity of the resolved-edge count.
With a known realization of all-zero factors, either endpoint resolves the edge,
giving an ordinary coverage function. The general objective mixes these effects.
Even independent positive-probability marker factors do not automatically give
adaptive diminishing returns: observing the first factor as positive increases
the chance that the second observation is decisive.

On an unrestricted abstract graph with one cell per physical field and all true
calls positive, k acquired fields resolve exactly the induced edges on k vertices.
Maximizing this utility contains densest-k-subgraph, and deciding whether k
vertices resolve k(k-1) directed edges contains CLIQUE. This establishes a
general-graph computational obstruction. It is not a hardness proof for the
specific radius graphs and 2-D camera geometry of the ORION experiment. Alternative
submodular surrogate goal functions may still exist.

## Development evidence and scientific gate

[Frozen development protocol](MULTICLAIM_PROTOCOL.md), [implementation](study_multiclaim.py),
[complete development results](results/multiclaim_study.json).

Three nested radius claims (40/60/80 native pixels), batch 16, joint alpha .05,
width .005 per claim, hypothetical camera side 1024; five repeats per stochastic
method. All 75 stochastic paths and three deterministic runs completed with all
targets covered at every checked boundary; this small experiment does not prove
coverage. The proof supplies coverage under the fixed-graph/exact-call assumptions.

| Patient | Envelope fields / rounds | Mixture | Focused reuse | Deterministic greedy |
|---|---:|---:|---:|---:|
| CRC02 | 145.0 / 10.6 | 170.2 / 12.8 | 164.8 / 12.2 | 160.0 / 10.0 |
| CRC11 | 148.0 / 10.8 | 173.2 / 12.6 | 169.0 / 12.4 | 160.0 / 10.0 |
| CRC19 | 62.4 / 5.8 | **59.8 / 5.8** | 60.4 / 5.6 | 64.0 / **4.0** |

Uniform and focused-without-statistical-reuse controls remain in the JSON. The
envelope proposal loses to mixture on CRC19; greedy can win with high feedback
cost. Range optimality has not become universal acquisition optimality.

The next gate is a frozen patient-disjoint confirmation against mixture,
focused reuse and greedy, followed by actual acquisition-cost measurement and
an independent cohort. Three radius variants of one marker pair are a starting
case, not a demonstration across proteins or clinically meaningful decisions.

## Relevant prior work

- [Multiple-importance proposal selection](https://arxiv.org/abs/1805.00829)
- [Optimal mixture weights](https://arxiv.org/abs/1411.3954)
- [Gambling confidence sets for bounded vectors](https://arxiv.org/abs/2402.03683)
- [Adaptive submodularity](https://arxiv.org/abs/1003.3967)
- [Stochastic Boolean function evaluation and submodular cover](https://arxiv.org/abs/1303.0726)

The contribution being investigated is cost-aware shared spatial certification,
with these limitations made explicit, rather than a new general confidence
sequence, a generic greedy approximation theorem, or a proprietary model evaluation.
