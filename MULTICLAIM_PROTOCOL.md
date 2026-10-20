# Shared spatial-claim screen

Exploratory development protocol. Retain all existing single-claim results. Use
three nested, geometry-fixed claims: positive PanCK/CD3e directed-pair fractions
within 40/60/80 native pixels. Shared measured cells identify pairs in all claims.
Patients CRC02/CRC11/CRC19, hypothetical camera side 1024, batches of 16 queries,
five random repeats, width .005 per claim, joint alpha .05 via alpha/3 per claim.
No confirmation patient is consumed.

Derive a common sampling distribution minimizing the largest predictable
normalized unresolved-group importance weight over all active claims. Compare:
envelope PPS, equal mixture of active-claim PPS distributions, uniform active
groups, and focused single-claim PPS with global observation caching. In the
focused control only the chosen claim receives a statistical update, but all
claims receive deterministic cache updates. Also compare a deterministic greedy
field control with batch size 16 and normalized unresolved incident-edge scores
summed across claims. All methods stop only after every claim reaches the width.

Freeze unresolved query subsets and the sampling distribution across each batch.
For all-claim updates use PRE-DRAW maximum importance weights as bounds. The
selected group's smaller bound is not predictable for a martingale whose mean
is established before drawing the group. Repeated draws retain repeated values.
Use simultaneous coverage, field exposures and feedback rounds as outcomes.

Prove the deterministic bounds are sharp under unrestricted unknown marker bits,
check whether resolution utility is submodular, and construct a small example
showing whether sharing can save measurement. Do not infer a general greedy
approximation guarantee from the geometry. Source calls are treated as exact;
actual instrument costs, clinical benefit and a new confidence-sequence theorem
are not established by this screen. Seed 20261109 is an arbitrary identifier.

Before outcomes, strengthen the focused control: add focused-reuse, which updates
every other claim whose remaining groups have full support under the chosen
claim's proposal, using each claim's own predictable maximum importance weight.
This prevents an apparent benefit caused merely by discarding usable samples.
Use a common total budget of 4096 randomized queries for every stochastic policy.
