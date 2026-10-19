# Residual acquisition stress protocol

Exploratory protocol written before running the new policy. The previously seen
CRC02, CRC11 and CRC19 panels are development data. No new confirmation is claimed.

The earlier fixed-scale experiment sampled with replacement. A stronger comparison
must avoid remeasuring known regions and share observations across nested scales.
Partition target edges by the 128-pixel anchor core. A larger measured core reveals
all its nested fine-core totals. Retain that information across scales; do not count
extra halo observations unless they are actually included by the implementation.

At each step, choose a scale using the history, then draw a region with probability
proportional to its remaining unknown edge count. Let K be the already known positive
count, U the remaining edge count, and M the fixed total. For selected region t,
Z = K/M + (U/M) * Y_unknown(t)/M_unknown(t) has conditional expectation T. Its
predictable support is [K/M,(K+U)/M]. Use bounded betting with those supports and
intersect with the exact known-data interval. This is an application of established
adaptive importance sampling and betting, not a new confidence-sequence theorem.

Compare fixed cores 128 and 2048 against a geometry-only adaptive policy choosing
among 128, 512 and 2048. The adaptive score is expected newly learned edge count
under residual PPS divided by rectangle area plus per-acquisition overhead. No
unmeasured labels enter scale selection. This is a heuristic, not an optimal policy.

Use 10 repeats, seed 20261027, width targets 0.04 and 0.005, overheads 0 and 500,000
native pixel-area equivalents per acquired rectangle, at most 2048 acquisitions.
Report exact union area, summed rectangle area, acquisition count, both cost models,
precision attainment, and time-uniform coverage at all inspected checkpoints.
Compare only attained precision. Cache sharing is by fully measured nested cores;
additional halo/cell-level reuse remains outside this experiment.

Prior work: Neufeld et al. ICML 2014, Adaptive Monte Carlo via Bandit Allocation;
Waudby-Smith and Ramdas, Confidence sequences for sampling without replacement
(NeurIPS 2020), and Estimating means of bounded random variables by betting (2024).
