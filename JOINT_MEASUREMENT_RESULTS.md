# Complete joint measurement plans: useful decisions, no algorithmic win

**Exact fixed-plan costs depend on both label and geometry uncertainty. A simple
redundancy-pruned greedy control matches the optimum in all six small cases.**
There is no supported optimization advantage to claim from this screen.

The [primary protocol](JOINT_MEASUREMENT_PROTOCOL.md) was committed at `521524a`
before evaluation. It enumerates every admissible categorical population and
every permitted fragile-edge toggle, rather than sampling a hypothesis bank.
These are four-cell synthetic problems with hypothetical costs and ideal exact
adjudications. They do not validate a clinical decision, scanner savings, or
real-world error budgets. The existing NuCLS reference disagreement continues
to rule out treating routine expert agreement as proof of truth.

## All cases and controls

Costs use routine cell review=1, exact cell adjudication=4 and exact edge
verification=2. E counts wrong distinct routine-reviewed cells; repeating a cell
does not provide another independent measurement. All plans are fixed before
observing responses. The target is the three-way covariance sign, including zero.

| Case | Complete worlds | E | Optimal | Original greedy | Pruned greedy | Full exact |
|---|---:|---:|---:|---:|---:|---:|
| Path | 48 | 0 | 2 | 4 | 2 | 20 |
| Path | 48 | 1 | 8 | 16 | 8 | 20 |
| Cycle | 48 | 0 | 8 | 8 | 8 | 20 |
| Cycle | 48 | 1 | 20 | 24 | 20 | 20 |
| Star | 39 | 0 | 6 | 7 | 6 | 20 |
| Star | 39 | 1 | 12 | 16 | 12 | 20 |

The stronger control was added after observing redundant actions in the original
greedy plans. Its deletion rule was committed at `6d76bb3` before evaluation;
see [the amendment](JOINT_MEASUREMENT_CONTROL_AMENDMENT.md). All original outputs
remain in the artifact. Its agreement with exact optimization prevents presenting
the original greedy gaps as evidence for a new algorithm.

For E=1, the path requires exact adjudication of cells 0 and 1 only; the star
requires those two adjudications and both edge verifications; the cycle requires
all four adjudications and both edge verifications. Label-only plans cannot
resolve the cycle or star even if all labels are adjudicated. Geometry-only
plans fail in every case. Routine-only plans fail everywhere at E=1, including
the path where routine reviews suffice when error-free.

## Exact decision-specific criterion

For each pair of admissible worlds with different signs, a plan must either
measure an exact quantity on which they disagree or review more than 2E distinct
cells on which their categorical identities disagree. This is necessary and
sufficient. When exact answers agree and the restricted Hamming distance is at
most 2E, split the differing positions between the two answers to construct a
transcript within E errors of both worlds. Conversely, the triangle inequality
rules out any shared transcript when their distance exceeds 2E.

Applying this to every opposite-decision pair is sufficient for a decoder:
all worlds compatible with any feasible transcript then have the same sign.
It does not require identifying the entire specimen. Enumerating all 1,024
action subsets certifies the minimum fixed-plan cost. Persistent-noise testing,
Hamming separation and pair-cover optimization are established foundations;
see [Jia et al.](https://jmlr.org/papers/v25/23-1484.html). This is an exact spatial
reference problem, not a novel general decision-tree theory or an adaptive optimum.

## Verification and next research gate

Run `python study_joint_measurement.py` and `python verify_joint_measurement.py`.
The [saved results](results/joint_measurement.json) include complete worlds,
integer targets, all selected actions, restrictions and costs. Independent replay
rebuilds the complete Cartesian population, reconstructs graph targets, rejects
every cheaper plan, and enumerates every feasible noisy transcript of each
reported plan. Unit tests additionally compare the pair criterion with transcript
enumeration for every action subset in another small problem.

The useful next problem is a scalable separation oracle for jointly uncertain
identities and geometries, paired with computable cost lower bounds. Require it
to beat a redundancy-pruned greedy control at equal verified decision validity,
not merely a deliberately expensive census. Realistic geometry constraints,
empirically supported references and independent patient data must precede an
instrument-saving claim. These six examples do not demonstrate that such an
algorithmic advantage exists; retain the option to stop or pivot if it fails.
