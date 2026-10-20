# Fallible review: when more observations cannot identify the answer

**Repeated review can resolve a bounded number of wrong response events, but
cannot guarantee resolution of systematic cell errors.** An exact finite model
separates these assumptions rather than treating a second response as truth.
A trusted exact reference action restores identifiability in the tiny example;
this reference is an explicit idealization, not a demonstrated perfect assay.

## Model and complete planning check

The [protocol](FINITE_REVIEW_PROTOCOL.md) fixes two cells, one acquired field per
cell, and two claims A0*B1 and (A0*B1+A1*B0)/2. Each field costs 1.1 and reveals
both recorded marker calls. At most one recorded call per marker differs from
its true binary label. Reviews reveal a response for one acquired cell/marker,
at cost .05 or .5, with up to one or three responses. No assay error rates or
actual cost units are inferred here.

The belief set enumerates all 16 binary true populations satisfying observed
raw-call and review-error budgets. Bounds are exact extrema over this set. The
adversary selects any next response that leaves at least one feasible population.
The exact Bellman value is

    V(s) = 0                                      if all target widths are zero
    V(s) = min_a [cost(a) + max_y V(update(s,a,y))] otherwise.

No remaining actions with unresolved targets gives infinity. The finite solver
uses feasible-policy upper bounds for pruning. An independent unpruned Bellman
test agrees in a smaller complete case. Controls greedily maximize immediate
worst-case width decrease per cost, or first acquire every field and then apply
that review rule. These are exact two-cell planning values, not scalable policy
performance or patient-level results.

## Costs and failures

Each entry is worst-case total hypothetical acquisition plus review cost.
Infinity means the required exact readout cannot be guaranteed.

| Review model | Review cost | Greedy | Census first | Optimal |
|---|---:|---:|---:|---:|
| Perfect, one response | .05 | 2.40 | 2.40 | 2.40 |
| One wrong response allowed per marker, one response/cell | .05 | infinity | infinity | infinity |
| One wrong response allowed per marker, up to three/cell | .05 | 2.80 | 2.75 | 2.70 |
| One systematically wrong cell allowed per marker | .05 | infinity | infinity | infinity |
| Systematic errors plus trusted reference at cost 1 | .05 | 6.35 | 6.40 | 6.20 |
| Perfect, one response | .50 | 4.20 | 4.20 | 4.20 |
| One wrong response allowed per marker, one response/cell | .50 | infinity | infinity | infinity |
| One wrong response allowed per marker, up to three/cell | .50 | 8.20 | 7.70 | 7.20 |
| One systematically wrong cell allowed per marker | .50 | infinity | infinity | infinity |
| Systematic errors plus trusted reference at cost 1 | .50 | 7.70 | 8.20 | 6.20 |

The event-budget optimum is 1.8%-6.5% below the stronger census-first control;
this is a small oracle gap under the specified model, not demonstrated savings.
The repeated-event problem explores 44,738 cached optimal-value states in each
cost configuration. A cell-error repeat is dominated in worst-case planning:
the adversary may repeat its previous response without changing the feasible
set, so spending again cannot improve the minimax value. Such actions are pruned;
this is not an assumption that different responses never occur. Full results
are in [`finite_review_study.json`](results/finite_review_study.json).

## Constructive systematic-error impossibility

Suppose every acquired call and every review response is one. Both following
true worlds remain feasible when errors are counted by wrong cells:

- All four labels are one: targets (1,1).
- A0 is zero, all other labels are one: targets (0,.5).

The second world has one wrong recorded A call and one wrongly reviewed A cell.
Every repeated A0 response can be wrong without exceeding the **cell** budget.
Any adaptive policy can therefore receive the same observation transcript in
both worlds, regardless of repetition count. It cannot certify narrower than
width one for the first claim. This is an information obstruction, independent
of computation or policy design. The witness also defeats one-response reviews
with an event budget of one; it does not defeat repeated reviews under that
budget, since each additional wrong response consumes another event.

This elementary indistinguishability argument is not a priority claim for a new
general inference theorem. It identifies an unjustified observation assumption
in the earlier ideal-review emulator. Repetition is useful only under a defensible
response-error model or additional evidence that distinguishes these worlds.

The trusted reference diagnostic was added after this impossibility result and
frozen before its evaluation. It directly reveals one true cell/marker label at
cost one. Four such reads plus two fields guarantee the 6.2 cost optimum here;
ordinary systematic reviews cannot beat it in worst case. Real references may
also be fallible, in which case their own error model must enter the belief set.

## Research consequence

Keep acquisition/review policies conditional on explicit observation assumptions.
The lab pilot should compare separate evidence sources and adjudication errors,
not count repeated classifier outputs as independent truth. Next test a reference
with a nonzero bounded error model, derive the remaining attainable precision,
and combine it with [coupled-marker bounds](COUPLED_MARKER_RESULTS.md). No clinical
validation, independent-cohort result or measured wall-clock saving is claimed.

```bash
python -m unittest discover -s tests -v
python study_finite_review.py
```
