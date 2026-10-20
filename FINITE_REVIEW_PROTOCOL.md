# Exact planning with two uncertain markers and fallible reviews

Tiny-population development/theory study, not biological validation. Two cells,
two directed-pair claims: A0*B1 and (A0*B1+A1*B0)/2. Two fields, one per cell.
Each field reveals both recorded calls. At most one recorded A call and one
recorded B call can be wrong. Review responses have externally assumed global
error budgets per marker over all review events. No error rate is estimated.

Compare perfect review with one review per cell/marker, fallible review with
one response per cell/marker, and fallible review with up to three responses
per cell/marker. Fallible budget is one per marker across all events. New field
cost 1.1; each review cost .05 or .5. These are hypothetical units. The target
width is zero (exact finite-population readout). Enumerate every binary true
population, retain all feasible completions, and solve the complete worst-case
finite decision problem by dynamic programming when computationally feasible.

Controls: greedy immediate width decrease per cost, and acquire-all before that
review rule. Report minimax costs, infeasible cases, solved state counts and
compute. Verify an exact optimum against direct Bellman evaluation on smaller
cases. Any reduction is only an oracle gap on this tiny problem, not a scalable
policy or an instrument-saving claim. A one-review failure must remain visible.

Theory-driven assumption check added before its outcomes: also count the review
budget by wrong cells rather than wrong response events, allowing the same cell
to be wrong on every repetition. Keep this separate from the original event-budget
experiment. It tests the observation assumption, not a newly tuned policy.

After the repetition impossibility result, a separate theory-driven diagnostic
adds a trusted exact reference action, cost 1.0 per cell/marker, to the systematic
cell-error model at both review costs. Freeze before evaluating this diagnostic.
This action is an ideal reference assumption, not a demonstrated perfect assay.
Its role is to distinguish missing information from a deficient search policy.
