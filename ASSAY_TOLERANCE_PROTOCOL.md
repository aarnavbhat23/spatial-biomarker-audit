# Assay-error tolerance: exploratory decision analysis

Added after the general-budget screen passed, before computing these outcomes.
On each development graph, find a bracket for the largest common per-marker
error budget compatible with width .005. Start at budget zero and use the
already failed ceil(n*.001) budget as an upper endpoint. At each midpoint,
accept only if the valid outer interval is narrow enough; reject only if two
attainable corruption witnesses are farther apart than .005. Stop with an
unresolved bracket if neither certificate decides. Never label a heuristic
witness as the exact optimum. Do not extrapolate these annotation-relative
requirements to actual assay error rates or whole-slide clinical reliability.

This analysis translates conditional precision into a requirement for a future
lab pilot. It is not a new-patient confirmation or error-rate estimation study.

Refinement amendment after the first bracket outcomes: binary search stopped
with an uninformative [0,56] bracket on one graph because its first midpoint was
unresolved. Before refined outcomes, fix an exhaustive scan of every integer
budget 0..ceil(n*.001) on all nine graphs. Retain original brackets unchanged.
Report the largest sufficiently narrow outer interval and smallest too-wide
attainable witness interval; all budgets below the former are sufficient and all
at/above the latter impossible by nesting, even if a heuristic score is not
monotone. This is deterministic optimization refinement, not fresh confirmation.
