# General coupled errors and fallible references: frozen development screen

Extend the one-A-error case to arbitrary A/B budgets. Use direct binary label
variables with per-edge product constraints and marker disagreement budgets.
Compare an established MILP formulation/HiGHS solver, its continuous relaxation,
and a fast graph-specific outer bound with explicit attainable corruption
witnesses. These ingredients are classical; test a useful restricted application,
not priority for integer optimization. Report unresolved primal/dual gaps.

Before outcomes: random geometric graphs of sizes 12/24/48, seeds 0/1/2,
radius .35, directed simple edges, random recorded A/B calls of probabilities
.3/.4. Budgets ceil(n*.1) and ceil(n*.25), both markers. Limit each solver
objective to five seconds. Exhaustively check all populations on smaller graphs.
Report bounds, attainable extrema, gaps, LP looseness and time; never equate an
LP optimum or a timed-out incumbent with an exact integer optimum.

Development panels CRC02/11/19, radii 40/60/80, both marker budgets
ceil(n*.001) and ceil(n*.005). Run fast outer bounds and feasible witnesses only;
no full-panel exactness claim. Compare the valid label-blind degree bound and
report certified endpoint gaps. Source labels are reference annotations.

Fallible-reference diagnostic: n=8 geometric graphs, seeds 0/1/2, raw budget two
errors per marker, source labels as reference truth for constructing the test.
Flip the first two cell labels per marker in raw calls. Review the first four
or all eight cells/marker with correct reference responses, while conservatively
allowing one response to be wrong per marker. Compare no reference and an ideal
zero-error reference control. This tests attainable precision conditional on
budget assumptions, not estimated assay errors or independent biological data.

Continue the algorithmic direction only if the fast bounds are valid by proof
and exhaustive tests, strictly tighten the degree outer interval in at least
half of full-panel cases, and attain <=.1 percentage-point endpoint gap in at
least half. Preserve failures. Even a pass requires new-patient confirmation,
calibrated error budgets and measured acquisition/reference costs.
