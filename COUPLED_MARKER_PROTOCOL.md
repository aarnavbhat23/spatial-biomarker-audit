# Coupled marker sensitivity: frozen development screen

Evaluate a derived sharp full-census bound for at most one wrong A call and
k wrong B calls on a fixed simple directed graph with unit edge weights. The
small-A-budget assumption is deliberately restrictive; it is not an estimated
assay error rate. Source calls are the reference, not biological truth.

Before reading these new outcomes, fix development patients CRC02/CRC11/CRC19,
radii 40/60/80 native pixels and B error fractions 0/.001/.005 (ceil(n*fraction)).
A budget is always one. Report attainable lower/upper bounds, width, wall time,
and the valid label-blind degree outer bound. Check simple-graph uniqueness for
every panel/radius, and record prepared-input SHA256 hashes. No parameter fitting,
patient selection or independent confirmation is part of this screen.

Validation: exhaust all recorded and true binary populations on a three-cell
example, including every B budget 0..3; test the joint-error counterexample in
which separate observed-label sensitivities both vanish. This establishes the
special-case algorithm, not priority over all combinatorial optimization work.
