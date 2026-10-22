# Frozen development audit: does a valid cell ontology resolve ambiguity?

The enrichment confirmation and its later refinement have already been read.
This audit reuses those same six patients and is explicitly development work.
It neither replaces their failed frozen gate nor introduces a fresh cohort.

Question: are opposite enrichment signs still feasible when each cell must be
exactly one of other, epithelial, or T cell? Independent binary-marker flips
previously permitted an epithelial/T-cell overlap. This audit checks whether that
modeling freedom was responsible for the reported ambiguity.

This categorical assumption concerns mutually exclusive *cell identities*, not
universal absence of protein coexpression or doublets. Verify that the source
mapping contains no overlap; otherwise stop. Keep the fixed measured kNN12 graph,
same cell order and same separate A/B Hamming budgets ceil(.005*n). No new error
budget calibration or geometry guarantee is claimed.

Construct feasible minimization/maximization witnesses from the raw labels and
the two existing witnesses after dropping any newly added overlapping marker.
Use deterministic best single-cell categorical moves, at most 50 moves per start,
with objective changes evaluated in integer arithmetic. Never call this a global
optimizer; a failed search is inconclusive. Save every correction index and check
all signs independently using integer contact counts and degree sums.

For any feasible opposite-sign pair x,y, a query reveals the full categorical
identity of a selected cell. If at most E reviewed cells have arbitrary incorrect
identities, a fixed query set Q can distinguish the pair for every answer vector
iff Hamming(x_Q,y_Q) > 2E. This is the standard disjoint-Hamming-balls criterion,
not a new coding/active-learning principle. Repeated reading of the same cell
with systematic error does not increase its distance. A shared transcript on an
adaptive path also cannot justify stopping while its restriction remains within
E of both worlds. The result is conditional on supplied E and label feasibility.

For E=0 and 1, report differing cells and check whether all available cells even
permit distinguishing the pair. Partition coordinates into the previous 4x4
normalized review fields, and solve the finite pair's cheapest unit-cost field
cover by sorting disagreement counts. This is an exact *pairwise lower bound*
on a worst-case nonadaptive sign-certifying acquisition plan, not a sufficient
plan for all feasible labelings, and not measured microscope cost. Do not infer
that this lower bound applies unchanged to every adaptive branch.

Validation: independent brute-force small categorical examples for objective
increments; exhaustive three-symbol answer vectors to check the 2E criterion;
exact small-field subset enumeration for the field bound; saved real-data
witnesses verified separately from the search implementation. Original frozen
binary bounds remain valid supersets; no new sign certificate is asserted.

Closest foundations: Wiener, Hanneke, El-Yaniv (2015), version-space compression,
https://jmlr.org/papers/v16/wiener15a.html ; classical error-correcting-code
minimum-distance decoding. The proposed research contribution would have to be
a useful joint measurement/acquisition guarantee on independently audited data,
not relabeling these standard foundations as a new learning theory.
