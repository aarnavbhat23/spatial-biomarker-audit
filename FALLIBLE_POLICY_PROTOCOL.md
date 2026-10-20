# Fallible two-marker review — frozen before new patient outcomes

Keep the general raw-label algorithm at 9d8c6da unchanged. Develop a separate
conditional bound using raw labels with <=B errors per marker and partial
reference responses with <=E systematic wrong cells per marker. Reviews are
not truth; repeated responses from the same cell cannot spend additional
independent error allowances. Source annotations remain proxy truth for emulation.

Derive an exact signed linear optimization oracle over the intersection of the
two Hamming constraints. With flips f relative to raw calls, partition cells
into reviewed disagreements D, reviewed agreements G and unreviewed U. Feasible
counts satisfy d+g+u<=B and |D|-d+g<=E. Sorting weights within groups and scanning
counts gives the exact linear optimum, including negative weights and forced
flips. This is classical cardinality optimization, not a new general optimizer.

Use the oracle to bound the two linear terms in a bilinear graph expansion
around the hybrid raw/reference labeling. Bound the signed interaction term by
feasible changed-cell counts and incident edges; construct feasible alternating
optimization witnesses. Validate bounds, oracle solutions, witnesses, empty-set
handling and reference floors by exhaustive enumeration on small graphs before
patient evaluation. Keep unresolved endpoint gaps visible. A conditional bound
is required for every transcript, including adaptively selected reviews.

Fresh confirmation data: the same Danenberg 2.1.0 source/mapping/QC as dd3041a,
but the seventh through twelfth eligible patients lexically, each with its
lowest eligible image and >=500 retained cells. The first six patients were
used for bound transfer and must not be reused as fresh policy confirmation.
No new target outcomes have been inspected. Use directed kNN k=12 and target
pair-fraction width .005. Partition each image's coordinate bounding box into
a fixed 4x4 grid of emulated fields; this is not calibrated scanner geometry.

For each patient use two seeds and two raw-corruption profiles: uniform random
flips or flips at highest source-target single-marker influence (an adversarial
generator unavailable to the policy). Both markers have exactly ceil(.005*n)
raw flips. The reference has one persistent wrong cell per marker, sampled
independently and fixed before the policy; E=1 is a hypothetical bound, not a
measured reference accuracy. All recorded calls and graph coordinates are known
initially; a field review reveals both reference markers for its cells. The
policy sees only raw calls, coordinates, budgets and already acquired responses.
This is review after raw-call census, not an end-to-end acquisition experiment.

Compare seeded uniform fields, degree-prioritized fields, static raw-target
influence, a bound-width lookahead control, and an endpoint-witness disagreement
policy. Degree/influence scores sum only over unreviewed cells. Lookahead ranks
fields by the worst future bound among three feasible-label proxy responses:
current feasible anchor and the two endpoint witnesses. These proxies guide
actions only, not validity; they are not a complete minimax response set.
Witness disagreement prioritizes fields exposing different labels in the two
attained endpoint populations, weighted by their affected claim edges. Both
adaptive rules divide their scores by quoted field-plus-cell cost. Ties break
by field index. Intersect bounds over history, stop on .005 width or full review,
and retain all failures. No policy may inspect unqueried reference responses.

Costs are hypothetical: one per field, .02 per reviewed cell (both markers),
plus feedback latency 0/1/8 per sequential field. Report fields, cells, rounds,
costs, validity, widths, endpoint gaps and time separately. The primary latency
is one. A policy continuation gate requires at least 10% mean cost reduction
over the strongest frozen control at equal attainment, with improvement in both
raw-error profiles. Compare the witness policy against uniform, degree, static
influence and lookahead; a failed gate stops any review-saving claim. Biological
units are six patients, not the 24 corruption scenarios or policy repetitions.
Physical savings and error calibration still require the lab pilot.
