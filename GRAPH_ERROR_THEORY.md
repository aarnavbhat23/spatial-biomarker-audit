# Propagating graph errors into a marked-pair fraction

Let an observed graph contain m edges and an unknown compatible labeling give
marked-edge total s in [L,U]. A corrected graph deletes at most D observed edges
and adds at most A new edges. Removed and added edge marks lie in [0,1]; shared
edges keep their marks. Correspondence is an assumption, not information handed
to the acquisition policy. Unmatched nuclei may contribute whole sets of added
or deleted edges. Separately wrong marker labels require the existing feasible
label set, not just a graph budget.

For m>D and any compatible nonempty corrected graph,

    lower = max(0,L-D) / (m-min(D,L)+A)
    upper = (U+A) / (m-min(D,m-U)+A).

Return [0,1] if m=0 or D>=m, because the uncertainty set then allows an empty
graph and its marked-edge fraction is undefined. Counts L,U may be fractional
(e.g. bounds on weights in [0,1]); D,A and m are nonnegative integers.

Proof: write the corrected total as s-r+t and its denominator as m-d+a,
where 0<=r<=min(d,s), 0<=t<=a, d<=D and a<=A. For the lower endpoint choose t=0
and add as many zero-mark edges as possible. Removing marked mass lowers the
fraction until the total is zero; removing zero-mark mass cannot lower it.
Hence delete min(D,s) marked mass and add A zero-mark mass. The resulting
lower expression increases with s, so substitute L. Symmetrically, for the
upper endpoint remove min(D,m-s) zero-mark mass and add A fully marked mass.
Its expression also increases with s, so substitute U. For binary integral
marks these endpoints are attainable in the unconstrained edge-edit relaxation.
Actual radius graphs and products of vertex labels impose further constraints;
this is an enclosure, not a sharp extremum for those graphs.

For the earlier directed graph, use directed edge counts and directed deletion/
addition budgets. An undirected edge can represent the average of its two
oriented marker products, a weight in [0,1]; use a consistent normalization.
Never combine an undirected error budget with a directed denominator unchanged.

## Coverage and adaptive review

On calibration image i use score_i=max(D_i,A_i)/max(1,m_i), from a frozen image-only
segmenter compared with manual annotation. With n calibration scores and rank
k=ceil((n+1)(1-alpha)), q is the kth ordered score (infinity if k>n). For a new
image use D=A=ceil(q*max(1,m)). Under exchangeability, standard split conformal
covers both edge budgets together with marginal probability at least 1-alpha.
The correspondence/matching procedure is fixed before calibration. This is an
image-level claim; cells and edges are not treated as independent samples.

Once an image's graph enclosure is valid, it encloses every labeling on the
shared vertices. Consequently it can be composed with simultaneous label bounds
at every adaptive review time without an additional graph union bound over
review times. If the label event fails with probability alpha_l and the graph
event with alpha_g, their intersection has probability at least
1-alpha_l-alpha_g, requiring no independence. Using 95% for each source gives
at most a 90% lower guarantee for their intersection, not 95%. Changing the
segmentation or graph after review is not covered by the original static-graph
event and requires a new joint argument/calibration design.

Conformal calibration does not establish exchangeability. BBBC039 official
splits share a screen and plate identifiers; they are not independent patients.
Nor does calibration imply useful precision: a high quantile can allow deletion
of every edge. The measured study records exactly that failure.

## Checks and precedents

Tests enumerate all binary marked-edge worlds on up to five observed edges,
all allowed deletion subsets and up to two additions. They also check composition
with label-count intervals, empty graphs, finite-sample quantile indexing,
split/merged nuclei and births/deaths in the offline correspondence.

The algebra is a basic edge-edit relaxation and split conformal is established
statistics. STARLING, downstream segmentation-error benchmarks and conformal
instance-segmentation sets precede this work (references in the protocol).
The contribution here is a precise integration and falsifiable empirical test
for this project's estimand, not a first-ever uncertainty method.
