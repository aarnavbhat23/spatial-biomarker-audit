# Measured graph-error calibration — frozen 2026-10-21

## Scope

Replace the privileged fixed-graph assumption with explicit edge deletion and
addition budgets. First establish whether those budgets can be calibrated on
manual annotations without becoming useless. This is a nuclear geometry study,
not validation of tumour marker accuracy, registration across modalities,
patient transport, or scanner feedback latency.

BBBC039v1 (Caicedo et al. 2018; Broad Bioimage Benchmark Collection,
https://bbbc.broadinstitute.org/BBBC039) releases 200 DNA-channel fluorescence
images, manual masks, and official train/validation/test partitions (100/50/50).
Images represent chemical perturbations from one screen, not independent human
patients. CC0; retain citation. Use official splits, not a random cell split.
Do not inspect validation/test outcomes before selecting the segmenter.

## Frozen measurement procedure

Decode source masks by connected components of the first PNG channel, following
the authors' published example. Use 8-connectivity as in skimage label default.
Prediction uses only the input image: Gaussian smoothing sigma=1, Otsu threshold,
remove objects smaller than 30 pixels, fill holes, Euclidean distance transform,
local maxima markers and watershed. Train-only grid: Otsu multiplier .8/1/1.2
and peak minimum distance 5/8/12, nine settings. Select best mean image-level
instance F1 under one-to-one IoU>=.5 matching on the 100 train images. Ties by
listed parameter order. No trained public weights or undocumented tuning.

Match predicted and reference nuclei by maximum-cardinality bipartite matching
of pairs with IoU>=.5 (deterministic order). Build undirected radius graphs at
40 pixels on centroids. Match IDs only for offline scoring; unmatched predicted
and reference nodes receive disjoint IDs. Thus false/missed/split/merged nuclei
and centroid errors contribute to deleted/added graph edges. Both graph vertices
and edges may be wrong. The policy never receives the reference graph.

For each image m=max(1, number of predicted edges), score=max(deleted,added)/m.
Calibrate the 95% split-conformal upper quantile on 50 validation images, rank
ceil((50+1)*.95)=49, with infinity if that rank is unavailable. On each test
image use D=A=ceil(q*m); the exact finite-sample statement is marginal under
exchangeability of calibration and future images with the frozen pipeline.
It is not a distribution-free patient/domain-shift guarantee or a confidence
statement about the realized 50-test coverage. Report 50 test scores, coverage,
exact binomial interval (descriptive only, if independent images), counts,
and a naive zero-error comparator. Official partitions may be nonrandom;
coverage is empirical, not proof of exchangeability.

## Functional bound and success gate

For m observed edges and marked-edge count s, at most D edges removed and A
added, the enclosing pair-frequency interval is
[max(0,s-D)/(m-min(D,s)+A), (s+A)/(m-min(D,m-s)+A)].
Clip to [0,1] and use [0,1] if an allowed empty target graph makes the fraction
undefined. For uncertain labels substitute the existing count bounds L,U in
the lower/upper formulas. Matching must preserve the marks of matched cells;
marker disagreements require the separate label uncertainty set. No biological
marker assignments exist here; evaluate the width at a hypothetical s=.25m
only as a precision diagnostic, never as an observed biomarker.

Frozen gate: at least 47/50 test images covered AND median diagnostic interval
width <=.10. Count manual calibration cost explicitly: all validation fields,
annotation count and pixels. Wall-clock computation is not acquisition time.
Failure means do not add useful measured-graph robustness to the headline;
retain the bound, calibration mechanism, negative result, and plan better
segmentation plus a paired H&E/multiplex annotation pilot.

## Prior-art difference

STARLING (Lee et al., 2025), https://doi.org/10.1038/s41467-024-55214-w , already
models segmentation error in phenotyping. Its released pellet annotations
include model-derived labels and mixture controls, not per-cell adjudicated
truth; do not use them to claim calibrated marker errors.
Bruhns et al. (2025), https://doi.org/10.1371/journal.pcbi.1013350 , already
quantifies downstream segmentation effects; its baseline mask is a curated
proxy for truth and released table is one patient.
Lu et al. (2026), https://proceedings.mlr.press/v337/lu26b.html , already constructs
conformal instance-segmentation sets. Split conformal and graph perturbation
bounds are standard ideas. This study's target is propagating deletion/addition
budgets into this project's marked-edge estimand, with explicit denominators;
no first-ever uncertainty or segmentation claim.
