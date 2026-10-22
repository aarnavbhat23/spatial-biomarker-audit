# Contrast-specific robustness screen — frozen 2026-10-21

Question: can jointly bounding an abundance-adjusted spatial contrast retain
biologically relevant sign information under supplied marker-error budgets
more efficiently than separate bounds on contacts and abundances?

This changes the estimand. It does not rescue the original contact fraction,
resolve BBBC039 graph failures, calibrate actual marker error, or demonstrate
H&E-only deployment. Graph and source phenotype annotations remain measured
proxies. No detector is tuned using already-evaluated BBBC039 images.

## Estimand and algorithm

For a fixed simple directed graph W, m>0, outgoing degrees d and incoming e,
set C=W/m - d*e'/m^2. Define theta=a'C b, the covariance of endpoint marker
bits under a uniformly drawn directed edge. Theta lies in [-.25,.25]. The
abundance correction is degree weighted; this is an established covariance/
modularity-style contrast, not a new biological interaction statistic.

Under <=B_A and <=B_B arbitrary binary call flips, expand about raw a,b:
Delta theta = (Cb)'Delta a + (C'a)'Delta b + Delta a'C Delta b.
Bound signed linear terms exactly by sorting signed flip gains, including the
zero-flip choice. Bound the interaction by the minimum of:
- min(B_A*B_B, top-B_A sum(d), top-B_B sum(e))/m
  + top-B_A sum(d)*top-B_B sum(e)/m^2;
- sqrt(B_A*B_B*||C||_1*||C||_infinity).
Compute absolute row/column sums using sparse edges plus the rank-one term,
without materializing a dense matrix. Include alternating feasible witnesses
(two starts, at most 20 passes); they do not establish global optimality.

Control: unchanged fast contact bounds plus exact degree-weighted abundance
bounds, propagated separately into contact minus abundance product. Primary
candidate is its intersection with the direct contrast bound. Record direct
unintersected intervals too. Tests must compare every two-marker binary world
on small graphs against the enclosure, including isolated vertices, empty
marker populations, complete directed graphs and zero budgets. Show witnesses
satisfy the budgets. Centering must not silently fix the true prevalences.

## Fresh patient screen

Use the existing Danenberg CSV and its pinned hash, QC and epithelial/selected-
T-cell phenotype mapping. Select the 13th–18th lexicographic eligible patients,
one image with >=500 retained cells per patient, exactly as earlier selection
rules except for skipping the first 12. Those patients have not been used in
this project's first/second six-patient studies. No selection by label prevalence,
contrast sign or bound width. Keep empty-marker cases and all failures.

Primary: directed kNN12 and B_A=B_B=ceil(.005*n). Secondary: kNN6/24 and
fractions .001/.005, all 36 combinations retained. Compare widths of bounds on
THE SAME enrichment estimand, not enrichment against raw-contact widths.
Success requires >=25% median width reduction against the separated control
and sign exclusion of zero in at least four of six primary patients. These are
conditional sensitivity results; there are no invented annotation errors or
claims that source phenotypes are biological truth. A positive result motivates
joint graph/label certification on a new annotated cohort, not a 9/10 grade.

## Prior art

Newman-style degree correction, network assortativity and modularity are known.
Binary metadata mixing constraints already have dedicated bounds:
https://arxiv.org/abs/1908.04588 . Robustness to graph/attribute perturbations
has substantial precedent:
https://papers.neurips.cc/paper_files/paper/2019/hash/e2f374c3418c50bc30d67d5f7454a5b4-Abstract.html .
Spatial K-statistics already address robustness and abundance effects:
https://arxiv.org/abs/2412.08498 . The possible contribution here is a sparse
conditional certificate that jointly propagates both marker-error budgets into
this contrast, exposing cancellations lost by independent interval propagation.
A broader priority review is required before a novelty claim.
