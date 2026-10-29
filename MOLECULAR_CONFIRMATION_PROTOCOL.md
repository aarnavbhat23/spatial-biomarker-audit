# Molecular association under bounded error: public-data protocol

Freeze 2026-10-29 before downloading/scoring the new coordinate table. Use the
filtered mouse VISp smFISH table linked by the original SSAM documentation:
https://s3.amazonaws.com/starfish.data.spacetx/spacetx-website/data/smFISH_Allen/s3_spot_table.csv
The documentation states that rotated_x/rotated_y are already micrometres:
https://ssam.readthedocs.io/en/develop/userguide/02-data.html
Save source checksum, schema and all selection counts. This is one prepared
public field, not a new independent patient cohort or assay-error reference.

Primary six genes: Gad2, Slc17a7, Pvalb, Sst, Vip, Cux2. Keep every missing gene
or sparse comparison as unevaluable; do not replace markers after seeing data.
Use all 15 pairs and bandwidths 5,10,20 micrometres. Split the field at the
midpoint of its observed x range. Left is development, right confirmation;
both are analyzed by the same frozen code, with no tuning on the left outcomes.
These regions are technical partitions, not independent subjects.

Within each half and gene, uniformly sample at most 200 measured molecules,
seed 20261031; retain comparisons with >=20 in both genes. Target only this
selected measured-molecule population. No assertion about unobserved molecules
or all cells follows. Conditional permutation null: labels are exchangeable
within fixed 2x2 coordinate strata in the selected pair, conditional on positions
and stratum counts. This adjusts coarse density but is not a universal null of
no interaction or a correction for within-cell coexpression.

For each pair, generate 1999 independent within-stratum label permutations,
preserving both channel counts. Reuse permutations across all three scales.
One-sided raw-kernel proximity statistic. Nominal Monte Carlo p=(1+#scores>=S0)/2000.
Robust upper p uses observed statistic lower bound and each permutation's upper
bound under location displacement rho and up to k A/B label swaps in the same
fixed strata. Pair-distance bounds from the prior study are conventional.
The swap allowance changes any bounded pair mean by <=k/nA+k/nB. Thus upper p
dominates each admissible latent-data permutation p under the stated null.

Evaluate (rho,k)=(0,0),(.5,0),(.5,1),(1,1). These are sensitivity assumptions,
not calibrated error rates. For the right-half primary family use Bonferroni45
at alpha=.05 even if some genes/pairs are unavailable. Also report left results
with the same correction, clearly labeled development. All scales/results stay.
No selection of a favorable radius, pair or direction after outcomes.

Continue a robust-association claim only if the confirmation region has at
least one adjusted rejection at (.5,1), and at least half the nominal adjusted
rejections survive (.5,0). If no nominal rejections exist the gate fails. This
is a usefulness screen for the stated exchangeability/sensitivity contract,
not new biology, causal communication, independent coverage or general novelty.
Missing molecules and changed channel membership are not silently absorbed
into the permutation theorem; those require a different observation model.

Verify the robust p bound by complete enumeration on small examples. Compare
global versus stratified nominal permutation results on the same selected
molecules, without using the global null to claim density-adjusted association.
