# Nature-inspired direction: conserve molecular observations before cell assignment

Date: 2026-10-29. [Protocol](CONSERVED_MASS_PROTOCOL.md) `1215d93` preceded the
screen. **The observable passes a synthetic feasibility gate; there is no new
algorithmic result.** This changes the scientific target from cell-type neighbor
contrast to cross-channel molecular proximity.

## Why this direction follows the failure

Cell segmentation can reassign molecules without changing their original
locations. A statistic computed directly from those locations is invariant to
that reassignment. It still depends on correctly detected molecular identities,
localization, field selection and missing signals. Conservation is an explicit
measurement assumption; it does not mean that all assay errors conserve mass.

For channels A and B, the target is

    S_h = mean_{i in A,j in B} exp(-||x_i-x_j||^2/(2*h^2)).

This familiar kernel statistic is raw proximity, not causal communication,
between-cell interaction, a cell-type label, or a calibrated enrichment test.
Within-cell and between-cell proximity are both counted. Increasing bandwidth
changes the target and can obscure biologically meaningful structure.

If each detected location is displaced by at most rho, every pair distance lies
between max(0,d-2rho) and d+2rho. Kernel monotonicity yields safe pairwise extrema;
averaging them is an outer interval even though a common configuration need not
attain all extrema simultaneously. The conventional global Lipschitz control
has error 2rho*exp(-1/2)/h. Neither construction is claimed as novel.

If at most one molecule per channel was missed among n observed, the difference
in each normalized channel measure has total variation at most 1/(n+1). Since
the kernel lies in [0,1], the extra score error is at most 2/(n+1). This requires
an actual bound on missing mass. Unknown or selective under-detection without a
bound is outside the certificate. Wrong gene calls also need separate treatment.

## Complete screen

- All six exact cases contain all 729 discrete location assignments each,
  checked by separate scalar-pair enumeration (4,374 evaluated configurations).
- Conventional pair-distance bounds distinguish shared and separated clusters
  in **12/12** scaling pairs; the gate required 9/12. The global Lipschitz
  control already distinguishes **9/12**. This is evidence for an observable's
  stability in this easy regime, not evidence of a new method outperforming
  established methods.
- The 12/12 result persists when allowing one missing molecule per channel.
- A declared out-of-budget example fails as expected: moving separated channels
  20 units into coincidence produces proximity one, outside their rho=1 interval.
- Permuting or regrouping the same molecules among cell assignments, and applying
  common rigid translations, leaves the score unchanged. Relative movement of
  the channels does change it; the statistic is not trivially constant.

All examples are synthetic and have intentionally large differences in cluster
location. Interval separation is not a sampling p-value, a test of subtle cell
communication, or independent biological confirmation. The old cell-contrast
failure remains unchanged because this is a different estimand.

## Raw-coordinate data feasibility

The [SSAM-lite public sample archive](https://zenodo.org/records/5517607) was
downloaded and its published MD5 verified. Size is 26,103,043 bytes; SHA-256 is
`da4862b4e11d860f081b8e74b90e2119f8b6eaf9bf35869f515da3be95cda17c`.
Only readmes, archive metadata and CSV headers were inspected. No marker-pair
scores or biological outcomes were examined.

The readmes describe mouse cortex osmFISH (1,802,589 spots, 33 genes) and human
pancreas ISS (461,078 spots, 138 genes). The coordinate table headers actually
contain gene identity and X/Y. These counts are source declarations, not counts
recomputed here. Source readmes do not establish physical units, sample/patient
replication, localization error or detection completeness. Those checks remain
mandatory before choosing physical bandwidths or claiming a held-out cohort.

The [feasibility manifest](results/molecular_data_feasibility.json) records
provenance and the actual headers. Raw data remain local and ignored by version
control. These are prepared teaching datasets from published studies, not
independent adjudicated ground truth.

## Selected research question and next gates

**Which molecular association conclusions remain valid across physically
justified localization and detection errors, and at what spatial scales?**
Nature supplies a conservation constraint and a scale hierarchy, not a novel
statistic. A potential contribution would be useful simultaneous uncertainty
certificates across marker pairs/scales, under a calibrated error model, with
an acquisition decision that beats strong simple controls at equal cost.

Before that claim: verify coordinate units and sample identities from original
releases; freeze a marker/scale family and spatial held-out design; account for
inhomogeneous tissue density and within-cell signal; add detection and gene-call
errors; calibrate against an independent reference; compare conventional kernel
analysis and existing correction pipelines on the same target. A multiscale
search must control selection instead of reporting only favorable bandwidths.
Do not equate pixels, technical fields, and independent patients.

Prior art directly limits broad novelty claims:

- [SSAM](https://www.nature.com/articles/s41467-021-23807-4) already performs
  segmentation-free analysis using molecular density fields.
- [spage2vec](https://pmc.ncbi.nlm.nih.gov/articles/PMC7983892/) already learns
  local segmentation-free molecular representations.
- [cellAdmix](https://www.nature.com/articles/s41588-025-02497-4) addresses
  misassigned molecular mixtures and downstream effects; its
  [implementation](https://github.com/kharchenkolab/cellAdmix) is a relevant control.

Reproduce with `python study_conserved_spatial_mass.py --verify` and the full
test suite. [All screen outcomes](results/conserved_spatial_mass.json) are saved.
The result justifies further research into calibrated molecular observables,
not a manuscript or a higher novelty rating today.
