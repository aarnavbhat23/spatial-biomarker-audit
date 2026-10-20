# Independent IMC cohort: frozen graph-bound confirmation

Source: Danenberg et al., Nature Genetics 2022, public processed release
10.5281/zenodo.7324285 (version 2.1.0, CC-BY-4.0). This cohort is independent of
the ORION colorectal development data. A CSV header, annotation dictionary and
initial-chunk phenotype vocabulary were inspected to establish eligibility and
mapping; no patient targets, interval widths or comparative outcomes were read.
The labels remain source phenotypes derived from measured channels, not
independently adjudicated biological truth. This tests tissue/platform transfer
of conditional graph sensitivity, not protein prediction or clinical validity.

Before target evaluation, freeze these rules:
- Read only patient/image/cell identifiers, QC flags, phenotype, epithelial flag
  and cell-centre coordinates. Retain tumour rows with is_hotAggregate=0,
  finite coordinates and nonmissing patient, phenotype and epithelial flag.
- Exclude the ambiguous source class `CD4^{+} T cells & APCs` entirely.
- A is the supplied is_epithelial flag. B is membership in exactly
  `CD4^{+} T cells`, `CD8^{+} T cells`, `T_{Reg} & T_{Ex}`. Other unambiguous source
  classes have B=0. This is an annotated epithelial-to-selected-T-phenotype
  claim, not the PanCK/CD3 positivity claim or every biological T cell.
- Sort metabric_id lexically; select the first six distinct patients with an
  eligible image containing >=500 retained cells. For each select its lowest
  eligible ImageNumber. No substitutions based on targets or outcomes.
- Within each image, reject duplicate (ImageNumber,ObjectNumber) identifiers
  and duplicate coordinates. Never connect different images or patients.
- Construct directed k-nearest-neighbour graphs for k=6/12/24 with deterministic
  distance/index tie-breaking. These topological proximity claims avoid assuming
  a physical calibration across platforms. They differ from radius-based claims.
- Apply the unchanged general fast-bound algorithm from commit 9d8c6da, both
  budgets ceil(n*.001) and ceil(n*.005), yielding 36 cases. Compare valid degree
  bounds and attained corruption witnesses. Keep every case, including absence
  of either class. No new tuning of witness starts/iterations or graph k.

Pass tightness transfer only if >=27/36 cases strictly tighten the degree outer
interval and >=27/36 have maximum certified endpoint gap <=.001 (.1 percentage
points). Report cases where width exceeds .005 and where attainable witnesses
prove it impossible. A pass concerns conditional algorithmic tightness on six
patients; it does not validate assay-error budgets, active acquisition savings,
independent annotation accuracy or new general optimization theory. Record source
checksums, retained/excluded counts, graph hashes and patient-level outcomes.
