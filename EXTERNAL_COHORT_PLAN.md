# Independent-cohort feasibility and frozen analysis requirements

Status as of 2026-10-19: metadata inspected; no external target outcomes evaluated.
No independent-cohort success is claimed. Public download attempts returned
HTTP 401 from the API and 403 from the listed file links in this environment.
Resolve legitimate access before treating the dataset as executable validation.

Candidate: [Hickey's multiplexed imaging release, used by STELLAR](https://doi.org/10.5061/dryad.g4f4qrfrc).
The release describes coordinates, region labels, normalized marker intensities
and expert cell-type annotations for donor B004. Its eight intestinal regions
belong to one donor. B005/B006 are unannotated, so inferred types cannot be
independent ground truth. A separate tonsil/esophagus table needs annotation
and donor verification. Release metadata was accessible at version 249582;
its listed CSV file records refer to version 188236. These are source facts,
not verified column-level data.

## Access and annotation gate

Retrieve only through the author's public distribution or another authorized
copy. Verify release SHA-256 checksums and preserve source terms/citation.
Inspect headers, units, donor/section/region identifiers and annotation
provenance before computing any target or policy outcome. Confirm that coordinate
origins are region-specific, exclude duplicates, and prohibit edges across
different regions or sections. Resolve physical pixel calibration from source
documentation; do not assume the ORION scale or convert normalized intensities
directly into validated positivity calls.

Candidate CSV sizes/hashes from the file API:

| File | Bytes | SHA-256 |
|---|---:|---|
| B004_training_dryad.csv | 165599687 | 487a74b1d2d61fa830b0c1ad9e05f0ff4e90cd913de69671fdabdac0235f6344 |
| B0056_unnanotated_dryad.csv | 387756050 | 54a48af81f5e428b123960d93510be81c6fd7040041d830d58b445612ae21393 |
| BE_Tonsil_l3_dryad.csv | 127298344 | 3b308d8cb11782cc628fc658a8cb9541643638a5c545a2b5b589639c6d456729 |

## Prospective transfer test after the gate

Before target evaluation, commit the exact cell-type-to-binary mapping, graph
radii/units, population boundaries, field geometry, policy settings, joint
error budget, precision, random seeds, repeat count and exclusion rules. Use
expert epithelial/T-cell annotations only if their definitions are documented
and suitable. That target would test algorithmic transfer to cell-type spatial
claims; it is not the same PanCK/CD3e assay target. If definitions are ambiguous,
stop and choose a better annotated cohort rather than tune thresholds to win.

Keep envelope, mixture, focused-reuse and deterministic greedy controls frozen.
Report every eligible region, all-time joint coverage relative to source
annotations, attainment, unique fields and feedback rounds over a prespecified
latency grid. Include the new error-sensitivity bounds with externally stated
budgets. Review failure cases, not only mean savings. Source annotations remain
proxy truth; this test cannot validate biological labeling accuracy.

One independent donor supports technical platform/tissue transfer, not a
patient-level replication claim. For a stronger confirmation, require a second
release with multiple independently annotated donors and preserve donor-level
reporting. New unlabeled data alone cannot establish valid coverage. If public
data cannot supply the necessary labels/units, prioritize the lab pilot and
record this gate as unmet; do not substitute predicted labels as truth.
