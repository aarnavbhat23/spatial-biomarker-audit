# Frozen external reference-choice screen

This screen tests a measurement assumption, not a prediction algorithm. The
independent NuCLS unbiased-control release contains matched expert annotations
without displayed algorithmic suggestions. The downloaded anchor table has only
been inspected for columns and grouping: 1,569 anchors, 54 fields, five TCGA
patients, unique anchor IDs. No class counts, disagreement rates or spatial
outcomes have been evaluated before this protocol commit.

Source: Amgad et al., NuCLS, GigaScience 2022,
https://doi.org/10.1093/gigascience/giac037 ; data CC0 1.0;
https://sites.google.com/view/nucls/multi-rater .
Direct file ID 1RTyfmkZW8um3qmWqtu2gbu3Sgql-DQQj,
v3.1_final_anchors_U-control_Ps_AreTruth.csv,
SHA256 07edc0db091a51f90e58fcc18332c347f5f62a9356320c9aba7377e666469ff6.
The anchor population is selected by the source's pathologist inference process;
its inferred P-truth is NOT independently adjudicated biological truth.

Use the three SP.1/SP.2/SP.3 annotations as alternative expert references, without
using EM/MV inferred class labels as truth. On each field retain only anchors
with finite bounding boxes and a detected, labeled biological class from all
three experts. Exclude undetected, unlabeled and missing annotations. Publish all
exclusions and observed class encodings. If encodings are unexpected, stop rather
than silently remap. Population selection and shared anchor geometry remain
limitations. This is a class-assignment comparison, not a detection-error study.

New target: tumor-nucleus to lymphocyte association on H&E, NOT epithelial-to-T
protein-marker association and NOT validation of virtual protein readouts.
A=exact class 'tumor'; B=exact class 'lymphocyte'; other detected biological
classes map to neither. Use shared anchor bounding-box midpoints, six directed
nearest neighbors, field size at least 20 and deterministic distance/index ties.
No cross-field edges. Duplicate coordinates trigger an explicit feasibility stop.

For each expert and field calculate the previous degree-adjusted edge covariance
with exact rational arithmetic. Patient target is the edge-count-weighted average
of within-field covariances over the same retained fields for all three experts.
Report every field and all five patients, including ineligible fields. Define
observed-reference sign stability as all three patient values strictly sharing
one sign; a zero value is not stable. Report sign reversals and full ranges even
if there are none. Also report pairwise categorical and A/B disagreement counts
on matched nuclei; these measure reference disagreement, NOT error versus truth.

A one-wrong-cell review assumption is compatible with two complete reference
vectors only if their categorical Hamming distance is <=2 (a common unknown
truth may be within one of each). Apply the analogous triangle bound to the
previous .5%-per-marker raw budgets, without interpreting a failure as an assay
error-rate estimate. No significance or clinical generalization claim is planned:
there are only five patients and all experts share the source acquisition.

This supplies external measured-reference evidence about assumptions. It does
not establish coverage against true protein calls, measured acquisition savings,
a new estimator, or a 9-level novelty claim. No gate will be selected from these
outcomes and no new patient subset will be chosen for a favorable result.
