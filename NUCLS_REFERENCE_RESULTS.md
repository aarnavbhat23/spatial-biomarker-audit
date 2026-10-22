# External expert annotations challenge the tiny-reference-error assumption

On the independent NuCLS unbiased-control release, **six of seven eligible
field comparisons have expert disagreements incompatible with both experts being
within one wrong target-category cell of a common truth**. This is measured
reference disagreement, not a measurement of which expert is wrong or a protein
assay error rate. The sign study itself is uninformative: one expert's selected
spatial contrast is zero in every eligible field.

The initial [three-expert protocol](NUCLS_REFERENCE_PROTOCOL.md), committed at
d666fa1, was unevaluable: SP.1 did not annotate any field in this release. An
[explicit sentinel note](NUCLS_INPUT_NOTE.md) preceded evaluation. Its artifact
is retained with zero eligible fields; zero comparisons is not zero disagreement.
A [separate two-expert extension](NUCLS_PAIR_PROTOCOL.md), committed at 5f30154,
selected SP.2/SP.3 based solely on common annotation coverage before any pairwise
spatial outcomes. It kept the same field size, graph, class mapping and five
patients. No outcome-selected patients, fields or changed target thresholds.

## What was actually compared

Source: Amgad et al., [NuCLS](https://doi.org/10.1093/gigascience/giac037),
[public multi-rater release](https://sites.google.com/view/nucls/multi-rater),
CC0 1.0. The downloaded table contains 1,569 anchors, 54 fields and five patients.
Only 622 anchors have detected, labeled annotations from both available senior
experts. Seven fields meet the frozen minimum of 20 common anchors, retaining
187 anchors from three patients. All fields/patients, including exclusions, are
reported. The source's pathologist-derived anchor population and shared bounding
boxes are selected proxies, not independently adjudicated geometry or truth.

The target is **tumor-nucleus to lymphocyte association on H&E**, not protein
marker positivity or the earlier epithelial-to-T-cell estimand. Six directed
neighbors define a fixed common graph within each field. Exact patient values
average within-field degree-adjusted covariances with edge-count weights.

| Patient | Eligible fields | Common anchors | SP.2 patient contrast | SP.3 patient contrast |
|---|---:|---:|---:|---:|
| TCGA-A1-A0SP | 3 | 75 | 0 | 0 |
| TCGA-A7-A0DA | 0 | 0 | unavailable | unavailable |
| TCGA-AR-A1AR | 2 | 45 | 0 | +.000322061 |
| TCGA-C8-A12V | 2 | 67 | 0 | -.009950249 |
| TCGA-E2-A158 | 0 | 0 | unavailable | unavailable |

No patient has strictly opposite reference signs; no patient has two strictly
same-sign values either, because SP.2 has no eligible field containing both tumor
and lymphocyte target labels. This composition limitation blocks meaningful sign
validation. Do not report 0/5 sign reversals as evidence of robust generalization.
There are only three evaluable patients, not five independent comparisons.

## The useful assumption check

Mapped three-category distances between SP.2/SP.3 across the seven fields are
4, 5, 5, 0, 4, 13 and 19 cells. If both references were within E=1 wrong cell of
the same underlying categorical truth, the triangle inequality would require
distance <=2. Six pairs violate that necessary condition. Analogously, six pairs
violate compatibility with both experts having <=ceil(.005*n) errors in each
binary target label. Here field sizes are 21–39, so ceiling makes that integer
allowance one error per marker (roughly 2.6%–4.8%, not an actual .5% rate).

This does NOT show the selected reference has more than one error, nor that
protein-marker error rates match H&E classification disagreement. It shows that
a uniformly tiny error allowance for *every available expert* cannot simply be
assumed in this population. The experts agree exactly on the tumor indicator
in these fields; lymphocyte assignments cause the mapped disagreements. FOV
absence, undetected cells and source anchor selection are excluded from this
class-only comparison and can add further uncertainty.

The next useful project contribution is an empirically calibrated reference
uncertainty model which preserves a scientifically useful decision on fields
containing both targets. It needs independently defined reference/adjudication
and patient-level holdouts. More synthetic acquisition savings at E=1 would
avoid the unresolved assumption. These external data add a real check on
measurement realism, not a new architecture or a 9-level novelty result.

Full source hash, classes, exclusions, integer field numerators, exact rational
patient values and comparisons: results/nucls_references.json (unevaluable
initial screen) and results/nucls_pair_references.json (coverage-driven extension).
Independent all-pairs graph construction and rational-arithmetic replay reproduce
all seven extension fields. This is independent computation, not third-party
experimental replication.

```bash
python study_nucls_references.py
python study_nucls_references.py --available-pair
python verify_nucls_references.py
```
