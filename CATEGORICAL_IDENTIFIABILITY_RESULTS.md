# Three sign ambiguities survive valid categorical identities

**Opposite signs remain feasible in three of the six previously evaluated
patients even when epithelial, T-cell and other identities are mutually
exclusive.** Thus the earlier ambiguity was not entirely an artifact of permitting
overlapping identities. This is a stronger conditional limitation, not a new
sign certificate or evidence of a new learning architecture.

The [development protocol](CATEGORICAL_IDENTIFIABILITY_PROTOCOL.md) was committed
at 217a54d before this audit. It reuses the enrichment patients, graphs, cell order
and supplied .5% per-marker Hamming budgets. No new cohort was selected. The
original frozen gate and results remain unchanged. These identities are a cell
ontology assumption, not a claim that protein coexpression or doublets cannot
occur. Source identities had no epithelial/T-cell overlap.

| Patient | Budget per marker | Feasible negative value | Feasible positive/max-search value | Categorical sign ambiguity proved? |
|---|---:|---:|---:|---|
| MB-0040 | 6 | -.00569243 | +.00137285 | yes |
| MB-0043 | 4 | -.00811530 | -.00030441 | inconclusive |
| MB-0045 | 10 | -.00895684 | -.00080050 | inconclusive |
| MB-0046 | 8 | -.01770480 | -.00550251 | no new result; earlier outer bound certified negative |
| MB-0050 | 6 | -.00713887 | +.00045592 | yes |
| MB-0056 | 7 | -.00726602 | +.00065747 | yes |

Values are attained feasible points, not global extrema or interval endpoints.
The deterministic local search uses three frozen starts and at most 50 accepted
single-cell moves per start. Its failure to find a positive world for MB-0043
or MB-0045 proves nothing about whether such a world exists. Existing outer
bounds remain valid supersets; categorical tightening was not substituted for a
full global solver. The previous independent binary-marker witnesses proved
ambiguity in four patients; this audit establishes the stronger categorical
claim in three. Most old positive witnesses contained overlaps, now recorded.

Each saved witness is a list of exact correction indices. Its sign is checked
by the integer numerator

    m * contact_count - (sum_i out_degree_i * a_i) * (sum_j in_degree_j * b_j)

with denominator m squared. Mutual exclusivity and both Hamming budgets are
checked separately. A source-to-artifact replay using independent all-pairs
nearest-neighbor construction reproduces the graphs and numerators.

## From ambiguity to necessary review information

For two feasible opposite-sign categorical worlds x,y and a fixed queried cell
set Q, arbitrary errors on at most E distinct reviewed cells can be tolerated
for this pair iff distance(x_Q,y_Q) > 2E. If distance <= 2E, partition their
positions of disagreement into two groups of at most E and construct an answer
vector matching x on one and y on the other. That vector is compatible with
both; a universally valid sign certificate cannot choose a sign. If distance
>2E, the triangle inequality makes the two answer balls disjoint. The same
shared transcript blocks an adaptive stopping decision while both worlds remain
compatible. This is classical minimum-distance/error-ball reasoning, not a new
active-learning principle. Systematic rereads of one cell do not add distance.

The three saved opposite-sign pairs differ on 18, 21 and 25 cells, respectively.
With E=1, at least three distinct disagreeing cell identities must be observed
to separate each pair robustly. In the existing hypothetical 4x4 field grid, one
carefully selected field suffices for each *pair*. That is a weak, exact pairwise
lower bound on a nonadaptive plan which must distinguish all feasible worlds;
it is not a complete acquisition policy, proof of optimal cost, or measured
instrument savings. Repeated same-cell answers cannot be treated as independent
truth measurements under systematic errors.

## Contribution boundary and next decisive work

The useful direction is to certify which conclusions survive realistic reference
ambiguity and what new measurements can resolve it. A worthwhile future result
would combine categorical identities, uncertain graph geometry and fallible
references in one feasible population, with acquisition plans compared to
certified cost lower bounds. No joint guarantee or global acquisition optimum
has been established here. Independent observed expert annotations are the next
assumption check, rather than more favorable simulations with invented error
budgets. The [NuCLS screen](NUCLS_REFERENCE_PROTOCOL.md) addresses reference
choice on a different H&E identity target; it cannot establish protein truth.

Classical foundations: [version-space compression](https://jmlr.org/papers/v16/wiener15a.html)
and disjoint Hamming balls. Neither the covariance estimand, local coordinate
search nor the pairwise distance criterion earns a standalone 9-level novelty
claim. Geometry/assay calibration and independent method validation remain open.

Data: Danenberg et al., DOI 10.5281/zenodo.7324285, CC BY 4.0. Raw public source
files remain outside version control; checksums and all six cases are retained
in results/categorical_identifiability.json.

```bash
python study_categorical_identifiability.py
python -m unittest discover -s tests -p 'test_categorical_identifiability.py' -v
python verify_categorical_witnesses.py
```
