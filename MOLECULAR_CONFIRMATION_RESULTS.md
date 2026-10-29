# Public molecular proximity under bounded measurement error

2026-10-29. **The frozen public-data usefulness gate fails.** There are no
Bonferroni-corrected positive-proximity rejections in either technical region,
even before measurement error. The result does not imply absence of biological
interaction. It stops this marker/scale/null combination as a useful validated
method claim.

## Frozen analysis and provenance

[Protocol](MOLECULAR_CONFIRMATION_PROTOCOL.md), committed as `a0c9526` before
download/scoring; [complete results](results/molecular_confirmation.json);
[source manifest](results/visp_source_manifest.json). The public VISp table
contains 530,783 measured molecules. Its [source documentation](https://ssam.readthedocs.io/en/develop/userguide/02-data.html)
explicitly gives coordinates in micrometres. SHA256:

    7121a8c65e0e493f995d84e764485759ad14b1eae2338b322ed2654b832ce9c6

Of six fixed markers, Gad2, Pvalb and Cux2 are available; Slc17a7, Sst and Vip
are missing. They were not replaced. Each half has three evaluable pairs at
three scales (9 tests), and 12 unevaluable pairs. Correction remains for the
original 45-test family. At most 200 molecules per marker per half are selected
with the frozen seed. The target is this selected measured population; the
study does not establish power for all molecules or population-level coverage.

| Outcome | Left development | Right confirmation |
|---|---:|---:|
| Evaluable pair/scale tests | 9 | 9 |
| Global label-null adjusted rejections | 0 | 0 |
| Stratified label-null adjusted rejections | 0 | 0 |
| Robust rejections, 0.5 micrometres, no swaps | 0 | 0 |
| Robust rejections, 0.5 micrometres, one swap | 0 | 0 |

Both halves use the same frozen analysis. They are regions of one prepared
mouse tissue field, not independent animals, patients or cohorts. Source cell
assignments are included only as a within-cell diagnostic, not ground truth.
Raw molecular proximity is not cell-type association or causal communication.

## Conditional robustness argument

For fixed recorded strata and a fixed measured population, assume true binary
gene labels are exchangeable within strata conditional on positions and counts.
Use the mean cross-channel Gaussian kernel S. For displacement at most rho per
point, each true distance lies in [max(0,d-2rho),d+2rho], giving kernel bounds.
Up to k A/B swaps within strata change a bounded pair mean by at most
e=k/nA+k/nB: couple uniform cross-pair draws and bound the probability that a
replaced member is sampled. This bound may be loose.

Let L0 be the observed score's lower bound and Ub each permuted score's upper
bound, with this swap allowance. Couple the same within-stratum permutation
actions on recorded and admissible true labels. They preserve the number of
label differences, so every admissible latent permutation p-value is at most

    p_upper = (1 + number of b with Ub >= L0)/(B+1).

Under the stated exchangeability null, uniformly sampled independent group
permutations yield the usual conservative Monte Carlo p-value. Pointwise
domination preserves validity. Bonferroni45 gives familywise control for each
fixed declared error budget. This is a composition of conventional bounds and
permutation inference, not an established novel statistical method.

The validity statement is **conditional on the error budget containing truth**.
The budgets have not been calibrated from instrument reference measurements.
Missing molecules, new channel members, swaps across the fixed strata and
gene-dependent selection outside the null are not covered. An unknown-error
confidence guarantee cannot be obtained by simply calling these sensitivity
assumptions realistic. Separate exact tests check 325 small position/label/budget
configurations under all eight allowed permutation actions.

## Why the next question needs to change

The following is a **post-outcome algebraic diagnostic**, not a replacement
test or a successful confirmation. Let nA=nB=m, K be the Gaussian Gram matrix
on the pooled 2m positions, T=sum(K), and z contain +1 for A and -1 for B. Then

    MMD_biased^2 = z^T K z / m^2 >= 0,
    S = T/(4m^2) - MMD_biased^2/4.

T is constant under label permutations. Therefore maximizing this proximity
score is equivalent to minimizing the conventional empirical kernel discrepancy.
For the unrestricted balanced label null,

    E_perm[S] = (T-2m)/(2m(2m-1)),
    S-E_perm[S] = (1-T/(4m^2))/(2m-1) - MMD_biased^2/4
               <= 1/(2m-1).

The identity also holds for each stratified assignment, but the displayed
expectation and ceiling are for global permutations only. This is elementary
algebra using established [MMD](https://www.jmlr.org/papers/v13/gretton12a.html),
not a priority claim. Random labeling addresses a different null from independent
spatial processes; see [the primary comparison of intensity tests](https://doi.org/10.1080/03610918.2021.1901118).
Shared tissue clustering can remain entirely compatible with label exchangeability.
The earlier shared-versus-separated synthetic comparison did not establish
detection under this public-data null.

For m=200 the global positive-excess ceiling is 1/399=.002506, while the generic
one-swap allowance is .01. This does not alone prove a failure for every dataset
or stratified test. In the actual saved 18 cases, all scores are below the
global permutation mean, and every lower bound at (rho=.5,k=1) is zero. Thus
those robust p upper bounds equal one. These diagnostics explain why further
bound variants or more permutations are not the next priority.

## Paper decision and next phase

Keep this negative result and defer the full manuscript. A promising research
question is **when can a spatial conclusion distinguish shared tissue anatomy
from cross-channel dependence under a measured error model?** First specify
the biological alternative and design a falsification benchmark: independent
channels with shared heterogeneous intensity, same-cell coexpression, distinct
domains, and genuine coupled signals. Require controlled false-positive rate
and power for the intended distinction before opening another tissue field.

Compare established conditional-label, intensity-adjusted and appropriate
spatial-process controls on that same estimand. A new guarantee must add useful
calibrated measurement robustness, rather than rename one of those nulls.
Then freeze marker eligibility using metadata, use fresh sample identities and
obtain localization/gene-call/missing-detection reference measurements. The
current field is now development data. Real review/scan costs still require the
previously planned lab pilot. None of these future stages is reported complete.

## Reproduction and verification

Download the table linked by the source documentation to
`data/visp/s3_spot_table.csv`; match the SHA256 above. The data remain untracked.

    python study_molecular_confirmation.py
    python verify_molecular_confirmation.py
    python study_molecular_null_diagnostic.py
    python -m unittest discover -s tests -p test_molecular_permutation.py -v

Separate source replay regenerates selection and permutations and recomputes all
72 sensitivity results with a different kernel-distance/statistic calculation.
Maximum score discrepancy is 7.81e-18. This is computational reproducibility,
not independent biological replication. Small exact checks and the algebraic
identity run in CI; the 55 MB public source replay is a separate local check.
