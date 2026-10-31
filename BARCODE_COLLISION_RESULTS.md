# Full barcode ambiguity and the value of an independent probe

Protocol `0315fe7` and estimator/tests `415bdf0` were committed before codebook
enumeration and photon simulations. This round completes the six-step research
workflow with an exact witness, distinguishing measurement, adversarial failure,
identification boundary, frozen comparison and equation/protocol prior-work audit.

## Source-design result

The unmodified public codebook has 140 codewords, including six blanks. Blanks
are retained in the design audit and excluded from the biological pair universe.

| Universe | Possible unordered pairs | Distinct summed codes | Ambiguous classes | Worlds in ambiguous classes |
|---|---:|---:|---:|---:|
| All codewords | 9,870 | 8,010 | 870 | 2,730 |
| Nonblank genes | 9,045 | 7,413 | **858** | **2,490** |

Thus 27.53% of the **enumerated equal-weight possible gene-pair worlds** belong
to ambiguous classes. This is not a collision rate in tissue: gene abundance,
spatial separation and brightness are not sampled from biological data.

The first recorded gene-only collision contains:

- Adamdec1 + Adgrg6
- Ccr2 + Pik3r1
- Eif2ak3 + Crtc1

All three produce the identical full 16-bit integer sum under the stated
equal-brightness, coincident-source model. Repeating or adaptively choosing
existing barcode rounds cannot distinguish these worlds. The proof is a boundary
case; it does not claim that normal separated spots are undecodable.

An exact independent probe of either gene in a selected world separates it from
all alternatives in its full-sum class. All 4,980 such gene-only separations are
verified. The [theory](BARCODE_COLLISION_THEORY.md) also explains why a negative
probe does not generally resolve the remaining alternatives.

## Photon screen and attempted defeat

The independent probe obeys K~Pois(n*(b+eta*z)). If b is only bounded by [bL,bU],
uniform identification requires eta_min>bU-bL. Otherwise a present-target and
absent-target world can have the same photon law at every exposure count.
Three such exact background witnesses are retained. Zero efficiency and
unconstrained cross-reactivity are additional failure modes.

Under the frozen bounded-background model, the standard calibrated Poisson
decision passes the technical screen. Across all null settings, maximum observed
rejection is **3.1%**, Wilson upper **3.95%**. For eta=8 and 64 exposures,
power is **100%** at both background settings (Wilson lower **99.81%**).
Weak-probe cases inside the overlap boundary remain unidentifiable uniformly.

Using nominal rather than bounded background can give 99.9% false positives
in the high-background stress. This is a deliberately misspecified comparator,
not evidence that our conventional test beats strong calibrated methods. Each
setting uses 2,000 replicates; finite Monte Carlo interval coverage varies around
its theoretical bound and is not empirical calibration of a microscope.

## Verification and limits

```sh
python study_barcode_collisions.py
python verify_barcode_collisions.py
python -m unittest discover -s tests -v
```

The 8 KB attributed source codebook is included; this screen needs no image
download. [All classes and outcomes](results/barcode_collisions.json),
[photon summary](results/barcode_probe_summary.csv) and original compressed count
arrays are retained. [Independent verification](results/barcode_collision_verification.json)
uses sorted integer sums and direct Poisson tail probabilities: 18,915 world
checks across the two overlapping universes and 192,000 method/replicate decisions
from 96,000 unique photon simulations. The full local suite has 95 passing tests.

A separately labeled post-screen deduction extends a binary-sum collision to
positive source separation d. If each readout has independent registration
uncertainty rho, the two image sets overlap exactly for d<=2*rho in the stated
two-source model. An independent Gaussian rendering test verifies the witness
and disjoint-center boundary. This does not use measured registration bounds;
shared fiducials/round transforms can invalidate its assumptions. It supplies
a concrete calibration requirement, not a new claimed universal resolution law.

## Research decision

The [prior-work audit](BARCODE_COLLISION_PRIOR_ART.md) finds that additive code
ambiguity, overlap-aware decoding and extra gene-specific measurements already
exist. The concrete source-codebook audit is useful, but it does not establish
a novel decoder or biological result. Keep novelty around 6.5/10 and research
quality around 7.5/10, with roughly one-point judgment uncertainty.

The next worthwhile question is a claim-specific **minimum extra measurement**
under near-coincident geometry, variable brightness and calibrated nuisance
uncertainty. Establish lower bounds and compare full decoding plus uncertainty,
full acquisition, and a simple independent probe. Independent calibration and
specimens remain necessary; submission is deferred.
