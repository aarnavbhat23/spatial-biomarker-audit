# A full-image ambiguity, and what it does not imply

2026-11-01. This extends a previously tested centroid construction to full
images. Equality of observation laws and binary decision lower bounds are
classical statistical tools. It is not a new decoding architecture.

Let c_a,c_b be binary codewords, each read out in known round/color slots.
In bit j, two equal-brightness objects at +/-d/2 have gene-image mean

    lambda_j(p) = b_j + A_j sum_{i in {a,b}} c_i[j] K_j(p-x_i-t_round(j)-s_color(j)).

The PSF K_j can be any common translated profile, not necessarily Gaussian.
Suppose another pair c_u,c_v satisfies c_a+c_b=c_u+c_v with every entry <=1.
There is then one emitter in every active bit. Let m0_j and m1_j be its location
in the two worlds before instrumental shifts. If

    m0-m1 = Ar*u + Ac*v,

choose world 0 shifts (-u/2,-v/2) and world 1 shifts (+u/2,+v/2).
The emitter locations after shifts coincide in every bit. With the same
background and active-bit brightness, the complete mean images coincide.
Independent Poisson pixel laws consequently coincide, as do any common
camera-noise kernels applied to them. This statement concerns distributions,
not only matching one noisy realization.

For any forced binary decision, write e0=P0(decide 1), e1=P1(decide 0).
Equality of laws implies e0+e1=1, so max(e0,e1)>=1/2. Any adaptive sequence of
repeated exposures to these same bit operators also preserves equality:
conditioned on any identical history, the next action is chosen by the same
policy and has the same conditional observation law. Induction gives equality
of the complete history. New probes or independently calibrated registration
can change the observation operator and are outside this impossibility claim.

The four saved Adamdec1 + Adgrg6 alternatives satisfy the hypotheses with the
actual source codebook and round/color layout. Independent verification checks
all 64 active-bit center identities at d=1 and d=3, including a non-Gaussian PSF.
The pixel benchmark checks rendered equality at both declared exposures.

## Why this is a limited, useful counterexample

* Equal per-object brightness is a least-favorable assumption, not an empirical
  statement about the tissue. Different brightness can supply distinguishing
  information if it is stable across bits. White et al. explicitly use a
  log-brightness-variance term in their decoding objective (Eq. 9).
* Known registration removes these particular aliases at nonzero separation.
  ISTDECO's perfect high-signal registered controls confirm that limitation.
* Unknown per-bit gain, background, dropout and probe cross-reactivity require
  additional nuisance modeling. The ideal specific-probe experiment does not
  establish a laboratory benefit without those calibrations.
* A lower bound at an admissible mechanism does not estimate how frequently
  ambiguous mechanisms occur in a biological population.
* The independent author reference table quantifies selected-feature alignment
  on one field. Its pixel errors cannot be transported to a different instrument
  or tissue as calibrated Gaussian variances.

The conditional Gaussian fiducial result in LOCALIZATION_INFORMATION_THEORY.md
addresses a different, explicitly restricted experiment. Do not combine its
numerical error floor with these Poisson images as if the covariance had been
measured. The common contribution is a way to ask what extra measurement is
needed; the current data do not establish a general-purpose practical guarantee.

Priority: [White et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/)
already identifies barcode recombination under imperfect registration and uses
strong codes, brightness and global optimization. [ISTDECO](https://github.com/axanderssonuu/istdeco)
already performs joint spatial/spectral deconvolution. The possible additional
contribution is a calibrated, claim-specific necessity result. Independent
priority feedback and independently annotated validation remain necessary.
