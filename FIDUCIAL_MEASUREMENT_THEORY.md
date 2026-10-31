# Minimum measurements conditional on calibration

The fixed claim is Adamdec1+Adgrg6 rather than Ccr2+Pik3r1 or Eif2ak3+Crtc1,
conditional on an exactly known two-molecule, binary codeword-sum class. This
is a restricted composition claim, not an assertion about cell interaction.

## Shared nuisance parameters change the ambiguity

Put the two equal-brightness molecules at fixed coordinates +/-d/2 on a fixed
axis. Every active bit in the selected sum has one emitter. Let h be the
difference of the unit-separation emitter centers in two proposed compositions.
Write an active bit's translation as t_round + c_color, with a shared translation
per round and a constant offset per color. A common PSF and background are assumed.
The two complete noiseless image collections coincide if and only if

    A_round * u + A_color * v = d*h,

where u and v are differences between the explanations' nuisance parameters.
Thus mean-center equality is also sufficient for image equality in this
single-emitter-per-active-bit model. It does not cover brighter unequal molecules,
nonbinary sums, variable PSFs, arbitrary local chromatic distortion, unknown
source geometry or finite-photon distinguishability.

Projecting h out of the nuisance column space is conventional linear algebra.
It reveals why independent bit shifts overstate ambiguity, and why perfectly
known round translations can remove the tested aliases. A nonzero residual gives
structural separation in the declared exact-data model, not practical decoder
accuracy. At d=0 all compositions in a sum class still coincide.

## Conditioning on the actual selected fiducials

Let v_ri be a matched bead displacement in a given round and coordinate. If its
error relative to one common round translation is bounded by ell, and matches
are correct, the feasible round translation interval is

    I_r(ell) = [max_i(v_ri)-ell, min_i(v_ri)+ell].

An empty interval contradicts that contract for the retained data. If nonempty,
the set of differences of two translations in I_r is exactly [-w_r,w_r], where
w_r is its full width. The interval center cancels because both explanations
must explain the SAME fiducial measurements. Round 1 defines the reference and
has relative translation zero. The maximum ambiguous separation is the LP

    maximize d >= 0
    subject to A_round*u + A_color*v = d*h,
               |u_r| <= w_r, u_1 = 0; v unrestricted.

Every feasible u can be realized by two translations at the interval midpoint
+/-u/2. Therefore this is an exact threshold within the stated scalar/axis
model, not merely a loose outer bound. Feasible separations include every value
from zero to the maximum by scaling the difference variables. Saved primal and
dual witnesses certify feasibility and optimality. Both horizontal and vertical
axes are evaluated; arbitrary orientations are not established.

The actual fiducials are recorded at 473 nm, whereas gene readouts use 561 and
650 nm. Their matched residuals do not measure gene-color error. Their detection
cap and <=2-pixel matching rule select/truncate residuals; the coarse starting
translation already used the full development image. Neither the smallest
compatible ell nor a residual percentile supplies a coverage guarantee. Missing
or incorrectly matched beads, deformation and varying chromatic shifts can
invalidate this model. An externally calibrated error contract is still needed.

## What one additional image could identify

An independently specific Adamdec1 probe has indicator z=1 for the first
composition and z=0 for either alternative. Suppose one added image contains
independent target and target-free reference ROIs, with reference area 16 times
larger. Counts are

    X ~ Poisson(E*(b*(1+s)+eta*z)),  Y ~ Poisson(E*16*b),
    b in [2,10], |s| <= delta, eta >= eta_min.

E is hypothetical photon dose. Probe specificity, efficiency and the
representativeness of the reference ROI must be independently justified. The
probe must cover both candidate source positions; the fixed model does not
include loss by target mislocalization. Raw camera intensities are not photon
counts without gain/noise calibration.

Even an infinitely bright reference does not solve arbitrary mismatch. At a
given b, the absence and minimum-presence target-intensity sets intersect if
eta_min <= 2*delta*b. For example b=10, delta=.25, eta_min=2 admits a positive
world with s=-.25 and an absent world with s=-.05, both with target rate 9.5 and
reference rate 160. Their entire joint count laws coincide at every dose. No
estimator can reliably settle this expanded ambiguity. Strict separation requires
eta_min > 2*delta*b; uniformly over the full b range it requires eta_min>20*delta.
This is a model identification boundary, not a universal imaging theorem.

The original conditional binomial test is a standard exact comparison of two
Poisson rates. It discards total intensity and uses a conservative b=10 bound.
The follow-up keeps absolute reference intensity through exact Poisson bounds.
With each one-sided tail at .0125, form intervals [LX,UX], [LB,UB]. Declare:

    presence if LX > (1+delta)*UB;
    absence if UX-(1-delta)*LB < eta_min;
    abstain if both or neither condition holds.

The rectangle has simultaneous coverage >=.95 by a union bound. Under absence,
a false presence requires either LX>lambda_X or UB<b. Under presence, a false
absence requires UX<lambda_X or LB>b. Each wrong declaration therefore has
probability <=.025 under its contract. This is an ordinary conservative interval
construction. It is not a new estimator, an optimal allocation or an equal-cost
win over a target-only image: it uses a larger reference ROI and its photons.

## Originality and the next falsification

Chromatic registration is established work: [Matsuda et al., Chromagnon](https://doi.org/10.1038/s41598-018-25922-7)
shows that multicolor alignment needs sample-relevant calibration. Bead residuals
and uncertainty at a target feature are also distinct in
[registration-error analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4487814/).
[Exact conditional Poisson-rate tests](https://www.statsmodels.org/v0.13.5/generated/statsmodels.stats.rates.test_poisson_2indep.html),
confidence rectangles and LP duality are conventional. The [preceding priority
audit](BARCODE_COLLISION_PRIOR_ART.md) already records overlap decoding and
measurement-design precedents. This combination is an application audit; priority
for a paper-level contribution has not been established.

The specific research opportunity is to certify when *this composition claim*
needs a gene probe rather than better geometric calibration, and when even that
probe cannot help. The next decisive test requires multicolor reference spots
or repeat captures with independent truth, plus measured target/background
mismatch and a lower probe-efficiency bound. Freeze those contracts on a
calibration field, then test coverage and decisiveness on a new field. A real
image outside the contract must permit abstention; do not tune the uncertainty
budget against its desired answer. Lab timing is needed before any cost claim.
