# When gene photons cannot compensate for uncertain registration

2026-11-01. Conditional theory for a known-variance Gaussian localization model.
The projection/test is classical. The application question is which measurement
must improve before a specified barcode-composition decision becomes reliable.
No claim of a new decoder, universal resolution theorem or biological validation.

## Experiment and exact decision risk

Compare two fixed two-molecule compositions with active-bit center vectors m0
and m1. Their difference is h=d*h_unit. The source codebook/readout map determines
round and color design matrices Ar, Ac. Assume

    G = m_w + Ar*t + Ac*c + epsilon, epsilon ~ N(0,vg*I),
    F = t + zeta,                    zeta ~ N(0,vf*I),

with independent errors and exactly known positive variances. Round translations
t and constant channel offsets c are unrestricted nuisance parameters. F denotes
an independently measured translation estimate in each round, not eight copies
of one shared-reference error. Real shared-reference correlations would require
the corresponding full covariance model.

Subtracting Ar*F gives Z with mean m_w+Ac*c and covariance

    C = vg*I + vf*Ar*Ar^T.

Let v minimize (h-Ac*v)^T C^{-1}(h-Ac*v), put r=h-Ac*v and

    D^2 = r^T*C^{-1}*r,  w=C^{-1}*r.

Then Ac^T*w=0. The threshold statistic w^T*(Z-(m0+m1)/2) has means -D^2/2
and +D^2/2 and common variance D^2. Hence the conventional invariant linear test
has error Phi(-D/2) in each world for every allowed nuisance value. D=0 gives
indistinguishable projected laws and error 1/2.

To show minimax optimality, whiten the full joint observation (G,F). Its two
mean sets are parallel affine spaces with common nuisance subspace. Their minimum
Euclidean distance is D. Choose the nuisance difference attaining that distance.
The equal-prior simple Gaussian testing problem between that pair has Bayes error
Phi(-D/2); every test has maximum error at least its average error on that pair.
The invariant statistic attains this lower bound uniformly. This proves the
minimax maximum-error value for THIS binary Gaussian experiment. It does not
imply every decoder has the same risk on real raw images or multiworld claims.

## Calibration ceiling

The same squared separation has the joint variational form

    D^2(vg,vf) = min_(u,v) ||h-Ar*u-Ac*v||^2/vg + ||u||^2/vf.

Keep the fiducial information fixed and let gene localization variance approach
zero. If Ar*u+Ac*v=h is feasible, then

    D_infinity^2 = min ||u||^2/vf subject to Ar*u+Ac*v=h.

The upper limit follows by any feasible pair. For the lower limit, bounded
objective controls u and forces the first residual to zero; projection onto the
finite-dimensional color subspace supplies a feasible limiting solution. Thus
the limit equals the minimum constrained energy. If the equation is infeasible,
the nonzero residual outside the nuisance span makes D^2 diverge instead.

For feasible aliases, the optimal error cannot fall below Phi(-D_infinity/2)
by increasing gene photons alone. This does not apply when new independent
fiducial information is acquired as well. If vf=tau^2/n, the ceiling information
scales with n/tau^2. To reach maximum error alpha, a necessary and sufficient
condition in this idealized experiment is D^2 >= 4*Phi^{-1}(1-alpha)^2.

An independent, calibrated scalar probe with known mean difference a and variance
s^2 adds a^2/s^2 to D^2. The SNR4 control therefore adds 16 by assumption. If an
unrestricted probe-background offset is included in the nuisance subspace, that
dimension can carry zero identifying information. Probe selectivity/background
calibration is a requirement, not a result supplied by these equations.

## Boundaries and prior work

The current application keeps one emitter per active bit, fixed positions,
known orientation/brightness assumptions and a Gaussian localization summary.
PSF fitting, unequal brightness, missing readouts, photon counting, nonlinear
chromatic distortion and uncertain variance are not modeled. The public
cross-channel benchmark tests registration behavior; its selected-match radius
does not estimate vg/vf or establish Gaussianity for the MERFISH target field.

Interference projection and invariant Gaussian detection are established:
[Scharf and Friedlander, 1994](https://doi.org/10.1109/78.301849).
The distinction between registration error and target-feature localization is
also established in [fluorescence registration analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4487814/).
Our closest raw-image comparator remains
[overlap-aware barcode decoding](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/).
The possible paper contribution is a useful claim-specific calibration/measurement
boundary and its validation. The general projection, Gaussian risk formula and
the existence of calibration-limited precision are not claimed as inventions.
Equation-level priority for the complete application result remains unresolved.
