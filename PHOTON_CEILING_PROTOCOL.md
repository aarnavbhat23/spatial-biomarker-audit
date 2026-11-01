# Finite localization information with uncertain registration — frozen 2026-11-01

Motivated by the failed reference-radius gate; this is a separate conditional
theory/simulation screen, not a rescue of that gate or new independent biology.
Before evaluation, retain all four oriented alternatives for the fixed
Adamdec1+Adgrg6 composition. Use the existing source codebook and round/color map.

Observation model: gene-center coordinates G=m_w+A_r*t+A_c*c+epsilon, and
independent fiducial translation estimates F=t+zeta, with epsilon~N(0,sigma_g^2 I)
and zeta~N(0,tau^2/n I). Nuisance t and constant channel offsets c are unbounded.
These are known, correctly specified variances and independent Gaussian
observations. The public selected-bead benchmark does NOT calibrate this model.

Subtract A_r F. The remaining covariance is C=sigma_g^2 I+(tau^2/n)A_r A_r^T.
For two compositions with mean difference d*h, project out A_c in the C^{-1}
metric. Their squared separation D^2 is this projected quadratic form. Check
against direct whitening/projection of the complete joint (G,F) experiment.
The conventional generalized least-squares test is the strong control itself.
Under this affine Gaussian family the minimax equal-error probability is
Phi(-D/2); derive the simple-pair lower bound and an attaining nuisance-invariant
test. Do not claim a new testing principle.

As gene photon multiplier E increases, sigma_g=.3/sqrt(E). At E=infinity,
compute D_infinity^2=min ||u||^2/(tau^2/n) subject to A_r*u+A_c*v=d*h. A finite
ceiling exists when this equation is feasible. This connects the previous exact
alias to finite calibration information. More gene photons cannot cross an
unmet error threshold at the ceiling. Extra independent reference observations
increase n; an independently calibrated gene probe of SNR4 adds 16 to D^2.
An unconstrained probe-background offset can eliminate that added information.

Grid: d=[.25,.5,1], tau=[.1,.3], n=[1,16], E=[1,16,256] for all four alternatives.
Also record the infinite-E limit. For the first original-orientation alternative
only, generate 2,000 replicates per truth and setting; seed2026110102. Generate
full independent gene/fiducial noise before projection. Save the projected
statistics and independent probe statistic for replay. The classifier uses equal
priors/thresholds, no abstention; it is not the previous count-based certificate.

Gate: direct joint-projection and covariance computations agree to 1e-8; all
empirical error deviations from the analytic binomial expectation are <=.05.
Report exact predicted errors, measured errors and every setting. A statistical
or model gate failure stays a failure. The practical threshold is <=5% error;
report how many settings each measurement design reaches rather than selecting
one favorable case. The ceiling result is limited to Gaussian localizations,
fixed two-source geometry, one emitter per active bit and known variances. It
does not validate unequal brightness, dropout, PSF fitting or an overlap decoder.
Those raw-image comparisons remain separate submission requirements.
