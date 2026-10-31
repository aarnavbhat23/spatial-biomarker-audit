# Measured fiducial-set follow-up — frozen 2026-10-30

This is a post-screen theory/development extension, not independent validation.
Keep previous fiducial selections and outcomes. Before evaluating this extension,
freeze coordinate-wise displacement allowances ell=[.1,.25,.5,1] pixels, both
horizontal and vertical fixed-source configurations and all four existing target
comparisons (two alternative compositions, both gene orientations).

For each of rounds 2..8 and each coordinate, condition on ALL previously saved
matched bead displacements v_i. Under an externally justified per-displacement
absolute error bound ell and correct bead matches, the common round translation
lies in [max(v_i)-ell, min(v_i)+ell]. Round 1 is the reference, so its relative
translation is exactly zero. If any interval is empty, report that the contract
is incompatible with the retained matches. Do not drop extreme matches.

For two explanations of the same observations, the allowed translation difference
in a round is the interval's width in either sign; the observed interval center
cancels. Constant gene-color offsets remain free. For each retained comparison,
maximize d>=0 subject to A_round*u+A_color*v=d*h and |u_r|<=interval_width_r,
with u_1=0. h is the already specified unit-separation mean-center difference.
Solve separately for x and y. The orthogonal coordinate can use zero nuisance
difference. Save LP witnesses and replay their equalities and inequalities.

This yields a sharp mean-center ambiguity threshold *conditional on the assumed
uniform displacement error, selected matches, shared translation, constant color
offsets and fixed source positions/orientation*. It is not a calibrated resolution
limit, a photon-level decoder or a bound on tissue error frequency. A small
required allowance is a consistency diagnostic, not a confidence guarantee.
Chromatic variation not captured by a constant color offset would enlarge the
model and can invalidate a distinguishability statement. Do not choose ell from
the observed percentiles or claim source-specific calibration was obtained.
