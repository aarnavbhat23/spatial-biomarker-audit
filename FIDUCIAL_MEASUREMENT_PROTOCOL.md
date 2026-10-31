# Shared fiducials and one additional image — frozen 2026-10-30

Scope: condition the previous barcode ambiguity on the actual source acquisition
layout, audit observed fiducial consistency, and test a new measurement contract.
No tissue molecule truth, new cohort or scanner savings are claimed.

## Fixed target and source

Retain the previously reported class: Adamdec1+Adgrg6 versus Ccr2+Pik3r1 and
Eif2ak3+Crtc1. The binary claim is the first composition within this exactly
known two-molecule sum class. A probe for Adamdec1 separates it if independently
selective. Do not substitute a more favorable target after results.

Use existing field 000 fiducials and source dataorganization.csv. Fiducials are
473 nm; readouts are 561/650 nm, with two colors per round. The source fiducial
array SHA-256 is 3219c1f7dc1b25ae4471ce6453f7119b22c17603053cca4ff6c766915827a247.
Original whole-image translations were already evaluated on this development
field. Neither those estimates nor tile disagreement are calibrated error bounds.

## Fiducial diagnostic, fixed before bead-level outcomes

Detect peaks independently in each of nine fiducial planes using positive
Gaussian(1)-Gaussian(4), 99.9th percentile threshold, minimum distance 8 pixels,
border exclusion 12, maximum 300 peaks ranked by filtered intensity. Refine
centroids within radius 3 using the positive filtered weights. Match to round 1
with the existing whole-image translation, mutual nearest neighbor, distance<=2
pixels. Save all detections and accept/reject counts. Matching truncates residuals
and uses a previously fitted shift: this is a consistency audit, not unbiased
error calibration. Do not retune if matches are scarce.

Fit translation from checkerboard-even reference tiles of side 128 pixels using
the median matched displacement; report residuals on odd tiles. Require >=10
training and >=10 held-out matches per evaluated round; otherwise flag that
round unevaluable. Report medians, 95th percentiles, maxima and the shift relative
to the original estimate. No percentile becomes a hard rho or cross-color bound.

## Structured geometry audit

For every nonblank binary-sum collision class, enumerate every pair of worlds
and both assignments of the alternative genes to fixed positions +/-d/2. Set
d=1 for normalized algebra. Compare allowed mean-center adjustments:
(i) independent per-bit shifts; (ii) shared round drift; (iii) shared round drift
plus constant channel offsets; (iv) exactly known round drift plus channel
offsets. Use the actual readout/round/color mapping. Report residual distance
from the relevant design-matrix column space and verify witnesses against
direct linear equations. Enumerate all cases, not only the target example.

For the fixed target, also minimize required per-world round-drift radius with
free constant color offsets by a linear program on the difference between two
worlds' nuisance parameters. Absolute observed drift is corrected in both worlds
and cancels; what matters is allowed residual error. Do not use the observed bead
residuals as that radius. Shared-color offsets across rounds are a model assumption.

## One-image probe/reference simulation

An added gene-specific image contains a target ROI and a preselected background
ROI of 16 times its area. They yield independent counts
X~Pois(E*(b*(1+s)+eta*z)), Y~Pois(E*16*b), with b in [2,10], |s|<=delta.
E is a hypothetical photon-dose multiplier for ONE image, not a measured
integration time. The background ROI must be target-free and representative;
its existence and calibration are not established in the public field.

Condition on N=X+Y. Under absence, p<=p0=(1+delta)/(17+delta). Under presence
with eta>=eta_min, p>=p1=((1-delta)*10+eta_min)/((17-delta)*10+eta_min).
Declare presence using Binomial upper-tail p<=.025 at p0; declare absence using
lower-tail p<=.025 at p1; if both or neither reject, abstain. This is standard
conditional Poisson-rate inference, not a new statistical test. Keep even
overlapping-hypothesis settings; do not silently force a binary answer.

Compare a target-only image using exact Poisson tails against the entire same
b range and efficiency bound. Both use the same target count; the paired method
uses additional background area/photons and is not an equal-information win.
The exact conditional test is the strong conventional method itself.

Grid: b=[2,10], delta=[0,.1,.25], eta_min=eta=[.5,2,8], E=[1,4,16,64],
z=[0,1], 2,000 replicates each; seed 2026103003. Use worst target-region offset
s=+delta under absence and -delta under presence. Retain all counts. Separately
stress absence with s=.5 while assuming delta=0, b=10, eta_min=8, E=64.

Technical gate: incorrect-decision Wilson upper bound<=.075 in every in-contract
setting; at eta=8,E=64, correct-decision Wilson lower bound>=.80 for both truths
at every b/delta. A pass is conditional simulated feasibility only. Failure does
not permit changing the gate. Interpret whole-case equality when eta<=2*delta*b:
background uncertainty between the two ROIs can still defeat the reference.

## Decision rule

Even if the simulation passes, the real answer is 'not yet established' unless
the source supplies representative cross-color/target localization calibration,
probe specificity/efficiency and a defensible background ROI. Save the minimum
missing measurements as an assay requirement. Check prior work on chromatic
registration, inverse-problem nuisance projections and conditional Poisson tests
before attributing any novelty.
