# Reference pilot needed before a real one-image decision

Prepared 2026-10-30; no collaborator, instrument booking or new acquisition exists.
The target remains the previously fixed two-molecule composition claim. The
purpose is to validate its observation contract before any new biological claim.

| Quantity needed | Reference measurement | Current gap |
|---|---|---|
| Gene-color and round geometry | Same features visible in 473/561/650 channels, repeat captures, independent feature identities; retain unmatched/saturated features | Current field supplies 473-nm fiducials only |
| Probe specificity | Independently identified Adamdec1-positive/negative controls, plus competing gene compositions where feasible | Proposed probe is not a validated specific readout |
| Minimum efficiency | Positive references across relevant depth/brightness conditions, calibrated gain/read noise and non-detections | eta_min=8 is simulated |
| Representative background | Preselected target-free regions near the target, repeated at the same exposure, with heterogeneity/bleed-through retained | Area 16 and mismatch delta are assumptions |
| Photon/noise model | Dark frames, gain/linearity checks and saturation limits | Source uint16 intensities are not identified photon counts |
| Acquisition cost | Probe/hybridization preparation, exposure duration, stage travel, usable reference area and feedback latency | One frame does not establish one cheap assay operation |

Keep calibration fields separate from final reference-test fields/specimens.
Never choose a target-free region because it makes the desired claim pass.
Preserve failed matching, probe-negative and out-of-contract cases. If constant
color shifts or common background are contradicted, expand the uncertainty model
and recheck identification before developing another estimator.

For a fixed claim, fixed dose and a 5% total false-decision target, one possible
error ledger allocates 2.5% to all calibration-contract failures jointly and
2.5% to the held-out count decision. The union bound gives <=5% only if the
calibration procedure actually covers the relevant joint contract with at least
97.5% probability and the new count model/independence assumptions hold.
The bead consistency audit does not supply that guarantee. Bounds for unseen
fields require explicit exchangeability/transfer assumptions; a hard maximum
cannot be learned from finite examples without further restrictions.

Before final capture, freeze the specimen/field split, calibrated uncertainty
construction, eligible ROI rule, dose, detector-noise/saturation treatment,
minimum useful claim-resolution rate and reference controls. Obtain sample size
from pilot effect/noise data and a declared precision target; an invented field
count is not a power calculation. Evaluate false decisions and abstentions by
independent specimen/field. Pixels and repeated photons are not biological
replicates. Sequential dose increases need a separate stopping rule; the current
fixed-time construction does not justify repeated peeking.

Advance to a minimum-measurement paper only after the contract is useful, the
independent reference experiment resolves the claim at valid error control and
the equation/protocol contribution survives comparison with overlap-aware
decoding and established experimental design.
