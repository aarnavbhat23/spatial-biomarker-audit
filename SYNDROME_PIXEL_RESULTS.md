# Closest published pipeline: controls pass; measured comparison is incomplete

2026-11-01. Frozen protocols a380e64 (initial run) and 4c10382 (convergence
follow-up). Execute the original SeqFISH_ADCG and SeqFISHSyndromeDecoding
packages at their pinned commits with GLPK; no upstream implementation was
rewritten or substituted. Exact dependencies and prepared inputs are in
`comparators/white_pixel`. Their academic/noncommercial licenses are preserved.

All 12 controlled stacks yield eight fitted spots, one in each active bit.
The decoder recovers the exact true pair in **8/8 noiseless, correctly registered
controls**. At the shared alias inputs it chooses one gene, with remaining dots
unassigned; a forced candidate-pair readout still has exactly 50% paired error.
This does not mean that the authors' decoder itself asserts either complete pair
with high confidence. It can leave observations unexplained.

An independent Python verifier enumerates every feasible set of up to two
disjoint codewords on each eight-dot network. Every published decoding output
attains the global optimum of its fixed spatial/brightness/free-dot objective,
within numerical tolerance. The actual fitter and decoder, rather than just
that objective, were executed. The controls establish implementation capability
under this constructed input; they do not establish population performance.

## Measured development crop and retained computational failures

The original fixed 48x48 crop gives 640 fitted spots, hitting max_iters=40 in
**all 16 bits**, and returns 39 decoded barcodes. Those counts are incomplete
under the fitting budget and have no independent molecular truth. They cannot
support an accuracy, FDR or biological claim.

The separately frozen follow-up increases only the measured max_iters to 200,
matching the authors' example cap. It did not yield a completed measured output
before the run was stopped. The log contains ten fitting-cap messages and an
upstream attempt to read optimizer primals with zero solutions available. The
precise solver termination status was not captured, so this is **not evidence
of mathematical infeasibility**. The convergence gate does not pass; no fitted
point table or decoded counts from that incomplete run are invented.

The sanitized log, unchanged-input comparison and failed outcome are saved and
verified. The runner now checkpoints fitted points after each bit, so future
decoder failures need not discard the expensive fit. This reporting/checkpoint
change does not retune the frozen fitter or decoder. The follow-up is a
post-outcome computational diagnostic, not a confirmation experiment.

## What remains before a fair empirical performance claim

The preprocessing, PSF/background model, fitting convergence and solver status
must be validated on a matched, independently labelled image benchmark. Raising
a cap alone cannot calibrate them. The 48x48 development crop is a transfer
diagnostic, not a faithful reproduction of the authors' original assay workflow.
Both the initial truncated result and the extended failure remain in the record.

For the narrower theory claim, the [full-image ambiguity](PIXEL_INFORMATION_THEORY.md)
and [noisy ISTDECO comparison](PIXEL_DECODER_RESULTS.md) remain valid in their
declared settings. Outside mathematical/priority feedback and gene-level
calibration remain open in the [review packet](INDEPENDENT_REVIEW_PACKET.md).
