# Review packet: when does a spatial claim need another measurement?

For a short introduction, start with the [research summary](RESEARCH_SUMMARY.md).
This file contains the more detailed assumptions and open questions.

Status, 2026-11-02: outside feedback has been requested; no technical assessment
has been received yet. A new decoder and a biologically calibrated decision
guarantee have not been established.

## Proposed focused contribution

Given a stated barcode-composition claim, find observationally equivalent
alternatives under explicit registration constraints, then identify what
independent reference or probe information can resolve that claim. The intended
paper is about measurement requirements. It does not require inventing another
image-processing architecture.

Read these core artifacts:

1. [Full-image theorem and limits](PIXEL_INFORMATION_THEORY.md).
2. [Actual published decoder comparison](PIXEL_DECODER_RESULTS.md), complete
   source hashes and saved mass readouts.
   The [closest published ADCG/syndrome pipeline](SYNDROME_PIXEL_RESULTS.md)
   passes its noiseless controls; its measured convergence follow-up fails to
   produce a completed decode. Both outcomes are retained.
3. [Conditional finite-localization theorem](LOCALIZATION_INFORMATION_THEORY.md)
   and [external calibration benchmark](SUBMISSION_RESEARCH_RESULTS.md).
4. `verify_pixel_decoder.py`, `verify_submission_research.py` and
   `verify_fiducial_measurement.py`, which independently check the recorded
   geometry, held-out predictions, decision calculations and LP certificates.

## Proof and priority questions

* Does the full-pixel equality construction preserve every claimed image law,
  and is its adaptive-exposure extension limited to the same observation
  operators? How should unequal, bit-dependent brightness restrict its use?
* In the Gaussian theorem, are the complete nuisance subspace, whitening and
  least-favorable-pair argument sufficient for the stated minimax error?
* Is the infinite-photon limit taken with fiducial information fixed, and are
  the fixed versus optimized nuisance parameters clearly distinguished?
* Which established identifiability, imaging or coding result already implies
  the specific calibrated measurement boundary? Does the application add enough
  beyond that result to warrant a methods paper?

[White et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/) already
addresses overlap under imperfect registration and uses spatial variance,
brightness variance and global integer programming. [ISTDECO](https://github.com/axanderssonuu/istdeco)
already combines spatial and spectral deconvolution. Matched-subspace testing
and additive-code uniqueness are established. Neither a generic 50% ambiguity
bound nor a new audit of one codebook alone justifies a broad novelty claim.

## Calibration: what is observed and what is missing

The author-released cross-channel table contains 680 selected matches from one
field. Keeping all appearances of a physical feature together gives 65 training,
32 calibration and 34 test groups in the frozen primary split. The calibrated
affine radius is 2.963 pixels, covering 34/34 test groups but **failing** the
<=1-pixel precision gate. The large calibration feature remains in the data.
No split was reassigned to obtain a narrower result.

These source-selected matches are not independently adjudicated molecule
identities. They omit unmatched detections and do not establish tissue gene
accuracy. Their residual distribution cannot supply calibrated Gaussian
variances on the different MERFISH instrument/field. Separate smFISH population
count tables do not resolve this molecule-level gap either.

For an empirical decision guarantee, collect or obtain held-out raw frames with
matched reference objects visible at the gene wavelengths, independent spot/
composition labels, detector background and localization characterization,
dropout/specificity controls, specimen identifiers, and a repeat or extra probe.
Split whole independent specimens before fitting/error calibration, preserve
unmatched detections, and keep the target extra-probe result blinded until the
decision rule is fixed. A collaborator must check that the reference objects
have the same relevant optical path as the target molecules.

The current specific-probe SNR is assumed, not measured. More gene photons may
help outside the exact alias model, and calibrated registration can remove its
particular ambiguities. The paper must retain these defeat conditions. If no
appropriate labeled cohort is available, restrict the manuscript to theory,
controlled pixels and selected-feature registration validation; omit claims of
biological identification, practical assay savings or general tissue accuracy.
