# Pixel-level ambiguity and a published decoder

2026-11-01. Freeze this protocol before running the decoder. This is a focused
measurement-identifiability benchmark, not a general decoder leaderboard.

## Question and exact construction

Use the four previously declared oriented alternatives to Adamdec1 + Adgrg6
in `study_photon_ceiling.cases()`, the source 140-word codebook (including all
six blanks), and the actual 8-round, 2-color acquisition layout. All these
worlds have the same binary additive signature, with at most one emitter in
each active bit. Positions are x=+/-d/2, y=0. Use the saved round/color alias
coefficients to place half the required displacement in each world's
registration error, with opposite signs. Independently render both worlds.

With a common translated PSF and equal brightness, matching the active emitter
centers makes the entire gene-image mean stack equal, not just its centroids.
Independent Poisson pixels with those means have the same observation law at
every exposure. Therefore any forced binary decision has worst-world error
at least 1/2. This is a restricted counterexample; arbitrary unequal brightness,
known registration and external references need not preserve it.

## Implementation and fixed inputs

Run the original, unchanged MIT-licensed ISTDECO implementation from
`axanderssonuu/istdeco`, commit `264c9383a2d73d5669561bff179fd4c359d2e103`.
Preserve its license and hash every executed source file. CPU, one thread,
float32 inputs, 100 multiplicative-update iterations, acc=1, no maximum
suppression before the claim score. No outcome-based parameter search.

Synthetic images: 25x25 pixels, Gaussian PSF sigma=1.2 pixels sampled at pixel
centers and normalized over the patch, background=2 photons/pixel, integrated
signal=200 photons per active emitter/bit. Exposure multiplies both means.
Separate d in {1,3}, exposure in {1,16}, and 12 fixed replicates per cell;
seed 2026110104. At each paired alias input, draw one stack and score the
identical output against both worlds; this pairing verifies the observation-law
argument rather than estimating performance on an uncontrolled population.

Controls receive the same full codebook and pixels:

* ISTDECO forced claim: compare total recovered intensity in the two candidate
  gene pairs, over the entire patch; ties select the first pair. This is a
  predefined claim readout, not the authors' spot-calling/FDR evaluation.
* Known-registration images, independently drawn for each world, with the same
  signal/noise, to test whether the decoder can resolve a recoverable version.
* A ninth round containing a specific Adamdec1 probe, channel 0, with the same
  signal/background. Append that measurement and its codebook bit for ISTDECO.
  Also use the conventional Poisson likelihood ratio for that known probe.
  Both methods receive the extra measurement; it is not exclusive information
  for a proposed method. Specificity and equal brightness are assumptions.
* The ambiguity certificate abstains on gene-only alias inputs. Abstention is
  a valid conclusion about this model, not a new decoding algorithm.

Gates: independently rendered means agree within 1e-9 in every alias cell;
paired forced error is exactly 1/2; decoder known-registration accuracy >=.8
on the d=3, E=16 controls; added-probe claim accuracy >=.95 for both methods
pooled over the full fixed grid. Retain every failed gate. No superiority claim
over ISTDECO follows from an information-theoretic counterexample.

## Measured image diagnostic

Use the already checksum-verified source TIFF-derived middle plane, field 000,
z=14, DOI 10.5281/zenodo.15678603. Fixed 48x48 crops with upper-left corners
(y,x)=(600,600),(600,1100),(1100,600),(1100,1100); all 16 gene bits, no replacement
of empty or difficult crops. Fit nominal fiducial shifts and each of four
previously saved fiducial-tile shift variants, using the SAME fixed crops.
Common preprocessing is max(Gaussian1-Gaussian3,0), then a common full-plane
99.9th-percentile scale per bit from the existing development screen. The
published method receives these prepared images; this is not unprocessed camera
data or a faithful reproduction of its original paper's full workflow.

Report per-gene recovered masses and changes in the dominant gene per crop
between registration variants. No molecule labels or gene identities are known
independently here. Real-image accuracy, FDR, acquisition savings and biological
calibration remain unmeasured. This development diagnostic has no pass gate and
cannot substitute for a held-out, independently annotated image benchmark.
