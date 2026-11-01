# Closest published spots-first pipeline: fixed image controls

2026-11-01. Additional comparator, frozen after the ISTDECO results and before
any ADCG/syndrome outputs. This is a development check, not fresh confirmation
of a newly selected method. The full noisy ISTDECO grid stays primary.

Execute unchanged author packages SeqFISH_ADCG.jl at
f0cb871d2f512bee0cb6dfb03e02089fec225241 and SeqFISHSyndromeDecoding.jl at
b122a9189ca543315e45b2a389f4857709e78eab, with free GLPK. Record Julia/dependency
manifest and package tree hashes. Honor the packages' upstream licenses.

Use the four existing oriented ambiguity pairs, d=3, E=16, sigma=1.2 and
background=32. Supply expected pixel means, without Poisson noise, for the
shared alias input and both known-registration controls: 12 full 16-bit stacks.
This isolates implementation correctness and information ambiguity; it cannot
estimate finite-noise performance. All 140 codewords remain available.

Fit each bit image using the package's fit_tile API: width bounds [1,1.4],
noise_mean=32, supplied tau=1e8 (the package replaces it internally), minimum
loss improvement=1, min_weight=1, max_iters=20, max_cd_iters=20. Preserve every
fit, including extra dots; do not replace fits by true locations. Coordinates
are converted from the package's 1-based to source 0-based pixel convention.

Represent each source bit as one binary coding block, pseudocolor 1. Derive a
parity-check matrix by exact GF(2) nullspace of the source codeword span and
verify H*c=0 for all 140 words. Configure zeros_probed=false, no drops,
xy_search_radius=4, z radius=0, lateral variance coefficient=7, log-brightness
variance coefficient=3.2, sigma/z variance coefficients=0, free-dot cost=1.
These are fixed API settings based on the authors' examples, not an optimized
reproduction of their complete experiment.

Report fitted point counts, decoded gene counts and candidate-pair claim scores
(sum first-pair counts minus sum alternative counts; ties choose first).
The same alias output is assessed against both worlds. Gates: registered
controls recover the exact two true genes with no extra genes in >=.8 of the
eight controls; paired alias forced error=1/2. Failed/empty decoder outputs are
reported explicitly, not silently replaced by known answers.

Also process the first fixed measured development crop (600,600), 48x48, nominal
registration, using the same bandpass/scales from PIXEL_DECODER_PROTOCOL. Multiply
the prepared floating intensities by 1000 before fitting. Use noise_mean=0,
min_weight=20, loss improvement=1, max_iters=40; other settings are unchanged.
This crop has no independent gene truth. Report output and point count only;
there is no accuracy gate, FDR claim or saved-acquisition claim.
