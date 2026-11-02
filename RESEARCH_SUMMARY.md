# Which extra measurement would settle an ambiguous gene call?

I'm studying a simple question in spatial transcriptomics: **when do the images
contain enough information to distinguish two molecular explanations, and what
should we measure when they don't?**

The aim is to help decide whether to improve alignment, collect more light, or
use an extra probe. So far, I have a mathematical example and controlled tests.
I haven't established a reliable decision rule for tissue experiments.

## A concrete example

In barcode imaging, a gene is identified by which imaging rounds and colors
light up. Two nearby molecules can make that pattern hard to interpret.

One example from the public codebook is:

| Explanation A | Explanation B |
|---|---|
| Adamdec1 + Adgrg6 | Ccr2 + Pik3r1 |

These pairs have the same combined on/off pattern. I constructed shifts between
rounds and colors that also make their complete pixel images indistinguishable.
The construction assumes equal molecule brightness, the same optical response
and background, and uncertainty about those shifts. Under the stated noise
model, the whole image distributions agree, rather than just one noisy image.

For these cases, collecting more light from the same barcode channels cannot
tell the pairs apart. A method forced to choose between equally likely
explanations has 50% average error. It can instead report that the call is
unresolved. Known alignment removes these particular ambiguities, and an ideal
extra probe specific to Adamdec1 distinguishes the two pairs.

## What I've tested

I ran the original ISTDECO code with the full 140-entry codebook,
including blanks. The fixed synthetic grid varies separation and exposure and
includes four alternative pair arrangements.

| Test | Result |
|---|---|
| Constructed indistinguishable pairs | 50% error when forced to choose a pair |
| Correctly aligned controls | 384/384 correct |
| Same cases with an ideal extra specific probe | 384/384 correct with both ISTDECO and a simple Poisson test |

These are controlled synthetic cases, not independent tissue measurements. The
score compares the total recovered signal for each candidate pair; it is not
ISTDECO's original spot-calling or false-discovery procedure.

I also ran the original ADCG/syndrome pipeline from White and colleagues. It
recovers the true pair in all eight noiseless, correctly aligned controls. On
the ambiguous inputs, it can leave spots unassigned, rather than confidently
asserting the wrong complete pair.

## What the real data show

Four fixed crops from one public MERFISH field were fairly stable under the
tested alignment changes: their dominant genes stayed the same, and recovered
signal changed by at most 2.19–2.57% across crops. This does not show a large
practical decoding failure. The ADCG/syndrome comparison on a measured crop is
incomplete because the extended fitting run did not finish successfully.

A separate registration study uses 680 author-selected cross-channel matches,
holding out whole reference features. Its uncertainty radius covers all 34 test
features, but the radius is 2.963 pixels and misses the planned one-pixel precision
target. These selected matches don't provide independent gene identities, and
their errors cannot be transferred to the different MERFISH field as calibration.

The extra probe's brightness and specificity are still assumptions. Different
molecule brightness or better alignment can break the constructed ambiguity.
I haven't measured how often these cases occur in tissue or whether an extra
probe would save acquisition time.

## Where I could use your feedback

[White et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/) already address
overlapping barcodes and imperfect registration.
[ISTDECO](https://github.com/axanderssonuu/istdeco) already combines spatial and
spectral deconvolution. Indistinguishable-distribution arguments are established
statistics. The possible contribution is a result that says which measurement
is necessary to resolve a particular gene-pair claim under measured alignment
uncertainty. That calibration step is still unfinished.

**Does this look like a useful research direction, and which prior result or
comparison would most decisively test it?**

For the details: [image-level proof](PIXEL_INFORMATION_THEORY.md),
[ISTDECO results and reproduction commands](PIXEL_DECODER_RESULTS.md),
[ADCG/syndrome results](SYNDROME_PIXEL_RESULTS.md), and
[registration study](SUBMISSION_RESEARCH_RESULTS.md).
