# When another imaging round cannot resolve a molecular claim

Date: 2026-10-30 (local). Frozen before the codebook enumeration and simulation.
This is a measurement-identifiability audit, not a new decoding or acquisition
heuristic. It is separate from the previous exact review-planning study.

## Smallest counterexample and proposed measurement

Four codewords 1100,0011,1010,0101. Two equal-brightness molecules at exactly
the same modeled position can be genes {0,1} or {2,3}. Both produce 1111 under
additive unsaturated imaging in every round. Their gene-pair co-presence differs.
This is co-presence within an unresolved volume, not physical interaction,
causal coupling or cell-type association. A gene-0-specific independent probe
distinguishes the two worlds if it has positive calibrated sensitivity and
controlled background. Repeating existing rounds cannot distinguish them.

Defeat the proposal with unknown probe efficiency (including zero), variable
background, cross-reactivity and changed geometry. Exact co-location is an
admissible mathematical boundary, not a claim that observed molecules coincide.
Near-coincident unequal-brightness images need separate resolution analysis.

## Fixed source-codebook audit

Use the existing development source codebook from Zenodo 15678603, SHA-256
416db756170566d6e46ed21a8272d41c75114d8c95e384d4f0289df92bcd8c34.
Enumerate every unordered pair with replacement over all 140 codewords, grouping
by its full 16-bit integer SUM, not logical OR. Report both the full codebook
and the restriction to the 134 nonblank genes: blanks are not actual molecules.
Do not choose a code subset to make the result look stronger.

Save every ambiguity class, collision size distribution and the first
lexicographic gene-only collision (if none exists, retain that negative result).
For each ambiguous class, verify an exact probe of either member of a specified
distinct-gene world separates it from all other worlds in that class. This is a
conditional algebraic result, not a claim that such a probe was manufactured.
Verify full integer sums independently of the grouping algorithm.

## Fixed photon-count screen

A target-specific probe has count per exposure Pois(b+eta*z), z in {0,1}.
Unknown b in [.8*b0,1.2*b0], b0 in [2,10]; eta in [.5,2,8]. At n in
[1,4,16,64] exposures, generate 2,000 independent replicates per setting,
seed 2026103002. Test worst positive background .8*b0 and worst null 1.2*b0.
Reject absent target if the exact two-sided 95% Poisson mean lower bound
exceeds 1.2*b0. Conventional comparator uses nominal b0 instead, illustrating
the cost and need for background calibration. Both receive the same photons.

When eta <= .4*b0, construct a positive and negative world with exactly equal
Poisson means. Repeated measurements cannot uniformly identify z over this
parameter class. Explicitly report eta=0 dropout as another impossibility.
The interval-separation threshold is structural, not a finite-sample guarantee.

Technical gate: every saved collision verifies exactly; robust null rejection
Wilson upper bound <=.075 in all settings; for eta=8 and n=64 the robust power
Wilson lower bound >=.80 at both backgrounds. An empty source collision set
fails the proposed source-specific ambiguity claim, without altering the toy
proof. No gate implies novel mathematics or empirical scanner savings.

## Originality and interpretation

Check overlapping-barcode decoding, uniquely decodable mixture codes, sparse
recovery and experimental encoding design at equation/protocol level. A general
overlap limitation, a target-specific probe or a Poisson interval is not a new
invention. The possible next contribution is a calibrated, claim-specific
measurement requirement with verified witnesses, beyond existing decoders.
Source codebook analysis is real design metadata; photon counts and molecule
configurations in this screen are synthetic. Existing raw images are not reused
as molecular ground truth, and no new tissue field is opened in this screen.
