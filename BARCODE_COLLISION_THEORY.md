# Identical full barcode images and what can separate them

## 1. Smallest explicit witness

Let four genes have binary codes c0=1100, c1=0011, c2=1010, c3=0101.
World W has molecules {0,1}; world V has {2,3}. Put both molecules at the same
modeled position, with the same brightness, point-spread function and background.
Then c0+c1=c2+c3=1111, using integer addition. Under additive unsaturated
imaging, every pixel in every barcode round has the same expected intensity.
With the same conditional noise model, the entire raw-data laws agree.

The claim 'genes 0 and 1 are both present in this unresolved volume' differs.
It is not a molecular-contact or tissue-interaction claim. Adaptive selection,
repetition or reordering of these existing readouts cannot distinguish W and V:
conditional transcript laws agree inductively. This does not rule out a new
probe, separation by expansion, or extra spatial resolution when positions differ.

More generally, a target T(w) is identified by an exact observation F(w) if and
only if T is constant on every feasible fiber {w:F(w)=y}. For linear additive
data F(w)=Cw, a differing-target pair in the null space is a witness of failure.
This is elementary inverse-problem theory, not a newly invented principle.

## 2. What one new probe identifies

In an ambiguity class of two-molecule sums from distinct binary codewords, two
different worlds cannot share a gene: c_i+c_j=c_i+c_k implies c_j=c_k. A
repeated-code pair 2c_i cannot collide with another pair of distinct binary
codes either. Thus within an exactly known full-sum class, an exact probe of
either member of a chosen world separates that world from every alternative.

This identifies that two-gene composition within the restricted class. It does
not recover arbitrary mixtures, unknown brightness or a larger unknown molecule
count. For resolving every member of a class of size m, one perfect single-gene
probe identifies its world on a positive response but leaves m-1 on a negative
response. With only such binary membership queries, worst-case m-1 probes are
necessary and sufficient. A pooled probe could bisect the candidate worlds, but
its signal and cross-reactivity contract would need separate justification.

## 3. Defeating the probe; a sharp background boundary

Suppose n independent exposures yield total photons K~Pois(n*(b+eta*z)), where
z is target presence. Calibration supplies b in [bL,bU] and eta>=etaL>0.
The negative mean set is [bL,bU]; the positive set starts at bL+etaL. They are
strictly separated exactly when etaL>bU-bL. If etaL<=bU-bL, choose z=1,b=bL,
eta=etaL and z=0,b=bL+etaL. Their complete count laws agree for every n.
Thus there is no uniformly consistent decision over this class. Known positive
efficiency alone is insufficient. Admissible zero efficiency or unconstrained
cross-reactivity also recreates the ambiguity. Independent background control
could change the class and must be counted as an additional measurement.

A standard exact Poisson mean interval supplies a fixed-time decision rule:
declare presence only if its lower endpoint exceeds bU. With a two-sided 95%
interval, the upper-tail false-positive probability is at most .025 under
the bounded-background null. This is a conventional test, not a new statistic.
An actual imaging protocol must address bleaching, exposure dependence and
background correlation rather than assuming repeated IID photon observations.

## 4. Scope and next research boundary

The source-codebook audit evaluates design metadata under exactly unresolved,
equal-brightness molecules. It does not estimate the prevalence of collisions
in the raw field. Distinct positions or gene-specific brightness can break exact
equality, while finite photons and uncertain PSFs can still make discrimination
hard. Quantifying that near-collision regime is necessary before imaging claims.

The candidate future contribution is a claim-specific measurement requirement
with calibrated photon/registration uncertainty and explicit counterexamples.
General mixture coding, overlap-aware decoding and experimental encoding design
already exist; see [the prior-work comparison](BARCODE_COLLISION_PRIOR_ART.md).
