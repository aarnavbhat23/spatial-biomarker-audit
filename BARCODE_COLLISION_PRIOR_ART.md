# Equation and experimental-protocol comparison — 2026-10-30

Targeted audit only; absence of an exact phrase in search results is not evidence
of priority. The source codebook is an existing measurement design, not ours.

| Closest work | Equation/protocol overlap | Difference and unresolved requirement |
|---|---|---|
| [Ahlswede and Balakirsky, binary adder channel, 1999](https://doi.org/10.1109/18.746834) | The receiver observes the componentwise integer sum of binary codewords; unique recovery requires distinct sums. This is exactly the algebra behind our two-molecule collision. | Unordered worlds within one gene codebook and claim-specific probes are an application, not invention of mixture-code identification. |
| [Untangling overlapping barcodes in spatial genomics data](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/) | Their forward image model adds PSFs; decoding uses spatial/brightness consistency and optimization. They explicitly compare overlap ambiguity across MERFISH and Reed–Solomon codes. | Our exact equal-position/equal-brightness witnesses and complete source-codebook inventory isolate a nonidentifiable boundary. They do not establish failure of that decoder on ordinary separated data or a new overlap problem. |
| [MERFISH with expansion microscopy, 2018](https://doi.org/10.1038/s41598-018-22297-7) | Physical separation and independent single-gene smFISH measurements already address high-density signal overlap. | Proposing an extra gene-specific measurement is established practice. We quantify a conditional background boundary, not experimental superiority. |
| [CISI, 2021](https://doi.org/10.1038/s41587-021-00883-x) | Composite measurements and reconstruction already reduce imaging requirements. | Claim certification would need a demonstrable advantage over reconstruction for the same observation and target, without silently adding a favorable prior. |
| [CIPHER](https://pmc.ncbi.nlm.nih.gov/articles/PMC13252831/) | Joint design of aggregate encodings and downstream inference already accounts for experimental measurability and noise. | A generic proposal to learn better codes is not novel. A distinct calibrated guarantee or lower bound is still required. |

The probe decision is the ordinary exact Poisson test against the largest
allowed null background. The mean-overlap impossibility is elementary. The
comparison with a nominal-background test demonstrates misspecification, not
superiority over strong calibrated inference.

## Candidate remaining contribution

Investigate the minimum additional measurement needed to settle a *specified
spatial claim* under calibrated PSF, brightness, background and missing-molecule
uncertainty. This must go beyond exact co-location: establish near-collision
information limits and test realistic raw-image controls. Compare against
overlap-aware decoding plus conventional uncertainty, full acquisition, and a
simple independent-probe strategy. All costs require measured acquisition and
calibration overhead. This remains an unvalidated direction; no priority claim
or 8–9 novelty rating is earned by the present collision enumeration.
