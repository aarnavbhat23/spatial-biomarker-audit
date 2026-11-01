# Spatial submission path — 2026-11-01

Progress: [the first execution round](SUBMISSION_RESEARCH_RESULTS.md) is complete.
The reference benchmark fails its primary radius gate; the conditional Gaussian
information result passes its computation checks. Raw-image/decoder comparison
and target-matched independent confirmation remain open. The assessment below
records the pre-round starting point; the new focused-project assessment is
approximately 7.5 overall, with substantial uncertainty and incomplete submission gates.

Decision: prioritize one focused spatial methods paper. The present project is
approximately 7/10 overall, with roughly one-point judgment uncertainty. A
defensible 7.5-8 target is plausible if the gates below pass; no rating increase
or publication probability follows from this plan.

Working question: **When can optical calibration settle a gene-composition claim,
and when is another molecular measurement necessary?**

Working title: **Calibration-dependent limits of barcode identification in spatial imaging**.

## Contribution and evidence ledger

| Paper claim | Available evidence | What must change before submission |
|---|---|---|
| Shared round/color calibration changes which compositions can explain an image | 5,274 audited oriented comparisons; actual acquisition layout; selected-bead conditional sets; 16 LP certificates | Establish useful finite-noise behavior under a defensible calibration procedure, rather than only exact mean-image overlap |
| Some claim ambiguities survive arbitrary additional photons | Exact counterexamples; target/background mismatch boundary | Extend beyond fixed equal-brightness two-source geometry or justify that restriction as the complete paper scope; obtain independent proof/priority review |
| An additional specific image can settle the fixed claim under a contract | Joint-count simulation passes strong-signal/error gates; original conditional test failure retained | Test independently calibrated geometry/efficiency/background, and retain abstentions and model failures |
| There is a practical advantage to claim-specific measurement | Not established | Compare a published overlap decoder plus uncertainty, a simple extra-probe strategy and full acquisition at matched information or explicit cost; a conventional control winning is a valid outcome |

Overlap decoding and additional smFISH references already exist. The candidate
contribution must be a useful calibration-conditioned requirement for a specified
claim. A confidence rectangle, another codebook inventory or a standard LP is
insufficient by itself. If priority review finds the same result/protocol, narrow
the work to a benchmark/resource only if it offers substantial independent reuse
value; changing the paper label does not create that value.

## A newly verified public-data lead

The authors of [Untangling overlapping barcodes](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/)
provide a [repository](https://github.com/CaiGroup/UntanglingBarcodes) with a
cross-channel reference-match table and separate smFISH gene-by-cell counts.
Metadata inspection confirms **680 rows** in
`figure_generating_notebooks/Fig2F-H/alignment/cross_channel_matched_reference_dots_pos0.csv`.
It includes raw reference/comparison coordinates, author-aligned coordinates and
channel fields. Four duplicated header names require explicit parsing. The
separate [fiducial-matching example](https://github.com/CaiGroup/seqfish_fm_match)
has **97 initial reference detections**. Exact tree/blob identifiers, available
tables and limitations are saved in [the metadata audit](results/submission_source_feasibility.json).

This is a usable feasibility lead, not completed independent calibration. Source
matches were already selected; author-aligned coordinates are fitted outcomes.
The smFISH count reference is not matched molecule-by-molecule truth. The
[documentation](https://caigroup.github.io/UntanglingBarcodes/input_files.html)
describes additional preprocessed images/fiducial fits, but their exact current
release paths must be checked against the live repository. The verified tables
are the basis for the next audit; do not assume every documented raw file is
present. Source software/data terms and attribution must be retained if reused.

## Ordered work, with explicit stop decisions

1. **Independent reference feasibility.** Resolve table semantics, coordinate
   units, source selection, reference identities, channels and repeated views.
   Inspect available positions/specimens. Split complete physical reference
   features and acquisitions, preserving all repeats in one split. If identities
   or acquisition groups cannot be recovered, report a limited consistency audit.
   A single position cannot become an independent-specimen cohort.
2. **Calibration benchmark.** Freeze the split, uncertainty construction, field
   eligibility and error/precision criteria before scoring residual outcomes.
   Fit translation and affine controls from raw coordinates. Refit without the
   held-out features; author-aligned coordinates cannot establish held-out
   accuracy. Compare calibration coverage AND width. Account for selection,
   local deformation and channel mismatch. If intervals are invalid or vacuous,
   stop the one-image guarantee claim for that assay.
3. **One substantial theory/measurement result.** State an identified-set or
   finite-photon requirement that couples calibration to the fixed claim, with
   both ambiguity witnesses and a useful sufficient design. Standard testing
   lower bounds may be tools, but the application result must distinguish itself
   at the equation/protocol level. Generalizing a proof without useful examples
   does not meet this gate.
4. **Frozen comparison.** Use realistic unequal brightness, background,
   localization, dropout and PSF variation. Include a published overlap-aware
   decoder, strong calibrated inference, the simple extra-probe control and full
   acquisition. Give reconstruction methods the same observations and an explicit
   decision/abstention rule. Evaluate false decisions, useful resolutions and
   photons/frames/chemistry costs separately. Freeze numerical tolerances from
   the target and calibration pilot, not the desired score.
5. **Independent confirmation.** Confirm on a new acquisition/study appropriate
   to the claim. Public multicolor beads can validate geometric components;
   a real one-probe composition guarantee additionally needs specific-probe
   positive/negative references, measured efficiency and background/noise controls.
   If those are unavailable, restrict the paper to conditional theory plus
   benchmark evidence and explicitly omit biological validation/cost claims.

The existing 99.15% strong-signal simulation result is development evidence,
not the independent confirmation. A fresh random seed in the same model does not
close this gap. Lab work is necessary for new experimental claims, but a rigorous
conditional methods/theory paper can be evaluated without claiming new biology.

## Submission gate and paper shape

Submit only after one contribution survives the closest-prior-work comparison,
the proof is checked independently, the frozen comparison answers a useful
question, and the manuscript's headline is supported by its actual validation.
Report all prespecified failures. A novel decoder is not required if a
substantive, original measurement-limit result is established.

Target a five-figure story: motivating ambiguity; calibration-dependent boundary;
finite-noise prediction; strong-control comparison; independent reference test.
Build the claim/figure skeleton now, but write the full manuscript after those
figures are supported. Keep unrelated graph, kernel and acquisition variants in
the research archive. The submission package needs a portable reproduction,
source/split manifest, immutable result tables, proof supplement, limitations and
all required source/tool disclosures for the chosen venue.

The numerical grade is a judgment, not an acceptance threshold. A conditional
paper with narrower claims can be strong; unsupported claims cannot be repaired
by stronger wording. Venue selection and author feedback follow the evidence
ledger, before submission. No manuscript is submitted by this plan.
