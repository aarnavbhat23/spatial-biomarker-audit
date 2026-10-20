# Instrument timing and feedback feasibility pilot

Status: planned, no collaborating lab or scanner logs available. No measured
time saving is claimed. This document and the blank
[timing table](pilot_timing_template.csv) are ready for a prospective partner.

## Question and first gate

Can a lab acquire another selected batch, reuse already acquired fields and
return reliable cell calls quickly enough for adaptive spatial certification?
The target is end-to-end time to a valid certificate, including feedback and
independent marker review. It is not nominal scanner throughput.

First request a short technical session using an existing nonclinical research
slide. Record instrument model/generation, software version, objective,
pixel size, channels, exposure, field dimensions and acquisition interface.
Check whether incremental ROIs can be added without a full rescan, whether
cached images/cell identifiers are stable, and whether per-batch output can
be processed before the next choice. Save actual logs and interface behavior.
Do not assume the simulated 1024-pixel square is a physical camera field.

The [manufacturer's workflow description](https://rarecyte.com/webinar-transcript-2026.06.02-Inside-Orion-HT.html)
supports remote ROI planning and scan jobs, but describes analysis of a
completed slide alongside scanning the next. It does not establish incremental
cell-table feedback or its latency. This newer instrument generation also
cannot be assumed identical to the instrument that generated the public cohort.

**Stop the online-hardware claim if these capabilities are absent.** The
fallback is offline certification or a fixed acquisition plan, evaluated as
such. A manual loop is feasible only if its operator time is explicitly charged.

## Timing design after the interface gate

Reserve a bounded first session, initially up to four hours. Use one stable
stained research slide and fixed settings, with two layouts: contiguous and
dispersed fields. Randomize job order within blocks. Test 1/4/16 unique-field
batches with two repeats and cold versus warm cache states; extend to 64 and
additional slides only if practical. A partial pilot is reported as partial,
not as a completed factorial design. Repeated scans may change fluorescence;
record scan order and signal drift separately from time.

For each job record monotonic start/end timestamps, request/first image/last
image/calls-ready/next-request times; unique new versus reused fields; stage
travel, focus, imaging, write, segmentation, phenotyping and transfer times
where separable; bytes, operator time, queue delays and failures. Overlapping
components must not be added twice. Keep staining, prescan, setup and full-slide
segmentation fixed costs separate from incremental imaging costs. If a graph
requires a complete measured prescan, charge that observation too.

Replay the frozen envelope, mixture, focused-reuse and deterministic greedy
policies against identical allowable fields and settings. Use paired layouts
and slides. An initial timing model C=F*t_field+R*t_feedback is only a candidate;
check held-out residuals for batch size, layout, caching and autofocus. If
nonadditive, replay job-level costs or fit a richer transparent model. Report
uncertainty across independent jobs/slides; fields are not biological replicates.

## Current break-even quantities, without invented minutes

At the hypothetical geometry already tested, envelope uses fewer fields than
greedy but more feedback rounds. It wins C=F+hR for h below approximately
4.69/26.92/20.62 field equivalents per round on CRC17/CRC18/CRC20, respectively.
These are panel-specific boundaries, not scanner measurements. Obtain h from
the actual feedback and incremental-field timing distributions. If the measured
range straddles a boundary, report an unresolved comparison. Include common
fixed costs when reporting percentage end-to-end savings.

## Marker-review pilot and follow-on decision

Keep timing feasibility separate from estimating biological error. Review a
blinded random cell sample plus a separately labeled influence-enriched sample,
using a qualified independent review procedure and resolving disagreements.
The enriched sample cannot estimate a global error fraction without its
sampling design. Record review time, failed calls and reviewer agreement.
Even under independent sampling and perfect adjudication, zero errors among
300 cells only gives an approximately 1% one-sided 95% error upper bound;
supporting .1% needs roughly 3,000 zero-error reviews. Finite-population and
stratified sampling require their appropriate bounds. Do not infer a tiny
error budget from a convenience sample or agreement between two related models.

Proceed to a frozen multi-slide efficiency study only when the interface works,
costs can be measured, and the chosen marker-error model permits the requested
precision. Prespecify endpoint (wall-clock time and total resource cost at valid
width), independent slide count, controls, failure handling and success margin
after the pilot, before evaluating new slides. A speedup without a valid
marker-aware certificate fails the scientific gate.
