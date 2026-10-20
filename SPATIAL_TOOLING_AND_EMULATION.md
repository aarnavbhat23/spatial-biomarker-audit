# Incremental spatial acquisition and review emulator

The executable prototype requests observations incrementally, returns only the
requested calls, and records a marker-error-aware interval after each action.
It runs entirely on public development tables and injected errors. Costs and
independent review are simulated. It does not operate an instrument, call a
proprietary predictor or establish real acquisition savings.

## Tools inspected online

| Tool | Useful role | Access and current integration |
|---|---|---|
| [Strand Python SDK](https://docs.strandai.com/sdks/python) | Vendor upload/inference/result interface, including AnnData and OME-Zarr outputs | Official SDK/API inspected. Actual prediction requires credentials and credits; no upload or paid job was started. Predictions remain unmeasured inputs, not a substitute for assay truth. |
| [QuPath](https://qupath.readthedocs.io/en/stable/docs/advanced/exporting_annotations.html) | Inspect tissue and exchange review/acquisition regions through GeoJSON | Source documentation inspected. Exported shapes use full-resolution pixel coordinates. Registration and origin must be checked before using our field geometry on an image. No live QuPath session was tested. |
| [SpatialData](https://spatialdata.scverse.org/en/stable/) | Preserve images, cell tables, coordinate transforms and region annotations together | Suitable future interchange layer; not installed as a new dependency for this prototype. |
| [SimPy 4.1.1](https://simpy.readthedocs.io/en/4.1.1/) | Resource/time-event emulator for scanner and review requests | Used in the observation service. Hypothetical costs are clearly separated from software wall time. |

The vendor also exposes a public prediction viewer on its website. That is useful
for interface exploration but, without paired measured labels and suitable data
access, cannot validate our certificate or the quality of its predictions.
Public SDK availability does not establish access to real model evaluation.

## Observation boundary and cost accounting

`acquisition_emulator.py` provides a local HTTP service with `GET /metadata`,
`POST /field` and `POST /review`. Metadata contains geometry-derived target
weights, field ownership and assumed budgets/costs. It contains no unobserved
A calls or reference A labels. A field request reveals recorded A calls for
that field; a review request is permitted only for an acquired cell and returns
an ideal independent reference A call. Repeated requests reuse the cache.
Requests are serialized; concurrent scanner/reviewer execution is not claimed.

The exact B calls and graph are assumed already acquired. Their full cost is a
common prerequisite excluded from incremental comparisons. Therefore even a
large incremental saving is not an end-to-end imaging saving. Under the stated
error budget and perfect-review assumptions, intervals are deterministic
identification bounds; no stochastic confidence calibration is needed.

## Results and decision

[Protocol](JOINT_REVIEW_PROTOCOL.md), [complete results](results/joint_review_study.json),
[interval and policy implementation](joint_review.py).
Thirty-six patient/corruption/cost/policy scenarios used the first eight occupied
fields from each of three development panels. All decision checkpoints covered
the source-reference targets under the injected-error budget.

| Scenario | Joint policy | Acquire-all then influence review | Incremental cost reduction |
|---|---|---|---:|
| CRC11 concentrated errors, review cost .05 / .5 | 2 fields, 4 reviews | 8 fields, 4 reviews | 73.3% / 61.1% |
| CRC19 concentrated errors, review cost .05 / .5 | 2 fields, 3 reviews | 8 fields, 3 reviews | 73.7% / 64.1% |
| CRC11/CRC19 random errors | Same fields/reviews as influence control | Same | 0% |
| Both CRC02 cases | No incremental observations needed | Same | Not defined |

CRC02's selected fields already have a narrow enough bound from the exact-B
prerequisite. Those zero-observation cases are retained and do not count as
policy wins. The policy improves only two of the six patient/corruption tables,
so it **fails the prespecified four-of-six continuation criterion**. Random
review is a much weaker control here; beating it is insufficient.

A CRC19 concentrated-error replay through localhost HTTP selected exactly the
same actions as the bulk in-process observation interface. Exhaustive small
binary populations verify that partial-observation bounds are sharp, and that
new-field worst-case ambiguity is exact. Tests also check hidden-label metadata,
unacquired-review rejection and zero-cost cache reuse. This is a correctness
and feasibility result for the software, not external biological confirmation.

The initial score implementation omitted field feedback from its ranking
denominator while still charging it in reported costs. The corrected replay
uses the full field cost and gives the same decisions/outcomes. The initial
results remain in `results/joint_review_initial_score_audit.json`.

## Theory behind the prototype

With B fixed, define w_ju as the B-positive neighbor count for A cell u in
claim j and M_j as its directed pair count. Let R be reviewed cells, S acquired
but unreviewed cells and U unseen cells. If e discrepancies have been found,
b'=b-e errors can remain. Let K_j be the weighted true calls on R plus weighted
recorded positives on S. The sharp interval is

    [K_j - top_b'(w_j on S with recorded A=1),
     K_j + top_b'(w_j on S with recorded A=0) + sum(w_j on U)] / M_j.

Unknown recorded calls can equal any chosen true calls without spending error
budget. Flipping the largest acquired weights attains each endpoint. This
extends the prior census diagnostic to partial observations and exact reviews.

For an unacquired field with weights F, the exact worst post-acquisition
ambiguity is the maximum over a,c in 0,...,b', a+c<=|F| of

    top_(b'-a)(existing positive weights)
    + top_(b'-c)(existing negative weights) + top_(a+c)(F).

Allocate the largest fresh weights between the two signs. This gives a closed
finite enumeration for field lookahead, with no hidden A labels. The policy
greedily chooses the largest worst-case width decrease per hypothetical cost.
A greedy lookahead policy and elementary weighted-selection identities are
not claimed as a new general inference framework.

## Reproduce and next research gate

Install `requirements.txt`, then run:

```bash
python -m unittest discover -s tests -v
python study_joint_review.py
python serve_acquisition.py --patient CRC19 --corruption adversarial
```

The last command prints a localhost URL for manual or scripted incremental
requests; stop it with Ctrl-C. No external platform account is needed.

Keep the interval/oracle infrastructure. Stop this particular greedy policy
as a broad efficiency claim. Next compare batched/multistep action planning
against exact small-population dynamic programs, with both markers uncertain
and fallible reviews. Only proceed to a fresh-cohort efficiency claim after
that gate, followed by measured scanner/reviewer costs. The existing lab pilot
remains necessary; online emulation cannot replace it.
