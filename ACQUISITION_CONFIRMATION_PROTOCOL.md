# Frozen region-scale confirmation

Freeze this protocol and executable estimator before downloading the four new
patients' marker tables. Development panels were CRC02, CRC11 and CRC19. They
selected the 128-native-pixel core; 256 and 512 are secondary, not alternatives
selected on new outcomes. The comparator is the original 2048-pixel core with
the **same** probability-proportional-to-pair-count (PPS) sampler, confidence
sequence, deterministic completion bounds and stopping rule.

## New patients and target

CRC05, CRC06, CRC07 and CRC09 are new to this analysis. Their mapping to P37_S33,
P37_S34, P37_S35 and P37_S37 is verified in the source authors' file_list.csv.
They are four distinct patients, not technical repeat sections. They are a
convenience subset with accessible processed tables, not a representative cohort.
Archive entries are pinned separately in acquisition_sources.json.

Use a coordinate-only random panel of 32 occupied 2048-pixel cores per patient,
with >=100 source cells, seed 20261023 + patient number. Exclude cores whose
80-pixel halo leaves the conservative coordinate bounding rectangle. Retain
every source cell in the panel's halos. Target all directed, nonself pairs with
an anchor in the selected cores and Euclidean distance <=80 native pixels.
Target readout is mean Pan-CK_pos(anchor) times CD3e_pos(neighbor).

## Observation and estimation

For each scale, partition anchor cells into fixed square cores. A core's measured
ROI is expanded by 80 pixels in all directions, ensuring both endpoints of every
owned pair can be measured. Let M_t be its known geometric pair count and Y_t its
marker-positive pair total. Independently draw region indices with replacement,
probability M_t/sum(M). Then X=Y_t/M_t is in [0,1] with expectation equal to the
target. Repeated draws reuse measurements but remain IID random indices from the
fixed population. This is classical PPS/Hansen-Hurwitz sampling.

Invert the two fixed mixtures of nonnegative linear betting martingales in
study_region_scale.py. Their union error is <=0.05 simultaneously over time.
Intersect with deterministic lower/upper bounds from already measured owned
totals. Inspect at draws 8,16,32,64,128,256,512,1024,2048 and stop when width <=0.04.
Otherwise stop at 2048 and report failure to meet precision. Coverage conditions
on the exact fixed graph/labels and does not address staining or segmentation
error, model predictions, or new-patient generalization.

## Primary result and costs

100 acquisition repetitions per patient, seed 20261024 + patient number. Primary:
128-pixel cores achieve width <=0.04 in >=95% of repetitions for every patient,
with lower mean physical union area than the 2048-pixel baseline for each patient.
Coverage-at-stopping is reported, not used to select a confidence method. Finite
simulation coverage is descriptive; validity comes from the stated martingale
argument, conditional on its assumptions.

Report both exact ROI-union area and summed acquired rectangle area. Report
distinct region counts. Compute area-equivalent per-region overhead at which
each method becomes cheaper: C = area + overhead * number_of_regions. This is a
sensitivity analysis, not a claim of real instrument savings. Both approaches
condition on the same fluorescence-derived graph. More overlapping observations
than the owned totals are not exploited.

## Novelty boundary

PPS, bounded-mean betting, optional stopping, geometric halos, and cost sensitivity
are established ideas. Confirmation success validates a candidate observation
protocol; it does not make this a new statistical method or a 9/10 novelty result.
The intended research opening is a theory of optimal acquisition geometry and
information limits for marker-defined spatial relationships. Existing graph
sampling and stereology results must be distinguished before any priority claim.
