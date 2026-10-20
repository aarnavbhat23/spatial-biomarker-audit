# Batched physical-field screen

Exploratory development screen; no new confirmation patients. Fixed patients
CRC02, CRC11 and CRC19, camera fields 512/1024/2048 native pixels, batch sizes
1/4/16/64, five repeats, target width .005, maximum 2048 randomized queries.
Seed 20261105 is an arbitrary reproducibility identifier, not a study date.

At a batch boundary freeze known positives K, unresolved edge count U and group
counts m_g. Draw b independent group indices with replacement using m_g/U.
Acquire the union of their physical fields. Each sampled value uses its group's
PRE-BATCH unresolved edges, including repeated draws. Only after constructing
all values may the global cell/edge cache update. Bet sizes/support are frozen
within the batch; inspect the interval only at completed batch boundaries.

Compare physical fields, feedback rounds and query count separately. Report
cost F + h R for h=0,1,10,100 field-equivalents per feedback round; these are
sensitivity scenarios, not instrument measurements. Use the same boundary
inspection rule for every batch size, including b=1.

Include a stronger nonstatistical control: greedily acquire the unobserved field
incident to the most unresolved edges and stop when sharp deterministic bounds
are within .005. It sees only geometry and already observed labels. A graph
certificate can invalidate the need for elaborate statistical acquisition.

Success means a useful cost frontier survives that control, without losing
validity. Batch testing and betting confidence sequences are existing methods;
combining them is not by itself a high-novelty contribution. No clinical,
whole-slide, staining-cost or physical-instrument claim follows from this screen.

Amendment after the first screen: the greedy control must also acquire fields in
batches of 1/4/16/64. Before revealing labels, select the highest incident-edge
scores for the whole batch, breaking ties by field index. Otherwise comparing
feedback-round costs to an obligatorily sequential control would be unfair.
Retain the statistical sampler's existing seeds and settings unchanged.
