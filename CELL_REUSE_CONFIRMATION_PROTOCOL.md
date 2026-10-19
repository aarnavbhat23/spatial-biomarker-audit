# Frozen confirmation: tight precision with cell-level reuse

Freeze the current cell-reuse implementation and this protocol before downloading
marker tables for CRC13, CRC14, CRC15 and CRC16. The author-hosted object names map
these patients to P37_S45, S46, S47 and S48 respectively. Record source hashes and
the protocol commit. A local/GitHub timestamp is not external preregistration.

Use coordinate-only selection of 32 occupied 2048-pixel cores per patient, as in
`confirm_region_scale.prepare_panel`. Retain all source cells within the 80-pixel
halos and verify 100 random neighborhoods against the full source table. Use
native pixels; do not infer new slide physical scales or clinical meaning.

Primary comparison: fixed 128-pixel versus fixed 2048-pixel cores, both using
`study_cell_reuse.cell_trial`, exact physical visibility, residual PPS, partial
negative-marker deductions and the same predictable-support betting interval.
The target is absolute interval width <=0.005 (0.5 percentage points), alpha .05,
maximum 2048 acquired regions. Run 20 independent acquisition randomizations per
patient and method, seed 20261102+1000*patient_number+repeat. Evaluate an additional
128-core nested-only reuse ablation using the same seeds. Do not tune after viewing
these results or replace a patient because of a poor result.

Pass the selected fine-cell candidate only if on every patient it (1) attains
precision in at least 95% of runs and (2) has lower mean union area than the coarse
cell-reuse comparator. Report coarse attainment as well; a failed comparator does
not support an equal-precision win. Report all coverage checks and empty sets.
Report summed area, acquisition counts, and the overhead at which each comparison
changes direction. No measured monetary or instrument-time saving is claimed.

A failed primary criterion remains a failed confirmation. Earlier broad-precision
confirmation is a separate experiment and cannot substitute for this test.
