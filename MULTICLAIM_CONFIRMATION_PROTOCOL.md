# Frozen multi-claim confirmation

Commit and push this protocol, its source mapping, implementation and tests
before downloading or inspecting the new patients' cell tables. This is a
repository-frozen protocol, not an externally registered trial.

## Patients and fixed preparation

Use CRC17, CRC18 and CRC20 from ORION. Exact author S3 CSV names verify the
P37_S49/S50/S52 mapping; source archive members and range offsets are in
multiclaim_sources.json. These patients have no existing cell/prepared files in
this workspace at protocol freezing. Metadata discovery does not inspect calls.

Reuse prepare_panel without modification: 32 coordinate-selected 2048-pixel
cores, 80-pixel halos, at least 100 cells per candidate, fixed coordinate seed
20261023+patient number, and 100 complete-source neighborhood checks. Keep the
conditional source graph and its radius-80 edges; retain native pixel units.

## Methods and outcomes

Three geometry-fixed target fractions for the same PanCK/CD3e marker pair at
radii 40/60/80. Hypothetical camera field 1024 native pixels, batch size 16,
width .005 for EVERY claim, joint alpha .05 (alpha/3 per claim), total randomized
query limit 4096. No tuning on confirmation outcomes.

Compare envelope, equal-mixture, focused-reuse and deterministic greedy methods
exactly as in study_multiclaim.py. Use 20 stochastic repeats per patient/method,
one greedy run because it is deterministic. Seed 20261110+1000*patient+repeat
is an arbitrary reproducibility identifier. Share the physical observation cache
and known-edge deductions for every method. Focused-reuse updates all supported
claims, preventing discarded usable samples from creating a false advantage.

The primary criterion is envelope attainment >=.95 in every patient and strictly
lower mean physical field count than BOTH stochastic controls in every patient,
with both controls also attaining >=.95. The combined criterion requires all
three patients to pass. This is a predeclared empirical criterion, not a
statistical significance test. Retain every patient's losses and individual paths.

Report joint coverage over all inspected batch boundaries, empty sets, all
target intervals, fields, feedback rounds, randomized queries, and cost sensitivity
fields+h*rounds for h=0/1/10/100. No universal superiority over greedy is expected;
compare the entire cost frontier. Monte Carlo paths quantify acquisition
randomness within a panel, not independent patient replication.

Preserve source CRC32/length and SHA256 checks, prepared-graph checksums, protocol
commit and implementation checksums in the results. No source image model,
clinical readout, actual field latency or independently adjudicated labels are
evaluated. This is patient-disjoint confirmation within the same cohort; a second
cohort and measurement-error robustness remain separate gates.
