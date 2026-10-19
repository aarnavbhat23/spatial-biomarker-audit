# Initial study protocol

Written before training or inspecting final readout results.

Question: can probability-sampled rectangular tissue regions correct a predicted protein-pair proximity score while accounting for the cost of observing both endpoint labels?

Use same-section paired H&E and measured multiplex tissue images from ORION. Start with a small, explicitly named subset; no claim of cohort-wide or clinical validation. Patient IDs define model train/validation/test splits. Frozen geometric graphs and train-defined marker transformations define the estimand. No target protein channel enters prediction features.

The initial model is a simple H&E-only baseline. A measured-image segmentation may define a conditional evaluation graph, which must be disclosed. The primary outcome is continuous Pan-CK/CD3e proximity, not a validated cell phenotype or functional interaction. Compare measured-only and prediction-assisted random tile sampling first. Treat halo pixels and unique acquired area as costs. Any learned allocation must be fitted on separate development data.

For an initial feasibility pilot, report all selected patients/regions, model accuracy, readout errors, sampling bias, interval coverage, and width. Repeat acquisition 1000 times per reported budget. Repeated acquisition is within-population uncertainty, not patient replication. Selection of convenient available image regions makes this a feasibility pilot, not a representative biological sample.

Use finite-population tile sampling with a core/halo observation rule and exact tile-total variance. Require at least two sampled tiles for an estimated variance. Report normal/Student intervals as approximations; do not infer finite-sample coverage from unbiasedness. Include exact small-population enumeration and disclose zero-variance or rare-event failures. A rigorous bounded interval is a conservative comparator when the statistic has a known range.

No methods will be called new unless the closest survey-sampling, spatial inference and virtual-staining precedents are distinguished. Any failure, changed protocol, or shortcut is recorded. A successful pilot is an executable pipeline and valid accounting, not a novel-method claim.

## Data-format amendment, before model training

The accessible processed cell tables contain the source authors' binary marker-positivity calls rather than continuous intensities. The pilot therefore uses their Pan-CK_pos/CD3e_pos label-pair fraction. This changes the initial continuous-signal endpoint before examining model/readout results. These are measured-channel-derived gates, not independently adjudicated tumor/T-cell identities.

Fixed pilot: CRC01/CRC03 train (24 core tiles per patient); CRC19 validation (12 tiles); CRC02/CRC11 test (32 tiles each). Select occupied 2048-pixel core tiles using coordinates only, with seed 1729, excluding slide-edge tiles. Geometric radius 80 native pixels; verify scale from each registered image. The graph uses the source fluorescence segmentation, so the evaluation is conditional on that cell detection. The inferential population is the fixed selected panel of core tiles, not the whole slide.

H&E features: color and optical-density statistics at fixed local scales, without protein labels. Compare logistic regression and small histogram gradient-boosted classifiers; choose by mean validation Brier score. Use no test-based tuning. Predict each marker separately and multiply endpoint probabilities for the pair-score baseline. Freeze a prediction coefficient in [0,1] from the validation tile totals; compare it with coefficients zero and one.

Acquire 2, 4, 8, 16, or all 32 core tiles by simple random sampling without replacement, repeat 1000 times, seed 20261018. The measurement halo is 80 pixels and unique union area is charged. Report Student-t approximate intervals and a conservative bounded-variable Hoeffding interval, including failures or uninformative widths. Never label these classical estimators new. No learned region allocation will be claimed in this pilot.
