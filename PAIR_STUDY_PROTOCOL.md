# Direct pair prediction: independent confirmation protocol

Written after the initial five-patient pilot and its confidence audit, before fetching or training on the two additional confirmation patients. The earlier patient results motivated this experiment and are development evidence only.

Hypothesis: predicting joint marker-pair positivity directly from both H&E neighborhoods improves the pair readout over multiplying independent single-marker probabilities, and over a strong geometry-only/pair-prevalence control.

Keep original marker classifiers, all feature transformations, graph radius, tile selection and training/validation patients fixed. Train a histogram gradient-boosted classifier on 100,000 uniformly sampled training edges from CRC01/CRC03, with max_iter=100, max_leaf_nodes=15 and l2_regularization=1.0. Input: the two endpoints' H&E descriptors, their absolute difference, and distance divided by 80 pixels. No measured marker value is an input.

Controls: independent-marker probability product; distance-only boosted classifier; constant training edge-positive prevalence (using the same sampled training edges). Compare all validation Brier scores on CRC19, with no hyperparameter search. Estimate a coefficient in [0,1] for each predictor from validation tile totals. Freeze all choices and model files before computing new-patient outcomes.

New confirmation patients: CRC04 and CRC12, 32 occupied core tiles each, selected by the already fixed coordinate-only rule. They were not used to select the approach or hyperparameters. This is two-patient internal confirmation in the same cohort, not independent external validation. Repeat uniform random tile acquisition 1000 times at 4, 8, 16 and all 32 tiles with seed 20261019. Report measured-only, each predictor with coefficient one, and each validation-tuned coefficient. Do not silently select the best coefficient on confirmation outcomes.

Primary comparison: at eight acquired tiles, compare direct-pair RMSE against independent-marker, distance-only and constant-prevalence controls, with all patient-specific results shown. A preliminary positive result requires lower average RMSE than each control and an improvement over the independent-marker predictor on both confirmation patients. If not met, reject the methodological-advantage claim for this version. Report test-time readout errors, edge-level Brier/AUC, approximate-interval coverage and classical rigorous-interval widths. A narrower interval that undercovers is not an efficiency win.

This is a targeted empirical hypothesis. Joint prediction, boosted trees, and sampling correction are established methods. Passing this small confirmation does not establish architectural novelty, a clinical biomarker, or a cost-saving deployment claim.

Before fitting: the primary comparison uses validation-tuned coefficients for every predictor. Coefficient-one variants are secondary and all are reported.
