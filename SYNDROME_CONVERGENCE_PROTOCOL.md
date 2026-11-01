# Resolve the measured fitting cap without tuning decoding thresholds

2026-11-01. A post-outcome computational follow-up, frozen before execution.
The original measured crop hit max_iters=40 in all 16 bits. Its 39 decoded
barcodes are incomplete outputs, not an accuracy result.

Repeat the same 13 stacks and every original parameter, changing only the
measured crop's max_iters from 40 to 200, matching the iteration cap in the
authors' ADCG example. No crop, bit, codeword, threshold, PSF bounds or decoding
cost is replaced. Save results separately; do not overwrite the original run.
The synthetic controls are rerun unchanged as a consistency check.

This is a convergence diagnostic. Pass only if no measured bit reaches 200
fitted points. Report the count/gene-composition change from the 40-cap run
regardless of direction. If the cap still binds, declare nonconvergence at this
budget and stop interpreting measured barcode counts. Even a pass does not
calibrate the PSF/background model, identify real molecules or measure FDR.
