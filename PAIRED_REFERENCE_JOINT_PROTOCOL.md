# Joint-count follow-up — frozen after the conditional-test failure

2026-10-30. The first paired-reference screen (protocol ee47e68) fails its
strong-signal decisiveness gate: at b=2, delta=.25, eta=8, E=64, true absence
is correctly resolved in 75.05% of 2,000 trials. Keep that failure. Conditional
ratio inference discards total intensity and bounds the presence ratio using
the worst b=10, even when the background ROI supports b near 2.

This development observation motivates a separate conventional joint-count
control. It is not a claim of a new statistical construction or fresh biology.
Freeze before evaluating this new estimator, using new independent RNG seed
2026103004 and the same full parameter grid, truth definitions, 2,000 replicates,
counterexample and gates as FIDUCIAL_MEASUREMENT_PROTOCOL.md. The seed is fresh;
the mechanistic family is deliberately reused and is not external validation.

For independent target X and background Y counts at dose E and background area
16, form exact Poisson mean bounds [LX,UX] for lambda_X=X/E and [LB,UB] for
b=Y/(16E). Each of the four one-sided tails has error probability <=.0125.
The rectangle has simultaneous coverage >=.95 by a union bound. Declare:

- presence if LX-(1+delta)*UB > 0;
- absence if UX-(1-delta)*LB < eta_min;
- abstain if both or neither declarations hold.

Under the common-background/mismatch and minimum-efficiency contract, a wrong
declaration implies failure of the relevant two one-sided bounds, hence has
probability <=.025. The rectangle is conservative and not claimed minimax.
Compare the original conditional rule, target-only rule and joint rectangle on
the exact SAME saved new X,Y draws. Do not compare only favorable parameter cells.

Primary gate unchanged: every in-contract wrong-decision Wilson upper<=.075;
for eta=8,E=64, every b/delta/truth correct-decision Wilson lower>=.80. Report
the result for each method, including failures. Record the out-of-contract
background mismatch stress separately. Preserve existing fiducial audit and
geometry outcomes; this follow-up does not rerun or retune image processing.
