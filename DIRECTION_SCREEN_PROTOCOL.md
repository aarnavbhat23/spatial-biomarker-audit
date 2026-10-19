# Independent acquisition direction screen

This project is separate from Extrusion Net. These are exploratory development
tests using CRC02, CRC11 and validation CRC19. Existing confirmation patients will
not be presented as untouched confirmation for the new method.

1. **Overlap-aware observation.** Count cell-coverage signatures for the existing
   physical regions. Compute the exact probability that both endpoints of each
   directed pair are observed under tile SRS, including endpoints in different
   selected regions. Compare expected observed pairs to core ownership at budgets
   4, 8 and 16. Continue if there is substantial recoverable information, rather
   than changing estimators for negligible overlap. This is established graph
   sampling mathematics; any novelty must go beyond inclusion weighting.
2. **Prediction-free acquisition at the pair scale.** Sample geometric directed
   pairs independently with replacement from the fixed panel and measure a small
   square covering their endpoints. Bernoulli outcomes permit exact binomial
   intervals conditional on the finite panel. Test widths and coverage, physical
   union-area proxies for 128/256/512-native-pixel squares, and overhead-cost
   break-even versus coarse tiles. Extra labels in a square are deliberately not
   used in the initial point estimate, to preserve a transparent comparison.
3. **Worst-case measurement requirements.** Construct two valid marker assignments
   identical outside a single high-mass tile but with different pair readouts.
   Calculate their probability of being indistinguishable under tile sampling.
   This is a classical information lower-bound construction, not automatically
   a new theorem. Use it to reject unjustified small-budget precision claims.

The graph and measurement coordinates remain derived from source fluorescence
segmentation. No deployment-level acquisition savings, instrument feasibility,
functional interaction, clinical validity, or new-patient guarantee is implied.
Costs are area-plus-per-region sensitivity analyses, not measured assay prices.

Selection favors a demonstrable cost/valid-precision improvement and an unresolved
research question. A simple acquisition baseline can be useful without being a
novel algorithm. A new confirmatory cohort and a frozen protocol are required for
any claim of replication. Sources, hashes, all failures and all screen results
will remain available.
