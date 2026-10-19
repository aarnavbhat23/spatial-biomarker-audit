# Cell-level observation reuse

Exploratory extension after the residual acquisition screen. Geometry-only adaptive
selection did not beat the best fixed scale consistently. The next experiment
changes how much information an acquisition retains, not its selection heuristic.

A measured rectangle reveals both marker labels for every fixed-frame cell inside
it, including halo cells. An edge product is known if both endpoints are measured,
or if its measured anchor is Pan-CK negative, or if its measured neighbor is CD3e
negative. These latter two cases certify a zero without measuring the other cell.
For independent unrestricted binary marker assignments, setting every unmeasured
bit to zero or one attains the lower or upper remaining-edge bound. These are
logical bounds for the stated labels, not biological coexpression constraints.

Select a fixed core (128 or 2048). Draw a region proportional to its currently
unresolved owned-edge count. Before expanding the known set, compute the residual
PPS observation from the selected region's previously unresolved owned edges.
Then reuse all newly identifiable edges anywhere in the panel. The conditional
expectation and betting proof from the residual study continue to apply.

Compare nested-core reuse with full cell reuse on CRC02, CRC11, CRC19; 10 repeats,
seed 20261031, widths 0.04 and 0.005, maximum 2048 acquisitions. Keep fixed scales
and the same confidence law. Report area, acquisitions, precision, coverage, and
cost sensitivity. No fresh patient confirmation is claimed. Extra reuse must be
independently tested against brute-force visibility and edge-logic calculations.

This applies standard observation reuse and logical propagation. The possible
research contribution is an acquisition-aware spatial assay with validated cost
and uncertainty, not claiming these mathematical ingredients were invented here.
