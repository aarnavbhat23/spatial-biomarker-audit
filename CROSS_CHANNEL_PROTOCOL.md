# Cross-channel reference calibration — frozen 2026-11-01

Freeze before fitting or scoring reference residuals. The metadata/identity audit
has already checked 680 source-selected rows, six directed channel pairs and
131 connected components (106 with three channels, 25 with two). No component
contains two raw coordinates from one channel. Duplicate reference-coordinate
columns agree. This audit used identities/headers only, not residual outcomes.

Use the exact source blob 024fe5233ece13c91db1116c9aa43bbbef791cf0;
SHA256 7ebb373a323bdc3254745c3083dd62ccaf559876596c704dacfdfc2fceaf32ab.
Preserve source license/attribution. Parse duplicate headers by fixed positions.
Nodes are (channel, exact raw x, exact raw y); union every source match's two
nodes. These author-assigned matches are presumed identities, not adjudication.
Save every row and component; flag repeated-channel/inconsistent components
instead of silently selecting a favorable subset. Require >=60 valid components
and >=10 training matches per directed pair for the primary benchmark.

Sort components by SHA256 of their sorted node triples encoded with float.hex.
First floor(N/2) are training, next floor(N/4) calibration, remaining test. All
directions/repeats of a feature stay in its component. Never use the author-aligned
coordinates to fit or claim held-out transform accuracy. Controls: median raw
translation and ordinary least-squares affine transforms, each directed pair
fitted separately using only training components. Transform input is comparison
raw xy; response is reference raw xy. No robust-fitting or feature-selection tuning.

For each component, use the maximum Euclidean residual over ALL its observed
directed matches as one score. A 95% split-conformal radius is calibration order
statistic k=ceil(.95*(n_cal+1)); if k>n_cal radius=infinity. Apply the same radius
to all observed directions for each held-out component. Report test coverage,
Wilson intervals, radius and individual/per-channel residuals. The exchangeability
guarantee is conditional on exchangeable independent feature groups after fitting;
one field, selected matches and spatial dependence do not establish that premise.
This is a within-field reference benchmark, not a biological validation cohort.

Primary technical usefulness gate for affine: test coverage Wilson lower>=.80,
radius<=1 source coordinate unit, and median test component residual at least
20% below translation. These are declared pilot tolerances, not clinical or
gene-composition error tolerances. Keep each gate separately and retain failures.
Report translation even if it wins. Source coordinate units must remain as such
unless pixel/physical scale is independently resolved from the authors' workflow.

Spatial-extrapolation stress: sort components by mean raw x, use first half for
training, next quarter calibration, last quarter test, with identical estimators.
Report without retuning. Keep unmatched/undetected-feature limitations explicit:
the source table does not permit evaluating them. Do not transfer this radius to
the existing MERFISH field or use it as an Adamdec1-probe efficiency calibration.
