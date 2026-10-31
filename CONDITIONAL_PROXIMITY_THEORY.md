# Conditional proximity under bounded position errors

This is a conservative construction from an established cross-nearest-neighbour
statistic, Monte Carlo inference and triangle inequalities. No priority claim is
made. The important constraint is the information required by the null model.

For finite A and B positions, define S(A,B;r) as the fraction of B positions
whose nearest A is within radius r. Suppose each recorded A is within rho_A
of its latent position and each recorded B is within rho_B. The nearest-distance
triangle inequality gives

    L = S(A_observed,B_observed;r-rho_A-rho_B) <= S(A_true,B_true;r).

A negative lower radius means L=0. The implementation rounds the lower radius
inward to avoid a floating-point equality error. For an exact simulated null B_j,

    S(A_true,B_j;r) <= U_j = S(A_observed,B_j;r+rho_A).

Round this upper comparison outward. These bounds apply in Euclidean space or
the periodic metric used by the simulation. They need not be jointly sharp.

Assume the true conditional null law Q for B is known given anatomy, counts and
other declared conditioning information, and B_true and B_1,...,B_m are independent
draws from Q conditional on that information and A_true. The ordinary conservative
Monte Carlo p value is

    p_true = (1 + sum_j 1[S(A_true,B_j;r) >= S(A_true,B_true;r)])/(m+1).

Every exceedance counted here is counted by

    p_upper = (1 + sum_j 1[U_j >= L])/(m+1).

Thus p_upper >= p_true for every admissible recorded-error configuration.
Conditional exchangeability gives the usual superuniformity of p_true; pointwise
domination preserves it. If an externally calibrated simultaneous error event
contains truth with probability at least 1-delta, the unconditional rejection
probability is at most alpha+delta. Independence of that event and the test is
not needed for this union bound. For a family, simultaneous containment plus
test levels summing to alpha-delta gives an FWER bound alpha when delta<alpha.

The measured dataset supplies neither Q nor an externally calibrated delta/rho.
An intensity estimate from the same tested points cannot simply be substituted
for Q while retaining this proof. Conditioning on broad tiles assumes the
remaining within-tile distribution is correct; the falsification study explicitly
demonstrates failure of that assumption. A fiducial residual does not calibrate
every molecular localization or chromatic error.

When rho_A+rho_B >= r and recorded cross-channel distances are positive, this
generic lower bound is zero and p_upper=1. Before that boundary, power can already
collapse. The simulation retains that failure instead of selecting a favorable
error radius after seeing the results. The next useful contribution would need
valid uncertainty in Q and measurement calibration with nonvacuous power; the
current oracle construction does not solve those problems.

The cross-nearest-neighbour statistic is established; see the
[spatstat cross-G documentation](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Gcross.inhom.html).
Shared intensity and interaction require distinct assumptions; see the
[primary comparison of spatial intensity tests](https://doi.org/10.1080/03610918.2021.1901118).
Barcode matching likewise has established precedents in
[MERFISH](https://pmc.ncbi.nlm.nih.gov/articles/PMC5047202/) and
[MERlin](https://github.com/emanuega/MERlin). The raw-image screen is a transparent
limited control, not a replacement for a full validated decoding pipeline.
