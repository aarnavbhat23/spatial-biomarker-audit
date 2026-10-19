# What the intervals do and do not guarantee

The population consists of K fixed core tiles with known geometric pair counts
M_i. Each tile owns pairs whose first endpoint is in its core. Its halo observes
every possible second endpoint within the fixed radius. Observing a tile reveals
its complete measured pair total Y_i. The target is sum(Y_i)/sum(M_i).

This conditions on the cell detector, marker gating, selected panel, registration
and prediction map. It does not address their measurement error or transfer to
another patient. The field of view is the selected panel, not the whole slide.

## Point estimate and variance

Freeze predicted totals F_i and coefficient c using development data. Define
R_i=Y_i-cF_i. For an SRS sample S of m distinct tiles,

```
estimate = [c sum(F_i) + K mean_{i in S}(R_i)] / sum(M_i).
Var(estimate) = (K/sum(M_i))^2 (1-m/K) S_R^2/m.
```

The sample variance of observed residuals estimates S_R^2 unbiasedly. None of
these facts implies that a Student-t interval has finite-sample coverage.
One positive tile among 100 is a counterexample: at m=2, a measured-only
Student interval is [0,0] in the 98% of samples that miss that tile. A small
empirical standard error is not a certificate when rare contributions can be missed.

## Conservative bounded-variable comparator

Known bounds give a=min(-cF_i), b=max(M_i-cF_i). Normalize X_i=(R_i-a)/(b-a)
to [0,1]. A standard without-replacement Hoeffding bound is valid but may be
vacuous at small m. It gives a useful correctness baseline, not a new method.

## Fixed-budget betting interval

For a candidate population mean mu of X, the conditional mean before draw j is

```
mu_j = [K mu - sum_{k<j} X_k] / (K-j+1).
```

For each fixed positive s, the factors

```
positive: (s + X_j)/(s + mu_j)
negative: (s + 1-X_j)/(s + 1-mu_j)
```

are nonnegative and have conditional expectation one at the true mu. Products
are therefore mean-one nonnegative martingales. A fixed equal-weight mixture
over s in {0.01,0.05,0.1,0.5,1.0} retains this property. Markov's inequality and
a union bound allow rejection when either one-sided wealth exceeds 2/alpha.
Invert those tests at the stated fixed budget and map the interval back to the
pair fraction. The positive wealth decreases with mu and negative wealth
increases with mu on the feasible range, so bisection is valid; outer endpoints
preserve coverage up to the stated floating-point numerical tolerance.

The feasible range is tightened by observed totals and the known maximum
contributions of all remaining tiles. When every tile is observed, it reduces
to the census value. An empty interval on a rejection event is recorded as
noncoverage, not silently repaired. The implementation's advertised guarantee
is fixed-budget; the experiments do not claim a simultaneous guarantee over
all displayed budgets. Different stabilizers, data-dependent betting choices,
or adaptive acquisition require their own validity argument.

This is an application of established betting/sampling inference. Relevant
primary references are [Waudby-Smith and Ramdas, NeurIPS 2020](https://arxiv.org/abs/2006.04347)
and [Waudby-Smith and Ramdas, JRSS B 2024](https://academic.oup.com/jrsssb/article/86/1/1/7043257).
Exact small-population tests independently check unbiasedness, variance,
martingale expectations and coverage; they are checks of an implementation,
not a substitute for the argument.

Rigorous coverage can come with extremely wide intervals. An acquisition
policy or a new predictor must be compared at matched coverage and actual
measurement cost. A narrow under-covering interval is not an improvement.
