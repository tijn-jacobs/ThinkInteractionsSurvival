# Bayesian Flexible Parametric AFT: Mathematical Notes
## Building Block 1: Constant Acceleration Factor

### The Model

We model on the log cumulative hazard scale:

```
log H(t|X) = s(u | gamma, k)
```

where `u = log(t * phi(X; beta)) = log(t) - X*beta` and `s()` is a restricted
cubic spline with knot vector `k` and coefficient vector `gamma`.

### Restricted Cubic Spline Basis

For K knots at positions k_1 < k_2 < ... < k_K, the basis functions are:

```
v_1(u) = 1                                          (intercept)
v_2(u) = u                                          (linear)
v_{j+2}(u) = (u - k_{j+1})^3_+ 
            - lambda_j * (u - k_1)^3_+ 
            - (1 - lambda_j) * (u - k_K)^3_+        for j = 1,...,K-2
```

where `lambda_j = (k_K - k_{j+1}) / (k_K - k_1)`.

Key property: The spline is LINEAR beyond the boundary knots. This means
extrapolation beyond the data range follows the local linear trend, which is
a sensible default for survival models.

### Degrees of Freedom

- K knots => K basis functions => K parameters (gamma_0, ..., gamma_{K-1})
- "degrees of freedom" in the paper = K - 1 (the number of interior + 1 parameters beyond intercept)
- df=1 (K=2 knots): Only intercept + linear => Weibull model
- df=3 (K=4 knots): 4 parameters => flexible cubic spline

### Hazard and Survival

From log H(t|X) = s(u), we derive:

```
H(t|X) = exp(s(u))
S(t|X) = exp(-H(t|X)) = exp(-exp(s(u)))
h(t|X) = dH/dt = exp(s(u)) * ds/dt
        = exp(s(u)) * s'(u) * du/dt
        = exp(s(u)) * s'(u) * (1/t)
```

where s'(u) = ds/du is the derivative of the spline function.

### Monotonicity Constraint

For a valid cumulative hazard, we need H(t) to be non-decreasing, which
requires s'(u) >= 0 for all u in the data range. This is NOT automatically
enforced by the spline parameterization.

In a Bayesian context, we handle this by:
1. Using priors that favor positive s'(u) (e.g., gamma[2] ~ Normal(1, 2))
2. Using smoothness priors on higher-order coefficients
3. The likelihood itself will penalize violations (log(s'(u)) becomes -Inf)

### Log-Likelihood

For observation i with event time y_i, event indicator delta_i, and 
possible delayed entry at t_{0i}:

```
l_i = delta_i * [s(u_i) + log(s'(u_i)) - log(y_i)]
    - exp(s(u_i))
    + exp(s(u_{0i}))     [only if t_{0i} > 0]
```

where u_i = log(y_i) - X_i * beta.

### Prior Choices

**beta (log acceleration factors):**
- Normal(0, 2): weakly informative, allows effects up to ~exp(4) ≈ 55x
- For clinical applications, Normal(0, 1) may be more appropriate

**gamma (baseline spline coefficients):**
- gamma[1] (intercept): Normal(0, 5) — controls overall scale of hazard
- gamma[2] (linear slope): Normal(1, 2) — centered positive for monotonicity
- gamma[3:K] (non-linear): Normal(0, 1) — smoothness/regularization prior

The prior on the non-linear gamma coefficients acts similarly to a roughness
penalty in penalized splines. Smaller variance => smoother baseline.

### Knot Placement

Standard approach: quantiles of uncensored log event times.
For df=d, use d+1 knots with boundary at min/max and d-1 internal knots
at equally-spaced centiles.

Note: The paper mentions knots should ideally be on log(t * phi), but since
phi depends on beta (unknown), we use log(t) as an approximation. This is
standard practice and works well in simulation.

### Connection to Standard Models

| df | K (knots) | Parameters | Equivalent to    |
|----|-----------|------------|-------------------|
| 1  | 2         | 2          | Weibull AFT       |
| 2  | 3         | 3          | ~                 |
| 3  | 4         | 4          | Flexible cubic    |
| 5  | 6         | 6          | Very flexible     |

### Sign Convention Warning

The paper uses: S(t|X) = S0(t * exp(-X*beta))

- beta > 0 AND X = 1: phi = exp(-beta) < 1, so time is "decelerated" 
  (stretched), meaning LONGER survival for treated
- beta < 0 AND X = 1: phi = exp(-beta) > 1, so time is "accelerated" 
  (compressed), meaning SHORTER survival for treated

The acceleration factor (AF) is exp(-beta). When reported:
- AF > 1: survival time extended (protective)
- AF < 1: survival time shortened (harmful)

## Next Steps: Time-Dependent Acceleration Factor

The extension replaces phi(X; beta) = exp(-X*beta) with the time-dependent:

```
phi(X, t; beta) = exp(-X*beta - sum_p x_p * s(log(t) | gamma_p, k_p))
```

This adds:
- New spline coefficient vectors gamma_p for each time-dependent covariate
- New knot vectors k_p (can differ from baseline knots)
- Modified hazard formula with additional derivative terms
- Additional constraint: the factor (1 - sum_p x_p * s'(log(t))) must be > 0

The key additional complexity in the hazard is:

```
h(t|X) = exp(s_0(...)) * s_0'(...) * (1/t) * (1 - sum_p x_p * s_p'(log(t)))
```

where s_0 is the baseline spline and s_p are the time-dependent effect splines.

This will be implemented in Step 2.
