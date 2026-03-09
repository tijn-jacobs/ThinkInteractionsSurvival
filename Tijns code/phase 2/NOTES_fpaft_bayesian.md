# Bayesian Flexible Parametric AFT: Complete Mathematical Notes

## 1. Overview

These notes document the mathematical framework and Bayesian implementation
of the flexible parametric accelerated failure time (FPAFT) model proposed
by Crowther, Royston & Clements (2023, *Biostatistics*). We build from the
constant acceleration factor case to the full time-dependent extension.

The core idea: model the log cumulative hazard as a restricted cubic spline
of log-transformed (accelerated) time, then extend to allow the acceleration
factor itself to vary over time using additional splines.

---

## 2. Restricted Cubic Splines

All model components use restricted cubic splines (RCS), also called natural
cubic splines. Understanding the basis is essential.

### 2.1 Basis Functions

Given K knots at positions k_1 < k_2 < ... < k_K, the RCS basis for a
scalar argument u consists of K functions:

    v_1(u) = 1                                              [intercept]
    v_2(u) = u                                              [linear]
    v_{j+2}(u) = (u - k_{j+1})^3_+
                - lambda_j (u - k_1)^3_+
                - (1 - lambda_j) (u - k_K)^3_+             [j = 1,...,K-2]

where (x)_+ = max(x, 0) and lambda_j = (k_K - k_{j+1}) / (k_K - k_1).

The spline function is then:

    s(u | gamma, k) = sum_{j=1}^{K} gamma_j * v_j(u)

### 2.2 Derivatives

The derivative s'(u) = ds/du is obtained by differentiating each basis function:

    v_1'(u) = 0
    v_2'(u) = 1
    v_{j+2}'(u) = 3(u - k_{j+1})^2_+ * I(u > k_{j+1})
                 - lambda_j * 3(u - k_1)^2_+ * I(u > k_1)
                 - (1 - lambda_j) * 3(u - k_K)^2_+ * I(u > k_K)

So s'(u) = sum_{j=1}^{K} gamma_j * v_j'(u).

### 2.3 Key Property

RCS are constrained to be **linear** beyond the boundary knots k_1 and k_K.
This means extrapolation follows the local linear trend, which is a sensible
default for survival models and avoids the wild behavior of unconstrained
polynomial splines in the tails.

### 2.4 Degrees of Freedom and Knot Count

    df = K - 1

    | df | K (knots) | Parameters | Equivalent to        |
    |----|-----------|------------|----------------------|
    | 1  | 2         | 2          | Weibull AFT          |
    | 2  | 3         | 3          | Moderate flexibility |
    | 3  | 4         | 4          | Flexible cubic       |
    | 5  | 6         | 6          | Very flexible        |
    | 7  | 8         | 8          | Highly flexible      |

### 2.5 Knot Placement

Standard approach: place knots at quantiles of uncensored log event times.
For df = d, use K = d+1 knots:
- Boundary knots at the minimum and maximum uncensored log(t)
- d-1 internal knots at equally-spaced centiles between boundary knots

Knots are **fixed** (not estimated). This is standard in both the frequentist
and our Bayesian implementation.

---

## 3. Constant Acceleration Factor Model

### 3.1 Model Specification

The AFT model with constant acceleration factor:

    S(t|X) = S0(t * phi(X; beta))

where phi(X; beta) = exp(-X beta) is the acceleration factor.

Working on the log cumulative hazard scale:

    log H(t|X) = s(u | gamma, k0)

where u = log(t * phi) = log(t) - X beta and k0 are the baseline knots.

### 3.2 Hazard Function

Differentiating H(t) = exp(s(u)):

    h(t|X) = dH/dt = exp(s(u)) * s'(u) * du/dt

Since u = log(t) - X*beta, we have du/dt = 1/t, giving:

    h(t|X) = exp(s(u)) * s'(u) * (1/t)

### 3.3 Survival Function

    S(t|X) = exp(-H(t|X)) = exp(-exp(s(u)))

### 3.4 Log-Likelihood

For observation i with event time y_i, event indicator delta_i, and
possible delayed entry (left truncation) at t_{0i}:

    l_i = delta_i * [s(u_i) + log(s'(u_i)) - log(y_i)]
        - exp(s(u_i))
        + exp(s(u_{0i}))                                [if t_{0i} > 0]

where u_i = log(y_i) - X_i beta and u_{0i} = log(t_{0i}) - X_i beta.

### 3.5 Monotonicity Constraint

For a valid cumulative hazard, H(t) must be non-decreasing, which requires:

    s'(u) > 0 for all u in the data range

This is **not** automatically enforced by the RCS parameterization. In the
likelihood, events contribute log(s'(u_i)), which becomes -infinity if
s'(u_i) <= 0, so the sampler will reject such proposals. But this can
make sampling inefficient. Priors help (see Section 5).

---

## 4. Time-Dependent Acceleration Factor Model

### 4.1 Motivation

The constant AF assumption (phi independent of t) is analogous to the
proportional hazards assumption. It can be violated in practice — the
breast cancer example in the paper shows an AF that changes substantially
over follow-up time.

### 4.2 Model Specification

The survival function for a time-dependent AFT:

    S(t) = S0( integral_0^t eta(X, u; beta) du )

where eta(X, u; beta) is the instantaneous time-varying acceleration
factor at time u.

**The key modeling insight:** Rather than modeling eta directly (which
would require numerical integration), we model on the cumulative scale:

    S0( integral_0^t eta du ) = S0( t * phi(X, t; beta) )

where phi(X, t; beta) is now time-dependent. The integral is absorbed
into the product t * phi(X, t).

### 4.3 Recovering the Instantaneous AF

Given the cumulative quantity t * phi(X, t), the instantaneous AF is:

    eta(X, t; beta) = d/dt [ t * phi(X, t; beta) ]
                    = phi(X, t) + t * dphi/dt              ... (*)

This is Equation 3.6 of the paper.

### 4.4 Parameterization of phi(X, t)

Using RCS of log(t) for the time-dependent part:

    phi(X, t; beta) = exp( -X beta - sum_{p=1}^{P} x_p * s_p(log(t) | gamma_p, k_p) )

Here:
- X beta is the constant (time-independent) part
- Each covariate x_p with a time-dependent effect gets its own spline
  function s_p of log(t), with its own knot vector k_p and coefficient
  vector gamma_p
- When all gamma_p = 0, this reduces to the constant AF model

### 4.5 Derivative of phi

    dphi/dt = phi(X, t) * d/dt [ -X beta - sum_p x_p * s_p(log(t)) ]
            = phi(X, t) * ( - sum_p x_p * s_p'(log(t)) * (1/t) )
            = phi(X, t) / t * ( - sum_p x_p * s_p'(log(t)) )

This is Equation 3.9 of the paper.

### 4.6 The "w" Variable

Define w = log(t * phi(X, t)):

    w = log(t) + log(phi(X, t))
      = log(t) - X beta - sum_p x_p * s_p(log(t))

Its time derivative:

    dw/dt = 1/t - sum_p x_p * s_p'(log(t)) * (1/t)
          = (1/t) * (1 - sum_p x_p * s_p'(log(t)))

Define the **correction factor**:

    C(t) = 1 - sum_p x_p * s_p'(log(t))

### 4.7 Hazard Function (Eq 3.10)

    h(t|X) = exp(s0(w)) * ds0/dt
           = exp(s0(w)) * s0'(w) * dw/dt
           = exp(s0(w)) * s0'(w) * (1/t) * C(t)

where s0 is the baseline spline and C(t) is the correction factor.

**Step-by-step derivation of the chain rule:**

Starting from H(t|X) = exp(s0(w(t))):

    h(t|X) = dH/dt = exp(s0(w)) * s0'(w) * dw/dt

We need dw/dt. Since w = log(t) - X*beta - sum_p x_p * s_p(log(t)):

    dw/dt = d/dt[log(t)] - d/dt[sum_p x_p * s_p(log(t))]
          = 1/t - sum_p x_p * s_p'(log(t)) * d(log(t))/dt
          = 1/t - sum_p x_p * s_p'(log(t)) * (1/t)
          = (1/t) * [1 - sum_p x_p * s_p'(log(t))]
          = C(t) / t

Therefore:

    h(t|X) = exp(s0(w)) * s0'(w) * C(t) / t

### 4.8 Survival Function (Eq 3.11)

    S(t|X) = exp( -exp(s0(w)) )

Identical in form to the constant case, but w now depends on the
time-dependent splines.

### 4.9 Log-Likelihood

For observation i:

    l_i = delta_i * [ s0(w_i) + log(s0'(w_i)) - log(y_i) + log(C(y_i)) ]
        - exp(s0(w_i))
        + exp(s0(w_{0i}))                               [if t_{0i} > 0]

where:
- w_i = log(y_i) - X_i beta - sum_p x_{ip} * s_p(log(y_i))
- C(y_i) = 1 - sum_p x_{ip} * s_p'(log(y_i))
- w_{0i} uses t_{0i} instead of y_i

Compared to the constant AF log-likelihood, the **only addition** for
events is the + log(C(y_i)) term. The survival contribution is unchanged
in form (just w is computed differently).

### 4.10 Computing the Instantaneous AF from Posterior Draws

From equation (*), for a specific covariate pattern (e.g., X=1 binary treatment):

    eta(t) = phi(t) + t * dphi/dt
           = phi(t) * (1 - sum_p x_p * s_p'(log(t)))
           = phi(t) * C(t)

where phi(t) = exp(-beta - sum_p x_p * s_p(log(t))).

For posterior inference on eta(t), compute this at each posterior draw
and each time point, then take quantiles for credible intervals.

---

## 5. Prior Specification

### 5.1 Covariate Effects (beta)

    beta ~ Normal(0, 2)

Weakly informative. Allows acceleration factors exp(-beta) ranging from
roughly exp(-4) ≈ 0.02 to exp(4) ≈ 55, which covers any plausible
clinical effect. For more informative settings, Normal(0, 1) may be used.

### 5.2 Baseline Spline Coefficients (gamma0)

    gamma0[1] ~ Normal(0, 5)       [intercept: overall hazard level]
    gamma0[2] ~ Normal(1, 2)       [linear slope: centered positive]
    gamma0[j] ~ Normal(0, 1)       [j >= 3: non-linear terms]

**Rationale for gamma0[2] ~ Normal(1, 2):** The linear coefficient
controls whether the log cumulative hazard is increasing in u. For a
valid model we need s'(u) > 0 in the data range, and centering gamma0[2]
at a positive value helps. The variance of 2 allows flexibility.

**Rationale for gamma0[j] ~ Normal(0, 1) (j >= 3):** These control
curvature. The Normal(0, 1) prior acts as a Bayesian roughness penalty,
analogous to P-splines. Smaller variance = smoother baseline. This
provides automatic regularization against overfitting.

### 5.3 Time-Dependent Effect Spline Coefficients (gamma_tde)

    gamma_tde[p][j] ~ Normal(0, 0.5)   [all terms, for each TDE covariate p]

**Rationale:** These represent departures from a constant acceleration
factor. Centering at 0 provides **shrinkage toward the constant AF model**.
This is a key advantage of the Bayesian approach: when the data do not
support time-dependence, the posterior for gamma_tde concentrates near 0.

The tighter variance (0.5 vs 1.0 for baseline) reflects the prior
belief that time-dependent effects should be modest departures from
constancy. This also helps satisfy the positivity constraint on C(t).

### 5.4 Hierarchical Smoothing (optional extension)

Replace the fixed-variance priors with:

    gamma0[j] ~ Normal(0, tau_0)        for j >= 3
    gamma_tde[p][j] ~ Normal(0, tau_p)
    tau_0 ~ Half-Cauchy(0, 1)
    tau_p ~ Half-Cauchy(0, 0.5)

This lets the data determine the appropriate smoothness level. Equivalent
to adaptive penalized smoothing in a frequentist context.

---

## 6. Constraints and Numerical Issues

### 6.1 Baseline Monotonicity: s0'(u) > 0

Required for a valid (non-decreasing) cumulative hazard. In the likelihood,
events contribute log(s0'(u)), which is -infinity when s0'(u) <= 0.

**Handling:** The prior on gamma0[2] centered at 1 helps, and the likelihood
automatically rejects violations. If sampling is inefficient (many rejected
proposals), options include:
- Tighter priors on non-linear gamma0 coefficients
- Reparameterizing gamma0[2] with a lower bound (e.g., using a lognormal)
- Using more warmup iterations

### 6.2 Hazard Positivity: C(t) > 0

The correction factor C(t) = 1 - sum_p x_p * s_p'(log(t)) must be
positive for the hazard to be valid. Since events contribute log(C(y_i)),
violations give -infinity log-likelihood.

This constraint is **more subtle** than baseline monotonicity because it
depends on the covariate values x_p. For binary covariates:
- When x_p = 0: C(t) = 1 (always positive, no constraint)
- When x_p = 1: need s_p'(log(t)) < 1 everywhere in the data range

**Handling:** The Normal(0, 0.5) prior on gamma_tde keeps the derivatives
small. For more complex TDE patterns, the prior variance may need tuning.

### 6.3 Identifiability Between beta and gamma_tde Intercept

The TDE spline s_p(log(t)) includes an intercept term (v_1 = 1). This
intercept is confounded with the constant effect beta_p: adding a constant
c to s_p and subtracting c from beta_p leaves phi unchanged.

In our implementation, we keep both free and let the prior handle
identifiability. The tight prior on gamma_tde shrinks the intercept
toward 0, effectively anchoring the decomposition. An alternative is to
fix the TDE intercept to 0 (forcing s_p(log(t)) to represent only
time-varying departures), but this requires modifying the basis.

In practice, the sum beta + gamma_tde[1] (the intercept of the TDE
spline) is well-identified even if the individual components have some
posterior correlation. Check for this by examining pairwise posterior
correlations between beta and gamma_tde[,1].

---

## 7. Sign Conventions

This is a common source of confusion. The paper uses:

    S(t|X) = S0(t * phi(X; beta))
    phi(X; beta) = exp(-X * beta)

Therefore:
- beta > 0, X = 1: phi = exp(-beta) < 1, time is "decelerated" (stretched)
  => LONGER survival for the exposed group
- beta < 0, X = 1: phi = exp(-beta) > 1, time is "accelerated" (compressed)
  => SHORTER survival for the exposed group

The **acceleration factor** (as commonly reported) is exp(-beta):
- AF = exp(-beta) > 1: survival time extended by a factor of AF (protective)
- AF = exp(-beta) < 1: survival time shortened (harmful)

In the time-dependent case, the instantaneous AF is eta(t), and:
- eta(t) > 1: at time t, survival time is being extended
- eta(t) < 1: at time t, survival time is being shortened

---

## 8. Connection to Collapsibility and Causal Inference

A key result from Section 2 of the paper: the AFT acceleration factor
is **collapsible**, meaning:

    Marginal causal effect = -beta_X

regardless of the distribution of omitted covariates Z (as long as Z is
independent of X). This is NOT true for proportional hazards models,
where the marginal hazard ratio is attenuated over time even when X and
Z are independent.

For the time-dependent case, collapsibility still holds:

    Marginal causal effect at time t = -beta_X(t)

This makes the AFT framework particularly attractive for causal inference.
The Bayesian implementation inherits this property. The posterior on beta
(or beta(t) via the TDE splines) can be interpreted causally under the
usual assumptions (no unmeasured confounding, positivity, consistency).

---

## 9. Implementation Summary

### Files

| File                   | Purpose                                      |
|------------------------|----------------------------------------------|
| `fpaft_constant.stan`  | Stan model, constant AF only                 |
| `fpaft_tde.stan`       | Stan model, constant + time-dependent AF     |
| `fit_fpaft_constant.R` | Simulate + fit constant AF model             |
| `fit_fpaft_tde.R`      | Simulate + fit TDE model with validation     |
| `compare_df.R`         | Compare df=3/5/7 for baseline flexibility    |

### Workflow

1. **Simulate or prepare data**: event times, event indicators, covariates,
   (optional) delayed entry times.
2. **Choose df for baseline** (df=3-5 is typical; use LOO-CV to compare).
3. **Compute knots** at quantiles of uncensored log event times.
4. **Fit constant AF model first** as a baseline.
5. **Decide which covariates need TDE**: examine residuals, clinical
   knowledge, or fit TDE and check if gamma_tde shrinks to ~0.
6. **Choose df for TDE** (df=1-2 is usually sufficient).
7. **Fit TDE model** via cmdstanr with adapt_delta >= 0.95.
8. **Diagnose**: check divergences, Rhat, ESS, posterior predictive checks.
9. **Compare models** via LOO-CV (using the `loo` package).
10. **Report**: posterior summaries for beta, plots of eta(t) with credible
    intervals, survival curves with uncertainty.

### Computational Notes

- The Stan implementation uses analytic spline basis and derivatives (no
  numerical differentiation), matching the paper's emphasis on efficiency.
- For large datasets, consider within-chain parallelization via
  `reduce_sum` in Stan.
- The model handles delayed entry (left truncation) and right censoring.
  Interval censoring would require extending the likelihood.
- Typical runtime: seconds to low minutes for N ~ 1000, df_baseline ~ 5,
  df_tde ~ 2, with 4 chains x 2000 post-warmup draws.

---

## 10. Validation Results

### Constant AF (Step 1)

Simulated from Scenario 1 of the paper (mixture Weibull baseline,
beta = -0.5, N = 1000):

- df=3: posterior mean beta ≈ -0.40, noticeable bias (~0.10). Consistent
  with paper's Table 2 showing this baseline needs higher df.
- df=5-7: bias expected to shrink substantially (verify with compare_df.R).
- No divergent transitions, Rhat = 1.00, ESS > 1500.

### TDE Model (Step 2)

- Part A (constant data + TDE model): TDE coefficients should shrink to
  ~0, confirming the prior provides appropriate regularization.
- Part B (time-dependent data): model should recover the known AF trajectory.

---

## 11. Potential Extensions

- **Multiple covariates with TDE**: Already supported via `n_tde` and
  `tde_idx` in `fpaft_tde.stan`.
- **Random effects / frailties**: Add log-normal frailty for clustered data.
- **Hierarchical smoothing**: Replace fixed prior variances with estimated
  smoothing parameters (Section 5.4).
- **Horseshoe prior on TDE**: For automatic variable selection of
  time-dependent effects.
- **Cure fraction models**: Extend baseline to allow S(t) -> c > 0
  as t -> infinity.
- **Interval censoring**: Modify likelihood to handle interval-censored
  observations (requires evaluating S at two time points per observation).
- **Competing risks**: Cause-specific FPAFT models, though the AFT
  interpretation becomes more complex with competing events.
