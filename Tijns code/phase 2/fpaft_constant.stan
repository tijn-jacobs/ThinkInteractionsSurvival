// Flexible Parametric Accelerated Failure Time Model (Constant AF)
// Based on Crowther, Royston & Clements (2023), Biostatistics
//
// Model:
//   log H(t|X) = s(log(t * phi(X; beta)) | gamma, k0)
//   phi(X; beta) = exp(-X * beta)
//
// where s() is a restricted cubic spline (natural spline) function.
//
// The baseline is modeled on the log cumulative hazard scale using
// restricted cubic splines of log(t * phi), giving enormous flexibility.
//
// Knots are fixed (pre-specified), typically at quantiles of
// uncensored log event times.

functions {
  // ---------------------------------------------------------------
  // Restricted cubic spline basis functions
  // ---------------------------------------------------------------
  // Given a scalar u and a vector of knot locations k (length K >= 3),
  // returns the basis vector of length K (intercept + K-1 basis functions).
  //
  // For K knots, we get K basis functions (including intercept):
  //   v1 = 1  (intercept)
  //   v2 = u  (linear)
  //   v_{j+2} = (u - k_j)^3_+ - lambda_j*(u - k_min)^3_+ 
  //             - (1 - lambda_j)*(u - k_max)^3_+
  //   for j = 1, ..., K-2
  //
  // where lambda_j = (k_max - k_{j+1}) / (k_max - k_min)
  //
  // This gives "degrees of freedom" = K - 1 for the spline 
  // (K basis functions total including intercept, so K parameters).
  // df=1 => K=2 knots => Weibull
  // df=3 => K=4 knots => cubic spline with 4 parameters
  
  // Positive part: max(x, 0)
  real pos_part(real x) {
    return fmax(x, 0.0);
  }
  
  // Compute the full basis vector for a single value u
  // knots: sorted knot locations, length n_knots (>= 2)
  // Returns vector of length n_knots (gamma_0, gamma_1, ..., gamma_{n_knots-1})
  vector rcs_basis(real u, vector knots) {
    int K = num_elements(knots);
    vector[K] basis;
    real kmin = knots[1];
    real kmax = knots[K];
    real denom = kmax - kmin;
    
    // Intercept
    basis[1] = 1.0;
    // Linear term
    basis[2] = u;
    
    // Non-linear terms (j = 1, ..., K-2)
    for (j in 1:(K - 2)) {
      real lambda_j = (kmax - knots[j + 1]) / denom;
      basis[j + 2] = pow(pos_part(u - knots[j + 1]), 3)
                     - lambda_j * pow(pos_part(u - kmin), 3)
                     - (1.0 - lambda_j) * pow(pos_part(u - kmax), 3);
    }
    
    return basis;
  }
  
  // Compute the derivative of the basis w.r.t. u
  // This is ds/du evaluated at u
  // Returns vector of length n_knots (to be dotted with gamma)
  vector rcs_basis_deriv(real u, vector knots) {
    int K = num_elements(knots);
    vector[K] dbasis;
    real kmin = knots[1];
    real kmax = knots[K];
    real denom = kmax - kmin;
    
    // d/du of intercept = 0
    dbasis[1] = 0.0;
    // d/du of u = 1
    dbasis[2] = 1.0;
    
    // d/du of non-linear terms
    for (j in 1:(K - 2)) {
      real lambda_j = (kmax - knots[j + 1]) / denom;
      dbasis[j + 2] = 3.0 * pow(pos_part(u - knots[j + 1]), 2) * (u > knots[j + 1] ? 1.0 : 0.0)
                     - lambda_j * 3.0 * pow(pos_part(u - kmin), 2) * (u > kmin ? 1.0 : 0.0)
                     - (1.0 - lambda_j) * 3.0 * pow(pos_part(u - kmax), 2) * (u > kmax ? 1.0 : 0.0);
    }
    
    return dbasis;
  }
  
  // ---------------------------------------------------------------
  // Log-likelihood contribution for a single observation
  // ---------------------------------------------------------------
  // For patient i with:
  //   y_i     = event/censoring time
  //   delta_i = event indicator (1 = event, 0 = censored)
  //   t0_i    = delayed entry time (0 if no delayed entry)
  //   xbeta_i = X_i * beta (linear predictor)
  //
  // log H(t|X) = s(log(t * exp(-xbeta)) | gamma, knots)
  //            = s(log(t) - xbeta | gamma, knots)
  //
  // Let u = log(t) - xbeta = log(t * phi)
  //
  // h(t|X) = exp(s(u)) * s'(u) / t
  // S(t|X) = exp(-exp(s(u)))
  //
  // log-lik_i = delta_i * [s(u_i) + log(s'(u_i)) - log(y_i)]
  //           - exp(s(u_i))
  //           + exp(s(u0_i))   [delayed entry correction]
  //
  // where u_i = log(y_i) - xbeta_i
  //       u0_i = log(t0_i) - xbeta_i  (only if t0_i > 0)
  
  real fpaft_log_lik(real y, int delta, real t0, real xbeta,
                     vector gamma, vector knots) {
    real u = log(y) - xbeta;
    vector[num_elements(knots)] basis = rcs_basis(u, knots);
    vector[num_elements(knots)] dbasis = rcs_basis_deriv(u, knots);
    
    real s_val = dot_product(gamma, basis);    // s(u)
    real ds_val = dot_product(gamma, dbasis);  // s'(u)
    
    real ll = 0.0;
    
    // Event contribution
    if (delta == 1) {
      ll += s_val + log(ds_val) - log(y);
      // NOTE: ds_val = s'(u) must be > 0 for a valid hazard
      // This is a monotonicity constraint on the cumulative hazard
    }
    
    // Survival contribution
    ll -= exp(s_val);
    
    // Delayed entry correction
    if (t0 > 0.0) {
      real u0 = log(t0) - xbeta;
      vector[num_elements(knots)] basis0 = rcs_basis(u0, knots);
      real s_val0 = dot_product(gamma, basis0);
      ll += exp(s_val0);
    }
    
    return ll;
  }
}

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> P;                    // number of covariates
  matrix[N, P] X;                    // covariate matrix
  vector<lower=0>[N] y;             // event/censoring times
  array[N] int<lower=0, upper=1> d; // event indicators
  vector<lower=0>[N] t0;            // delayed entry times (0 = none)
  
  int<lower=2> n_knots;             // number of knots for baseline spline
  vector[n_knots] knots;            // knot locations (on log-time scale)
}

transformed data {
  int n_gamma = n_knots;  // number of baseline spline coefficients
}

parameters {
  vector[P] beta;              // log acceleration factors
  vector[n_gamma] gamma;       // baseline spline coefficients
}

model {
  // --- Priors ---
  // Weakly informative priors on covariate effects
  beta ~ normal(0, 2);
  
  // Prior on baseline spline coefficients
  // gamma[1] is the intercept (log-scale location)
  // gamma[2] is the linear slope (must be > 0 for monotonicity in simple cases)
  // Remaining gammas control curvature
  gamma[1] ~ normal(0, 5);
  gamma[2] ~ normal(1, 2);  // centered near positive value (Weibull-like)
  if (n_gamma > 2) {
    // Smoothness prior: penalize large non-linear spline coefficients
    // This acts as a roughness penalty similar to P-splines
    for (j in 3:n_gamma) {
      gamma[j] ~ normal(0, 1);
    }
  }
  
  // --- Likelihood ---
  {
    vector[N] xbeta = X * beta;
    for (i in 1:N) {
      target += fpaft_log_lik(y[i], d[i], t0[i], xbeta[i], gamma, knots);
    }
  }
}

generated quantities {
  // Log-likelihood for LOO-CV / WAIC
  vector[N] log_lik;
  {
    vector[N] xbeta = X * beta;
    for (i in 1:N) {
      log_lik[i] = fpaft_log_lik(y[i], d[i], t0[i], xbeta[i], gamma, knots);
    }
  }
}
