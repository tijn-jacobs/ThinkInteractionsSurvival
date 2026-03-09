// Flexible Parametric AFT Model with Time-Dependent Acceleration Factors
// Based on Crowther, Royston & Clements (2023), Biostatistics, Eqs 3.8-3.11
//
// Model:
//   log H(t|X) = s0(log(t * phi(X,t; beta)) | gamma0, k0)
//
// where the cumulative acceleration factor is now time-dependent:
//   phi(X,t; beta) = exp(-X*beta - sum_p x_p * s_p(log(t) | gamma_p, k_p))
//
// The instantaneous (interpretable) acceleration factor is recovered via:
//   eta(X,t) = d/dt [t * phi(X,t)] = phi + t * dphi/dt
//
// Hazard (Eq 3.10):
//   h(t|X) = exp(s0(w)) * s0'(w) * (1/t) * (1 - sum_p x_p * s_p'(log(t)))
//
// where w = log(t * phi(X,t))
//
// Survival (Eq 3.11):
//   S(t|X) = exp(-exp(s0(w)))
//
// IMPORTANT CONSTRAINT: The factor (1 - sum_p x_p * s_p'(log(t))) must be > 0
// for the hazard to be valid. This constrains the time-dependent spline coefficients.

functions {
  // ---------------------------------------------------------------
  // Restricted cubic spline basis and derivative (same as constant model)
  // ---------------------------------------------------------------
  real pos_part(real x) {
    return fmax(x, 0.0);
  }
  
  vector rcs_basis(real u, vector knots) {
    int K = num_elements(knots);
    vector[K] basis;
    real kmin = knots[1];
    real kmax = knots[K];
    real denom = kmax - kmin;
    
    basis[1] = 1.0;
    basis[2] = u;
    
    for (j in 1:(K - 2)) {
      real lambda_j = (kmax - knots[j + 1]) / denom;
      basis[j + 2] = pow(pos_part(u - knots[j + 1]), 3)
                     - lambda_j * pow(pos_part(u - kmin), 3)
                     - (1.0 - lambda_j) * pow(pos_part(u - kmax), 3);
    }
    
    return basis;
  }
  
  vector rcs_basis_deriv(real u, vector knots) {
    int K = num_elements(knots);
    vector[K] dbasis;
    real kmin = knots[1];
    real kmax = knots[K];
    real denom = kmax - kmin;
    
    dbasis[1] = 0.0;
    dbasis[2] = 1.0;
    
    for (j in 1:(K - 2)) {
      real lambda_j = (kmax - knots[j + 1]) / denom;
      dbasis[j + 2] = 3.0 * pow(pos_part(u - knots[j + 1]), 2) * (u > knots[j + 1] ? 1.0 : 0.0)
                     - lambda_j * 3.0 * pow(pos_part(u - kmin), 2) * (u > kmin ? 1.0 : 0.0)
                     - (1.0 - lambda_j) * 3.0 * pow(pos_part(u - kmax), 2) * (u > kmax ? 1.0 : 0.0);
    }
    
    return dbasis;
  }
  
  // ---------------------------------------------------------------
  // Spline function value: s(u) = gamma' * basis(u)
  // ---------------------------------------------------------------
  real rcs_eval(real u, vector gamma, vector knots) {
    return dot_product(gamma, rcs_basis(u, knots));
  }
  
  // Spline derivative value: s'(u) = gamma' * dbasis(u) 
  real rcs_deriv(real u, vector gamma, vector knots) {
    return dot_product(gamma, rcs_basis_deriv(u, knots));
  }
  
  // ---------------------------------------------------------------
  // Time-dependent FPAFT log-likelihood for one observation
  // ---------------------------------------------------------------
  // 
  // Parameters:
  //   y       : event/censoring time
  //   delta   : event indicator
  //   t0      : delayed entry time (0 if none)
  //   x       : covariate vector (length P)
  //   beta    : constant AF coefficients (length P)
  //   gamma0  : baseline spline coefficients (length n_knots0)
  //   knots0  : baseline knot locations
  //   
  //   For time-dependent effects (n_tde covariates):
  //   tde_idx     : indices into x for covariates with time-dependent effects
  //   gamma_tde   : concatenated TDE spline coefficients
  //   knots_tde   : TDE knot locations (assumed same for all TDE covariates for simplicity)
  //   n_knots_tde : number of TDE knots
  //   n_tde       : number of time-dependent covariates
  //
  // The cumulative AF is:
  //   phi(x,t) = exp(-x'beta - sum_p x[tde_idx[p]] * s_p(log(t)))
  //
  // Let w = log(t * phi(x,t)) = log(t) + log(phi(x,t))
  //       = log(t) - x'beta - sum_p x[tde_idx[p]] * s_p(log(t))
  //
  // Hazard:
  //   h(t) = exp(s0(w)) * s0'(w) * (1/t) * (1 - sum_p x[tde_idx[p]] * s_p'(log(t)))
  //
  // Note: the factor (1 - sum_p ...) comes from d/dt[log(t*phi)] via chain rule
  
  real fpaft_tde_log_lik(real y, int delta, real t0, vector x,
                          vector beta, 
                          vector gamma0, vector knots0,
                          array[] int tde_idx, int n_tde,
                          array[] vector gamma_tde,
                          vector knots_tde) {
    int P = num_elements(beta);
    real xbeta = dot_product(x, beta);
    real log_y = log(y);
    
    // Compute the time-dependent part of log(phi) at time y
    // log(phi(x,y)) = -xbeta - sum_p x[tde_p] * s_p(log(y))
    real tde_sum = 0.0;      // sum_p x_p * s_p(log(y))
    real tde_deriv_sum = 0.0; // sum_p x_p * s_p'(log(y))
    
    for (p in 1:n_tde) {
      real sp_val = rcs_eval(log_y, gamma_tde[p], knots_tde);
      real sp_deriv = rcs_deriv(log_y, gamma_tde[p], knots_tde);
      tde_sum += x[tde_idx[p]] * sp_val;
      tde_deriv_sum += x[tde_idx[p]] * sp_deriv;
    }
    
    // w = log(t * phi) = log(t) - xbeta - tde_sum
    real w = log_y - xbeta - tde_sum;
    
    // Baseline spline at w
    real s0_val = rcs_eval(w, gamma0, knots0);
    real s0_deriv = rcs_deriv(w, gamma0, knots0);
    
    // The hazard correction factor from Eq 3.10:
    // d/dt[log(t*phi)] = (1/t) * (1 - sum_p x_p * s_p'(log(t)))
    // So h(t) = exp(s0(w)) * s0'(w) * (1/t) * correction
    real correction = 1.0 - tde_deriv_sum;
    
    real ll = 0.0;
    
    if (delta == 1) {
      // log h(y) = s0(w) + log(s0'(w)) - log(y) + log(correction)
      // Both s0_deriv > 0 and correction > 0 required
      ll += s0_val + log(s0_deriv) - log_y + log(correction);
    }
    
    // Survival: log S(y) = -exp(s0(w))
    ll -= exp(s0_val);
    
    // Delayed entry
    if (t0 > 0.0) {
      real log_t0 = log(t0);
      real tde_sum_0 = 0.0;
      for (p in 1:n_tde) {
        tde_sum_0 += x[tde_idx[p]] * rcs_eval(log_t0, gamma_tde[p], knots_tde);
      }
      real w0 = log_t0 - xbeta - tde_sum_0;
      real s0_val_0 = rcs_eval(w0, gamma0, knots0);
      ll += exp(s0_val_0);
    }
    
    return ll;
  }
  
  // ---------------------------------------------------------------
  // Constant AF version (n_tde = 0) - for convenience / comparison
  // ---------------------------------------------------------------
  real fpaft_const_log_lik(real y, int delta, real t0, real xbeta,
                            vector gamma0, vector knots0) {
    real log_y = log(y);
    real w = log_y - xbeta;
    
    real s0_val = rcs_eval(w, gamma0, knots0);
    real s0_deriv = rcs_deriv(w, gamma0, knots0);
    
    real ll = 0.0;
    
    if (delta == 1) {
      ll += s0_val + log(s0_deriv) - log_y;
    }
    
    ll -= exp(s0_val);
    
    if (t0 > 0.0) {
      real w0 = log(t0) - xbeta;
      ll += exp(rcs_eval(w0, gamma0, knots0));
    }
    
    return ll;
  }
}

data {
  int<lower=1> N;                         // observations
  int<lower=1> P;                         // covariates
  matrix[N, P] X;                         // covariate matrix
  vector<lower=0>[N] y;                  // event/censoring times
  array[N] int<lower=0, upper=1> d;      // event indicators
  vector<lower=0>[N] t0;                 // delayed entry times
  
  // Baseline spline specification
  int<lower=2> n_knots0;                  // baseline knots
  vector[n_knots0] knots0;               // baseline knot locations
  
  // Time-dependent effect specification
  int<lower=0> n_tde;                     // number of covariates with TDE
  array[n_tde > 0 ? n_tde : 1] int<lower=1, upper=P> tde_idx;  // which covariates
  int<lower=2> n_knots_tde;              // TDE knots (same for all TDE covariates)
  vector[n_knots_tde] knots_tde;         // TDE knot locations
}

transformed data {
  int n_gamma0 = n_knots0;       // baseline spline parameters
  int n_gamma_tde = n_knots_tde; // TDE spline parameters per covariate
}

parameters {
  vector[P] beta;                                    // constant AF coefficients
  vector[n_gamma0] gamma0;                           // baseline spline coefficients
  array[n_tde > 0 ? n_tde : 0] vector[n_gamma_tde] gamma_tde;  // TDE spline coefficients
}

model {
  // --- Priors ---
  
  // Covariate effects
  beta ~ normal(0, 2);
  
  // Baseline spline
  gamma0[1] ~ normal(0, 5);        // intercept
  gamma0[2] ~ normal(1, 2);        // linear (favor positive = monotone)
  if (n_gamma0 > 2) {
    for (j in 3:n_gamma0) {
      gamma0[j] ~ normal(0, 1);    // smoothness prior
    }
  }
  
  // TDE spline coefficients
  // These represent DEPARTURES from the constant AF
  // Prior centered at 0 = no time-dependence (shrinkage toward constant AF)
  if (n_tde > 0) {
    for (p in 1:n_tde) {
      // The intercept and linear term of the TDE spline
      // are absorbed into beta, so we constrain:
      //   gamma_tde[p][1] = 0 (intercept, would be redundant with beta)
      //   gamma_tde[p][2] = free (linear time-dependence)
      //   gamma_tde[p][3:] = free (non-linear time-dependence)
      //
      // Actually, for identifiability, the TDE spline should have
      // NO intercept term (it's already in beta). We handle this by
      // zeroing the intercept contribution. But in our parameterization,
      // having the intercept free adds a constant to log(phi) which 
      // shifts the overall AF level — this IS distinct from beta because
      // it applies only to the TDE covariates.
      //
      // The cleaner approach: let gamma_tde include all terms but use
      // a prior centered at 0 for shrinkage. The intercept of the TDE
      // spline effectively modifies the main effect beta.
      
      gamma_tde[p] ~ normal(0, 0.5);  // tighter prior: shrink toward constant AF
    }
  }
  
  // --- Likelihood ---
  if (n_tde > 0) {
    for (i in 1:N) {
      target += fpaft_tde_log_lik(y[i], d[i], t0[i], to_vector(X[i]),
                                   beta, gamma0, knots0,
                                   tde_idx, n_tde, gamma_tde, knots_tde);
    }
  } else {
    // Fall back to constant AF (more efficient, no TDE overhead)
    vector[N] xbeta = X * beta;
    for (i in 1:N) {
      target += fpaft_const_log_lik(y[i], d[i], t0[i], xbeta[i], gamma0, knots0);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  
  if (n_tde > 0) {
    for (i in 1:N) {
      log_lik[i] = fpaft_tde_log_lik(y[i], d[i], t0[i], to_vector(X[i]),
                                      beta, gamma0, knots0,
                                      tde_idx, n_tde, gamma_tde, knots_tde);
    }
  } else {
    vector[N] xbeta = X * beta;
    for (i in 1:N) {
      log_lik[i] = fpaft_const_log_lik(y[i], d[i], t0[i], xbeta[i], gamma0, knots0);
    }
  }
}
