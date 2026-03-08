data {
  int<lower=1> N;                     // number of observations
  int<lower=1> K;                     // number of covariates
  matrix[N, K] X;                     // design matrix
  vector<lower=0>[N] y;               // observed times (event or censoring)
  array[N] int<lower=0, upper=1> delta; // 1 = event, 0 = censored
  real<lower=0> t0;                   // reference time (e.g. 1.0)
}

parameters {
  vector[K] beta;                     // AFT baseline effect
  vector[K] gamma;                    // time-varying effect (linear in log t)
  real<lower=0> alpha;                // Weibull shape
  real<lower=0> lambda;               // Weibull scale
}

transformed parameters {
  real log_t0 = log(t0);
}

model {
  // Priors (you can tighten/loosen these)
  beta   ~ normal(0, 1);
  gamma  ~ normal(0, 0.5);            // keep time-variation modest
  alpha  ~ lognormal(0, 0.5);
  lambda ~ lognormal(0, 1);

  // Likelihood
  for (n in 1:N) {
    real log_t = log(y[n]);
    real xb_beta  = dot_product(X[n], beta);
    real xb_gamma = dot_product(X[n], gamma);

    // log phi(X,t) = -X beta - (log t - log t0) X gamma
    real log_phi = -xb_beta - (log_t - log_t0) * xb_gamma;
    real phi = exp(log_phi);

    // instantaneous acceleration factor eta = phi * (1 - X gamma)
    real eta = phi * (1 - xb_gamma);

    // cumulative hazard H(t|X) = H0(t * phi)
    // using Weibull H0(u) = (u / lambda)^alpha
    real u = y[n] * phi;
    real log_z = log(u) - log(lambda); // z = u / lambda
    real H = exp(alpha * log_z);       // z^alpha

    if (delta[n] == 1) {
      // event: log f = log h - H
      // h0(u) = (alpha / lambda) * z^(alpha - 1)
      // log h0 = log(alpha) - log(lambda) + (alpha - 1) * log_z
      real log_h0 = log(alpha) - log(lambda) + (alpha - 1) * log_z;
      real log_h  = log_h0 + log(eta);  // h = h0 * eta
      target += log_h - H;
    } else {
      // censored: log S = -H
      target += -H;
    }
  }
}
