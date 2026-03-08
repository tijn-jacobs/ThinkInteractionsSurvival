functions {
  real pos(real x) {
    return x > 0 ? x : 0;
  }

  real cube(real x) {
    return x * x * x;
  }

  // Harrell-style restricted cubic spline basis component h_j(u)
  // knots must be sorted ascending, length Kk >= 3
  real rcs_h(real u, vector knots, int j) {
    int Kk = num_elements(knots);
    real kj = knots[j];
    real kKm1 = knots[Kk - 1];
    real kK   = knots[Kk];
    real denom = kK - kKm1;

    return cube(pos(u - kj))
      - cube(pos(u - kKm1)) * (kK - kj) / denom
      + cube(pos(u - kK))   * (kKm1 - kj) / denom;
  }

  // derivative wrt u
  real rcs_h_prime(real u, vector knots, int j) {
    int Kk = num_elements(knots);
    real kj = knots[j];
    real kKm1 = knots[Kk - 1];
    real kK   = knots[Kk];
    real denom = kK - kKm1;

    return 3 * square(pos(u - kj))
      - 3 * square(pos(u - kKm1)) * (kK - kj) / denom
      + 3 * square(pos(u - kK))   * (kKm1 - kj) / denom;
  }

  // Evaluate s(u) and s'(u) for restricted cubic spline:
  // s(u) = gamma[1] + gamma[2]*u + sum_{j=1}^{Kk-2} gamma[j+2]*h_j(u)
  // gamma length must be Kk: intercept + linear + (Kk-2) nonlinear terms
  real s_rcs(real u, vector knots, vector gamma) {
    int Kk = num_elements(knots);
    real out = gamma[1] + gamma[2] * u;

    for (j in 1:(Kk - 2)) {
      out += gamma[j + 2] * rcs_h(u, knots, j);
    }
    return out;
  }

  real s_rcs_prime(real u, vector knots, vector gamma) {
    int Kk = num_elements(knots);
    real out = gamma[2];

    for (j in 1:(Kk - 2)) {
      out += gamma[j + 2] * rcs_h_prime(u, knots, j);
    }
    return out;
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector<lower=0>[N] y;                       // observed time
  array[N] int<lower=0, upper=1> delta;       // 1=event, 0=censored

  int<lower=3> Kk;                            // number of knots
  vector[Kk] knots;                           // knots on u-scale (sorted)
}

parameters {
  vector[K] beta;                             // AFT coefficients
  vector[Kk] gamma;                           // spline coeffs for log H
}

model {
  // Simple priors (start here; extend later)
  beta  ~ normal(0, 2);

  // gamma[1] is intercept on log cumulative hazard scale
  gamma[1] ~ normal(0, 5);

  // remaining spline terms: modest regularization
  gamma[2:Kk] ~ normal(0, 2);

  for (n in 1:N) {
    real xb = dot_product(X[n], beta);
    real u  = log(y[n]) - xb;                 // u = log(t*phi) with phi=exp(-xb)

    real s      = s_rcs(u, knots, gamma);
    real s_prime = s_rcs_prime(u, knots, gamma);

    // Cumulative hazard H = exp(s)
    real H = exp(s);

    // Enforce positive hazard: need s'(u) > 0
    if (s_prime <= 0) {
      target += negative_infinity();
    } else {
      if (delta[n] == 1) {
        // log h = s + log(s') - log(t)
        target += (s + log(s_prime) - log(y[n])) - H;
      } else {
        target += -H;
      }
    }
  }
}
