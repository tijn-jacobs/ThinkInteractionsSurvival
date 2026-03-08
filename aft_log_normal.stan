data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] log_y;                      // log of Y = min(T, C)
  array[N] int<lower=0, upper=1> delta; // 1 = event, 0 = censored

  int<lower=1, upper=4> err_dist;
  // 1 = normal, 2 = logistic, 3 = student_t, 4 = gumbel (Weibull on T)
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;                  // scale for normal/logistic/gumbel
  real<lower=2> nu;                     // df for Student-t (if used)
}

model {
  // Priors
  beta  ~ normal(0, 5);
  sigma ~ student_t(3, 0, 2.5);         // half-Student-t on scale
  nu    ~ gamma(2, 0.1);                // only really used when err_dist == 3

  for (n in 1:N) {
    if (err_dist == 1) {
      // Normal errors  -> log-normal AFT
      if (delta[n] == 1)
        target += normal_lpdf(log_y[n] | X[n] * beta, sigma);
      else
        target += normal_lccdf(log_y[n] | X[n] * beta, sigma);

    } else if (err_dist == 2) {
      // Logistic errors -> log-logistic AFT
      if (delta[n] == 1)
        target += logistic_lpdf(log_y[n] | X[n] * beta, sigma);
      else
        target += logistic_lccdf(log_y[n] | X[n] * beta, sigma);

    } else if (err_dist == 3) {
      // Student-t errors -> robust AFT
      if (delta[n] == 1)
        target += student_t_lpdf(log_y[n] | nu, X[n] * beta, sigma);
      else
        target += student_t_lccdf(log_y[n] | nu, X[n] * beta, sigma);

    } else if (err_dist == 4) {
      // Gumbel errors on log T -> Weibull AFT on T
      if (delta[n] == 1)
        target += gumbel_lpdf(log_y[n] | X[n] * beta, sigma);
      else
        target += gumbel_lccdf(log_y[n] | X[n] * beta, sigma);
    }
  }
}
