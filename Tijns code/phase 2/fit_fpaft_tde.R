# =============================================================================
# Bayesian FPAFT Model - Step 2: Time-Dependent Acceleration Factor
# =============================================================================
#
# We test the TDE model in two ways:
#   A) Fit to data with NO time-dependence (should recover constant AF, 
#      TDE coefficients should shrink to ~0)
#   B) Fit to data WITH time-dependence (should recover the time-varying AF)
#
# For (B), we simulate data where the log-AF changes over time.

library(cmdstanr)
library(survival)

# ============================================================
# Helper functions (reused from constant model)
# ============================================================

simulate_mixture_weibull_aft <- function(n, beta, 
                                          lambda1 = 0.1, gamma1 = 3,
                                          lambda2 = 0.1, gamma2 = 1.6,
                                          p = 0.8,
                                          max_time = 5, seed = 42) {
  set.seed(seed)
  X <- rbinom(n, 1, 0.5)
  phi <- exp(-X * beta)
  U <- runif(n)
  
  surv_func <- function(t, phi_i) {
    t_star <- t * phi_i
    p * exp(-lambda1 * t_star^gamma1) + (1 - p) * exp(-lambda2 * t_star^gamma2)
  }
  
  event_times <- numeric(n)
  for (i in 1:n) {
    lo <- 1e-10; hi <- 100
    for (iter in 1:200) {
      mid <- (lo + hi) / 2
      if (surv_func(mid, phi[i]) > U[i]) lo <- mid else hi <- mid
    }
    event_times[i] <- (lo + hi) / 2
  }
  
  obs_time <- pmin(event_times, max_time)
  event <- as.integer(event_times <= max_time)
  data.frame(time = obs_time, event = event, X = X)
}

# Simulate data with a TIME-DEPENDENT acceleration factor
# 
# Here we simulate from the general framework:
#   S(t|X) = S0(t * phi(X, t))
#   phi(X, t) = exp(-X * beta(t))
#   beta(t) = beta_const + beta_td * f(t)
#
# For simplicity, we use a Weibull baseline (so we can use the 
# closed-form inverse CDF approach) and a linear time-dependent AF:
#   beta(t) = beta0 + beta1 * log(t)
#
# This means:
#   phi(X, t) = exp(-X * (beta0 + beta1 * log(t)))
#             = exp(-X * beta0) * t^(-X * beta1)
#
# And S(t|X) = S0(t * phi(X,t))
#            = exp(-lambda * (t * exp(-X*beta0) * t^(-X*beta1))^gamma)
#            = exp(-lambda * t^(gamma*(1-X*beta1)) * exp(-X*beta0*gamma))
#
# This is still a Weibull-type model but with modified shape by group.

simulate_tde_weibull_aft <- function(n, beta0, beta1,
                                      lambda = 0.3, gamma_shape = 1.5,
                                      max_time = 5, seed = 42) {
  set.seed(seed)
  X <- rbinom(n, 1, 0.5)
  U <- runif(n)
  
  # For X=0: S(t) = exp(-lambda * t^gamma_shape) => standard Weibull
  # For X=1: S(t) = exp(-lambda * exp(-beta0*gamma_shape) * t^(gamma_shape*(1-beta1)))
  #
  # Inverse CDF: t = (-log(U) / (lambda * c))^(1/g)
  # where for X=1: c = exp(-beta0*gamma_shape), g = gamma_shape*(1-beta1)
  
  event_times <- numeric(n)
  for (i in 1:n) {
    if (X[i] == 0) {
      event_times[i] <- (-log(U[i]) / lambda)^(1/gamma_shape)
    } else {
      c_i <- exp(-beta0 * gamma_shape)
      g_i <- gamma_shape * (1 - beta1)
      event_times[i] <- (-log(U[i]) / (lambda * c_i))^(1/g_i)
    }
  }
  
  obs_time <- pmin(event_times, max_time)
  event <- as.integer(event_times <= max_time)
  
  data.frame(time = obs_time, event = event, X = X)
}

compute_knots <- function(times, events, df = 3) {
  log_t_uncensored <- log(times[events == 1])
  n_knots <- df + 1
  n_internal <- df - 1
  if (n_internal == 0) {
    knot_locations <- quantile(log_t_uncensored, c(0, 1))
  } else {
    internal_centiles <- seq(0, 1, length.out = n_internal + 2)
    internal_centiles <- internal_centiles[-c(1, length(internal_centiles))]
    boundary <- quantile(log_t_uncensored, c(0, 1))
    internal <- quantile(log_t_uncensored, internal_centiles)
    knot_locations <- c(boundary[1], internal, boundary[2])
  }
  names(knot_locations) <- NULL
  knot_locations
}

# ============================================================
# PART A: Fit TDE model to constant-AF data (sanity check)
# ============================================================

cat("===========================================\n")
cat("PART A: TDE model on constant-AF data\n")
cat("(TDE coefficients should shrink to ~0)\n")
cat("===========================================\n\n")

true_beta <- -0.5
dat_const <- simulate_mixture_weibull_aft(n = 1000, beta = true_beta, seed = 123)

df_baseline <- 5  # use df=5 for better baseline fit
df_tde <- 2       # simple TDE: 3 knots = intercept + linear + 1 nonlinear

knots0 <- compute_knots(dat_const$time, dat_const$event, df = df_baseline)
knots_tde <- compute_knots(dat_const$time, dat_const$event, df = df_tde)

cat(sprintf("Baseline: df=%d, %d knots\n", df_baseline, length(knots0)))
cat(sprintf("TDE:      df=%d, %d knots\n", df_tde, length(knots_tde)))

mod_tde <- cmdstan_model("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival/Tijns code/phase 2/fpaft_tde.stan")

stan_data_a <- list(
  N = nrow(dat_const),
  P = 1,
  X = matrix(dat_const$X, ncol = 1),
  y = dat_const$time,
  d = dat_const$event,
  t0 = rep(0, nrow(dat_const)),
  n_knots0 = length(knots0),
  knots0 = knots0,
  n_tde = 1,
  tde_idx = array(1, dim = 1),  # X (column 1) has TDE
  n_knots_tde = length(knots_tde),
  knots_tde = knots_tde
)

fit_a <- mod_tde$sample(
  data = stan_data_a,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 500
)

cat("\n--- Part A Results ---\n")
print(fit_a$summary(variables = c("beta", "gamma0", "gamma_tde")))

beta_a <- fit_a$draws("beta[1]", format = "matrix")
cat(sprintf("\nTrue beta (constant): %.3f\n", true_beta))
cat(sprintf("Posterior mean beta:  %.3f (SD: %.3f)\n", mean(beta_a), sd(beta_a)))

# TDE coefficients should be near zero
gamma_tde_a <- fit_a$draws("gamma_tde", format = "matrix")
cat(sprintf("\nTDE coefficients (should be ~0):\n"))
for (j in 1:ncol(gamma_tde_a)) {
  cat(sprintf("  gamma_tde[1,%d]: mean=%.4f, sd=%.4f\n", 
              j, mean(gamma_tde_a[,j]), sd(gamma_tde_a[,j])))
}

cat(sprintf("\nDivergences: %d\n", sum(fit_a$diagnostic_summary()$num_divergent)))

# ============================================================
# PART B: Fit TDE model to time-dependent AF data
# ============================================================

cat("\n\n===========================================\n")
cat("PART B: TDE model on time-dependent AF data\n")
cat("===========================================\n\n")

# Simulate with time-dependent AF
# beta(t) = beta0 + beta1 * log(t)
# beta0 = -0.3, beta1 = 0.2
# So early on (t small, log(t) negative): AF more negative (faster event)
# Later (t large, log(t) positive): AF less negative (attenuation)
# This mirrors the breast cancer example in the paper (Fig 5)

true_beta0 <- -0.3
true_beta1 <- 0.2

dat_tde <- simulate_tde_weibull_aft(n = 1500, beta0 = true_beta0, beta1 = true_beta1,
                                     lambda = 0.3, gamma_shape = 1.5, seed = 456)

cat(sprintf("Data: N=%d, Events=%d (%.1f%%)\n", 
            nrow(dat_tde), sum(dat_tde$event), 100*mean(dat_tde$event)))

# Knots
knots0_b <- compute_knots(dat_tde$time, dat_tde$event, df = 5)
knots_tde_b <- compute_knots(dat_tde$time, dat_tde$event, df = 2)

cat(sprintf("Baseline: df=5, %d knots\n", length(knots0_b)))
cat(sprintf("TDE:      df=2, %d knots\n", length(knots_tde_b)))

stan_data_b <- list(
  N = nrow(dat_tde),
  P = 1,
  X = matrix(dat_tde$X, ncol = 1),
  y = dat_tde$time,
  d = dat_tde$event,
  t0 = rep(0, nrow(dat_tde)),
  n_knots0 = length(knots0_b),
  knots0 = knots0_b,
  n_tde = 1,
  tde_idx = array(1, dim = 1),
  n_knots_tde = length(knots_tde_b),
  knots_tde = knots_tde_b
)

fit_b <- mod_tde$sample(
  data = stan_data_b,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 500
)

cat("\n--- Part B Results ---\n")
print(fit_b$summary(variables = c("beta", "gamma0", "gamma_tde")))

beta_b <- fit_b$draws("beta[1]", format = "matrix")
gamma_tde_b <- fit_b$draws("gamma_tde", format = "matrix")

cat(sprintf("\nTrue beta0 (constant part): %.3f\n", true_beta0))
cat(sprintf("True beta1 (log(t) slope):  %.3f\n", true_beta1))
cat(sprintf("Posterior mean beta: %.3f (SD: %.3f)\n", mean(beta_b), sd(beta_b)))

cat(sprintf("\nTDE coefficients:\n"))
for (j in 1:ncol(gamma_tde_b)) {
  cat(sprintf("  gamma_tde[1,%d]: mean=%.4f, sd=%.4f\n", 
              j, mean(gamma_tde_b[,j]), sd(gamma_tde_b[,j])))
}

cat(sprintf("\nDivergences: %d\n", sum(fit_b$diagnostic_summary()$num_divergent)))

# ============================================================
# PART C: Compute and plot the estimated time-dependent AF
# ============================================================

# The instantaneous AF at time t is:
#   eta(X=1, t) = phi(t) + t * dphi/dt
# where
#   phi(t) = exp(-beta - s_tde(log(t)))
#   dphi/dt = phi(t)/t * (-s_tde'(log(t)))
# So:
#   eta(t) = phi(t) * (1 - s_tde'(log(t)))
#
# The acceleration factor (as reported) is exp(-beta(t)):
#   AF(t) = eta(t) for X=1 vs X=0
#
# For the true model: AF(t) = exp(-(beta0 + beta1*log(t))) 
# which in terms of eta:
#   True cumulative: t * phi(t) = t * exp(-beta0) * t^(-beta1) = exp(-beta0) * t^(1-beta1)
#   True eta(t) = d/dt[exp(-beta0) * t^(1-beta1)] = exp(-beta0) * (1-beta1) * t^(-beta1)

# R functions for RCS
rcs_basis_r <- function(u, knots) {
  K <- length(knots)
  kmin <- knots[1]
  kmax <- knots[K]
  denom <- kmax - kmin
  basis <- numeric(K)
  basis[1] <- 1
  basis[2] <- u
  for (j in 1:(K-2)) {
    lj <- (kmax - knots[j+1]) / denom
    basis[j+2] <- max(u - knots[j+1], 0)^3 - 
                  lj * max(u - kmin, 0)^3 - 
                  (1-lj) * max(u - kmax, 0)^3
  }
  basis
}

rcs_deriv_r <- function(u, knots) {
  K <- length(knots)
  kmin <- knots[1]
  kmax <- knots[K]
  denom <- kmax - kmin
  dbasis <- numeric(K)
  dbasis[1] <- 0
  dbasis[2] <- 1
  for (j in 1:(K-2)) {
    lj <- (kmax - knots[j+1]) / denom
    dbasis[j+2] <- 3 * max(u - knots[j+1], 0)^2 * (u > knots[j+1]) -
                   lj * 3 * max(u - kmin, 0)^2 * (u > kmin) -
                   (1-lj) * 3 * max(u - kmax, 0)^2 * (u > kmax)
  }
  dbasis
}

# Posterior means
beta_post <- mean(beta_b)
gamma_tde_post <- colMeans(gamma_tde_b)

# Compute estimated AF on a time grid
t_grid <- seq(0.1, 5, length.out = 200)

# Estimated instantaneous AF: eta(t) for X=1
# phi(t) = exp(-beta - s_tde(log(t)))
# eta(t) = phi(t) * (1 - s_tde'(log(t)))
eta_est <- sapply(t_grid, function(t) {
  log_t <- log(t)
  s_val <- sum(gamma_tde_post * rcs_basis_r(log_t, knots_tde_b))
  s_deriv <- sum(gamma_tde_post * rcs_deriv_r(log_t, knots_tde_b))
  phi <- exp(-beta_post - s_val)
  phi * (1 - s_deriv)
})

# True instantaneous AF
# eta_true(t) = exp(-beta0) * (1 - beta1) * t^(-beta1)
eta_true <- exp(-true_beta0) * (1 - true_beta1) * t_grid^(-true_beta1)

# Also compute posterior credible intervals for eta
n_draws <- nrow(beta_b)
n_sample <- min(500, n_draws)  # subsample for speed
idx <- sample(n_draws, n_sample)

eta_draws <- matrix(NA, n_sample, length(t_grid))
for (s in 1:n_sample) {
  b <- beta_b[idx[s], 1]
  g <- gamma_tde_b[idx[s], ]
  for (j in seq_along(t_grid)) {
    log_t <- log(t_grid[j])
    s_val <- sum(g * rcs_basis_r(log_t, knots_tde_b))
    s_deriv <- sum(g * rcs_deriv_r(log_t, knots_tde_b))
    phi <- exp(-b - s_val)
    eta_draws[s, j] <- phi * (1 - s_deriv)
  }
}

eta_lo <- apply(eta_draws, 2, quantile, 0.025)
eta_hi <- apply(eta_draws, 2, quantile, 0.975)

# Plot
pdf("fpaft_tde_af_plot.pdf", width = 9, height = 6)

plot(t_grid, eta_est, type = "l", col = "blue", lwd = 2,
     xlab = "Follow-up time", ylab = "Acceleration factor (eta)",
     main = "Time-Dependent Acceleration Factor\nBayesian FPAFT vs Truth",
     ylim = range(c(eta_lo, eta_hi, eta_true)))

polygon(c(t_grid, rev(t_grid)), c(eta_lo, rev(eta_hi)),
        col = rgb(0, 0, 1, 0.15), border = NA)

lines(t_grid, eta_true, col = "red", lwd = 2, lty = 2)
abline(h = 1, lty = 3, col = "gray50")

legend("topright",
       c("Estimated AF (posterior mean)", "95% credible interval", "True AF"),
       col = c("blue", rgb(0,0,1,0.3), "red"),
       lty = c(1, NA, 2), lwd = c(2, NA, 2),
       fill = c(NA, rgb(0,0,1,0.15), NA),
       border = c(NA, NA, NA))

dev.off()

cat("\nTime-dependent AF plot saved to fpaft_tde_af_plot.pdf\n")

# ============================================================
# PART D: Model comparison via LOO
# ============================================================

cat("\n--- Model Comparison (LOO-CV) ---\n")

tryCatch({
  library(loo)
  
  # Also fit constant AF model for comparison
  stan_data_const <- list(
    N = nrow(dat_tde),
    P = 1,
    X = matrix(dat_tde$X, ncol = 1),
    y = dat_tde$time,
    d = dat_tde$event,
    t0 = rep(0, nrow(dat_tde)),
    n_knots0 = length(knots0_b),
    knots0 = knots0_b,
    n_tde = 0,
    tde_idx = array(1, dim = 1),  # dummy, won't be used
    n_knots_tde = length(knots_tde_b),
    knots_tde = knots_tde_b
  )
  
  fit_const <- mod_tde$sample(
    data = stan_data_const,
    seed = 42,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 0
  )
  
  loo_tde <- loo(fit_b$draws("log_lik"))
  loo_const <- loo(fit_const$draws("log_lik"))
  
  comp <- loo_compare(loo_const, loo_tde)
  print(comp)
  
}, error = function(e) {
  cat("loo package not available. Install with: install.packages('loo')\n")
  cat("Error:", conditionMessage(e), "\n")
})

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
