# =============================================================================
# Bayesian Flexible Parametric AFT Model - Step 1: Constant Acceleration Factor
# =============================================================================
# 
# This script:
# 1. Simulates survival data from a known mixture-Weibull AFT model
# 2. Computes knots for the restricted cubic spline
# 3. Fits the Bayesian FPAFT model using cmdstanr
# 4. Validates against the known truth
#
# Dependencies: cmdstanr, survival, rstpm2 (for comparison)

library(cmdstanr)
library(survival)
library(ggplot2)

# ============================================================
# PART 1: Data Simulation
# ============================================================

# We simulate from a mixture Weibull AFT, following the paper's approach.
# 
# Baseline survival (mixture Weibull):
#   S0(t) = p * exp(-lambda1 * t^gamma1) + (1-p) * exp(-lambda2 * t^gamma2)
#
# Under AFT with binary treatment X and log-AF = beta:
#   S(t|X) = S0(t * exp(-X*beta))
#
# We use the paper's Scenario 1 parameters for a complex baseline:
#   lambda1=0.1, gamma1=3, lambda2=0.1, gamma2=1.6, p=0.8

simulate_mixture_weibull_aft <- function(n, beta, 
                                          lambda1 = 0.1, gamma1 = 3,
                                          lambda2 = 0.1, gamma2 = 1.6,
                                          p = 0.8,
                                          max_time = 5, seed = 42) {
  set.seed(seed)
  
  # Binary treatment
  X <- rbinom(n, 1, 0.5)
  
  # Acceleration factor
  phi <- exp(-X * beta)
  
  # Simulate event times using inverse CDF via bisection
  # S(t|X) = p * exp(-lambda1 * (t*phi)^gamma1) + (1-p) * exp(-lambda2 * (t*phi)^gamma2)
  U <- runif(n)
  
  surv_func <- function(t, phi_i) {
    t_star <- t * phi_i
    p * exp(-lambda1 * t_star^gamma1) + (1 - p) * exp(-lambda2 * t_star^gamma2)
  }
  
  event_times <- numeric(n)
  for (i in 1:n) {
    # Bisection to find t such that S(t|X_i) = U_i
    lo <- 1e-10
    hi <- 100
    for (iter in 1:200) {
      mid <- (lo + hi) / 2
      if (surv_func(mid, phi[i]) > U[i]) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    event_times[i] <- (lo + hi) / 2
  }
  
  # Administrative censoring at max_time
  obs_time <- pmin(event_times, max_time)
  event <- as.integer(event_times <= max_time)
  
  data.frame(
    time = obs_time,
    event = event,
    X = X,
    true_time = event_times
  )
}

# Simulate data
# Paper uses beta = 0.5 (log acceleration factor)
# This means AF = exp(0.5) ≈ 1.65 for X=1 vs X=0
# In AFT: phi = exp(-X*beta) = exp(-0.5) for treated
# So treated have time scaled by exp(-0.5) ≈ 0.607 (shorter survival)
# But note the sign convention: the paper's beta = -0.5 means 
# log(T) = ... - beta*X + eps, so negative beta => shorter time.
# Let's use beta = -0.5 (treated live shorter, AF = exp(-0.5))
# Actually, let's be very careful with signs.
#
# The paper defines: S(t|X) = S0(t * phi(X;beta)) with phi = exp(-X*beta)
# So if beta = -0.5 and X=1: phi = exp(0.5) > 1, meaning time is "sped up" 
# => shorter survival for treated.
# If beta = 0.5 and X=1: phi = exp(-0.5) < 1, meaning time is "slowed down"
# => longer survival for treated.
#
# We'll use beta = 0.5 (positive, treated have longer survival)
# and also beta = -0.5 (treated have shorter survival) as in the paper.

true_beta <- -0.5
n <- 200

dat <- simulate_mixture_weibull_aft(n = n, beta = true_beta, seed = 123)

cat("Data summary:\n")
cat(sprintf("  N = %d, Events = %d (%.1f%%)\n", 
            nrow(dat), sum(dat$event), 100*mean(dat$event)))
cat(sprintf("  True log-AF (beta) = %.2f\n", true_beta))
cat(sprintf("  Median time (X=0): %.3f\n", median(dat$time[dat$X == 0])))
cat(sprintf("  Median time (X=1): %.3f\n", median(dat$time[dat$X == 1])))

# ============================================================
# PART 2: Knot Placement
# ============================================================

# Following the paper: knots are placed at quantiles of log(t) 
# among uncensored observations.
# 
# For df = d, we use d+1 knots:
#   - Boundary knots at min and max of uncensored log(t)
#   - d-1 internal knots at equally-spaced centiles
#
# Standard choices: df=3 => 4 knots, df=5 => 6 knots

compute_knots <- function(times, events, df = 3) {
  # Use uncensored log-times for knot placement
  log_t_uncensored <- log(times[events == 1])
  
  n_knots <- df + 1
  n_internal <- df - 1
  
  if (n_internal == 0) {
    # df=1: only boundary knots (Weibull equivalent)
    knot_locations <- quantile(log_t_uncensored, c(0, 1))
  } else {
    # Centile positions for internal knots
    internal_centiles <- seq(0, 1, length.out = n_internal + 2)
    internal_centiles <- internal_centiles[-c(1, length(internal_centiles))]
    
    boundary <- quantile(log_t_uncensored, c(0, 1))
    internal <- quantile(log_t_uncensored, internal_centiles)
    
    knot_locations <- c(boundary[1], internal, boundary[2])
  }
  
  names(knot_locations) <- NULL
  knot_locations
}

# Note: The paper places knots on log(t * phi), but since phi depends
# on the unknown beta, they iterate or just use log(t) for uncensored obs.
# For initial knot placement, using log(t) is standard practice.

df_baseline <- 3  # degrees of freedom for baseline
knots <- compute_knots(dat$time, dat$event, df = df_baseline)
n_knots <- length(knots)

cat(sprintf("\nBaseline spline: df = %d, n_knots = %d\n", df_baseline, n_knots))
cat("Knot locations (log-time scale):", round(knots, 3), "\n")

# ============================================================
# PART 3: Prepare Data and Fit Stan Model
# ============================================================

stan_data <- list(
  N = nrow(dat),
  P = 1,
  X = matrix(dat$X, ncol = 1),
  y = dat$time,
  d = dat$event,
  t0 = rep(0, nrow(dat)),  # no delayed entry
  n_knots = n_knots,
  knots = knots
)

# Compile model
mod <- cmdstan_model("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival/Tijns code/phase 2/fpaft_constant.stan")

# Fit
fit <- mod$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 500
)

# ============================================================
# PART 4: Results
# ============================================================

cat("\n========================================\n")
cat("POSTERIOR SUMMARY\n")
cat("========================================\n")

# Beta (log acceleration factor)
fit$summary(variables = c("beta", "gamma"))

# Print key results
beta_draws <- fit$draws("beta[1]", format = "matrix")
cat(sprintf("\nTrue beta:      %.3f\n", true_beta))
cat(sprintf("Posterior mean: %.3f\n", mean(beta_draws)))
cat(sprintf("Posterior sd:   %.3f\n", sd(beta_draws)))
cat(sprintf("95%% CI:        [%.3f, %.3f]\n", 
            quantile(beta_draws, 0.025), quantile(beta_draws, 0.975)))

# Check: does the 95% CI contain the truth?
ci <- quantile(beta_draws, c(0.025, 0.975))
cat(sprintf("Truth in 95%% CI: %s\n", 
            ifelse(true_beta >= ci[1] & true_beta <= ci[2], "YES", "NO")))

# ============================================================
# PART 5: Diagnostics
# ============================================================

cat("\n========================================\n")
cat("DIAGNOSTICS\n")
cat("========================================\n")

diag_summary <- fit$summary(variables = c("beta", "gamma"))
print(diag_summary[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")])

# Check for divergences
cat(sprintf("\nDivergent transitions: %d\n", sum(fit$diagnostic_summary()$num_divergent)))

# ============================================================
# PART 6: Compare with frequentist (if rstpm2 available)
# ============================================================

# Try rstpm2 for comparison
tryCatch({
  library(rstpm2)
  
  surv_obj <- Surv(dat$time, dat$event)
  
  # Fit flexible parametric AFT
  freq_fit <- aft(surv_obj ~ X, data = dat, df = df_baseline)
  
  cat("\n========================================\n")
  cat("FREQUENTIST COMPARISON (rstpm2::aft)\n")
  cat("========================================\n")
  cat(sprintf("rstpm2 beta estimate: %.3f (SE: %.3f)\n",
              coef(freq_fit)["X"], sqrt(vcov(freq_fit)["X", "X"])))
  
}, error = function(e) {
  cat("\nrstpm2 not available for comparison. Install with: install.packages('rstpm2')\n")
})

# ============================================================
# PART 7: Posterior Predictive Checks
# ============================================================

# Compute predicted survival curves and compare to Kaplan-Meier

# Function to compute survival at time t given parameters
compute_survival <- function(t, xbeta, gamma, knots) {
  u <- log(t) - xbeta
  
  K <- length(knots)
  kmin <- knots[1]
  kmax <- knots[K]
  denom <- kmax - kmin
  
  # Compute basis
  basis <- numeric(K)
  basis[1] <- 1
  basis[2] <- u
  for (j in 1:(K-2)) {
    lambda_j <- (kmax - knots[j+1]) / denom
    basis[j+2] <- max(u - knots[j+1], 0)^3 - 
                  lambda_j * max(u - kmin, 0)^3 - 
                  (1 - lambda_j) * max(u - kmax, 0)^3
  }
  
  s_val <- sum(gamma * basis)
  exp(-exp(s_val))
}

# Posterior mean parameters
gamma_post <- colMeans(fit$draws("gamma", format = "matrix"))
beta_post <- mean(beta_draws)

# Compute predicted survival on a grid
t_grid <- seq(0.01, 5, length.out = 200)

surv_pred_x0 <- sapply(t_grid, function(t) {
  compute_survival(t, xbeta = 0 * beta_post, gamma = gamma_post, knots = knots)
})

surv_pred_x1 <- sapply(t_grid, function(t) {
  compute_survival(t, xbeta = 1 * beta_post, gamma = gamma_post, knots = knots)
})

# Kaplan-Meier
km_fit <- survfit(Surv(time, event) ~ X, data = dat)

# Save results for plotting
results <- list(
  t_grid = t_grid,
  surv_pred_x0 = surv_pred_x0,
  surv_pred_x1 = surv_pred_x1,
  km_fit = km_fit,
  fit = fit,
  stan_data = stan_data,
  true_beta = true_beta
)

save(results, file = "fpaft_constant_results.RData")

cat("\n========================================\n")
cat("DONE. Results saved to fpaft_constant_results.RData\n")
cat("========================================\n")

# ============================================================
# PART 8: Quick survival plot
# ============================================================

# Create a simple comparison plot
pdf("fpaft_constant_survplot.pdf", width = 8, height = 6)

plot(km_fit, col = c("blue", "red"), lwd = 2, 
     xlab = "Time", ylab = "Survival probability",
     main = sprintf("Bayesian FPAFT (df=%d) vs Kaplan-Meier\nTrue beta = %.2f, Posterior mean = %.3f",
                    df_baseline, true_beta, beta_post))

lines(t_grid, surv_pred_x0, col = "blue", lty = 2, lwd = 2)
lines(t_grid, surv_pred_x1, col = "red", lty = 2, lwd = 2)

legend("bottomleft", 
       c("KM (X=0)", "KM (X=1)", "FPAFT (X=0)", "FPAFT (X=1)"),
       col = c("blue", "red", "blue", "red"),
       lty = c(1, 1, 2, 2), lwd = 2)

dev.off()

cat("Survival plot saved to fpaft_constant_survplot.pdf\n")
