# =============================================================================
# Quick df comparison for Bayesian FPAFT (constant AF)
# Run after fit_fpaft_constant.R to reuse simulation + model compilation
# =============================================================================

library(cmdstanr)
library(survival)

# --- Reuse data simulation from the main script ---
source_env <- new.env()

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

true_beta <- -0.5
n <- 1000
dat <- simulate_mixture_weibull_aft(n = n, beta = true_beta, seed = 123)

# Compile model once
mod <- cmdstan_model("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival/Tijns code/phase 2/fpaft_constant.stan")

# --- Fit for df = 3, 5, 7 ---
results_df <- list()

for (df_val in c(3, 5, 7)) {
  cat(sprintf("\n===== Fitting df = %d =====\n", df_val))
  
  knots <- compute_knots(dat$time, dat$event, df = df_val)
  cat("Knots:", round(knots, 3), "\n")
  
  stan_data <- list(
    N = nrow(dat),
    P = 1,
    X = matrix(dat$X, ncol = 1),
    y = dat$time,
    d = dat$event,
    t0 = rep(0, nrow(dat)),
    n_knots = length(knots),
    knots = knots
  )
  
  fit <- mod$sample(
    data = stan_data,
    seed = 42,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 0
  )
  
  beta_draws <- fit$draws("beta[1]", format = "matrix")
  
  cat(sprintf("  Posterior mean: %.4f\n", mean(beta_draws)))
  cat(sprintf("  Posterior sd:   %.4f\n", sd(beta_draws)))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
              quantile(beta_draws, 0.025), quantile(beta_draws, 0.975)))
  cat(sprintf("  Bias: %.4f\n", mean(beta_draws) - true_beta))
  cat(sprintf("  Divergences: %d\n", sum(fit$diagnostic_summary()$num_divergent)))
  
  results_df[[as.character(df_val)]] <- list(
    df = df_val,
    mean = mean(beta_draws),
    sd = sd(beta_draws),
    ci = quantile(beta_draws, c(0.025, 0.975)),
    bias = mean(beta_draws) - true_beta,
    fit = fit
  )
}

# --- Summary table ---
cat("\n\n==========================================\n")
cat("SUMMARY: Effect of df on beta estimation\n")
cat("==========================================\n")
cat(sprintf("True beta = %.3f\n\n", true_beta))
cat(sprintf("%-6s %-10s %-10s %-10s %-22s\n", "df", "Mean", "SD", "Bias", "95% CI"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (nm in names(results_df)) {
  r <- results_df[[nm]]
  cat(sprintf("%-6d %-10.4f %-10.4f %-10.4f [%-8.4f, %-8.4f]\n",
              r$df, r$mean, r$sd, r$bias, r$ci[1], r$ci[2]))
}

