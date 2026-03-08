# =============================================================
# 1. Load libraries
# =============================================================
library(cmdstanr)
library(posterior)
library(dplyr)

set.seed(123)

# =============================================================
# 2. Simulate data from time-dependent Weibull AFT
# =============================================================

# Sample size
N <- 500

# Covariates: intercept + treatment
x1 <- rbinom(N, 1, 0.5)
X  <- cbind(1, x1)   # design matrix

K <- ncol(X)

# True parameters
beta_true  <- c(0.3, -0.5)   # time-constant AFT effect
gamma_true <- c(0.0,  0.4)   # time-varying effect (linear in log time)

alpha_true  <- 1.2           # Weibull shape
lambda_true <- 2.0           # Weibull scale

# Reference time
t0 <- 1.0

# -----------------------------------------
# Define phi(x,t) and eta(x,t)
# -----------------------------------------

phi_fun <- function(Xrow, t) {
  xb_beta  <- sum(Xrow * beta_true)
  xb_gamma <- sum(Xrow * gamma_true)
  log_phi  <- -xb_beta - (log(t) - log(t0)) * xb_gamma
  exp(log_phi)
}

eta_fun <- function(Xrow, t) {
  xb_gamma <- sum(Xrow * gamma_true)
  phi_fun(Xrow, t) * (1 - xb_gamma)
}

# -----------------------------------------
# Generate event times T
# -----------------------------------------

T <- numeric(N)

for (i in 1:N) {
  
  # Draw uniform for inverse survival sampling
  u <- runif(1)
  
  # Target cumulative hazard
  target <- -log(u)
  
  # Function whose root gives T
  f <- function(t) {
    phi <- phi_fun(X[i,], t)
    z   <- (t * phi) / lambda_true
    H   <- z^alpha_true
    H - target
  }
  
  # Solve for T numerically
  T[i] <- uniroot(f, interval = c(1e-6, 1000))$root
}

# -----------------------------------------
# Add censoring
# -----------------------------------------

C <- rexp(N, rate = 0.1)   # censoring distribution

y     <- pmin(T, C)
delta <- as.integer(T <= C)

table(delta)  # check proportion censored

# -----------------------------------------
# Pack data for Stan
# -----------------------------------------

stan_data <- list(
  N = N,
  K = K,
  X = X,
  y = y,
  delta = delta,
  t0 = t0
)

# =============================================================
# 3. Compile the Stan model
# =============================================================

mod_td <- cmdstan_model("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival/Tijns code/aft_td_weibull.stan")

# =============================================================
# 4. Fit the model
# =============================================================

fit_td <- mod_td$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4
)

# =============================================================
# 5. Summaries
# =============================================================

print(
  fit_td$summary(c("beta","gamma","alpha","lambda")),
  digits = 3
)

# Extract posterior draws for further analysis
post <- as_draws_df(fit_td)

head(post)







library(tidybayes)
library(ggplot2)
library(dplyr)

draws_long <- fit_td$draws() %>% 
  spread_draws(beta[k], gamma[j], alpha, lambda)

# ---- Beta ----
draws_long %>%
  ggplot(aes(x = beta, fill = factor(k))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior of Beta",
    x = "beta",
    fill = "Coefficient index"
  ) +
  theme_minimal()

# ---- Gamma ----
draws_long %>%
  ggplot(aes(x = gamma, fill = factor(j))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior of Gamma (time-varying effects)",
    x = "gamma",
    fill = "Coefficient index"
  ) +
  theme_minimal()

# ---- Alpha and Lambda ----
draws_long %>%
  pivot_longer(cols = c(alpha, lambda), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior of Weibull parameters",
    x = "",
    fill = "Parameter"
  ) +
  theme_minimal()
