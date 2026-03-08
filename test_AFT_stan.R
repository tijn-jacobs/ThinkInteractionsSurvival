## Testing AFT models in stan

library(cmdstanr)
library(dplyr)
library(tidybayes)
library(ggplot2)


setwd("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival")

# -----------------------------
# Simulate AFT data
# -----------------------------
N <- 100
K <- 5

X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))  # intercept + covariates
beta_true <- rep(0, K)
beta_true[1:3] <- c(0.5, -1, 0.7)
sigma_true <- 0.5

# latent event times
logT <- X %*% beta_true + rnorm(N, 0, sigma_true)
T <- exp(logT)

# censoring times
C <- rexp(N, rate = 1/4)

Y <- pmin(T, C)
delta <- as.integer(T <= C)  # 1 = event, 0 = censored
mean(1 - delta)

# Stan input
stan_data <- list(
  N = N,
  K = K,
  X = X,
  log_y = as.numeric(log(Y)),
  delta = delta,
  err_dist = 4
)

# -----------------------------
# Compile the Stan model
# -----------------------------
mod <- cmdstan_model("aft_log_normal.stan")

# -----------------------------
# Sample
# -----------------------------
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2500,
  iter_sampling = 2500
)

print(fit$summary())

# Extract posterior draws
posterior <- fit$draws()





# Plot the posteriors

posterior_long <- fit$draws() %>%
  spread_draws(beta[k], sigma)

# Plot betas
posterior_long %>%
  ggplot(aes(x = beta, fill = factor(k))) +
  geom_density(alpha = 0.5) +
  labs(x = "beta", fill = "Coefficient") +
  theme_minimal()

# Plot sigma
posterior_long %>%
  ggplot(aes(x = sigma)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  labs(x = "sigma") +
  theme_minimal()
