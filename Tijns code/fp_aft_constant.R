library(cmdstanr)

# Inputs you must already have:
# N, K, X, y, delta

# Choose knots on the u-scale.
# A practical default: knots at quantiles of log(y).
# (In the constant AFT model, u is log(y) - X beta; beta unknown,
# so start with log(y) quantiles.)
Kk <- 5
knots <- as.numeric(quantile(log(y), probs = seq(0.05, 0.95, length.out = Kk)))

stan_data <- list(
  N = N,
  K = K,
  X = X,
  y = as.numeric(y),
  delta = as.integer(delta),
  Kk = Kk,
  knots = knots
)

mod <- cmdstan_model("~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ThinkInteractionsSurvival/Tijns code/fp_aft_constant.stan")

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(fit$summary(c("beta", "gamma")))











library(posterior)
library(tidyr)

library(dplyr)
library(ggplot2)

draws <- as_draws_df(fit)

beta_long <- draws %>%
  select(starts_with("beta")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value")

ggplot(beta_long, aes(x = value, fill = param)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior of AFT coefficients β",
    x = "Value",
    fill = "Parameter"
  ) +
  theme_minimal()

# Restricted cubic spline basis (same as Stan)
pos <- function(x) pmax(x, 0)

rcs_h <- function(u, knots, j) {
  Kk <- length(knots)
  kj <- knots[j]
  kKm1 <- knots[Kk - 1]
  kK <- knots[Kk]
  denom <- kK - kKm1
  
  pos(u - kj)^3 -
    pos(u - kKm1)^3 * (kK - kj) / denom +
    pos(u - kK)^3   * (kKm1 - kj) / denom
}

s_rcs <- function(u, knots, gamma) {
  out <- gamma[1] + gamma[2] * u
  for (j in seq_len(length(knots) - 2)) {
    out <- out + gamma[j + 2] * rcs_h(u, knots, j)
  }
  out
}

t_grid <- seq(min(y), max(y), length.out = 100)
u_grid <- log(t_grid)

set.seed(1)
idx <- sample(seq_len(nrow(draws)), 200)

H_post <- lapply(idx, function(i) {
  gamma_i <- as.numeric(draws[i, grep("^gamma", names(draws))])
  s_vals <- s_rcs(u_grid, knots, gamma_i)
  exp(s_vals)
})

H_mat <- do.call(cbind, H_post)

df_H <- data.frame(
  t = rep(t_grid, ncol(H_mat)),
  H = as.vector(H_mat),
  draw = rep(seq_len(ncol(H_mat)), each = length(t_grid))
)

ggplot(df_H, aes(x = t, y = H, group = draw)) +
  geom_line(alpha = 0.05, color = "steelblue") +
  scale_x_log10() +
  labs(
    title = "Posterior baseline cumulative hazard H₀(t)",
    x = "Time (log scale)",
    y = "Cumulative hazard"
  ) +
  theme_minimal()

X0 <- c(1, 0)
X1 <- c(1, 1)

beta_draws <- as.matrix(draws[, grep("^beta", names(draws))])

phi0 <- exp(-beta_draws %*% X0)
phi1 <- exp(-beta_draws %*% X1)

df_phi <- data.frame(
  phi = c(phi0, phi1),
  group = rep(c("Control", "Treatment"), each = length(phi0))
)

ggplot(df_phi, aes(x = phi, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior distribution of time-acceleration factor φ(X)",
    x = "φ(X)",
    fill = ""
  ) +
  theme_minimal()

