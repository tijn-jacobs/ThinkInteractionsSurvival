# Simulation summary and visualization for linear and non-linear treatment 
# effect scenarios

library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)
library(forcats)
library(ggh4x)

# Set base path
base_dir <- "~/GitHub/ShrinkageTrees/simulations/Main sim"

# Load data
linear <- readRDS(file.path(base_dir, "main_linear_output.rds"))
nonlinear <- readRDS(file.path(base_dir, "main_nonlinear_output.rds"))

results_linear <- linear$results_linear
results_nonlinear <- nonlinear$results_nonlinear

# Add ATE_RMSE as squared bias
add_rmse <- function(df) {
  df %>% mutate(ATE_RMSE = ATE_bias^2)
}

results_linear <- add_rmse(results_linear)
results_nonlinear <- add_rmse(results_nonlinear)

# Define metrics
metrics <- c("sigma", "ATE", "ATE_bias", "ATE_RMSE", "ATE_coverage",
             "ATE_CI_Length", "CATE_RPEHE", "CATE_coverage",
             "CATE_CI_Length", "RMSE", "C_Index")

# Summary function
summarize_df <- function(df) {
  df %>%
    group_by(Method, p_feat) %>%
    summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(ATE_RMSE = sqrt(ATE_RMSE))
}

sum_linear <- summarize_df(results_linear)
sum_nonlinear <- summarize_df(results_nonlinear)

# LaTeX table function
make_latex_table <- function(df, scenario_label = "linear") {
  method_order <- c("CHF (k=0.1)", "CHF (CV)", "IndivAFT")
  method_display <- c(
    "CHF (k=0.1)" = "Causal Horseshoe Forest",
    "CHF (CV)"    = "Causal Horseshoe Forest (CV)",
    "IndivAFT"    = "AFT-BART"
  )
  
  p_vals <- sort(unique(df$p_feat))
  
  cat("\\begin{table}[ht]\n\\centering\n")
  cat("\\begin{tabular}{cl|ccc|ccc}\n")
  cat("\\toprule\n")
  cat("$p$ & Method & \\multicolumn{3}{c|}{CATE} & ")
  cat("\\multicolumn{3}{c}{ATE} \\\\\n")
  cat("\\midrule\n\\midrule\n")
  cat(" & & RMSE & Cover. & Len. & RMSE & Cover. & Len. \\\\\n")
  cat("\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}\n")
  
  for (p in p_vals) {
    sub_df <- df[df$p_feat == p, ]
    sub_df <- sub_df[match(method_order, sub_df$Method), ]
    
    for (i in seq_len(nrow(sub_df))) {
      row <- sub_df[i, ]
      line <- sprintf(
        "%s & %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
        if (i == 1) as.character(p) else " ",
        method_display[row$Method],
        round(row$CATE_RPEHE, 3),
        round(row$CATE_coverage, 3),
        round(row$CATE_CI_Length, 3),
        round(row$ATE_RMSE, 3),
        round(row$ATE_coverage, 3),
        round(row$ATE_CI_Length, 3)
      )
      cat(line, "\n")
    }
    cat("\\midrule\n")
  }
  
  cat("\\bottomrule\n\\end{tabular}\n")
  cat(sprintf(
    "\\caption{Simulation results for the {%s} treatment effect scenario, ",
    scenario_label))
  cat("summarising RMSE, coverage, and average credible interval length ")
  cat("for CATE and ATE across methods and dimensions.}\n")
  cat(sprintf("\\label{tab:main_%s}\n", scenario_label))
  cat("\\end{table}\n\n")
}

make_latex_table(sum_linear, "linear")
make_latex_table(sum_nonlinear, "nonlinear")

# Add setting labels
df_l <- sum_linear %>% mutate(Setting = "Linear")
df_n <- sum_nonlinear %>% mutate(Setting = "Non-linear")
df <- rbind(df_l, df_n)

# Clean and reshape for plotting
df_clean <- df %>%
  rename(
    ATE_RMSE       = ATE_RMSE,
    ATE_Coverage   = ATE_coverage,
    ATE_CILength   = ATE_CI_Length,
    CATE_RMSE      = CATE_RPEHE,
    CATE_Coverage  = CATE_coverage,
    CATE_CILength  = CATE_CI_Length
  )

df_long <- df_clean %>%
  pivot_longer(
    cols = c("ATE_RMSE", "ATE_Coverage", "ATE_CILength",
             "CATE_RMSE", "CATE_Coverage", "CATE_CILength"),
    names_to = c("Estimand", "Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    Setting  = factor(Setting, levels = c("Linear", "Non-linear")),
    Estimand = factor(Estimand, levels = c("ATE", "CATE")),
    Metric   = recode(Metric,
                      "RMSE"     = "RMSE",
                      "Coverage" = "Coverage",
                      "CILength" = "Length"),
    Metric   = factor(Metric, levels = c("RMSE", "Coverage", "Length")),
    p_feat   = as.factor(p_feat),
    Method   = factor(Method, levels = c("CHF (k=0.1)", "CHF (CV)", "IndivAFT")),
    Method   = fct_recode(Method,
                          "Causal Horseshoe Forest"     = "CHF (k=0.1)",
                          "Causal Horseshoe Forest (CV)"= "CHF (CV)",
                          "AFT-BART"                    = "IndivAFT")
  )

# Plot
p <- ggplot(df_long, aes(x = p_feat, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(
    data = df_long %>% filter(Metric == "Coverage") %>% distinct(Metric),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black"
  ) +
  ggh4x::facet_nested(Metric ~ Setting + Estimand,
                      scales = "free_y", switch = "y") +
  labs(x = "Number of covariates", y = NULL, fill = "Method", title = NULL) +
  scale_fill_manual(values = c(
    "Causal Horseshoe Forest"      = "#c7e9c0",
    "Causal Horseshoe Forest (CV)" = "#74c476",
    "AFT-BART"                     = "#238b45"
  )) +
  theme_bw(base_size = 19) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(face = "bold"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

ggsave(
  filename = file.path(base_dir, "main_sim_plot.png"),
  plot = p,
  width = 1.2 * 250,
  height = 1.3 * 150,
  units = "mm",
  dpi = 320
)

