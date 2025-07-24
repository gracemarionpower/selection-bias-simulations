
# ------------------------------------------------------------------------------
# Title: Selection bias simulations
#        Childhood body size and breast cancer risk
# Authors: Grace M. Power, Gibran Hemani
# Date: 24 July 2025
# Purpose: Assess whether selection bias can reproduce the observed protective 
#          MR effect (OR = 0.59; log(OR) ≈ -0.527) of childhood adiposity on 
#          breast cancer risk, assuming no true causal effect.
# ------------------------------------------------------------------------------

# Load required libraries
library(dplyr)
library(mvtnorm)
library(ivreg)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)

# ------------------------------------------------------------------------------
# Combined data-generating model: additive + interaction selection
# ------------------------------------------------------------------------------
dgm_combined <- function(n, rg, prs_beta_child, prs_beta_adult, cancer_prev,
                         bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_sel) {
  prs <- mvrnorm(n, mu = c(0, 0),
                 Sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))

  bodysize_child_latent <- prs[, 1] * prs_beta_child + rnorm(n, sd = sqrt(1 - prs_beta_child^2))
  bodysize_adult_latent <- prs[, 2] * prs_beta_adult + rnorm(n, sd = sqrt(1 - prs_beta_adult^2))
  cancer <- rbinom(n, 1, cancer_prev)

# UK Biobank empirical proportions:
# thinner  = 174048 / 522653 ≈ 0.333
# plumper  =  83032 / 522653 ≈ 0.159
# average  = 265573 / 522653 ≈ 0.508
# Cutoffs: P(thinner) = 0.333 → 33.3rd percentile
#          P(plumper) = 1 - 0.159 = 0.841 → 84.1st percentile

cut_child <- quantile(bodysize_child_latent, probs = c(0.333, 0.841))
cut_adult <- quantile(bodysize_adult_latent, probs = c(0.333, 0.841))

bodysize_child <- as.numeric(cut(
  bodysize_child_latent,
  breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
  labels = c(0, 1, 2),  # 0 = thinner, 1 = average, 2 = plumper
  right = TRUE
))

bodysize_adult <- as.numeric(cut(
  bodysize_adult_latent,
  breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
  labels = c(0, 1, 2),  # 0 = thinner, 1 = average, 2 = plumper
  right = TRUE
))

  selection_liability <-
    bodysize_child_latent * bodysize_sel_child +
    bodysize_adult_latent * bodysize_sel_adult +
    cancer * cancer_sel +
    (bodysize_child_latent * cancer) * interaction_sel

  selection_prob <- plogis(selection_liability)
  selection <- rbinom(n, 1, selection_prob)

  tibble(
    prs_child = prs[, 1],
    prs_adult = prs[, 2],
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# ------------------------------------------------------------------------------
# IV regression function
# ------------------------------------------------------------------------------
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)

  iv_fit <- summary(ivreg(
    cancer ~ bodysize_child + bodysize_adult |
      prs_child + prs_adult, data = dat_selected
  ))

  tibble(
    beta  = iv_fit$coefficients["bodysize_child", 1],
    se    = iv_fit$coefficients["bodysize_child", 2],
    lower = iv_fit$coefficients["bodysize_child", 1] - 1.96 * iv_fit$coefficients["bodysize_child", 2],
    upper = iv_fit$coefficients["bodysize_child", 1] + 1.96 * iv_fit$coefficients["bodysize_child", 2]
  )
}

# ------------------------------------------------------------------------------
# Simulation across grid of selection parameters
# ------------------------------------------------------------------------------
simulate_joint_selection <- function(child_vals, cancer_vals, interaction_vals, n = 1e6) {
  grid <- expand.grid(
    bodysize_sel_child = child_vals,
    cancer_sel = cancer_vals,
    interaction_sel = interaction_vals
  )

  results <- lapply(1:nrow(grid), function(i) {
    dat <- dgm_combined(
      n = n,
      rg = 0.67,
      prs_beta_child = 0.1,
      prs_beta_adult = 0.1,
      cancer_prev = 1/7,
      bodysize_sel_child = grid$bodysize_sel_child[i],
      bodysize_sel_adult = 0,
      cancer_sel = grid$cancer_sel[i],
      interaction_sel = grid$interaction_sel[i]
    )

    est <- estimation(dat)
    if (is.null(est)) return(NULL)

    est %>% mutate(
      bodysize_sel_child = grid$bodysize_sel_child[i],
      cancer_sel = grid$cancer_sel[i],
      interaction_sel = grid$interaction_sel[i]
    )
  })

  bind_rows(results)
}

# ------------------------------------------------------------------------------
# Run simulation and plot
# ------------------------------------------------------------------------------
bodysize_range <- c(0, -2, -4)
cancer_range <- c(0, -2, -4)
interaction_range <- c(0, -2, -4)

sim_results <- simulate_joint_selection(bodysize_range, cancer_range, interaction_range)

sim_results <- sim_results %>%
  mutate(interaction_label = factor(case_when(
    interaction_sel == 0  ~ "No interaction",
    interaction_sel == -2 ~ "Some interaction",
    interaction_sel == -4 ~ "Strong interaction",
    TRUE ~ paste("Interaction =", interaction_sel)
  ), levels = c(
    "No interaction",
    "Some interaction",
    "Strong interaction"
  )))

sim_results <- sim_results %>%
  mutate(
    cancer_label = factor(case_when(
      cancer_sel == 0   ~ "No cancer selection",
      cancer_sel == -2  ~ "Some cancer under-selection",
      cancer_sel == -4  ~ "Strong cancer under-selection",
      TRUE ~ paste("Cancer selection =", cancer_sel)
    ), levels = c(
      "No cancer selection",
      "Some cancer under-selection",
      "Strong cancer under-selection"
    ))
  )

sim_results <- sim_results %>%
  mutate(bodysize_label = factor(case_when(
    bodysize_sel_child == 0   ~ "No body size selection",
    bodysize_sel_child == -2  ~ "Some selection favoring thin children",
    bodysize_sel_child == -4  ~ "Strong selection favoring thin children",
    TRUE ~ paste("Body size selection =", bodysize_sel_child)
  ), levels = c(
    "No body size selection",
    "Some selection favoring thin children",
    "Strong selection favoring thin children"
  )))

# Plot results
ggplot(sim_results, aes(x = bodysize_sel_child, y = beta, color = cancer_label)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cancer_label), alpha = 0.2, color = NA) +
  facet_wrap(~ interaction_label) +
  geom_hline(yintercept = log(0.59), linetype = "dashed", color = "red") +
annotate(
  "text",
  x = 0,
  y = log(0.59),
  label = "Observed MR log(OR) ≈ -0.527",
  color = "red",
  size = 3,
  vjust = -1,
  hjust = 0
  ) +
 scale_x_reverse(
  limits = c(0.2, -4.2),  # 
  breaks = c(0, -2, -4),
  labels = c(
    "No body size selection",
    "Some selection\nfavouring thin children",
    "Strong selection\nfavouring thin children"
  )
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "Selection on childhood body size",
    y = "Estimated log(OR) per category increase in childhood body size",
    color = "Cancer selection bias",
    fill = "Cancer selection bias"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(2, "lines"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

