# ------------------------------------------------------------------------------
# Title: Selection bias simulations
#        Childhood body size and breast cancer risk
# Authors: Grace M. Power, Gibran Hemani
# Date: 30 June 2025
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

# ------------------------------------------------------------------------------
# Data-generating model: ordinal body size + selection
# ------------------------------------------------------------------------------
dgm <- function(n, rg, prs_beta_child, prs_beta_adult, cancer_prev,
                bodysize_sel_child, bodysize_sel_adult, cancer_sel) {
  prs <- rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  # Latent continuous traits
  bodysize_child_latent <- prs[, 1] * prs_beta_child + rnorm(n, sd = sqrt(1 - prs_beta_child^2))
  bodysize_adult_latent <- prs[, 2] * prs_beta_adult + rnorm(n, sd = sqrt(1 - prs_beta_adult^2))
  cancer <- rbinom(n, 1, cancer_prev)
  
  # Convert latent body size to ordinal categories (0 = thinner, 1 = average, 2 = plumper)
  cut_child <- quantile(bodysize_child_latent, probs = c(0.23, 0.77))
  cut_adult <- quantile(bodysize_adult_latent, probs = c(0.23, 0.77))
  
  bodysize_child <- as.numeric(cut(
    bodysize_child_latent,
    breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
    labels = c(0, 1, 2), right = TRUE
  ))
  
  bodysize_adult <- as.numeric(cut(
    bodysize_adult_latent,
    breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
    labels = c(0, 1, 2), right = TRUE
  ))
  
  # Selection model uses latent traits
  selection_liability <- bodysize_child_latent * bodysize_sel_child +
    bodysize_adult_latent * bodysize_sel_adult +
    cancer * cancer_sel
  
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
# IV regression: effect of childhood body size on cancer, adjusted for adult body size
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
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

# ------------------------------------------------------------------------------
# Vary selection on childhood body size and cancer
# ------------------------------------------------------------------------------
simulate_lines <- function(bodysize_vals, cancer_vals, n = 1e6) {
  grid <- expand.grid(
    bodysize_sel_child = bodysize_vals,
    cancer_sel = cancer_vals
  )
  
  results <- lapply(1:nrow(grid), function(i) {
    dat <- dgm(
      n = n,
      rg = 0.67,
      prs_beta_child = 0.1,
      prs_beta_adult = 0.1,
      cancer_prev = 1/7,
      bodysize_sel_child = grid$bodysize_sel_child[i],
      bodysize_sel_adult = 0,
      cancer_sel = grid$cancer_sel[i]
    )
    
    est <- estimation(dat)
    if (is.null(est)) return(NULL)
    
    est %>% mutate(
      bodysize_sel_child = grid$bodysize_sel_child[i],
      cancer_sel = grid$cancer_sel[i]
    )
  })
  
  bind_rows(results)
}

# ------------------------------------------------------------------------------
# Plot additive selection bias simulations
# ------------------------------------------------------------------------------
bodysize_range <- seq(-5, 0, by = 0.2)
cancer_levels <- seq(0, -5, by = -1)
line_results <- simulate_lines(bodysize_range, cancer_levels)

ggplot(line_results, aes(x = bodysize_sel_child, y = beta, color = as.factor(cancer_sel))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(cancer_sel)), alpha = 0.2, color = NA) +
  geom_hline(yintercept = log(0.59), linetype = "dashed", color = "red") +
  labs(
    title = "Simulated bias in IV estimates of childhood body size on breast cancer under additive selection (no true effect)",
    x = "Selection on childhood body size (negative = thinner more likely selected)",
    y = "Estimated log(OR) per category increase",
    color = "Cancer selection",
    fill = "Cancer selection"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Interaction-only selection model
# ------------------------------------------------------------------------------
dgm_interaction <- function(n, rg, prs_beta_child, prs_beta_adult, cancer_prev, interaction_sel) {
  prs <- rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  bodysize_child_latent <- prs[, 1] * prs_beta_child + rnorm(n, sd = sqrt(1 - prs_beta_child^2))
  bodysize_adult_latent <- prs[, 2] * prs_beta_adult + rnorm(n, sd = sqrt(1 - prs_beta_adult^2))
  cancer <- rbinom(n, 1, cancer_prev)
  
  cut_child <- quantile(bodysize_child_latent, probs = c(0.23, 0.77))
  cut_adult <- quantile(bodysize_adult_latent, probs = c(0.23, 0.77))
  
  bodysize_child <- as.numeric(cut(
    bodysize_child_latent,
    breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
    labels = c(0, 1, 2), right = TRUE
  ))
  
  bodysize_adult <- as.numeric(cut(
    bodysize_adult_latent,
    breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
    labels = c(0, 1, 2), right = TRUE
  ))
  
  selection_liability <- (bodysize_child_latent * cancer) * interaction_sel
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
# Simulate interaction-only selection bias
# ------------------------------------------------------------------------------
simulate_interaction_only <- function(interaction_vals, n = 1e6) {
  results <- lapply(interaction_vals, function(inter_sel) {
    dat <- dgm_interaction(
      n = n,
      rg = 0.67,
      prs_beta_child = 0.1,
      prs_beta_adult = 0.1,
      cancer_prev = 1/7,
      interaction_sel = inter_sel
    )
    
    est <- estimation(dat)
    if (is.null(est)) return(NULL)
    
    est %>% mutate(interaction_sel = inter_sel)
  })
  
  bind_rows(results)
}

# ------------------------------------------------------------------------------
# Plot interaction-only selection bias
# ------------------------------------------------------------------------------
interaction_range <- seq(0, -5, by = -1)
sim_interaction_only <- simulate_interaction_only(interaction_range)

ggplot(sim_interaction_only, aes(x = interaction_sel, y = beta)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = log(0.59), linetype = "dashed", color = "red") +
  labs(
    title = "Simulated bias in IV estimates of childhood body size on breast cancer under interaction-dependent selection (no true effect)",
    
    x = "Interaction: childhood body size × cancer",
    y = "Estimated log(OR) per category increase"
  ) +
  theme_minimal()
