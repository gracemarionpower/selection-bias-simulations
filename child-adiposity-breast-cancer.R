# ------------------------------------------------------------------------------
# Title: Selection bias simulations
#        Childhood body size and breast cancer risk
# Author: Grace M. Power
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
# Purpose: Generate data with selection on childhood/adult body size and cancer
# ------------------------------------------------------------------------------
dgm <- function(n, rg, prs_beta_child, prs_beta_adult, cancer_prev,
                bodysize_sel_child, bodysize_sel_adult, cancer_sel) {
  prs <- rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  bodysize_child <- prs[, 1] * prs_beta_child + 
    rnorm(n, sd = sqrt(1 - prs_beta_child^2))
  bodysize_adult <- prs[, 2] * prs_beta_adult + 
    rnorm(n, sd = sqrt(1 - prs_beta_adult^2))
  cancer <- rbinom(n, 1, cancer_prev)
  
  selection_liability <- bodysize_child * bodysize_sel_child +
    bodysize_adult * bodysize_sel_adult +
    cancer * cancer_sel
  selection_prob <- plogis(selection_liability)
  selection <- rbinom(n, 1, selection_prob)
  
  tibble(
    prs_child = prs[,1],
    prs_adult = prs[,2],
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# ------------------------------------------------------------------------------
# Purpose: IV regression adjusting for adulthood body size
# ------------------------------------------------------------------------------
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)
  
  iv_fit <- summary(ivreg(cancer ~ bodysize_child + bodysize_adult |
                            prs_child + prs_adult, data = dat_selected))
  
  tibble(
    beta  = iv_fit$coefficients["bodysize_child", 1],
    se    = iv_fit$coefficients["bodysize_child", 2],
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

# ------------------------------------------------------------------------------
# Purpose: Vary selection on childhood body size and cancer
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
      prs_beta_child = 0.3,
      prs_beta_adult = 0.3,
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

# Run additive selection simulations
bodysize_range <- seq(-5, 0, by = 0.2)
cancer_levels <- seq(0, -5, by = -1)
line_results <- simulate_lines(bodysize_range, cancer_levels)

# Plot: Additive Selection
ggplot(line_results, aes(x = bodysize_sel_child, y = beta, color = as.factor(cancer_sel))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(cancer_sel)), alpha = 0.2, color = NA) +
  geom_hline(yintercept = log(0.59), linetype = "dashed", color = "red") +
  labs(
    title = "Simulated bias in IV estimates of childhood body size on breast cancer under additive selection (no true effect)",
    x = "Selection on childhood body size (negative = thinner more likely selected)",
    y = "Estimated log(OR)",
    color = "Cancer selection",
    fill = "Cancer selection"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Interaction-only selection model
# ------------------------------------------------------------------------------

# Function: dgm_interaction()
dgm_interaction <- function(n, rg, prs_beta_child, prs_beta_adult, cancer_prev, interaction_sel) {
  prs <- rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  bodysize_child <- prs[,1] * prs_beta_child + rnorm(n, sd = sqrt(1 - prs_beta_child^2))
  bodysize_adult <- prs[,2] * prs_beta_adult + rnorm(n, sd = sqrt(1 - prs_beta_adult^2))
  cancer <- rbinom(n, 1, cancer_prev)
  
  selection_liability <- (bodysize_child * cancer) * interaction_sel
  selection_prob <- plogis(selection_liability)
  selection <- rbinom(n, 1, selection_prob)
  
  tibble(
    prs_child = prs[,1],
    prs_adult = prs[,2],
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# Simulate interaction-only bias
simulate_interaction_only <- function(interaction_vals, n = 1e6) {
  results <- lapply(interaction_vals, function(inter_sel) {
    dat <- dgm_interaction(
      n = n,
      rg = 0.67,
      prs_beta_child = 0.3,
      prs_beta_adult = 0.3,
      cancer_prev = 1/7,
      interaction_sel = inter_sel
    )
    
    est <- estimation(dat)
    if (is.null(est)) return(NULL)
    
    est %>% mutate(interaction_sel = inter_sel)
  })
  
  bind_rows(results)
}

interaction_range <- seq(0, -5, by = -1)
sim_interaction_only <- simulate_interaction_only(interaction_range)

# Plot: Interaction-Dependent Selection
ggplot(sim_interaction_only, aes(x = interaction_sel, y = beta)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = log(0.59), linetype = "dashed", color = "red") +
  labs(
    title = "Simulated bias in IV estimates of childhood body size on breast cancer under interaction-dependent selection (no true effect)",
    x = "Interaction term: childhood body size × breast cancer",
    y = "Estimated log(OR)"
  ) +
  theme_minimal()
