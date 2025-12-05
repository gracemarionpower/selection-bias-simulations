# Load required libraries
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)

dgm_combined <- function(
    n, rg, cancer_prev,
    bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction__child_sel, interaction_adult_sel, interaction_child_adult_sel,
    phi_track = 0.35,
    h2_child = 0.10, h2_adult_target = 0.10
) {
  # Correlated GRSs
  prs <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, rg, rg, 1), 2))
  
  # Build latents
  beta_child   <- sqrt(h2_child)
  sd_child_eps <- sqrt(1 - h2_child)
  beta_adult   <- sqrt(h2_adult_target * (1 + phi_track^2))
  sd_adult_eps <- sqrt(1 - beta_adult^2)
  
  child_latent_raw <- prs[,1] * beta_child + rnorm(n, sd = sd_child_eps)
  adult_noise      <- prs[,2] * beta_adult + rnorm(n, sd = sd_adult_eps)
  adult_latent_raw <- phi_track * child_latent_raw + adult_noise
  
  # Standardise so UKB cutpoints work
  child_latent <- as.numeric(scale(child_latent_raw))
  adult_latent <- as.numeric(scale(adult_latent_raw))
  
  # Outcome (population prevalence)
  cancer <- rbinom(n, 1, cancer_prev)
  
  # UK Biobank empirical proportions:
  # thinner  = 174048 / 522653 ≈ 0.333
  # plumper  =  83032 / 522653 ≈ 0.159
  # average  = 265573 / 522653 ≈ 0.508
  # Cutoffs: P(thinner) = 0.333 → 33.3rd percentile
  #          P(plumper) = 1 - 0.159 = 0.841 → 84.1st percentile
  #
  # UK Biobank empirical proportions:
  # thinner ≈ 0.333, average ≈ 0.508, plumper ≈ 0.159
  cut_child <- quantile(child_latent, probs = c(0.333, 0.841))
  cut_adult <- quantile(adult_latent, probs = c(0.333, 0.841))
  
  # Return numeric 0/1/2 (0 = thinner, 1 = average, 2 = plumper)
  bodysize_child <- cut(child_latent, breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
                        labels = FALSE, right = TRUE) - 1
  bodysize_adult <- cut(adult_latent, breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
                        labels = FALSE, right = TRUE) - 1
  
  # Selection on mean-centered categories (0/1/2 centered)
  child_c <- scale(as.numeric(bodysize_child), center = TRUE, scale = FALSE)
  adult_c <- scale(as.numeric(bodysize_adult), center = TRUE, scale = FALSE)
  
  selection_liability <-
    child_c  * bodysize_sel_child +   # 0, -0.25, -0.5 (logit per +1 cat)
    adult_c  * bodysize_sel_adult +   # set to 0 in the grid
    cancer   * cancer_sel +           # 0, -0.25, -0.5 (logit per +1 cat)
    (child_c * cancer) * interaction_child_sel +   # 0, -0.25, -0.5 (logit per +1 cat)
    (adult_c * cancer) * interaction_adult_sel +
    (adult_c * child_c * cancer) * interaction_child_adult_sel
  
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

estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)
  
  # First stages (linear) for ordered 0/1/2 on both GRS
  fs_child <- lm(bodysize_child ~ prs_child + prs_adult, data = dat_selected)
  fs_adult <- lm(bodysize_adult ~ prs_child + prs_adult, data = dat_selected)
  
  dat_selected$pred_child <- fitted(fs_child)
  dat_selected$pred_adult <- fitted(fs_adult)
  dat_selected$r_child    <- resid(fs_child)
  dat_selected$r_adult    <- resid(fs_adult)
  
  # Second stage: logistic with control-function residuals (2SRI)
  m2 <- glm(
    cancer ~ pred_child + pred_adult + r_child + r_adult,
    family = binomial(),
    data   = dat_selected
  )
  
  co <- summary(m2)$coefficients
  beta <- co["pred_child", "Estimate"]
  se   <- co["pred_child", "Std. Error"]
  
  tibble(
    beta  = beta,              # log(OR) per +1 category
    se    = se,
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

simulate_joint_selection <- function(child_vals, cancer_vals, interaction_vals, n = 246511) {
  grid <- expand.grid(
    bodysize_sel_child = child_vals,
    cancer_sel         = cancer_vals,
    interaction_sel    = interaction_vals
  )
  
  results <- lapply(1:nrow(grid), function(i) {
    dat <- dgm_combined(
      n = n,
      rg = 0.67,
      cancer_prev = 1/7,
      bodysize_sel_child = grid$bodysize_sel_child[i],
      bodysize_sel_adult = 0,
      cancer_sel         = grid$cancer_sel[i],
      interaction_sel    = grid$interaction_sel[i],
      phi_track          = 0.35,
      h2_child           = 0.10,
      h2_adult_target    = 0.10
    )
    
    est <- estimation(dat)
    if (is.null(est)) return(NULL)
    
    est %>% mutate(
      bodysize_sel_child = grid$bodysize_sel_child[i],
      cancer_sel         = grid$cancer_sel[i],
      interaction_sel    = grid$interaction_sel[i]
    )
  })
  
  bind_rows(results)
}