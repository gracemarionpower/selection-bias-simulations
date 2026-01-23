# ------------------------------------------------------------------------------
# Title: Simulations with a common confounder U
#        Childhood body size and breast cancer risk
# Author: Grace M. Power
# Date:   30 Oct 2025
# Purpose: Add U -> {X1, X2, Y} so that conditioning on S can open
#          Z1 -> X1 <- U -> Y in addition to Z1 -> X1 -> S <- Y
# ------------------------------------------------------------------------------

rm(list = ls())

# -- libraries
library(dplyr)
library(MASS)
library(ggplot2)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)
library(grid)
library(furrr)
library(here)
library(fastglm)

# ------------------------------ Data generating model -------------------------
dgm_confounded <- function(
  n,
  rg,
  cancer_prev,
  # selection parameters (same interpretation as your main sims)
  bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_sel,
  # tracking & heritability
  phi_track = 0.35,
  h2_child = 0.10, h2_adult_target = 0.10,
  # NEW: confounding strengths U -> X1, U -> X2, U -> Y
  gamma_u_x1 = 0.30,
  gamma_u_x2 = 0.20,
  gamma_u_y  = 0.50
) {
  # Instruments (GRSs)
  grs <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, rg, rg, 1), 2))

  # Unobserved confounder U
  U <- rnorm(n)

  # Child & adult latent body size (add U)
  beta_child   <- sqrt(h2_child)
  sd_child_eps <- sqrt(1 - h2_child)

  beta_adult   <- sqrt(h2_adult_target * (1 + phi_track^2))
  sd_adult_eps <- sqrt(1 - beta_adult^2)

  child_latent_raw <- grs[,1] * beta_child + gamma_u_x1 * U + rnorm(n, sd = sd_child_eps)
  adult_noise      <- grs[,2] * beta_adult + gamma_u_x2 * U + rnorm(n, sd = sd_adult_eps)
  adult_latent_raw <- phi_track * child_latent_raw + adult_noise

  # Standardise latent scores
  child_latent <- as.numeric(scale(child_latent_raw))
  adult_latent <- as.numeric(scale(adult_latent_raw))

  # Ordinal body size categories (as in your main sims)
  cut_child <- quantile(child_latent, probs = c(0.333, 0.841))
  cut_adult <- quantile(adult_latent, probs = c(0.333, 0.841))

  bodysize_child <- cut(child_latent, breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
                        labels = FALSE, right = TRUE) - 1
  bodysize_adult <- cut(adult_latent, breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
                        labels = FALSE, right = TRUE) - 1

  # Cancer depends on U only (no causal effect of X1/X2 on Y)
  # Intercept set to hit desired prevalence approximately
  alpha_y <- qlogis(cancer_prev)
  p_y     <- plogis(alpha_y + gamma_u_y * U)
  cancer  <- rbinom(n, 1, p_y)

  # Selection S depends on X1, X2, Y and X1×Y (same structure as your paper)
  child_c <- scale(as.numeric(bodysize_child), center = TRUE, scale = FALSE)
  adult_c <- scale(as.numeric(bodysize_adult), center = TRUE, scale = FALSE)

  selection_liability <-
    child_c  * bodysize_sel_child +
    adult_c  * bodysize_sel_adult +
    cancer   * cancer_sel +
    (child_c * cancer) * interaction_sel

  selection_prob <- plogis(selection_liability)
  selection      <- rbinom(n, 1, selection_prob)

  tibble(
    grs_child = grs[,1],
    grs_adult = grs[,2],
    U         = U,
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# ------------------------------- Estimation (2SRI) ----------------------------
# identical to your estimator; we do NOT (and cannot) adjust for U
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)

  fs_child <- lm(bodysize_child ~ grs_child + grs_adult, data = dat_selected)
  fs_adult <- lm(bodysize_adult ~ grs_child + grs_adult, data = dat_selected)

  dat_selected$pred_child <- fitted(fs_child)
  dat_selected$pred_adult <- fitted(fs_adult)
  dat_selected$r_child    <- resid(fs_child)
  dat_selected$r_adult    <- resid(fs_adult)

  x <- data.frame(intercept=rep(1, nrow(dat_selected)),
             pred_child=dat_selected$pred_child,
             pred_adult=dat_selected$pred_adult,
             r_child=dat_selected$r_child,
             r_adult=dat_selected$r_adult) %>% as.matrix()

  m2 <- fastglm(
    y = dat_selected$cancer, x = x,
    family = binomial(),
    data   = dat_selected
  )

  co <- summary(m2)$coefficients
  grab <- function(k) if (k %in% rownames(co)) c(co[k, "Estimate"], co[k, "Std. Error"]) else c(NA, NA)
  ch <- grab("pred_child"); ad <- grab("pred_adult")

  tibble(
    term  = c("child","adult"),
    beta  = c(ch[1], ad[1]),
    se    = c(ch[2], ad[2])
  )
}

# ----------------------------- Simulation wrapper -----------------------------
simulate_confounded <- function(
  child_vals, cancer_vals, interaction_vals,
  gamma_sets = list( # (γ1, γ2, γy) triplets for U->X1, U->X2, U->Y
    none   = c(0,    0,    0),
    mild   = c(0.3,  0.2,  0.5),
    strong = c(0.6,  0.4,  0.9)
  ),
  n = 246511
) {
  grid <- expand.grid(
    bodysize_sel_child = child_vals,
    cancer_sel         = cancer_vals,
    interaction_sel    = interaction_vals,
    confounding_label  = names(gamma_sets),
    stringsAsFactors   = FALSE
  )

  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- gamma_sets[[ grid$confounding_label[i] ]]
    dat <- dgm_confounded(
      n = n,
      rg = 0.67,
      cancer_prev = 1/7,
      bodysize_sel_child = grid$bodysize_sel_child[i],
      bodysize_sel_adult = 0,
      cancer_sel         = grid$cancer_sel[i],
      interaction_sel    = grid$interaction_sel[i],
      phi_track          = 0.35,
      h2_child           = 0.10,
      h2_adult_target    = 0.10,
      gamma_u_x1 = g[1], gamma_u_x2 = g[2], gamma_u_y = g[3]
    )
    est <- estimation(dat)
    if (is.null(est)) return(NULL)
    est %>%
      mutate(
        bodysize_sel_child = grid$bodysize_sel_child[i],
        cancer_sel         = grid$cancer_sel[i],
        interaction_sel    = grid$interaction_sel[i],
        confounding_label  = grid$confounding_label[i]
      )
  })
  bind_rows(res)
}

# ---------------------------------- Run sims ----------------------------------
bodysize_range    <- c(0, -0.25, -0.5)
cancer_range      <- c(0, -0.25, -0.5)
interaction_range <- c(0, -0.25, -0.5)

simulate_confounded(
  child_vals = bodysize_range,
  cancer_vals = cancer_range,
  interaction_vals = interaction_range,
  n = 250000
)


plan(multicore, workers = 150)
set.seed(915)
sim_results_U <- furrr::future_map_dfr(1:300, function(i) {
  simulate_confounded(
    child_vals = bodysize_range,
    cancer_vals = cancer_range,
    interaction_vals = interaction_range,
    # adjust n or the gamma_sets above as needed
    n = 246511
  ) %>% mutate(replicate = i)
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

saveRDS(sim_results_U, file=here("sims/sim_results_U.rds"))
