# Load required libraries
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)
library(furrr)
library(fastglm)

dgm <- function(
  n,
  rg,
  cancer_prev,
  # selection parameters (same interpretation as your main sims)
  bodysize_sel_child, bodysize_sel_adult, cancer_sel, 
  interaction_child_sel, interaction_adult_sel, interaction_ca_sel, interaction_child_adult_sel, confounding_sel,
  # tracking & heritability
  phi_track = 0.35,
  b_adult = 0,
  h2_child = 0.10, h2_adult_target = 0.10,
  # NEW: confounding strengths U -> X1, U -> X2, U -> Y
  gamma_u_x1 = 0.30,
  gamma_u_x2 = 0.20,
  gamma_u_y  = 0.50,
  sim_id = 1,
  simrep = 1
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
  p_y     <- plogis(alpha_y + gamma_u_y * U + b_adult * adult_latent)
  cancer  <- rbinom(n, 1, p_y)

  # Selection S depends on X1, X2, Y and X1×Y (same structure as your paper)
  child_c <- scale(as.numeric(bodysize_child), center = TRUE, scale = FALSE)
  adult_c <- scale(as.numeric(bodysize_adult), center = TRUE, scale = FALSE)

  selection_liability <-
    child_c  * bodysize_sel_child +
    adult_c  * bodysize_sel_adult +
    cancer   * cancer_sel +
    U * confounding_sel +
    (child_c * cancer) * interaction_child_sel +   # 0, -0.25, -0.5 (logit per +1 cat)
    (adult_c * cancer) * interaction_adult_sel +
  (adult_c * child_c) * interaction_ca_sel +
    (adult_c * child_c * cancer) * interaction_child_adult_sel

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
  
  # First-stage regressions
  fs_child <- lm(bodysize_child ~ grs_child + grs_adult, data = dat_selected)
  fs_adult <- lm(bodysize_adult ~ grs_child + grs_adult, data = dat_selected)
  
  dat_selected$pred_child <- fitted(fs_child)
  dat_selected$pred_adult <- fitted(fs_adult)
  dat_selected$r_child    <- resid(fs_child)
  dat_selected$r_adult    <- resid(fs_adult)
  
  # Univariable IV models
  m_child_uni <- fastglm(
    x = model.matrix(~ pred_child + r_child, data = dat_selected),
    y = dat_selected$cancer,
    family = binomial()
  )
  
  m_adult_uni <- fastglm(
    x = model.matrix(~ pred_adult + r_adult, data = dat_selected),
    y = dat_selected$cancer,
    family = binomial()
  )
  
  # Multivariable IV model
  m_multi <- fastglm(
    x = model.matrix(~ pred_child + pred_adult + r_child + r_adult, data = dat_selected),
    y = dat_selected$cancer,
    family = binomial()
  )
  
  # Extract coefficients
  grab <- function(model, k) {
    co <- summary(model)$coefficients
    if (k %in% rownames(co)) c(co[k, "Estimate"], co[k, "Std. Error"]) else c(NA, NA)
  }
  
  ch_uni <- grab(m_child_uni, "pred_child")
  ad_uni <- grab(m_adult_uni, "pred_adult")
  ch_multi <- grab(m_multi, "pred_child")
  ad_multi <- grab(m_multi, "pred_adult")
  
  tibble(
    timepoint = rep(c("child", "adult"), each = 2),
    model     = rep(c("univariable", "multivariable"), times = 2),
    term  = c("child_univariable", "adult_univariable", "child_multivariable", "adult_multivariable"),
    beta  = c(ch_uni[1], ad_uni[1], ch_multi[1], ad_multi[1]),
    se    = c(ch_uni[2], ad_uni[2], ch_multi[2], ad_multi[2]),
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

p_overlap <- function(b1, b2, se1, se2) {
  z_diff <- abs(b1 - b2) / sqrt(se1^2 + se2^2)
  pval <- 2 * pnorm(-z_diff)
  return(pval)
}

cochran_q <- function(betas, ses) {
  w <- 1 / (ses^2)
  b_fixed <- sum(betas * w) / sum(w)
  q <- sum(w * (betas - b_fixed)^2)
  pval <- pchisq(q, df = length(betas) - 1, lower.tail = FALSE)
  return(pval)
}

confint_overlap <- function(upper1, lower1, upper2, lower2) {
  X_lower <- max(lower1, lower2)
  X_upper <- min(upper1, upper2)
  as.numeric(X_lower < X_upper)
}

calculate_overlap <- function(target_vals, est) {
  results <- target_vals %>%
    left_join(est, by = c("timepoint", "model"), suffix = c("_target", "_est")) %>%
    rowwise() %>%
    mutate(
      poverlap = p_overlap(beta_target, beta_est, se_target, se_est)
      # overlap = confint_overlap(upper_target, lower_target, upper_est, lower_est)
    ) %>%
    ungroup() %>%
    dplyr::select(timepoint, model, poverlap, everything())
  results
}

target_vals <- tibble(
    beta = c(log(0.59), log(1.08), log(0.63), log(0.82)),
    lower = log(c(0.5, 0.93, 0.55, 0.73)),
    upper = log(c(0.71, 1.27, 0.73, 0.92)),
    se = (upper - lower) / (2 * 1.96),
    timepoint = c("child", "adult", "child", "adult"),
    model = c("multivariable", "multivariable", "univariable", "univariable")
)

gr <- expand.grid(
    n = 250000,
    rg = 0.67,
    cancer_prev = 1/7,
    phi_track = 0.35,
    h2_child = 0.10,
    h2_adult_target = 0.10,
    confounding_sel = seq(0.5, -0.5, by = -0.25),
    bodysize_sel_child = seq(0.5, -0.5, by = -0.25),
    bodysize_sel_adult = seq(0.5, -0.5, by = -0.25),
    cancer_sel         = seq(0.5, -0.5, by = -0.25),
    interaction_child_sel = seq(0.5, -0.5, by = -0.25),
    interaction_adult_sel = seq(0.5, -0.5, by = -0.25),
    interaction_ca_sel = seq(0.5, -0.5, by = -0.25),
    interaction_child_adult_sel = seq(0.5, -0.5, by = -0.25),
    b_adult = c(0, -0.2)
) %>% mutate(sim_id = row_number())
gr <- lapply(1, function(x) gr %>% mutate(simrep=x)) %>% bind_rows()
dim(gr)
str(gr)

dat <- do.call(dgm, as.list(gr[1,]))
estimation(dat)

sim <- function(...) {
  params <- list(...)
  dat <- do.call(dgm, params)
  if (is.null(dat)) return(NULL)
  est <- estimation(dat)
  if (is.null(est)) return(NULL)
  
  ov <- calculate_overlap(target_vals, est)
  
  params_df <- as_tibble(params)
  bind_cols(params_df, ov)
}


plan(multicore, workers = 200)

results <- future_pmap(gr, sim, .options=furrr_options(seed=TRUE), .progress=TRUE) %>%
  bind_rows()

saveRDS(results, file="results_upd.rds")
dim(results)
