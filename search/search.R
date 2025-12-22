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

dgm_combined <- function(
    n, rg, cancer_prev,
    bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_child_sel, interaction_adult_sel, interaction_ca_sel, interaction_child_adult_sel, 
    confounding, confounding_sel,
    phi_track = 0.35,
    h2_child = 0.10, h2_adult_target = 0.10, b_adult=0, b_child=0, sim_id = NA
) {
  # Correlated GRSs
  prs <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, rg, rg, 1), 2))
  
  # Build latents
  beta_child   <- sqrt(h2_child)
  sd_child_eps <- sqrt(1 - h2_child)
  beta_adult   <- sqrt(h2_adult_target * (1 + phi_track^2))
  sd_adult_eps <- sqrt(1 - beta_adult^2)

  conf <- rnorm(n)
  
  child_latent_raw <- prs[,1] * beta_child + conf * confounding + rnorm(n, sd = sd_child_eps)
  adult_noise      <- prs[,2] * beta_adult + conf * confounding + rnorm(n, sd = sd_adult_eps)
  adult_latent_raw <- phi_track * child_latent_raw + adult_noise
  
  # Standardise so UKB cutpoints work
  child_latent <- as.numeric(scale(child_latent_raw))
  adult_latent <- as.numeric(scale(adult_latent_raw))
  
  
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

  liability <- bodysize_child * b_child + bodysize_adult * b_adult + conf * confounding
  
  determine_mean <- function(prev, s) {
    f <- function(m) {
      mean(plogis(rnorm(10000, mean = m, sd = s))) - prev
    }
    uniroot(f, c(-20, 20))$root
  }

  m <- determine_mean(cancer_prev, sd(liability))
  liability <- m + liability

  # Outcome (population prevalence)
  cancer <- rbinom(n, 1, plogis(liability))

  selection_liability <-
    child_c  * bodysize_sel_child +   # 0, -0.25, -0.5 (logit per +1 cat)
    adult_c  * bodysize_sel_adult +   # set to 0 in the grid
    cancer   * cancer_sel +           # 0, -0.25, -0.5 (logit per +1 cat)
    conf * confounding_sel +
    (child_c * cancer) * interaction_child_sel +   # 0, -0.25, -0.5 (logit per +1 cat)
    (adult_c * cancer) * interaction_adult_sel +
    (adult_c * child_c) * interaction_ca_sel +
    (adult_c * child_c * cancer) * interaction_child_adult_sel
  
  selection_prob <- plogis(selection_liability)
  selection <- rbinom(n, 1, selection_prob)
  
  tibble(
    prs_child = prs[,1],
    prs_adult = prs[,2],
    child_latent,
    adult_latent,
    bodysize_child,
    bodysize_adult,
    liability,
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
  
  m3 <- glm( cancer ~ pred_child + r_child, family = binomial(), data = dat_selected)

  m4 <- glm( cancer ~ pred_adult + r_adult, family = binomial(), data = dat_selected)  

  co <- summary(m2)$coefficients
  beta <- co["pred_child", "Estimate"]
  se   <- co["pred_child", "Std. Error"]
  
  t1 <- tibble(
    beta  = beta,              # log(OR) per +1 category
    se    = se,
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    timepoint = "child",
    model = "multivariable"
  )

  beta <- co["pred_adult", "Estimate"]
  se   <- co["pred_adult", "Std. Error"]
  
  t2 <- tibble(
    beta  = beta,              # log(OR) per +1 category
    se    = se,
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    timepoint = "adult",
    model = "multivariable"
  )

  co <- summary(m3)$coefficients
    beta <- co["pred_child", "Estimate"]
    se   <- co["pred_child", "Std. Error"]
    t3 <- tibble(
      beta  = beta,              # log(OR) per +1 category
      se    = se,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se,
      timepoint = "child",
      model = "univariable"
    )

  co <- summary(m4)$coefficients
    beta <- co["pred_adult", "Estimate"]
    se   <- co["pred_adult", "Std. Error"]
    t4 <- tibble(
      beta  = beta,              # log(OR) per +1 category
      se    = se,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se,
      timepoint = "adult",
      model = "univariable"
    )


  return(bind_rows(
    t1, t2, t3, t4
  ))

}

estimation2 <- function(dat) {
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
  # Second stage: logistic with control-function residuals (2SRI)
  X_full <- cbind(1, dat_selected$pred_child, dat_selected$pred_adult, 
                  dat_selected$r_child, dat_selected$r_adult)
  y <- dat_selected$cancer
  m2 <- fastglm(x = X_full, y = y, family = binomial())
  
  X_full <- cbind(1, dat_selected$pred_child, dat_selected$r_child)
  m3 <- fastglm(x = X_full, y = y, family = binomial())

  X_full <- cbind(1, dat_selected$pred_adult, dat_selected$r_adult)
  m4 <- fastglm(x = X_full, y = y, family = binomial())

  co <- summary(m2)$coefficients
  beta <- co[2, "Estimate"]
  se   <- co[2, "Std. Error"]
  
  t1 <- tibble(
    beta  = beta,              # log(OR) per +1 category
    se    = se,
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    timepoint = "child",
    model = "multivariable"
  )

  beta <- co[3, "Estimate"]
  se   <- co[3, "Std. Error"]
  
  t2 <- tibble(
    beta  = beta,              # log(OR) per +1 category
    se    = se,
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    timepoint = "adult",
    model = "multivariable"
  )

  co <- summary(m3)$coefficients
    beta <- co[2, "Estimate"]
    se   <- co[2, "Std. Error"]
    t3 <- tibble(
      beta  = beta,              # log(OR) per +1 category
      se    = se,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se,
      timepoint = "child",
      model = "univariable"
    )

  co <- summary(m4)$coefficients
    beta <- co[2, "Estimate"]
    se   <- co[2, "Std. Error"]
    t4 <- tibble(
      beta  = beta,              # log(OR) per +1 category
      se    = se,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se,
      timepoint = "adult",
      model = "univariable"
    )


  return(bind_rows(
    t1, t2, t3, t4
  ))

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
      Qpval = cochran_q(c(beta_target, beta_est), c(se_target, se_est)),
      overlap = confint_overlap(upper_target, lower_target, upper_est, lower_est)
    ) %>%
    ungroup() %>%
    dplyr::select(timepoint, model, Qpval, overlap, everything())
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
    n = 500000,
    rg = 0.67,
    cancer_prev = 1/7,
    phi_track = 0.35,
    h2_child = 0.10,
    h2_adult_target = 0.10,
    confounding = c(0, 0.05, 0.1),
    confounding_sel = c(0, 0.05, 0.1),
    bodysize_sel_child = seq(0, -0.6, by = -0.2),
    bodysize_sel_adult = seq(0, -0.6, by = -0.2),
    cancer_sel         = seq(0, -0.6, by = -0.2),
    interaction_child_sel = seq(0, -0.6, by = -0.2),
    interaction_adult_sel = seq(0, -0.6, by = -0.2),
    interaction_ca_sel = seq(0, -0.6, by = -0.2),
    interaction_child_adult_sel = seq(0, -0.6, by = -0.2),
    b_adult = c(0, -0.2)
) %>% mutate(sim_id = row_number())
dim(gr)
str(gr)

dat <- dat <- do.call(dgm_combined, as.list(gr[1,]))

sim <- function(...) {
  params <- list(...)
  dat <- do.call(dgm_combined, params)
  if (is.null(dat)) return(NULL)
  est <- estimation(dat)
  if (is.null(est)) return(NULL)
  
  ov <- calculate_overlap(target_vals, est)
  
  params_df <- as_tibble(params)
  bind_cols(params_df, ov)
}

sim2 <- function(...) {
  params <- list(...)
  dat <- do.call(dgm_combined, params)
  if (is.null(dat)) return(NULL)
  est <- estimation2(dat)
  if (is.null(est)) return(NULL)
  
  ov <- calculate_overlap(target_vals, est)
  
  params_df <- as_tibble(params)
  bind_cols(params_df, ov)
}

plan(multicore, workers = 200)

results <- future_pmap(gr, sim2, .options=furrr_options(seed=TRUE), .progress=TRUE) %>%
  bind_rows()


saveRDS(results, file="results.rds")
dim(results)




# set.seed(1)
# do.call(sim, as.list(gr[1,]))$beta_est
# set.seed(1)
# do.call(sim, as.list(gr[1,]))$beta_est


# microbenchmark(
#   sim1 = do.call(sim, as.list(gr[1,])),
#   sim2 = do.call(sim2, as.list(gr[1,])),
#   times = 3
# )


# a <- rbinom(500000, 1, 0.3)
# b <- rnorm(500000)
# bmat <- cbind(1, b)

# library(microbenchmark)

# microbenchmark(
#   fastglm = summary(fastglm(x=bmat, y=a, family=binomial())),
#   glm = {mtemp <- glm(a ~ b, family=binomial()); vcov(mtemp)},
#   sglm = summary(glm(a ~ b, family=binomial())),
#   times = 5
# )


# mtemp <- glm.fit(bmat, a, family=binomial())
# mtemp <- fastglm(x=bmat, y=a, family=binomial())
# mtemp <- glm(a ~ b, family=binomial())

# summary(mtemp)

# vcov(mtemp)

# summary(fastglm(x=bmat, y=a, family=binomial()))
