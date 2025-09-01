# ------------------------------------------------------------------------------
# Title: Selection bias simulations
#        Childhood body size and breast cancer risk
# Authors: Grace M. Power, Gibran Hemani
# Date: 29 August 2025
# Purpose: Assess whether selection bias can reproduce the observed protective 
#          MR effect (OR = 0.59; log(OR) ≈ -0.527) of childhood adiposity on 
#          breast cancer risk, assuming no true causal effect.
# ------------------------------------------------------------------------------

rm(list = ls())

# Load required libraries
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)

# ------------------------------------------------------------------------------
# Combined data-generating model: additive + interaction selection
# ------------------------------------------------------------------------------

# Data-generating model
#   - GRS correlation rg
#   - Adult body size "tracks" childhood (phi_track)
#   - Selection acts on centered categories (0/1/2 mean-centered), not latents
dgm_combined <- function(
    n, rg, cancer_prev,
    bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_sel,
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
    (child_c * cancer) * interaction_sel   # 0, -0.25, -0.5 (logit per +1 cat)
  
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
# Estimation: 2SRI logistic (returns log(OR) per +1 category in childhood size)
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Simulation wrapper
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Run sims
# ------------------------------------------------------------------------------
bodysize_range    <- c(0, -0.25, -0.5)
cancer_range      <- c(0, -0.25, -0.5)
interaction_range <- c(0, -0.25, -0.5)

set.seed(1407)
sim_results <- purrr::map_dfr(1:500, function(i) {
  simulate_joint_selection(
    bodysize_range, cancer_range, interaction_range, n = 246511
  ) %>% mutate(replicate = i)
})

# ------------------------------------------------------------------------------
# Labels
# ------------------------------------------------------------------------------
sim_results <- sim_results %>%
  mutate(
    interaction_label = factor(case_when(
      interaction_sel == 0     ~ "No interaction",
      interaction_sel == -0.25 ~ "Some interaction",
      interaction_sel == -0.5  ~ "Strong interaction"
    ), levels = c("No interaction","Some interaction","Strong interaction")),
    cancer_label = factor(case_when(
      cancer_sel == 0     ~ "No cancer selection",
      cancer_sel == -0.25 ~ "Some cancer under-selection",
      cancer_sel == -0.5  ~ "Strong cancer under-selection"
    ), levels = c("No cancer selection","Some cancer under-selection","Strong cancer under-selection")),
    bodysize_label = factor(case_when(
      bodysize_sel_child == 0     ~ "No body size selection",
      bodysize_sel_child == -0.25 ~ "Some selection favouring thin children",
      bodysize_sel_child == -0.5  ~ "Strong selection favouring thin children"
    ), levels = c("No body size selection","Some selection favouring thin children","Strong selection favouring thin children"))
  )

# ------------------------------------------------------------------------------
# Summarise across replicates (mean logOR and 2.5–97.5% quantiles)
# ------------------------------------------------------------------------------
summary_results <- sim_results %>%
  group_by(bodysize_label, cancer_label, interaction_label) %>%
  summarise(
    beta  = mean(beta, na.rm = TRUE),
    lower = quantile(beta, 0.025, na.rm = TRUE),
    upper = quantile(beta, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# Mean logOR, SD, and mean SE across replicates (by parameter cell)
# ------------------------------------------------------------------------------
summary_logOR <- sim_results %>%
  group_by(bodysize_sel_child, cancer_sel, interaction_sel) %>%
  summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta, na.rm = TRUE),
    mean_SE    = mean(se, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Plot: log(OR) per category; dashed MR line at log(0.59)
# ------------------------------------------------------------------------------
mr <- log(0.59)

df <- summary_logOR %>%
  mutate(
    bodysize_label = factor(bodysize_sel_child,
                            levels = c(0, -0.25, -0.50),
                            labels = c("No body size selection",
                                       "Some selection favouring thin children",
                                       "Strong selection favouring thin children")),
    cancer_label = factor(cancer_sel,
                          levels = c(0, -0.25, -0.50),
                          labels = c("No cancer selection",
                                     "Some cancer under-selection",
                                     "Strong cancer under-selection")),
    interaction_label = factor(interaction_sel,
                               levels = c(0, -0.25, -0.50),
                               labels = c("No interaction","Some interaction","Strong interaction")),
    env_lo = mean_logOR - 1.96 * sd_logOR,
    env_hi = mean_logOR + 1.96 * sd_logOR
  )

pd <- position_dodge(width = 0.35)

yr   <- range(df$env_lo, df$env_hi, mr)    # values to show
pad  <- diff(yr) * 0.12
ymin <- min(yr[1] - pad, mr - 0.02)
ymax <- yr[2] + pad
off  <- 0.03 * (ymax - ymin)               # vertical gap below the line

ggplot(df, aes(bodysize_label, mean_logOR,
               color = cancer_label, group = cancer_label)) +
  geom_errorbar(aes(ymin = env_lo, ymax = env_hi),
                position = position_dodge(0.35), width = 0.15) +
  geom_point(position = position_dodge(0.35), size = 2) +
  geom_line(position = position_dodge(0.35), linewidth = 0.7) +
  facet_wrap(~ interaction_label) +
  geom_hline(yintercept = mr, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = mr - off,                 # <— below the line
           label = "Observed MR log(OR) \u2248 -0.53",
           color = "red", size = 3, hjust = 0, vjust = 1.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(x = "Selection on childhood body size",
       y = "Mean log(OR) across replicates",
       color = "Breast cancer selection") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text = element_text(face = "bold"))




# ================================================================
# Selection R² decomposition (McFadden, Shapley-averaged)
# Uses  existing: dgm_combined(), summary_logOR
# Produces: r2_selection_table
# ================================================================

# --- libs (namespace dplyr verbs to avoid MASS::select masking) ---
suppressWarnings(suppressMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
}))

# ---------------------------
# Helper 1: simulate selection
# ---------------------------
.gen_selection_data <- function(
    n, bodysize_sel_child, cancer_sel, interaction_sel,
    rg = 0.67, cancer_prev = 1/7,
    phi_track = 0.35, h2_child = 0.10, h2_adult_target = 0.10
) {
  dat <- dgm_combined(
    n = n, rg = rg, cancer_prev = cancer_prev,
    bodysize_sel_child = bodysize_sel_child,
    bodysize_sel_adult = 0,
    cancer_sel = cancer_sel,
    interaction_sel = interaction_sel,
    phi_track = phi_track, h2_child = h2_child, h2_adult_target = h2_adult_target
  )
  # centered child body-size score used in selection
  dat$child_c <- as.numeric(scale(as.numeric(dat$bodysize_child), center = TRUE, scale = FALSE))
  dat
}

# ------------------------------------------------------------
# Helper 2: Shapley shares using McFadden R² (safe caching)
# ------------------------------------------------------------
.shapley_mcfadden <- function(dat) {
  terms <- c("child_c","cancer","child_c:cancer")
  
  # tiny cache for log-likelihoods (list; explicit key for null)
  ll_cache <- list()
  fit_ll <- function(v) {
    if (length(v) == 0) {
      k <- "<NULL>"
      if (!is.null(ll_cache[[k]])) return(ll_cache[[k]])
      fit <- glm(selection ~ 1, family = binomial(), data = dat)
      ll  <- as.numeric(logLik(fit))
      ll_cache[[k]] <- ll
      return(ll)
    } else {
      k <- paste(sort(v), collapse = "+")
      if (!is.null(ll_cache[[k]])) return(ll_cache[[k]])
      fit <- glm(reformulate(v, response = "selection"), family = binomial(), data = dat)
      ll  <- as.numeric(logLik(fit))
      ll_cache[[k]] <- ll
      return(ll)
    }
  }
  
  ll_null <- fit_ll(character(0))
  ll_full <- fit_ll(terms)
  denom   <- ll_full - ll_null
  
  # if selection has no variation or degenerate fit
  if (!is.finite(denom) || denom <= 0) {
    return(tibble(
      R2_total   = NA_real_,
      share_child= NA_real_,
      share_cancer=NA_real_,
      share_inter= NA_real_
    ))
  }
  
  # all 3! term orders
  orders <- list(
    c("child_c","cancer","child_c:cancer"),
    c("child_c","child_c:cancer","cancer"),
    c("cancer","child_c","child_c:cancer"),
    c("cancer","child_c:cancer","child_c"),
    c("child_c:cancer","child_c","cancer"),
    c("child_c:cancer","cancer","child_c")
  )
  
  contribs <- lapply(orders, function(ord) {
    base <- character(0)
    contrib <- setNames(c(0,0,0), terms)
    for (t in ord) {
      ll_base <- fit_ll(base)
      ll_add  <- fit_ll(c(base, t))
      dR2 <- (ll_add - ll_base) / denom
      if (!is.finite(dR2) || dR2 < 0) dR2 <- 0   # clamp numerical wiggles
      contrib[t] <- contrib[t] + dR2
      base <- c(base, t)
    }
    contrib
  })
  
  avg <- Reduce(`+`, contribs) / length(contribs)
  tot_R2 <- 1 - (ll_full / ll_null)  # McFadden's R² overall
  
  tibble(
    R2_total    = as.numeric(tot_R2),
    share_child = as.numeric(avg["child_c"]),
    share_cancer= as.numeric(avg["cancer"]),
    share_inter = as.numeric(avg["child_c:cancer"])
  )
}

# ----------------------------------------
# Build grid from summary_logOR object
# ----------------------------------------
if (!exists("summary_logOR")) {
  stop("summary_logOR not found. Run simulation code first so summary_logOR exists.")
}

grid_cells <- summary_logOR %>%
  dplyr::select(bodysize_sel_child, cancer_sel, interaction_sel) %>%
  dplyr::distinct()

# ----------------------------------------
# Compute per-parameter-cell Shapley R² shares
# ----------------------------------------
set.seed(20250901)
n_for_r2 <- 246511  # per your note

r2_raw <- grid_cells %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  split(.$row_id) %>%
  purrr::map_dfr(function(row) {
    row <- as.list(row)
    dat <- .gen_selection_data(
      n = n_for_r2,
      bodysize_sel_child = row$bodysize_sel_child,
      cancer_sel         = row$cancer_sel,
      interaction_sel    = row$interaction_sel
    )
    
    # if selection has no variance (rare), return NA row
    if (length(unique(dat$selection)) < 2) {
      return(tibble(
        bodysize_sel_child = row$bodysize_sel_child,
        cancer_sel         = row$cancer_sel,
        interaction_sel    = row$interaction_sel,
        R2_total   = NA_real_,
        share_child= NA_real_,
        share_cancer=NA_real_,
        share_inter= NA_real_
      ))
    }
    
    out <- .shapley_mcfadden(dat)
    tibble(
      bodysize_sel_child = row$bodysize_sel_child,
      cancer_sel         = row$cancer_sel,
      interaction_sel    = row$interaction_sel,
      R2_total    = out$R2_total,
      share_child = out$share_child,
      share_cancer= out$share_cancer,
      share_inter = out$share_inter
    )
  })

# ----------------------------------------
# Format
# ----------------------------------------
fmt_pct <- function(x) ifelse(is.na(x), NA_character_, paste0(round(100 * x), "%"))

r2_selection_table <- r2_raw %>%
  mutate(
    `β1 - child body size selection` = bodysize_sel_child,
    `β3 - breast cancer selection`   = cancer_sel,
    `β4 - interaction selection`     = interaction_sel,
    `R² of selection - child body size` = fmt_pct(share_child),
    `R² of selection - breast cancer`   = fmt_pct(share_cancer),
    `R² of selection - interaction`     = fmt_pct(share_inter)
  ) %>%
  # all-zero selection betas → R² columns should be NA (to match your example)
  mutate(across(
    c(`R² of selection - child body size`,
      `R² of selection - breast cancer`,
      `R² of selection - interaction`),
    ~ ifelse(`β1 - child body size selection` == 0 &
               `β3 - breast cancer selection`   == 0 &
               `β4 - interaction selection`     == 0, NA, .)
  )) %>%
  dplyr::select(
    `β1 - child body size selection`,
    `R² of selection - child body size`,
    `β3 - breast cancer selection`,
    `R² of selection - breast cancer`,
    `β4 - interaction selection`,
    `R² of selection - interaction`
  ) %>%
  arrange(`β1 - child body size selection`,
          `β3 - breast cancer selection`,
          `β4 - interaction selection`)

# show table
r2_selection_table
