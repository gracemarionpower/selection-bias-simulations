# ------------------------------------------------------------------------------
# Title: Selection bias simulations
#        Childhood body size and breast cancer risk
# Authors: Grace M. Power, Gibran Hemani
# Date: 4 September 2025
# Purpose: Assess whether selection bias can reproduce the observed protective 
#          MR effect (OR = 0.59; log(OR) ≈ -0.53) of childhood adiposity on 
#          breast cancer risk, assuming no true causal effect.
# ------------------------------------------------------------------------------

rm(list = ls())

setwd("/Users/sd20930/Library/CloudStorage/OneDrive-UniversityofBristol/2. Projects at submission stage/Applied_GWAS_BMI_BC/Manuscript/Simulation paper")

# -- libraries
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)
library(tidyr)
library(stringr)

# -- data-generating model
dgm_combined <- function(
    n, rg, cancer_prev,
    bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_sel,
    phi_track = 0.35,
    h2_child = 0.10, h2_adult_target = 0.10
) {
  grs <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, rg, rg, 1), 2))
  beta_child   <- sqrt(h2_child)
  sd_child_eps <- sqrt(1 - h2_child)
  beta_adult   <- sqrt(h2_adult_target * (1 + phi_track^2))
  sd_adult_eps <- sqrt(1 - beta_adult^2)
  child_latent_raw <- grs[,1] * beta_child + rnorm(n, sd = sd_child_eps)
  adult_noise      <- grs[,2] * beta_adult + rnorm(n, sd = sd_adult_eps)
  adult_latent_raw <- phi_track * child_latent_raw + adult_noise
  child_latent <- as.numeric(scale(child_latent_raw))
  adult_latent <- as.numeric(scale(adult_latent_raw))
  cancer <- rbinom(n, 1, cancer_prev)
  cut_child <- quantile(child_latent, probs = c(0.333, 0.841))
  cut_adult <- quantile(adult_latent, probs = c(0.333, 0.841))
  bodysize_child <- cut(child_latent, breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
                        labels = FALSE, right = TRUE) - 1
  bodysize_adult <- cut(adult_latent, breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
                        labels = FALSE, right = TRUE) - 1
  child_c <- scale(as.numeric(bodysize_child), center = TRUE, scale = FALSE)
  adult_c <- scale(as.numeric(bodysize_adult), center = TRUE, scale = FALSE)
  selection_liability <-
    child_c  * bodysize_sel_child +
    adult_c  * bodysize_sel_adult +
    cancer   * cancer_sel +
    (child_c * cancer) * interaction_sel
  selection_prob <- plogis(selection_liability)
  selection <- rbinom(n, 1, selection_prob)
  tibble(
    grs_child = grs[,1],
    grs_adult = grs[,2],
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# -- 2SRI estimation (returns child & adult effects)
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)
  fs_child <- lm(bodysize_child ~ grs_child + grs_adult, data = dat_selected)
  fs_adult <- lm(bodysize_adult ~ grs_child + grs_adult, data = dat_selected)
  dat_selected$pred_child <- fitted(fs_child)
  dat_selected$pred_adult <- fitted(fs_adult)
  dat_selected$r_child    <- resid(fs_child)
  dat_selected$r_adult    <- resid(fs_adult)
  m2 <- glm(
    cancer ~ pred_child + pred_adult + r_child + r_adult,
    family = binomial(),
    data   = dat_selected
  )
  co <- summary(m2)$coefficients
  get_row <- function(name) {
    if (name %in% rownames(co)) {
      beta <- co[name, "Estimate"]; se <- co[name, "Std. Error"]
      c(beta = beta, se = se, lower = beta - 1.96*se, upper = beta + 1.96*se)
    } else c(beta = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  }
  ch <- get_row("pred_child"); ad <- get_row("pred_adult")
  tibble(
    term  = c("child", "adult"),
    beta  = c(ch["beta"],  ad["beta"]),
    se    = c(ch["se"],    ad["se"]),
    lower = c(ch["lower"], ad["lower"]),
    upper = c(ch["upper"], ad["upper"])
  )
}

# -- simulation wrapper (adult selection fixed at 0)
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

# -- run simulations
bodysize_range    <- c(0, -0.25, -0.5)
cancer_range      <- c(0, -0.25, -0.5)
interaction_range <- c(0, -0.25, -0.5)

set.seed(1407)
sim_results <- purrr::map_dfr(1:500, function(i) {
  simulate_joint_selection(bodysize_range, cancer_range, interaction_range, n = 246511) %>%
    mutate(replicate = i)
})

# -- label factors
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

# -- summarise across replicates
summary_results <- sim_results %>%
  group_by(term, bodysize_label, cancer_label, interaction_label) %>%
  summarise(
    beta  = mean(beta, na.rm = TRUE),
    lower = quantile(beta, 0.025, na.rm = TRUE),
    upper = quantile(beta, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

summary_logOR <- sim_results %>%
  group_by(term, bodysize_sel_child, cancer_sel, interaction_sel) %>%
  summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta, na.rm = TRUE),
    mean_SE    = mean(se, na.rm = TRUE),
    .groups = "drop"
  )

# -- plot (child solid; adult short-dashed; MR label once)
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
    env_hi = mean_logOR + 1.96 * sd_logOR,
    term_label = dplyr::recode(term, child = "Child effect", adult = "Adult effect")
  )

pd <- position_dodge(width = 0.35)
yr   <- range(df$env_lo, df$env_hi, mr, na.rm = TRUE)
pad  <- diff(yr) * 0.12
ymin <- min(yr[1] - pad, mr - 0.02)
ymax <- yr[2] + pad
off  <- 0.03 * (ymax - ymin)

first_facet <- levels(df$interaction_label)[1]
mr_label_df <- data.frame(
  interaction_label = factor(first_facet, levels = levels(df$interaction_label)),
  x = 1, y = mr - off,
  label = "Observed MR log(OR) \u2248 -0.53"
)

ggplot(df,
       aes(bodysize_label, mean_logOR,
           color = cancer_label,
           linetype = term_label,
           group = interaction(cancer_label, term_label))) +
  geom_errorbar(aes(ymin = env_lo, ymax = env_hi),
                position = pd, width = 0.15, show.legend = FALSE) +
  geom_point(position = pd, size = 2, show.legend = FALSE) +
  geom_line(position = pd, linewidth = 0.9, lineend = "butt", show.legend = TRUE) +
  facet_wrap(~ interaction_label) +
  geom_hline(yintercept = mr, linetype = "dashed", color = "red") +
  geom_text(data = mr_label_df, inherit.aes = FALSE,
            aes(x = x, y = y, label = label),
            color = "red", size = 3, hjust = 0, vjust = 1.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(x = "Selection on childhood body size",
       y = "Mean log(OR) across replicates",
       color = "Breast cancer selection",
       linetype = "Effect") +
  scale_linetype_manual(values = c("Child effect" = "solid",
                                   "Adult effect" = "33")) +
  guides(
    linetype = guide_legend(override.aes = list(shape = NA, linewidth = 1.3, lineend = "butt")),
    color    = guide_legend(override.aes = list(linetype = "solid", linewidth = 1.2))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.key.width = grid::unit(2, "cm"))



# ==== R² decomposition + final table  ====

# Build parameter grid from simulation summary
grid_cells <- summary_logOR %>%
  dplyr::select(bodysize_sel_child, cancer_sel, interaction_sel) %>%
  dplyr::distinct()

# Compute per-parameter-cell Shapley R² shares
set.seed(20250901)
n_for_r2 <- 246511

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
    if (length(unique(dat$selection)) < 2) {
      return(tibble::tibble(
        bodysize_sel_child = row$bodysize_sel_child,
        cancer_sel         = row$cancer_sel,
        interaction_sel    = row$interaction_sel,
        R2_total = NA_real_, share_child = NA_real_, share_cancer = NA_real_, share_inter = NA_real_
      ))
    }
    out <- .shapley_mcfadden(dat)
    tibble::tibble(
      bodysize_sel_child = row$bodysize_sel_child,
      cancer_sel         = row$cancer_sel,
      interaction_sel    = row$interaction_sel,
      R2_total    = out$R2_total,
      share_child = out$share_child,
      share_cancer= out$share_cancer,
      share_inter = out$share_inter
    )
  })

# Format R² table (percent strings for display)
fmt_pct <- function(x) ifelse(is.na(x), NA_character_, paste0(round(100 * x), "%"))

r2_selection_table <- r2_raw %>%
  dplyr::mutate(
    `β1 - child body size selection` = bodysize_sel_child,
    `β3 - breast cancer selection`   = cancer_sel,
    `β4 - interaction selection`     = interaction_sel,
    `R² of selection - child body size` = fmt_pct(share_child),
    `R² of selection - breast cancer`   = fmt_pct(share_cancer),
    `R² of selection - interaction`     = fmt_pct(share_inter)
  ) %>%
  dplyr::mutate(across(
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
  dplyr::arrange(`β1 - child body size selection`,
                 `β3 - breast cancer selection`,
                 `β4 - interaction selection`)

# Numeric version for merging
r2_prepared <- r2_selection_table %>%
  dplyr::transmute(
    bodysize_sel_child = as.numeric(`β1 - child body size selection`),
    cancer_sel         = as.numeric(`β3 - breast cancer selection`),
    interaction_sel    = as.numeric(`β4 - interaction selection`),
    R2_child       = as.numeric(stringr::str_remove(`R² of selection - child body size`, "%")),
    R2_cancer      = as.numeric(stringr::str_remove(`R² of selection - breast cancer`, "%")),
    R2_interaction = as.numeric(stringr::str_remove(`R² of selection - interaction`, "%"))
  )

# Child & adult estimates side-by-side from simulation summary
effects_wide <- df %>%
  dplyr::select(bodysize_sel_child, cancer_sel, interaction_sel, term,
                mean_logOR, sd_logOR, mean_SE) %>%
  dplyr::mutate(term = dplyr::recode(term, child = "Child", adult = "Adult")) %>%
  tidyr::pivot_wider(
    names_from  = term,
    values_from = c(mean_logOR, sd_logOR, mean_SE),
    names_glue  = "{.value} ({term})"
  )

# Final table for manuscript
final_table <- effects_wide %>%
  dplyr::left_join(r2_prepared, by = c("bodysize_sel_child","cancer_sel","interaction_sel")) %>%
  dplyr::transmute(
    `β1 - child body size selection`      = bodysize_sel_child,
    `R² of selection - child body size %` = ifelse(is.na(R2_child), NA, paste0(R2_child, "%")),
    `β3 - breast cancer selection`        = cancer_sel,
    `R² of selection - breast cancer %`   = ifelse(is.na(R2_cancer), NA, paste0(R2_cancer, "%")),
    `β4 - interaction selection`          = interaction_sel,
    `R² of selection - interaction %`     = ifelse(is.na(R2_interaction), NA, paste0(R2_interaction, "%")),
    `Mean logOR (Child)` = round(`mean_logOR (Child)`, 3),
    `SD logOR (Child)`   = round(`sd_logOR (Child)`, 3),
    `Mean SE (Child)`    = round(`mean_SE (Child)`, 3),
    `Mean logOR (Adult)` = round(`mean_logOR (Adult)`, 3),
    `SD logOR (Adult)`   = round(`sd_logOR (Adult)`, 3),
    `Mean SE (Adult)`    = round(`mean_SE (Adult)`, 3)
  ) %>%
  dplyr::arrange(`β1 - child body size selection`,
                 `β3 - breast cancer selection`,
                 `β4 - interaction selection`)

final_table

