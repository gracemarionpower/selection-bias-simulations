# ------------------------------------------------------------------------------
# Title: Sensitivity sims – selection on adult body size
#        Childhood & adulthood body size and breast cancer risk
# Author: Grace M. Power
# Date:   23 Jan 2026
# ------------------------------------------------------------------------------

rm(list = ls())

# -- libraries
library(dplyr)
library(MASS)
library(purrr)
library(tidyr)
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
  # selection parameters
  bodysize_sel_child, bodysize_sel_adult, cancer_sel, interaction_sel,
  # tracking & heritability
  phi_track = 0.35,
  h2_child = 0.10, h2_adult_target = 0.10,
  # confounding strengths U -> X1, U -> X2, U -> Y
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

  child_latent_raw <- grs[,1] * beta_child +
    gamma_u_x1 * U +
    rnorm(n, sd = sd_child_eps)

  beta_adult  <- sqrt(h2_adult_target)
  adult_noise <- grs[,2] * beta_adult +
    gamma_u_x2 * U +
    rnorm(n, sd = 1)

  adult_latent_raw <- phi_track * child_latent_raw + adult_noise

  # Standardise latent scores
  child_latent <- as.numeric(scale(child_latent_raw))
  adult_latent <- as.numeric(scale(adult_latent_raw))

  # Ordinal body size categories
  cut_child <- quantile(child_latent, probs = c(0.333, 0.841))
  cut_adult <- quantile(adult_latent, probs = c(0.333, 0.841))

  bodysize_child <- cut(child_latent,
    breaks = c(-Inf, cut_child[1], cut_child[2], Inf),
    labels = FALSE, right = TRUE) - 1

  bodysize_adult <- cut(adult_latent,
    breaks = c(-Inf, cut_adult[1], cut_adult[2], Inf),
    labels = FALSE, right = TRUE) - 1

  # Cancer depends on U only (sharp null for X -> Y)
  alpha_y <- qlogis(cancer_prev)
  p_y     <- plogis(alpha_y + gamma_u_y * U)
  cancer  <- rbinom(n, 1, p_y)

  # Selection depends on ADULT body size, cancer and interaction
  child_c <- scale(as.numeric(bodysize_child), center = TRUE, scale = FALSE)
  adult_c <- scale(as.numeric(bodysize_adult), center = TRUE, scale = FALSE)

  selection_liability <-
    child_c * bodysize_sel_child +
    adult_c * bodysize_sel_adult +
    cancer  * cancer_sel +
    (adult_c * cancer) * interaction_sel

  selection_prob <- plogis(selection_liability)
  selection      <- rbinom(n, 1, selection_prob)

  tibble(
    grs_child = grs[,1],
    grs_adult = grs[,2],
    bodysize_child,
    bodysize_adult,
    cancer,
    selection
  )
}

# ------------------------------- Estimation (2SRI) ----------------------------
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)

  fs_child <- lm(bodysize_child ~ grs_child + grs_adult, data = dat_selected)
  fs_adult <- lm(bodysize_adult ~ grs_child + grs_adult, data = dat_selected)

  dat_selected$pred_child <- fitted(fs_child)
  dat_selected$pred_adult <- fitted(fs_adult)
  dat_selected$r_child    <- resid(fs_child)
  dat_selected$r_adult    <- resid(fs_adult)

  x <- data.frame(
    intercept = rep(1, nrow(dat_selected)),
    pred_child = dat_selected$pred_child,
    pred_adult = dat_selected$pred_adult,
    r_child = dat_selected$r_child,
    r_adult = dat_selected$r_adult
  ) %>% as.matrix()

  m2 <- fastglm(
    y = dat_selected$cancer, x = x,
    family = binomial()
  )

  co <- summary(m2)$coefficients
  grab <- function(k) if (k %in% rownames(co))
    c(co[k, "Estimate"], co[k, "Std. Error"]) else c(NA, NA)

  ch <- grab("pred_child")
  ad <- grab("pred_adult")

  tibble(
    term = c("child","adult"),
    beta = c(ch[1], ad[1]),
    se   = c(ch[2], ad[2])
  )
}

# ----------------------------- Simulation wrapper -----------------------------
simulate_adult_selection <- function(
  adult_vals, cancer_vals, interaction_vals,
  gamma_sets = list(
    none   = c(0,    0,    0),
    mild   = c(0.3,  0.2,  0.5),
    strong = c(0.6,  0.4,  0.9)
  ),
  n = 246511
) {
  grid <- expand.grid(
    bodysize_sel_adult = adult_vals,
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
      bodysize_sel_child = 0,  # <-- key: no child-driven selection
      bodysize_sel_adult = grid$bodysize_sel_adult[i],
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
        bodysize_sel_adult = grid$bodysize_sel_adult[i],
        cancer_sel         = grid$cancer_sel[i],
        interaction_sel    = grid$interaction_sel[i],
        confounding_label  = grid$confounding_label[i]
      )
  })

  bind_rows(res)
}

# ---------------------------------- Run sims ----------------------------------
adult_range       <- c(0, -0.25, -0.5)
cancer_range      <- c(0, -0.25, -0.5)
interaction_range <- c(0, -0.25, -0.5)

plan(multicore, workers = 150)
set.seed(915)

sim_results_adult_sel <- furrr::future_map_dfr(1:50, function(i) {
  simulate_adult_selection(
    adult_vals = adult_range,
    cancer_vals = cancer_range,
    interaction_vals = interaction_range,
    n = 246511
  ) %>% mutate(replicate = i)
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

saveRDS(
  sim_results_adult_sel,
  file = here("sims/sim_results_adult_selection.rds")
)


# -------------------- Summarise (adult selection) --------------------
plot_df <- sim_results_adult_sel %>%
  dplyr::group_by(confounding_label,
                  bodysize_sel_adult,
                  cancer_sel,
                  interaction_sel,
                  term) %>%
  dplyr::summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # force confounding order
    confounding_label = factor(confounding_label,
                               levels = c("none","mild","strong"),
                               labels = c("No confounding",
                                          "Moderate confounding",
                                          "Strong confounding")),

    bodysize_label = factor(bodysize_sel_adult,
                            levels = c(0, -0.25, -0.5),
                            labels = c("No body size selection",
                                       "Some selection favouring thin adults",
                                       "Strong selection favouring thin adults")),

    cancer_label = factor(cancer_sel,
                          levels = c(0, -0.25, -0.5),
                          labels = c("No cancer selection",
                                     "Some cancer under-selection",
                                     "Strong cancer under-selection")),

    interaction_label = factor(interaction_sel,
                               levels = c(0, -0.25, -0.5),
                               labels = c("No interaction",
                                          "Some interaction",
                                          "Strong interaction")),

    term_label = ifelse(term == "child", "Child effect", "Adult effect"),

    env_lo = mean_logOR - 1.96 * sd_logOR,
    env_hi = mean_logOR + 1.96 * sd_logOR
  )

# -------------------- Plot --------------------
mr <- log(0.59)
pd <- position_dodge(width = 0.35)

ggplot(plot_df,
       aes(bodysize_label, mean_logOR,
           colour = cancer_label,
           linetype = term_label,
           group = interaction(confounding_label, cancer_label, term_label))) +
  geom_errorbar(aes(ymin = env_lo, ymax = env_hi),
                position = pd, width = 0.15, show.legend = FALSE) +
  geom_point(position = pd, size = 2, show.legend = FALSE) +
  geom_line(position = pd, linewidth = 0.9) +
  facet_grid(confounding_label ~ interaction_label) +
  geom_hline(yintercept = mr, linetype = "dashed", colour = "red") +
  coord_cartesian(ylim = c(-0.75, 0.2)) +
  labs(
    x = "Selection on adulthood body size",
    y = "Mean log(OR) across replicates",
    colour = "Breast cancer selection",
    linetype = "Effect"
  ) +
  scale_linetype_manual(values = c("Child effect" = "solid",
                                   "Adult effect" = "33")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.key.width = unit(2, "cm")
  )
