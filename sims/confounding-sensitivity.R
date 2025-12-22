# ------------------------------------------------------------------------------
# Title: Sensitivity sims with a common confounder U
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

  m2 <- glm(
    cancer ~ pred_child + pred_adult + r_child + r_adult,
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

library(furrr)
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

# ----------------------- Summaries and plotting objects -----------------------

summary_logOR_U <- sim_results_U %>%
  filter(confounding_label %in% c("none", "mild","strong")) %>%     # <<< keep only mild/strong
  group_by(term, bodysize_sel_child, cancer_sel, interaction_sel, confounding_label) %>%
  summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta, na.rm = TRUE),
    mean_SE    = mean(se,  na.rm = TRUE),
    .groups = "drop"
  )

# nice, ordered table to print/save
table_U <- summary_logOR_U %>%
  mutate(
    confounding_label = factor(confounding_label, levels = c("none", "mild","strong"),
                               labels = c("No confounding", "Mild confounding","Strong confounding")),
    term_label = dplyr::recode(term, child = "Child effect", adult = "Adult effect"),
    bodysize_label = factor(bodysize_sel_child, levels = c(0, -0.25, -0.50),
                            labels = c("No body size selection",
                                       "Some selection favouring thin children",
                                       "Strong selection favouring thin children")),
    cancer_label = factor(cancer_sel, levels = c(0, -0.25, -0.50),
                          labels = c("No cancer selection",
                                     "Some cancer under-selection",
                                     "Strong cancer under-selection")),
    interaction_label = factor(interaction_sel, levels = c(0, -0.25, -0.50),
                               labels = c("No interaction","Some interaction","Strong interaction"))
  ) %>%
  arrange(confounding_label, interaction_label, cancer_label, bodysize_label, term_label)

print(table_U)

# ----------------------- plot -----------------------

mr <- log(0.59)

dfU <- summary_logOR_U %>%                                
  mutate(
    bodysize_label = factor(bodysize_sel_child, levels = c(0, -0.25, -0.50),
                            labels = c("No body size selection",
                                       "Some selection favouring thin children",
                                       "Strong selection favouring thin children")),
    cancer_label = factor(cancer_sel, levels = c(0, -0.25, -0.50),
                          labels = c("No cancer selection",
                                     "Some cancer under-selection",
                                     "Strong cancer under-selection")),
    interaction_label = factor(interaction_sel, levels = c(0, -0.25, -0.50),
                               labels = c("No interaction","Some interaction","Strong interaction")),
    confounding_label = factor(confounding_label, levels = c("none", "mild","strong"),
                               labels = c("No confounding", "Mild confounding","Strong confounding")),
    env_lo = mean_logOR - 1.96 * sd_logOR,
    env_hi = mean_logOR + 1.96 * sd_logOR,
    term_label = dplyr::recode(term, child = "Child effect", adult = "Adult effect")
  )

pd  <- position_dodge(width = 0.35)
yr  <- range(dfU$env_lo, dfU$env_hi, mr, na.rm = TRUE)
pad <- diff(yr) * 0.12
ymin <- min(yr[1] - pad, mr - 0.02); ymax <- yr[2] + pad
off  <- 0.03 * (ymax - ymin)

# put the MR label on the first column ("No interaction") and both rows (mild/strong)
mr_label_df <- expand.grid(
  interaction_label = factor("No interaction", levels = levels(dfU$interaction_label)),
  confounding_label = factor(levels(dfU$confounding_label), levels = levels(dfU$confounding_label))
) |>
  transform(x = 1, y = mr - off, label = "Observed MR log(OR) ≈ -0.53")

pU <- ggplot(dfU,
             aes(bodysize_label, mean_logOR,
                 color = cancer_label,
                 linetype = term_label,
                 group = interaction(cancer_label, term_label))) +
  geom_errorbar(aes(ymin = env_lo, ymax = env_hi),
                position = pd, width = 0.15, show.legend = FALSE) +
  geom_point(position = pd, size = 2, show.legend = FALSE) +
  geom_line(position = pd, linewidth = 0.9, lineend = "butt", show.legend = TRUE) +
  facet_grid(confounding_label ~ interaction_label) +
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.key.width = grid::unit(2, "cm"), legend.position = "bottom", legend.direction="vertical")
ggsave(pU, file="sensitivity-plot.pdf", width = 10, height = 10)
print(pU)

