
library(dplyr)
library(glue)

# 0) Reference MR estimate
ref_logOR <- log(0.59)  # ≈ -0.53

# 1) Coverage helper
coverage_pct_vs <- function(df, ref) {
  if (all(c("lower","upper") %in% names(df))) {
    mean(ref >= pmin(df$lower, df$upper) & ref <= pmax(df$lower, df$upper), na.rm = TRUE) * 100
  } else {
    lo <- df$beta - 1.96 * df$se
    hi <- df$beta + 1.96 * df$se
    mean(ref >= pmin(lo, hi) & ref <= pmax(lo, hi), na.rm = TRUE) * 100
  }
}

# 2) Summarise each (β1,β3,β4) combo across replicates
summary_by_combo <- sim_results %>%
  group_by(beta1 = bodysize_sel_child,
           beta3 = cancer_sel,
           beta4 = interaction_sel) %>%
  summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta,   na.rm = TRUE),
    mean_se    = mean(se,   na.rm = TRUE),
    coverage_obsMR = coverage_pct_vs(cur_data_all(), ref_logOR),
    n_reps     = n(),
    .groups    = "drop"
  ) %>%
  mutate(
    abs_bias_from_null = abs(mean_logOR),               # |bias| from 0
    is_strong_betas    = (beta1 == -4 & beta3 == -4 & beta4 == -4)
  )

# 3) Identify scenarios programmatically
most_biased   <- summary_by_combo %>% arrange(desc(abs_bias_from_null)) %>% slice(1)
strong_betas  <- summary_by_combo %>% filter(is_strong_betas)

top_cov <- summary_by_combo %>%
  arrange(desc(coverage_obsMR), beta1, beta3, beta4) %>%
  slice_head(n = 2)

num_zero_cov <- sum(summary_by_combo$coverage_obsMR == 0, na.rm = TRUE)

