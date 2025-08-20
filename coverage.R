# Reference MR estimate from Richardson et al.
ref_logOR <- log(0.59)  # ≈ -0.53

# Summarise across replicates
summary_results <- sim_results %>%
  group_by(bodysize_sel_child, cancer_sel, interaction_sel) %>%
  summarise(
    mean_logOR = mean(beta, na.rm = TRUE),
    sd_logOR   = sd(beta, na.rm = TRUE),
    mean_se    = mean(se, na.rm = TRUE),
    coverage   = mean(
      (ref_logOR >= (beta - 1.96 * se)) &
      (ref_logOR <= (beta + 1.96 * se)),
      na.rm = TRUE
    ) * 100,  # as %
    n_reps     = n(),
    .groups = "drop"
  ) %>%
  arrange(bodysize_sel_child, cancer_sel, interaction_sel)

# Nice labels for table readability
summary_results <- summary_results %>%
  mutate(
    bodysize_label = case_when(
      bodysize_sel_child == 0   ~ "No body size selection",
      bodysize_sel_child == -2  ~ "Moderate body size selection",
      bodysize_sel_child == -4  ~ "Strong body size selection"
    ),
    cancer_label = case_when(
      cancer_sel == 0   ~ "No cancer selection",
      cancer_sel == -2  ~ "Moderate cancer under-selection",
      cancer_sel == -4  ~ "Strong cancer under-selection"
    ),
    interaction_label = case_when(
      interaction_sel == 0   ~ "No interaction",
      interaction_sel == -2  ~ "Moderate interaction",
      interaction_sel == -4  ~ "Strong interaction"
    )
  ) %>%
  select(
    bodysize_label, cancer_label, interaction_label,
    mean_logOR, sd_logOR, mean_se, coverage, n_reps
  )

# Preview
print(summary_results, n = 27)
