# ------------------------------------------------------------------------------
# Variance decomposition of selection probability
# ------------------------------------------------------------------------------

# summarise the 27 scenarios from  finished sim_results
summary_by_combo <- sim_results %>%
    group_by(beta1 = bodysize_sel_child,
             beta3 = cancer_sel,
             beta4 = interaction_sel) %>%
    summarise(
        mean_logOR = mean(beta, na.rm = TRUE),
        sd_logOR   = sd(beta,   na.rm = TRUE),
        mean_se    = mean(se,   na.rm = TRUE),
        coverage_obsMR = mean(
            log(0.59) >= (beta - 1.96*se) & log(0.59) <= (beta + 1.96*se),
            na.rm = TRUE
        ) * 100,
        n_reps = n(),
        .groups = "drop"
    )

# merge R^2 contributions
summary_with_r2 <- summary_by_combo %>%
    left_join(r2_table, by = c("beta1","beta3","beta4")) %>%
    # optional: express as percentages with 1 decimal
    mutate(
        r2_child_pct       = 100 * r2_child,
        r2_cancer_pct      = 100 * r2_cancer,
        r2_interaction_pct = 100 * r2_interaction,
        r2_cross_pct       = 100 * r2_cross
    )

# final table columns in the order you described (first three = r2 columns)
final_table <- summary_with_r2 %>%
    select(beta1, beta3, beta4,
           r2_child_pct, r2_cancer_pct, r2_interaction_pct,
           mean_logOR, sd_logOR, mean_se, coverage_obsMR, n_reps)

print(final_table, n = 27)
