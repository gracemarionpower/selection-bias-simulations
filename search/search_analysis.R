library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)
library(furrr)
library(UpSetR)


results <- readRDS("results_upd.rds")
prob_same <- function(b1, b2, se1, se2) {
  diff <- b1 - b2
  se_diff <- sqrt(se1^2 + se2^2)
  pval <- 2 * pnorm(-abs(diff / se_diff))
  pval
}

results <- results %>%
  rowwise() %>%
  mutate(
    p_same = prob_same(beta_target, beta_est, se_target, se_est)
  )


results_sum <- group_by(results, sim_id, simrep) %>%
  summarise(
    n = n(),
    total_overlap = sum(poverlap >0.1),
    any_overlap = as.numeric(any(poverlap == 1)),
    p_same_at_least_one = 1 - prod(1 - p_same),
    p_same_all = prod(p_same),
    bodysize_sel_adult = first(bodysize_sel_adult),
    bodysize_sel_child = first(bodysize_sel_child),
    cancer_sel = first(cancer_sel),
    interaction_child_sel = first(interaction_child_sel),
    interaction_adult_sel = first(interaction_adult_sel),
    interaction_child_adult_sel = first(interaction_child_adult_sel),
    b_adult = first(b_adult),
    confounding_sel = first(confounding_sel)
)
results_sum$all_overlap <- as.numeric(results_sum$total_overlap == 4)
dim(results_sum)

summary(glm(total_overlap ~ bodysize_sel_adult + bodysize_sel_child + cancer_sel + interaction_child_sel + interaction_adult_sel + interaction_child_adult_sel + b_adult + confounding_sel, data=results_sum), family="binomial")


prop.table(table(results_sum$total_overlap))

summary(results_sum$prob_het)

table(results_sum$prob_het > 0.1)
table(results_sum$p_same_all > 0.1) %>% prop.table() %>% {. * 100}
table(results_sum$p_same_at_least_one > 0.1) %>% prop.table() %>% {. * 100}


a <- subset(results_sum, p_same > 0.1)
b <- subset(results, sim_id %in% a$sim_id) %>% ungroup() %>% dplyr::select(contains("_sel"), b_adult, beta_target, beta_est, timepoint, model, p_same, lower_target, upper_target, lower_est, upper_est, sim_id)

summary(b)

a <- subset(results_sum, p_same > 0.1 & interaction_child_adult_sel == 0)
b <- subset(results, sim_id %in% a$sim_id) %>% ungroup() %>% dplyr::select(contains("_sel"), b_adult, beta_target, beta_est, timepoint, model, p_same, lower_target, upper_target, lower_est, upper_est, sim_id)

summary(b)

## Upset plot


temp <- subset(results, poverlap > 0.1) %>%
  mutate(name = paste0(timepoint, "_", model), sim = paste0(sim_id, "_", simrep))

l <- list(
  child_multivariable = subset(temp, name == "child_multivariable")$sim %>% unique(),
  adult_multivariable = subset(temp, name == "adult_multivariable")$sim %>% unique(),
  child_univariable = subset(temp, name == "child_univariable")$sim %>% unique(),
  adult_univariable = subset(temp, name == "adult_univariable")$sim %>% unique()
)

pdf("upset_plot.pdf", width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

