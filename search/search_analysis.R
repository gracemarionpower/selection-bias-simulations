library(dplyr)
library(mvtnorm)
library(ggplot2)
library(purrr)
library(tibble)
library(MASS)
library(grid)
library(furrr)
library(UpSetR)


results <- readRDS("results.rds")
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


results_sum <- group_by(results, sim_id) %>%
  summarise(
    n = n(),
    total_overlap = sum(overlap),
    any_overlap = as.numeric(any(overlap == 1)),
    prob_het = prod(Qpval),
    mean_het = mean(Qpval),
    p_same_at_least_one = 1 - prod(1 - p_same),
    p_same_all = prod(p_same),
    interaction_child_sel = first(interaction_child_sel),
    interaction_adult_sel = first(interaction_adult_sel),
    interaction_child_adult_sel = first(interaction_child_adult_sel),
    b_adult = first(b_adult),
    confounding = first(confounding),
    confounding_sel = first(confounding_sel)
)

prop.table(table(results_sum$total_overlap))

summary(results_sum$prob_het)

table(results_sum$prob_het > 0.1)
table(results_sum$p_same_all > 0.1) %>% prop.table() %>% {. * 100}
table(results_sum$p_same_at_least_one > 0.1) %>% prop.table() %>% {. * 100}


group_by(results_sum, b_adult, confounding, confounding_sel) %>%
  summarise(
    n = n(),
    prop = mean(p_same_all > 0.1)
  ) %>%
  ungroup()





a <- subset(results_sum, p_same > 0.1)
b <- subset(results, sim_id %in% a$sim_id) %>% ungroup() %>% dplyr::select(contains("_sel"), b_adult, beta_target, beta_est, timepoint, model, p_same, lower_target, upper_target, lower_est, upper_est, sim_id)

summary(b)

a <- subset(results_sum, p_same > 0.1 & interaction_child_adult_sel == 0)
b <- subset(results, sim_id %in% a$sim_id) %>% ungroup() %>% dplyr::select(contains("_sel"), b_adult, beta_target, beta_est, timepoint, model, p_same, lower_target, upper_target, lower_est, upper_est, sim_id)

summary(b)

## Upset plot


temp <- subset(results, overlap == 1) %>%
  mutate(name = paste0(timepoint, "_", model))

l <- list(
  child_multivariable = subset(temp, name == "child_multivariable")$sim_id %>% unique(),
  adult_multivariable = subset(temp, name == "adult_multivariable")$sim_id %>% unique(),
  child_univariable = subset(temp, name == "child_univariable")$sim_id %>% unique(),
  adult_univariable = subset(temp, name == "adult_univariable")$sim_id %>% unique()
)

pdf("upset_plot.pdf", width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()


eg <- expand.grid(
  b_adult = c(-0.2, 0),
  confounding = c(0, 0.05, 0.1)
) 

params <- as.list(eg[1,])

plot_function <- function(params) {
  print(params)
  temp2 <- subset(temp, b_adult == params$b_adult & confounding == params$confounding)

  l <- list(
    child_multivariable = subset(temp2, name == "child_multivariable")$sim_id %>% unique(),
    adult_multivariable = subset(temp2, name == "adult_multivariable")$sim_id %>% unique(),
    child_univariable = subset(temp2, name == "child_univariable")$sim_id %>% unique(),
    adult_univariable = subset(temp2, name == "adult_univariable")$sim_id %>% unique()
  )
  return(l)
}

l <- plot_function(as.list(eg[1,]))
pdf(paste0("upset_badult_", eg$b_adult[1], "_confounding_", eg$confounding[1], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

l <- plot_function(as.list(eg[2,]))
pdf(paste0("upset_badult_", eg$b_adult[2], "_confounding_", eg$confounding[2], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

l <- plot_function(as.list(eg[3,]))
pdf(paste0("upset_badult_", eg$b_adult[3], "_confounding_", eg$confounding[3], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

l <- plot_function(as.list(eg[4,]))
pdf(paste0("upset_badult_", eg$b_adult[4], "_confounding_", eg$confounding[4], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

l <- plot_function(as.list(eg[5,]))
pdf(paste0("upset_badult_", eg$b_adult[5], "_confounding_", eg$confounding[5], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

l <- plot_function(as.list(eg[6,]))
pdf(paste0("upset_badult_", eg$b_adult[6], "_confounding_", eg$confounding[6], ".pdf"), width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()


eg %>% pmap(., function(...) {
  params <- list(...)
  print(params)
  temp2 <- subset(temp, b_adult == params$b_adult & confounding == params$confounding)

  l <- list(
    child_multivariable = subset(temp2, name == "child_multivariable")$sim_id %>% unique(),
    adult_multivariable = subset(temp2, name == "adult_multivariable")$sim_id %>% unique(),
    child_univariable = subset(temp2, name == "child_univariable")$sim_id %>% unique(),
    adult_univariable = subset(temp2, name == "adult_univariable")$sim_id %>% unique()
  )

  pdf(paste0("upset_badult_", params$b_adult, "_confounding_", params$confounding, ".pdf"), width=6, height=4)
  upset(fromList(l), order.by = "freq")
  dev.off()
  Sys.sleep(1)
})

l <- list(
  child_multivariable = subset(temp, name == "child_multivariable")$sim_id %>% unique(),
  adult_multivariable = subset(temp, name == "adult_multivariable")$sim_id %>% unique(),
  child_univariable = subset(temp, name == "child_univariable")$sim_id %>% unique(),
  adult_univariable = subset(temp, name == "adult_univariable")$sim_id %>% unique()
)

pdf("upset_plot.pdf", width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()



all4 <- Reduce(intersect, l)
length(all4)/nrow(results_sum) * 100

all4_res <- subset(results_sum, sim_id %in% all4)
all4_res %>% summary()

pdf("all4_params.pdf", width=8, height=6)
par(mfrow=c(2,2))
hist(all4_res$interaction_child_sel, main="Child body size selection", xlab="Log-odds per +1 category")
hist(all4_res$interaction_adult_sel, main="Adult body size selection", xlab="Log-odds per +1 category")
hist(all4_res$interaction_child_adult_sel, main="Child x Adult interaction selection", xlab="Log-odds per +1 category")
dev.off()

a <- subset(results_sum, p_same > 0.1)
b <- subset(results, sim_id %in% a$sim_id) %>% ungroup() %>% dplyr::select(contains("_sel"), b_adult, beta_target, beta_est, timepoint, model, p_same, lower_target, upper_target, lower_est, upper_est, sim_id)
b <- subset(b, !duplicated(sim_id))

summary(b)

pdf("all_params.pdf", width=8, height=8)
par(mfrow=c(3,2))
hist(b$bodysize_sel_child, main="Child body size selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
hist(b$bodysize_sel_adult, main="Adult body size selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
hist(b$cancer_sel, main="Cancer selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
hist(b$interaction_child_sel, main="Child body size x Cancer interaction selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
hist(b$interaction_adult_sel, main="Adult body size x Cancer interaction selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
hist(b$interaction_child_adult_sel, main="Child x Adult body size x Cancer interaction selection", xlab="Log-odds per +1 category", xlim=c(-0.6, 0))
dev.off()

table(results$p_same > 0.1) %>% prop.table() %>% {. * 100}

temp <- subset(results, p_same > 0.05) %>%
  mutate(name = paste0(timepoint, "_", model)) 

l <- list(
  child_multivariable = subset(temp, name == "child_multivariable")$sim_id %>% unique(),
  adult_multivariable = subset(temp, name == "adult_multivariable")$sim_id %>% unique(),
  child_univariable = subset(temp, name == "child_univariable")$sim_id %>% unique(),
  adult_univariable = subset(temp, name == "adult_univariable")$sim_id %>% unique()
)

pdf("upset_plot_prob.pdf", width=6, height=4)
upset(fromList(l), order.by = "freq")
dev.off()

all4_prob <- Reduce(intersect, l)
length(all4_prob)/nrow(results_sum) * 100


## Test


dat <- dgm_combined(
  n = 246511,
  rg = 0.67,
  cancer_prev = 1/7,
  bodysize_sel_child = -0.25,
  bodysize_sel_adult = -0.25,
  cancer_sel         = -0.25,
  interaction_child_sel = -0.5,
  interaction_adult_sel = -0.5,
  interaction_child_adult_sel = -0.5,
  phi_track          = 0.35,
  h2_child           = 0.10,
  h2_adult_target    = 0.10,
  b_adult = -0.3,
  b_child = -0.2
)
dat
summary(lm(liability ~ prs_child, dat))
summary(lm(liability ~ prs_adult, dat))
summary(lm(liability ~ prs_adult + prs_child, dat))

glm(cancer ~ bodysize_child + bodysize_adult, family=binomial(), dat) %>% summary()
glm(cancer ~ bodysize_child + bodysize_adult, family=binomial(), dat, subset=dat$selection==1) %>% summary()
estimation(dat)

dat <- dgm_combined(
  n = 246511,
  rg = 0.67,
  cancer_prev = 1/7,
  bodysize_sel_child = -0.25,
  bodysize_sel_adult = -0.25,
  cancer_sel         = -0.25,
  interaction_child_sel = -0.5,
  interaction_adult_sel = -0.5,
  interaction_child_adult_sel = -0.5,
  phi_track          = 0.35,
  h2_child           = 0.10,
  h2_adult_target    = 0.10,
  beta_adult = -0.5,
  beta_child = 0
)
estimation(dat)

dat <- dgm_combined(
  n = 500000,
  rg = 0.67,
  cancer_prev = 1/7,
  bodysize_sel_child = 0,
  bodysize_sel_adult = 0,
  cancer_sel         = 0,
  interaction_child_sel = 0,
  interaction_adult_sel = 0,
  interaction_child_adult_sel = 0,
  phi_track          = 0.35,
  h2_child           = 0.10,
  h2_adult_target    = 0.10,
  beta_adult = 0,
  beta_child = 0
)
estimation(dat)

calculate_overlap(target_vals, estimation(dat))
