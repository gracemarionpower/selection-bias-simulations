# ------------------------------------------------------------------------------
# Variance decomposition of selection probability
# ------------------------------------------------------------------------------

variance_decomp <- function(n = 246511, beta1, beta3, beta4) {
  
  # Simulate latent predictors
  Xchild <- rnorm(n)
  Xadult <- rnorm(n)  # not used (beta2=0)
  Ycancer <- rbinom(n, 1, 1/7)
  
  # Linear predictor for selection
  lp <- beta1 * Xchild + beta3 * Ycancer + beta4 * (Xchild * Ycancer)
  
  # Logistic link
  p <- plogis(lp)
  
  # Decompose variance of selection liability
  var_total <- var(lp)
  var_child <- var(beta1 * Xchild)
  var_cancer <- var(beta3 * Ycancer)
  var_interaction <- var(beta4 * (Xchild * Ycancer))
  var_residual <- var_total - (var_child + var_cancer + var_interaction)
  
  tibble(
    beta1 = beta1,
    beta3 = beta3,
    beta4 = beta4,
    child_prop = 100 * var_child / var_total,
    cancer_prop = 100 * var_cancer / var_total,
    interaction_prop = 100 * var_interaction / var_total,
    residual_prop = 100 * var_residual / var_total
  )
}

# ------------------------------------------------------------------------------
# Run across all combinations of selection parameters
# ------------------------------------------------------------------------------
grid <- expand.grid(
  beta1 = c(0, -2, -4),
  beta3 = c(0, -2, -4),
  beta4 = c(0, -2, -4)
)

decomp_table <- purrr::map_dfr(1:nrow(grid), function(i) {
  variance_decomp(
    n = 246511,
    beta1 = grid$beta1[i],
    beta3 = grid$beta3[i],
    beta4 = grid$beta4[i]
  )
})

# Round nicely for manuscript
decomp_table <- decomp_table %>%
  mutate(across(child_prop:residual_prop, ~ round(.x, 1)))

# Print clean results
print(decomp_table, n = Inf)
