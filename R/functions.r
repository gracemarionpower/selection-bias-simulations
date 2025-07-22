#' Data Generating Mechanism for Mendelian Randomization with Selection Bias
#'
#' Simulates a dataset with polygenic risk scores (PRS), exposures x1 and x2,
#' outcome y, and selection indicators. The function models selection bias where
#' participation depends on the exposures and outcome.
#'
#' @param n Integer. Sample size for the simulation.
#' @param rg Numeric. Genetic correlation between x1 and x2 PRS (0-1).
#' @param prs_beta_x1 Numeric. Effect size of x1 PRS on exposure x1.
#' @param prs_beta_x2 Numeric. Effect size of x2 PRS on exposure x2.
#' @param y_prev Numeric. Prevalence of outcome y in the population (0-1).
#' @param x1_sel Numeric. Selection coefficient for exposure x1.
#'   Positive values indicate higher x1 increases selection probability.
#' @param x2_sel Numeric. Selection coefficient for exposure x2.
#'   Positive values indicate higher x2 increases selection probability.
#' @param y_sel Numeric. Selection coefficient for outcome y.
#'   Positive values indicate outcome y increases selection probability.
#'
#' @return A tibble with the following columns:
#'   \describe{
#'     \item{prs_x1}{Polygenic risk score for exposure x1}
#'     \item{prs_x2}{Polygenic risk score for exposure x2}
#'     \item{x1}{Simulated exposure x1}
#'     \item{x2}{Simulated exposure x2}
#'     \item{y}{Binary outcome y (0/1)}
#'     \item{selection}{Binary selection indicator (0/1)}
#'   }
#'
#' @details
#' The function generates correlated polygenic risk scores from a bivariate normal
#' distribution, then simulates phenotypes with residual variance. Selection probability
#' is modeled using a logistic function of the linear combination of exposures and
#' outcome, weighted by their respective selection coefficients.
#'
#' @examples
#' \dontrun{
#' # Basic simulation with moderate selection bias
#' data <- dgm(
#'   n = 10000,
#'   rg = 0.67,
#'   prs_beta_x1 = 0.3,
#'   prs_beta_x2 = 0.3,
#'   y_prev = 1/7,
#'   x1_sel = 0.1,
#'   x2_sel = 0,
#'   y_sel = 0.2
#' )
#' }
#'
#' @seealso \code{\link{dgm_interaction}} for interaction-based selection bias
#' @export
dgm <- function(n, rg, prs_beta_x1, prs_beta_x2, y_prev,
                x1_sel, x2_sel, y_sel) {
  prs <- mvtnorm::rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  x1 <- prs[, 1] * prs_beta_x1 + 
    stats::rnorm(n, sd = sqrt(1 - prs_beta_x1^2))
  x2 <- prs[, 2] * prs_beta_x2 + 
    stats::rnorm(n, sd = sqrt(1 - prs_beta_x2^2))
  y <- stats::rbinom(n, 1, y_prev)
  
  selection_liability <- x1 * x1_sel +
    x2 * x2_sel +
    y * y_sel
  selection_prob <- stats::plogis(selection_liability)
  selection <- stats::rbinom(n, 1, selection_prob)
  
  tibble::tibble(
    prs_x1 = prs[,1],
    prs_x2 = prs[,2],
    x1,
    x2,
    y,
    selection
  )
}

#' Instrumental Variable Estimation with Selection Adjustment
#'
#' Performs instrumental variable regression to estimate the causal effect of exposure x1 
#' on outcome y, adjusting for exposure x2, using only selected participants. The function
#' uses polygenic risk scores as instruments for the exposures.
#'
#' @param dat A tibble containing the simulated data with columns:
#'   - `selection`: Binary selection indicator (0/1)
#'   - `y`: Binary outcome variable (0/1)
#'   - `x1`: Continuous exposure variable 1
#'   - `x2`: Continuous exposure variable 2 (confounder/mediator)
#'   - `prs_x1`: Polygenic risk score instrument for x1
#'   - `prs_x2`: Polygenic risk score instrument for x2
#'
#' @return A tibble with estimation results, or NULL if insufficient data:
#'   \describe{
#'     \item{beta}{Point estimate of the causal effect of x1 on y}
#'     \item{se}{Standard error of the estimate}
#'     \item{lower}{Lower bound of 95% confidence interval}
#'     \item{upper}{Upper bound of 95% confidence interval}
#'   }
#'
#' @details
#' The function first filters to selected participants only (selection == 1).
#' If fewer than 50 selected participants remain, returns NULL. Otherwise,
#' performs two-stage least squares regression using `ivreg()` with PRS as
#' instruments. The model regresses y on x1 and x2, instrumented by prs_x1
#' and prs_x2 respectively.
#'
#' @examples
#' \dontrun{
#' # Generate data and estimate causal effect
#' data <- dgm(n = 10000, rg = 0.67, prs_beta_x1 = 0.3, prs_beta_x2 = 0.3,
#'             y_prev = 1/7, x1_sel = 0.1, x2_sel = 0, y_sel = 0.2)
#' results <- estimation(data)
#' print(results)
#' }
#'
#' @seealso \code{\link{dgm}} for data generation, \code{\link{simulate_lines}} for simulation studies
#' @export
estimation <- function(dat) {
  dat_selected <- dat[dat$selection == 1, ]
  if (nrow(dat_selected) < 50) return(NULL)
  
  iv_fit <- summary(ivreg::ivreg(y ~ x1 + x2 |
                            prs_x1 + prs_x2, data = dat_selected))
  
  tibble::tibble(
    beta  = iv_fit$coefficients["x1", 1],
    se    = iv_fit$coefficients["x1", 2],
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

#' Simulate Selection Bias Across Parameter Grid
#'
#' Performs simulation study across a grid of selection coefficients for exposure x1
#' and outcome y. For each combination, generates data and estimates the causal effect
#' to study how selection bias affects inference.
#'
#' @param x1_vals Numeric vector of selection coefficients for exposure x1 to test.
#' @param y_vals Numeric vector of selection coefficients for outcome y to test.
#' @param n Integer. Sample size for each simulation (default: 1e6).
#'
#' @return A tibble with estimation results for each parameter combination:
#'   \describe{
#'     \item{beta}{Point estimate of the causal effect of x1 on y}
#'     \item{se}{Standard error of the estimate}
#'     \item{lower}{Lower bound of 95% confidence interval}
#'     \item{upper}{Upper bound of 95% confidence interval}
#'     \item{x1_sel}{Selection coefficient for x1 used in this simulation}
#'     \item{y_sel}{Selection coefficient for y used in this simulation}
#'   }
#'
#' @details
#' The function creates a grid of all combinations of `x1_vals` and `y_vals`.
#' For each combination, it generates data using `dgm()` with fixed parameters
#' (rg=0.67, prs_beta_x1=0.3, prs_beta_x2=0.3, y_prev=1/7, x2_sel=0) and
#' the varying selection coefficients. It then estimates the causal effect using
#' `estimation()`. Results from failed estimations (NULL returns) are excluded.
#'
#' @examples
#' \dontrun{
#' # Test selection bias across different coefficient values
#' x1_values <- seq(-0.2, 0.2, 0.1)
#' y_values <- seq(-0.2, 0.2, 0.1)
#' results <- simulate_lines(x1_values, y_values, n = 1e5)
#' 
#' # Visualize results
#' library(ggplot2)
#' ggplot(results, aes(x1_sel, y_sel, fill = beta)) +
#'   geom_tile() +
#'   scale_fill_gradient2()
#' }
#'
#' @seealso \code{\link{dgm}}, \code{\link{estimation}}, \code{\link{simulate_interaction_only}}
#' @export
simulate_lines <- function(x1_vals, y_vals, n = 1e6) {
  grid <- base::expand.grid(
    x1_sel = x1_vals,
    y_sel = y_vals
  )
  
  results <- base::lapply(1:nrow(grid), function(i) {
    dat <- dgm(
      n = n,
      rg = 0.67,
      prs_beta_x1 = 0.3,
      prs_beta_x2 = 0.3,
      y_prev = 1/7,
      x1_sel = grid$x1_sel[i],
      x2_sel = 0,
      y_sel = grid$y_sel[i]
    )
    
    est <- estimation(dat)
    if (base::is.null(est)) return(NULL)
    
    est %>% dplyr::mutate(
      x1_sel = grid$x1_sel[i],
      y_sel = grid$y_sel[i]
    )
  })
  
  dplyr::bind_rows(results)
}

#' Data Generating Mechanism with Interaction-Based Selection
#'
#' Simulates a dataset with polygenic risk scores (PRS), exposures x1 and x2,
#' outcome y, and selection indicators. Selection bias is based on the interaction
#' between exposure x1 and outcome y, rather than their main effects.
#'
#' @param n Integer. Sample size for the simulation.
#' @param rg Numeric. Genetic correlation between x1 and x2 PRS (0-1).
#' @param prs_beta_x1 Numeric. Effect size of x1 PRS on exposure x1.
#' @param prs_beta_x2 Numeric. Effect size of x2 PRS on exposure x2.
#' @param y_prev Numeric. Prevalence of outcome y in the population (0-1).
#' @param interaction_sel Numeric. Selection coefficient for the x1 * y interaction.
#'   Positive values indicate the interaction increases selection probability.
#'
#' @return A tibble with the following columns:
#'   \describe{
#'     \item{prs_x1}{Polygenic risk score for exposure x1}
#'     \item{prs_x2}{Polygenic risk score for exposure x2}
#'     \item{x1}{Simulated exposure x1}
#'     \item{x2}{Simulated exposure x2}
#'     \item{y}{Binary outcome y (0/1)}
#'     \item{selection}{Binary selection indicator (0/1)}
#'   }
#'
#' @details
#' Unlike the main `dgm()` function which uses main effects for selection,
#' this function models selection based solely on the interaction between x1 and y.
#' The selection liability is calculated as (x1 * y) * interaction_sel, then
#' converted to probability using the logistic function. This creates a specific
#' type of selection bias where selection depends on the joint values of exposure
#' and outcome.
#'
#' @examples
#' \dontrun{
#' # Simulate data with interaction-based selection
#' data <- dgm_interaction(
#'   n = 10000,
#'   rg = 0.67,
#'   prs_beta_x1 = 0.3,
#'   prs_beta_x2 = 0.3,
#'   y_prev = 1/7,
#'   interaction_sel = 0.5
#' )
#' }
#'
#' @seealso \code{\link{dgm}} for main effects selection, \code{\link{simulate_interaction_only}}
#' @export
dgm_interaction <- function(n, rg, prs_beta_x1, prs_beta_x2, y_prev, interaction_sel) {
  prs <- mvtnorm::rmvnorm(n, mean = c(0, 0),
                 sigma = matrix(c(1, sqrt(rg), sqrt(rg), 1), 2))
  
  x1 <- prs[,1] * prs_beta_x1 + stats::rnorm(n, sd = sqrt(1 - prs_beta_x1^2))
  x2 <- prs[,2] * prs_beta_x2 + stats::rnorm(n, sd = sqrt(1 - prs_beta_x2^2))
  y <- stats::rbinom(n, 1, y_prev)
  
  selection_liability <- (x1 * y) * interaction_sel
  selection_prob <- stats::plogis(selection_liability)
  selection <- stats::rbinom(n, 1, selection_prob)
  
  tibble::tibble(
    prs_x1 = prs[,1],
    prs_x2 = prs[,2],
    x1,
    x2,
    y,
    selection
  )
}

#' Simulate Interaction-Only Selection Bias
#'
#' Performs simulation study to examine how interaction-based selection bias
#' affects causal inference. Tests different interaction selection coefficients
#' while keeping main effect selection at zero.
#'
#' @param interaction_vals Numeric vector of interaction selection coefficients to test.
#' @param n Integer. Sample size for each simulation (default: 1e6).
#'
#' @return A tibble with estimation results for each interaction coefficient:
#'   \describe{
#'     \item{beta}{Point estimate of the causal effect of x1 on y}
#'     \item{se}{Standard error of the estimate}
#'     \item{lower}{Lower bound of 95% confidence interval}
#'     \item{upper}{Upper bound of 95% confidence interval}
#'     \item{interaction_sel}{Interaction selection coefficient used in this simulation}
#'   }
#'
#' @details
#' For each interaction coefficient in `interaction_vals`, the function generates
#' data using `dgm_interaction()` with fixed parameters (rg=0.67, prs_beta_x1=0.3,
#' prs_beta_x2=0.3, y_prev=1/7) and estimates the causal effect using `estimation()`.
#' This isolates the impact of interaction-based selection bias on causal inference,
#' separate from main effect selection bias.
#'
#' @examples
#' \dontrun{
#' # Test different interaction selection coefficients
#' interaction_values <- seq(-1, 1, 0.2)
#' results <- simulate_interaction_only(interaction_values, n = 1e5)
#' 
#' # Plot bias as function of interaction selection
#' library(ggplot2)
#' ggplot(results, aes(interaction_sel, beta)) +
#'   geom_line() +
#'   geom_hline(yintercept = 0, linetype = "dashed") +
#'   labs(x = "Interaction Selection Coefficient", y = "Estimated Effect")
#' }
#'
#' @seealso \code{\link{dgm_interaction}}, \code{\link{simulate_lines}}, \code{\link{estimation}}
#' @export
# Simulate interaction-only bias
simulate_interaction_only <- function(interaction_vals, n = 1e6) {
  results <- base::lapply(interaction_vals, function(inter_sel) {
    dat <- dgm_interaction(
      n = n,
      rg = 0.67,
      prs_beta_x1 = 0.3,
      prs_beta_x2 = 0.3,
      y_prev = 1/7,
      interaction_sel = inter_sel
    )
    
    est <- estimation(dat)
    if (base::is.null(est)) return(NULL)
    
    est %>% dplyr::mutate(interaction_sel = inter_sel)
  })
  
  dplyr::bind_rows(results)
}