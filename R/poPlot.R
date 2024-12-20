#' Plot of non-ordered coefficient values
#'
#' @name poPlot
#' @param model A fitted model object from the `orderedReg` function.
#' @param n_boot An integer specifying the number of bootstrap samples. Default is 1000.
#' @param seed An integer specifying the random seed for reproducibility. Default is 123.
#' @param iterlim Maximum number of iterations for the optimization algorithm. Default is 1000.
#'
#' @return A plot of the non-ordered coefficient values
#'
#' @details
#'  This function generates a plot with confidence intervals (using
#' resampling) for coefficients that vary across the levels of the outcome
#' to help identify which coefficients are non-ordered. It only checks the
#' non-ordered coefficients. For coefficients that are not sigificantly
#' different across levels, the proportional odds assumption is met. If they are
#' significantly different across level, the coefficients do not meet the
#' proportional odds assumption. This provides a graphical check for the
#' proportional odds assumption.
#'
#' @examples
#' library(orderedReg)
#'
#' set.seed(10101)
#'
#' data <- data.frame(
#'   x1 = rnorm(100),
#'   x2 = rbinom(100, 1, 0.5),
#'   z1 = rnorm(100),
#'   z2 = rbinom(100, 1, 0.5)
#' )
#'
#' y.hat.1 <- -1.5 + 0.5*data$x1 + 0.5*data$x2
#' y.hat.2 <- -0.5 + 0.5*data$x1 + 0.65*data$x2
#'
#' # Compute probabilities
#' p2 <- plogis(y.hat.2, lower.tail = FALSE)
#' p1 <- plogis(y.hat.1, lower.tail = FALSE) - p2
#' p0 <- 1 - p1 - p2
#'
#' # Combine probabilities into a matrix
#' prob_matrix <- cbind(p0, p1, p2)
#'
#' # Generate the ordinal outcome for each observation
#' y_values <- apply(prob_matrix, 1, function(p) {
#'   outcome <- rmultinom(1, size = 1, prob = p)
#'   which(outcome == 1)
#' })
#'
#' # Create the ordered factor
#' data$y <- factor(y_values, levels = c(1,2,3), ordered = TRUE)
#'
#' # Partial proportional odds model with a single formula
#' result2 <- orderedReg(
#'   formula = y ~ x1,
#'   partial_formula = ~ z1 + z2,
#'   data = data,
#'   family = "logit"
#' )
#'
#' # Plot the non-ordered coefficients
#' poPlot(result2)
#'
#' @import ggplot2 broom dplyr patchwork rsample purrr
#'
#' @export
bootstrap_coefs <- function(model, n_boot = 1000, seed = 123, iterlim = 3000) {
  # Set seed for reproducibility
  set.seed(seed)

  # Extract parameters from the fitted model's attributes
  data <- attr(model, "data")
  formula <- attr(model, "formula")
  family <- attr(model, "family")
  weights <- attr(model, "weights")
  partial_formula <- attr(model, "partial_formula")
  est_method <- attr(model, "est_method")

  # Generate bootstrap samples with rsample
  resample_data <- rsample::bootstraps(data, times = n_boot)

  # Function to fit the model for each bootstrap sample
  fit_ordinal <- function(split) {
    tryCatch({
      orderedReg(
        formula = formula,
        partial_formula = partial_formula,
        data = rsample::analysis(split),
        family = family,
        method = est_method,
        weights = weights,
        iterlim = iterlim
      )
    }, error = function(e) {
      # Return NULL if model fitting fails
      NULL
    })
  }

  # Fit models to each bootstrap sample and extract coefficients
  boot_models <- resample_data %>%
    dplyr::mutate(
      model = purrr::map(splits, fit_ordinal),
      # Only extract coefficients from successful fits
      coef_info = purrr::map(model, ~if(!is.null(.x)) {
        data.frame(
          term = names(coef(.x)),
          estimate = coef(.x),
          stringsAsFactors = FALSE
        )
      } else NULL)
    )

  # Remove failed fits and unnest coefficients
  boot_coefs <- boot_models %>%
    dplyr::filter(!purrr::map_lgl(coef_info, is.null)) %>%
    tidyr::unnest(coef_info)

  # Calculate percentile intervals
  percentile_intervals <- boot_coefs %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      lower = quantile(estimate, 0.025, na.rm = TRUE),
      upper = quantile(estimate, 0.975, na.rm = TRUE)
    )

  return(percentile_intervals)
}

#' @export
poPlot <- function(model, n_boot = 1000, seed = 123, iterlim = 1000) {
  # Get bootstrap confidence intervals
  bootstrap_ci <- bootstrap_coefs(model, n_boot, seed, iterlim)
  print(cbind(tidy(model),bootstrap_ci))
  # Extract and filter model coefficients
  coefs <- tidy(model)  %>% left_join(bootstrap_ci, by = "term") %>%
    filter(!grepl("Threshold|Intercept", term)) %>%
    filter(grepl(":", term)) %>%                    # Keep only terms containing ":"
    mutate(
      group = sub(":.*", "", term)
    )

  print(coefs)

  # Create a list of plots, one for each group
  plots <- coefs %>%
    group_by(group) %>%
    group_split() %>%
    map(function(data) {
      ggplot(data, aes(x = term, y = estimate, ymin = lower, ymax = upper)) +
        geom_errorbar(width = 0.2) +
        geom_point(size = 1.5) +
        coord_flip() +
        theme_minimal() +
        labs(
          title = paste("Coefficients for Group:", unique(data$group)),
          x = "Coefficient Term",
          y = "Estimate"
        )
    })

  # Combine the plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 1)
  return(combined_plot)
}





