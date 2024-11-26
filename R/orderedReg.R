#' Estimate Ordinal and Partial Proportional Odds Models with Transformed Thresholds
#'
#' @name orderedReg
#' @param formula A formula specifying the outcome and variables for proportional odds
#' @param partial_formula A formula specifying variables for non-proportional odds (optional)
#' @param data The data frame
#' @param method The estimation method from \link[maxLik]{maxLik}. Options include: "NR", "CG", "BFGS", "BFGS-R", "NM", and "SANN"
#' @param family The family of the model ("logit" or "probit")
#' @param verbose Logical indicating whether to print output (default: FALSE)
#' @param weights Name of variable in the data with sampling weights (Optional)
#' @param se_type Standard error type. Options include "default", "bootstrap"
#' @param bootstrap_reps Number of bootstrap replications (default: 1000)
#' @return A list with model results and fit metrics
#'
#' @details
#' This function estimates ordered and generalized ordered regression models
#' with transformed threshold parameters to ensure they are strictly increasing.
#'
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom rsample bootstraps
#' @importFrom dplyr mutate map unnest group_by summarise
#' @useDynLib orderedReg
#' @export
orderedReg <- function(
    formula,
    partial_formula = NULL,
    data,
    method = "BHHH",
    family = "logit",
    verbose = FALSE,
    weights = NULL,
    se_type = "default",
    bootstrap_reps = 1000
) {
  # Check inputs
  if (!family %in% c("logit", "probit")) stop("Family must be 'logit' or 'probit'.")
  if (!se_type %in% c("default", "bootstrap")) stop("se_type must be 'default' or 'bootstrap'.")

  # get thresholds
  get_thresholds <- function(params, num_thresholds) {
    inits <- params[1:num_thresholds]
    thresholds <- inits
    for (i in 2:num_thresholds) {
      thresholds[i] <- thresholds[i-1] + abs(inits[i])
    }
    return(thresholds)
  }


  # Default regression function
  run_regression <- function(data, formula, partial_formula, method, family, verbose, weights) {
    formula <- update(formula, . ~ . - 1)
    if (!is.null(partial_formula)) {
      partial_formula <- update(partial_formula, . ~ . - 1)
    }

    # Prepare the design matrix for the proportional odds model
    mf <- model.frame(formula, data)
    X <- as.matrix(modelr::model_matrix(data, formula))
    X1 <- X
    y <- stats::model.response(mf)
    categories <- sort(unique(y))
    num_thresholds <- length(categories) - 1

    # Variable names for parameters
    x_names_betas <- colnames(X)
    x_names_intercepts <- paste0("Threshold ", head(categories, -1), ":", tail(categories, -1))
    x_names <- c(x_names_intercepts, x_names_betas)

    # Prepare the design matrix for the non-proportional odds model if specified
    if (!is.null(partial_formula)) {
      parterms <- terms(partial_formula)
      parpreds <- attributes(parterms)$term.labels
      partial_formula <- reformulate(parpreds, response = NULL, intercept = FALSE)
      Z <- model.matrix(partial_formula, data = data)
      z_names <- colnames(Z)
      for (i in seq_along(head(categories, -1))) {
        x_names <- c(x_names, paste0(z_names, ":", categories[i + 1]))
      }
    } else {
      Z <- NULL
    }

    # Method of moments to determine initial values
    get_initial_values <- function(y, X, Z, family) {
      thresholds <- numeric(length(categories) - 1)
      for (i in seq_along(thresholds)) {
        current_category <- categories[i]
        mean_y <- mean(y <= current_category)
        thresholds[i] <- if (family == "logit") qlogis(mean_y) else qnorm(mean_y)
        thresholds_copy <- thresholds
      }
      for (i in seq_along(thresholds)) {
        if (i > 1) {
          thresholds[i] <- thresholds_copy[i] - thresholds_copy[i-1]
        }
      }


      beta <- rep(0, ncol(X))
      params <- c(thresholds, beta)

      if (!is.null(Z)) {
        Ncoef <- ncol(Z) * (length(categories) - 1)
        gamma <- rep(0, Ncoef)
        params <- c(params, gamma)
      }

      return(params)
    }

    initial_params <- get_initial_values(y, X1, Z, family)
    names(initial_params) <- x_names

    # Handle weights
    if (is.null(weights)) {
      weights <- rep(1, length(y))
    } else {
      weights <- data[[weights]]
      if (length(weights) != length(y)) stop("Length of weights must match the number of observations.")
    }

    # Define the log-likelihood function that calls the C++ function
    logLikFunction <- function(params) {
      # Transform thresholds to ensure increasing values
      thresholds <- get_thresholds(params, num_thresholds)
      names(thresholds) <- x_names_intercepts[1:num_thresholds]
      transformed_params <- c(thresholds, params[(num_thresholds + 1):length(params)])

      logLikValue <- logLikFunctionCpp(
        params = transformed_params,
        categories = categories,
        X1 = X1,
        y = y,
        family = family,
        Z_ = Z,
        partial = !is.null(partial_formula),
        weights_ = weights
      )
      # Return the sum of weighted log-likelihood contributions
      return(logLikValue)
    }

    # Use maxLik to estimate parameters
    if (verbose) printLevel <- 2 else printLevel <- 0
    result <- maxLik::maxLik(
      logLik = logLikFunction,
      start = initial_params,
      method = method,
      printLevel = printLevel
    )

    # Modify the result to show actual thresholds
    original_estimate <- result$estimate
    thresholds <- get_thresholds(original_estimate, num_thresholds)
    names(thresholds) <- x_names_intercepts[1:num_thresholds]
    result$estimate[1:num_thresholds] <- thresholds

    class(result) <- c("orderedReg", class(result))

    # Add additional attributes that will be useful for methods
    attr(result, "formula") <- formula
    attr(result, "partial_formula") <- partial_formula
    attr(result, "data") <- data
    attr(result, "family") <- family
    attr(result, "categories") <- categories
    attr(result, "N_thresholds") <- num_thresholds
    attr(result, "X") <- X1
    attr(result, "y") <- y

    # If partial formula exists, add Z matrix
    if (!is.null(partial_formula)) {
      attr(result, "Z") <- Z
    }

    return(result)

  }

  # Default estimation
  if (se_type == "default") {
    return(run_regression(
      data,
      formula,
      partial_formula,
      method,
      family,
      verbose,
      weights
    ))
  }

  # Bootstrap standard errors
  if (se_type == "bootstrap") {
    # Create bootstrap resamples
    bootstrap_samples <- rsample::bootstraps(data, times = bootstrap_reps)

    # Function to fit model on each bootstrap sample
    boot_model <- function(split) {
      boot_data <- rsample::analysis(split)
      result <- run_regression(
        boot_data,
        formula,
        partial_formula,
        method,
        family,
        verbose,
        weights
      )

      # Extract coefficients
      coefs <- result$estimate
      names(coefs) <- names(result$estimate)

      return(tibble::tibble(
        term = names(coefs),
        estimate = unname(coefs)
      ))
    }

    # Perform bootstrap
    bootstrap_results <- bootstrap_samples %>%
      dplyr::mutate(results = purrr::map(splits, boot_model)) %>%
      tidyr::unnest(results)

    # Compute bootstrapped standard errors
    bootstrap_se <- bootstrap_results %>%
      dplyr::group_by(term) %>%
      dplyr::summarise(
        bootstrap_estimate = mean(estimate),
        bootstrap_std_error = sd(estimate),
        boostrap_t_value = bootstrap_estimate / bootstrap_std_error,
        p_value = 2 * stats::pnorm(-1*abs(boostrap_t_value)),
        bootstrap_lower_ci = quantile(estimate, 0.025),
        bootstrap_upper_ci = quantile(estimate, 0.975)
      )

    # Run the original model to combine with bootstrap results
    original_model <- run_regression(
      data,
      formula,
      partial_formula,
      method,
      family,
      verbose,
      weights
    )

    # Combine original model with bootstrap results
    original_model$bootstrap_se <- bootstrap_se

    return(original_model)
  }
}
