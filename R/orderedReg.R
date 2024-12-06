#' Estimate Ordinal and Partial Proportional Odds Models with Transformed Thresholds
#'
#' @name orderedReg
#' @param formula A formula specifying the outcome and variables for the proportional odds model.
#' @param partial_formula A formula or list of formulas specifying variables for the non-proportional odds part of the model.
#' If a list is provided, each element corresponds to a specific outcome level (excluding the first level).
#' @param data A data frame containing the variables in the model.
#' @param method The estimation method to use, from \link[maxLik]{maxLik}. Options include: "NR", "CG", "BFGS", "BFGS-R", "NM", and "SANN".
#' @param family The family of the model, either "logit" or "probit". Default is "logit".
#' @param verbose Logical indicating whether to print optimization details during fitting. Default is `FALSE`.
#' @param weights Name of the variable in the data containing sampling weights (optional).
#' @param se_type The type of standard error to compute. Options include "default" or "bootstrap". Default is "default".
#' @param bootstrap_reps Number of bootstrap replications for calculating standard errors if `se_type = "bootstrap"`. Default is `1000`.
#' @return A list containing model results, including:
#' \describe{
#'   \item{`estimate`}{Parameter estimates.}
#'   \item{`logLik`}{Log-likelihood value at convergence.}
#'   \item{`gradient`}{Gradient vector at convergence.}
#'   \item{`hessian`}{Hessian matrix at convergence.}
#'   \item{`bootstrap_se`}{If bootstrap is used, contains a data frame with bootstrap standard errors and confidence intervals.}
#' }
#' Additional attributes, such as the design matrices, categories, and formula information, are included for further analysis.
#'
#' @details
#' This function estimates ordered and generalized ordered regression models
#' with transformed threshold parameters to ensure they are strictly increasing.
#'
#' For non-proportional odds, users can specify a single formula for `partial_formula` to apply the same variables to all thresholds or a list of formulas where each corresponds to a specific outcome level (excluding the first level). This flexibility allows for level-specific covariates in the model.
#'
#' @examples
#' # Example 1: Proportional odds model
#' set.seed(123)
#' data <- data.frame(
#'   y = factor(sample(1:3, 100, replace = TRUE), ordered = TRUE),
#'   x1 = rnorm(100),
#'   x2 = rbinom(100, 1, 0.5),
#'   z1 = rnorm(100),
#'   z2 = rbinom(100, 1, 0.5)
#' )
#'
#' result1 <- orderedReg(
#'   formula = y ~ x1 + x2,
#'   data = data,
#'   family = "logit"
#' )
#' summary(result1)
#'
#' # Example 2: Partial proportional odds model with a single formula
#' result2 <- orderedReg(
#'   formula = y ~ x1,
#'   partial_formula = ~ z1 + z2,
#'   data = data,
#'   family = "logit"
#' )
#' summary(result2)
#'
#' # Example 3: Partial proportional odds model with a list of formulas
#' partial_formulas <- list(
#'   ~ z1,         # Variables for the second level
#'   ~ z1 + z2     # Variables for the third level
#' )
#' result3 <- orderedReg(
#'   formula = y ~ x1,
#'   partial_formula = partial_formulas,
#'   data = data,
#'   family = "logit"
#' )
#' summary(result3)
#'
#' # Example 4: Weighted proportional odds model
#' data$weights <- runif(100, 0.5, 2)
#' result4 <- orderedReg(
#'   formula = y ~ x1 + x2,
#'   data = data,
#'   family = "logit",
#'   weights = "weights"
#' )
#' summary(result4)
#'
#' # Example 5: Bootstrap standard errors
#' result5 <- orderedReg(
#'   formula = y ~ x1 + x2,
#'   partial_formula = ~ z1 + z2,
#'   data = data,
#'   family = "logit",
#'   se_type = "bootstrap",
#'   bootstrap_reps = 500
#' )
#' result5$bootstrap_se
#'
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm quantile reformulate sd terms update
#' @importFrom Rcpp sourceCpp
#' @importFrom rsample bootstraps
#' @importFrom dplyr mutate group_by summarise %>%
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom utils head tail
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

  # Prepare design matrices for non-proportional odds
  # Prepare design matrices for non-proportional odds
  Z_list <- list()
  if (!is.null(partial_formula)) {
    if (is.list(partial_formula)) {
      if (length(partial_formula) != num_thresholds) {
        stop("If partial_formula is a list, it must have one formula for each threshold (excluding the first level).")
      }
      for (i in seq_along(partial_formula)) {
        Z <- model.matrix(partial_formula[[i]], data = data)
        Z <- Z[, -1, drop = FALSE]  # Drop intercept if present
        Z_list[[i]] <- Z
      }
    } else {
      Z <- model.matrix(partial_formula, data = data)
      Z <- Z[, -1, drop = FALSE]  # Drop intercept if present
      Z_list <- replicate(num_thresholds, Z, simplify = FALSE)
    }
  }


  # get thresholds
  get_thresholds <- function(params, num_thresholds) {
    if (length(params) < num_thresholds) {
      stop("Insufficient parameters to compute thresholds.")
    }
    inits <- params[1:num_thresholds]
    thresholds <- inits
    for (i in 2:num_thresholds) {
      thresholds[i] <- thresholds[i - 1] + abs(inits[i])
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
    # Method of moments to determine initial values
    get_initial_values <- function(y, X, Z_list, family) {
      thresholds <- numeric(length(categories) - 1)
      for (i in seq_along(thresholds)) {
        thresholds[i] <- if (family == "logit") qlogis(mean(y <= categories[i])) else qnorm(mean(y <= categories[i]))
      }
      beta <- rep(0, ncol(X))
      params <- c(thresholds, beta)

      # Check and handle Z_list
      if (!is.null(Z_list) && length(Z_list) > 0) {
        gamma <- unlist(lapply(Z_list, function(Z) rep(0, ncol(Z))))
        params <- c(params, gamma)
        expected_gamma_size <- sum(sapply(Z_list, ncol, simplify = TRUE))
        if (length(params) != (length(thresholds) + ncol(X) + expected_gamma_size)) {
          stop("Mismatch in parameter vector length.")
        }
      } else if (!is.null(Z_list) && length(Z_list) == 0) {
        # If Z_list is an empty list, do not add gamma parameters
        expected_gamma_size <- 0
      } else {
        # If Z_list is NULL, skip gamma
        expected_gamma_size <- 0
      }

      return(params)
    }


    initial_params <- get_initial_values(y, X1, Z_list, family)
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
      thresholds <- get_thresholds(params, num_thresholds)
      names(thresholds) <- x_names_intercepts[1:num_thresholds]
      transformed_params <- c(thresholds, params[(num_thresholds + 1):length(params)])

      logLikValue <- logLikFunctionCpp(
        params = transformed_params,
        categories = categories,
        X1 = X1,
        y = y,
        family = family,
        Z_list = if (length(Z_list) > 0) Z_list else list(),
        weights_ = weights
      )
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
      dplyr::mutate(results = map(splits, boot_model)) %>%
      unnest(results)

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
