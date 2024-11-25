#' Estimate Ordinal and Partial Proportional Odds Models
#'
#' @name orderedReg
#' @param formula A formula specifying the outcome and variables for proportional odds
#' @param partial_formula A formula specifying variables for non-proportional odds (optional)
#' @param data The data frame
#' @param method The estimation method from \link[maxLik]{maxLik}. Options include: "NR", "CG", "BFGS", "BFGS-R", "NM", and "SANN"
#' @param family The family of the model ("logit" or "probit")
#' @param verbose Logical indicating whether to print output (default: FALSE)
#' @param weights Name of variable in the data with sampling weights (Optional)
#' @return A list with model results and fit metrics
#'
#' @details
#' This function estimates ordered and generalized ordered (aka, partial proportional odds) regression models. This includes logit and probit ordinal regression models.
#'
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm
#' @importFrom Rcpp sourceCpp
#' @useDynLib orderedReg
#' @export
orderedReg <- function(formula, partial_formula = NULL, data, method = "BHHH", family = "logit", verbose = FALSE, weights = NULL) {
  # Check family
  if (!family %in% c("logit", "probit")) stop("Family must be 'logit' or 'probit'.")

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
    logLikValue <- logLikFunctionCpp(
      params = params,
      categories = categories,
      X1 = X1,
      y = y,
      family = family,
      Z_ = Z,
      partial = !is.null(partial_formula),
      weights_ = weights
    )
    # Return the sum of weighted log-likelihood contributions
    return(sum(logLikValue))
  }

  # Use maxLik to estimate parameters
  if (verbose) printLevel <- 2 else printLevel <- 0
  result <- maxLik::maxLik(logLik = logLikFunction, start = initial_params, method = method, printLevel = printLevel)
  return(result)
}
