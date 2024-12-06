#' Estimate Ordinal and Partial Proportional Odds Models with Transformed Thresholds
#'
#' @name orderedReg
#' @param formula A formula specifying the outcome and variables for the proportional odds model.
#' @param partial_formula A formula or list of formulas specifying variables for the non-proportional odds part of the model.
#' If a list is provided, each element corresponds to a specific outcome level (excluding the first level).
#' @param data A data frame containing the variables in the model.
#' @param method The estimation method to use, from \link[maxLik]{maxLik}. Options include: "NR", "BFGS", "BFGS-R", "CG", "NM", and "SANN".
#' @param family The family of the model, either "logit" or "probit". Default is "logit".
#' @param verbose Logical indicating whether to print optimization details during fitting. Default is `FALSE`.
#' @param weights Name of the variable in the data containing sampling weights (optional).
#' @return A list containing model results, including parameter estimates, log-likelihood, gradient, Hessian, etc.
#' @details
#' This function estimates ordered and generalized ordered regression models with strictly increasing thresholds. For non-proportional odds, either a single partial formula or a list of partial formulas may be provided.
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   y = factor(sample(1:3, 100, replace = TRUE), ordered = TRUE),
#'   x1 = rnorm(100),
#'   x2 = rbinom(100, 1, 0.5),
#'   z1 = rnorm(100),
#'   z2 = rbinom(100, 1, 0.5)
#' )
#'
#' # Proportional odds model
#' result1 <- orderedReg(
#'   formula = y ~ x1 + x2,
#'   data = data,
#'   family = "logit"
#' )
#' summary(result1)
#'
#' # Partial proportional odds model with a single formula
#' result2 <- orderedReg(
#'   formula = y ~ x1,
#'   partial_formula = ~ z1 + z2,
#'   data = data,
#'   family = "logit"
#' )
#' summary(result2)
#'
#' # Partial proportional odds model with a list of formulas
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
#' # Weighted model and using the verbose option
#' data$weights <- runif(100, 0.5, 2)
#' result4 <- orderedReg(
#'   formula = y ~ x1 + x2,
#'   data = data,
#'   family = "logit",
#'   weights = "weights",
#'   verbose=TRUE
#' )
#' summary(result4)
#'
#'
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm quantile reformulate sd terms update
#' @importFrom Rcpp sourceCpp
#' @importFrom rsample bootstraps
#' @importFrom dplyr mutate group_by summarise %>%
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom utils head tail globalVariables
#' @useDynLib orderedReg
#' @export
orderedReg <- function(formula, partial_formula = NULL, data, method = "BFGS", family = "logit", weights = NULL, verbose=FALSE) {
  # Check family
  if (!family %in% c("logit", "probit")) stop("Family must be 'logit' or 'probit'")

  # Extract outcome and ensure it is ordered
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  if (!is.ordered(y)) {
    stop("The response variable must be an ordered factor.")
  }

  categories <- sort(unique(y))
  N_thresholds <- length(categories) - 1

  # Construct X matrix without intercept (thresholds act as intercepts)
  X <- model.matrix(update(formula, . ~ . - 1), data = data)

  # Process weights
  if (is.null(weights)) {
    w <- rep(1, length(y))
  } else {
    w <- data[[weights]]
    if (length(w) != length(y)) stop("Length of weights must match number of observations.")
  }

  # Process partial formulas
  # If partial_formula is provided, we create a Z_list of design matrices.
  Z_list <- NULL
  if (!is.null(partial_formula)) {
    if (is.list(partial_formula)) {
      # Multiple formulas, one per threshold
      if (length(partial_formula) != N_thresholds) {
        stop("If partial_formula is a list, it must have one formula per threshold.")
      }
      Z_list <- lapply(partial_formula, function(pf) {
        # Construct a formula without intercept using reformulate to avoid '.'
        pfterms <- terms(pf)
        predictors <- attr(pfterms, "term.labels")
        # reformulate ensures no '.' is used; intercept = FALSE removes intercept
        pf_no_intercept <- reformulate(predictors, response = NULL, intercept = FALSE)
        model.matrix(pf_no_intercept, data = data)
      })
    } else {
      # Single partial formula applies to all thresholds
      pfterms <- terms(partial_formula)
      predictors <- attr(pfterms, "term.labels")
      pf_no_intercept <- reformulate(predictors, response = NULL, intercept = FALSE)
      Z_mat <- model.matrix(pf_no_intercept, data = data)
      Z_list <- replicate(N_thresholds, Z_mat, simplify = FALSE)
    }
  }

  # Function to get initial thresholds based on cumulative proportions
  get_initial_thresholds <- function(y, categories, family) {
    thresholds <- numeric(length(categories) - 1)
    for (i in seq_along(thresholds)) {
      p <- mean(y <= categories[i])
      thresholds[i] <- if (family == "logit") qlogis(p) else qnorm(p)
    }
    thresholds
  }

  init_thresholds <- get_initial_thresholds(y, categories, family)
  init_betas <- rep(0, ncol(X))
  if (!is.null(Z_list)) {
    init_gamma <- unlist(lapply(Z_list, function(Z) rep(0, ncol(Z))))
    initial_params <- c(init_thresholds, init_betas, init_gamma)
  } else {
    initial_params <- c(init_thresholds, init_betas)
  }

  # Transform raw thresholds to strictly increasing
  get_strictly_increasing_thresholds <- function(params, num_thresholds) {
    raw_thresholds <- params[1:num_thresholds]
    out <- numeric(num_thresholds)
    out[1] <- raw_thresholds[1]
    if (num_thresholds > 1) {
      for (i in 2:num_thresholds) {
        out[i] <- out[i - 1] + abs(raw_thresholds[i])
      }
    }
    out
  }

  # Log-likelihood function calling the C++ code
  logLikFunction <- function(params) {
    thr <- get_strictly_increasing_thresholds(params, N_thresholds)
    transformed_params <- c(thr, params[(N_thresholds + 1):length(params)])
    val <- logLikFunctionCpp(
      transformed_params,
      categories,
      X,
      y,
      family,
      Z_list,
      w
    )
    sum(val)
  }

  # Fit model
  if (isTRUE(method == "BHHH")) method <- "BFGS"
  fit <- maxLik::maxLik(logLik = logLikFunction, start = initial_params, method = method, printLevel=ifelse(verbose,2,0))

  # Replace raw thresholds with final strictly increasing thresholds
  est_params <- fit$estimate
  final_thresholds <- get_strictly_increasing_thresholds(est_params, N_thresholds)
  fit$estimate[1:N_thresholds] <- final_thresholds

  class(fit) <- c("orderedReg", class(fit))
  attr(fit, "formula") <- formula
  attr(fit, "partial_formula") <- partial_formula
  attr(fit, "data") <- data
  attr(fit, "family") <- family
  attr(fit, "categories") <- categories
  attr(fit, "N_thresholds") <- N_thresholds
  attr(fit, "X") <- X
  attr(fit, "y") <- y
  attr(fit, "weights") <- w
  if (!is.null(Z_list)) attr(fit, "Z_list") <- Z_list

  return(fit)
}
