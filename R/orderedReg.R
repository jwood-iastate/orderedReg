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
#' @param delta Step size for numerical derivatives. Default is 1e-4.
#' @return A list containing model results, including parameter estimates, log-likelihood, gradient, Hessian, etc.
#' @details
#' This function estimates ordered and generalized ordered regression models with strictly increasing thresholds.
#' If partial formulas are provided, additional parameters for these non-proportional odds are estimated.
#' The gradient and Hessian are approximated numerically via finite differences, using the Rcpp functions.
#'
#'
#' @examples
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
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm quantile reformulate sd terms
#' @importFrom Rcpp sourceCpp
#' @importFrom rsample bootstraps
#' @importFrom dplyr mutate group_by summarise %>%
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom utils head tail globalVariables
#' @useDynLib orderedReg
#' @export
orderedReg <- function(
    formula,
    partial_formula = NULL,
    data,
    family = "logit",
    weights = NULL,
    verbose=FALSE,
    method="BFGS",
    delta = 1e-4  # step size for numeric derivatives
) {
  if (!family %in% c("probit","logit")) {
    stop("This code currently supports family='probit' or 'logit'.")
  }

  mf <- model.frame(formula, data)
  y <- model.response(mf)
  if (!is.ordered(y)) stop("Response must be an ordered factor.")

  categories <- sort(unique(y))
  J <- length(categories)
  if (J<3) stop("Need at least 3 categories.")

  X <- model.matrix(reformulate(attr(terms(formula), "term.labels"), intercept = FALSE),data = data)

  if (is.null(weights)) {
    w <- rep(1, length(y))
  } else {
    w <- data[[weights]]
    if (length(w) != length(y)) stop("Weights length must match number of observations.")
  }

  # Handle partial formula
  k <- J - 1
  Z_list <- NULL
  if (!is.null(partial_formula)) {
    if (is.list(partial_formula)) {
      if (length(partial_formula) != k) {
        stop("If partial_formula is a list, it must have one formula per threshold.")
      }
      Z_list <- lapply(partial_formula, function(pf) {
        pf_no_int <- reformulate(attr(terms(pf), "term.labels"), intercept = FALSE)
        model.matrix(pf_no_int, data=data)
      })
    } else {
      pf_no_int <- reformulate(attr(terms(partial_formula), "term.labels"), intercept = FALSE)
      Z_mat <- model.matrix(pf_no_int, data=data)
      Z_list <- replicate(k, Z_mat, simplify=FALSE)
    }
  }

  # Parameter initialization
  linpred <- as.vector(X %*% rep(0, ncol(X)))
  cuts_init <- quantile(linpred, probs=seq(1/J, (J-1)/J, length.out=J-1))
  diffs <- diff(c(0, cuts_init))
  diffs[diffs<=0] <- 0.1
  alpha_init <- log(diffs)
  beta_init <- rep(0, ncol(X))
  if (!is.null(Z_list)) {
    gamma_init <- unlist(lapply(Z_list, function(Z) rep(0, ncol(Z))))
    initial_params <- c(alpha_init, beta_init, gamma_init)
  } else {
    initial_params <- c(alpha_init, beta_init)
  }

  # Define log-likelihood function
  logLikFunction <- function(par) {
    logLik <- logLikFunctionCpp(par, categories, X, y, family, Z_list, w)
    return(logLik)
  }

  # Define gradient and Hessian functions via numeric approximation
  gradFunction <- function(par) {
    # Use numeric gradient approximation
    gradients <- gradApproxCpp(logLikFunction, par, delta = delta)
    return(gradients)
  }

  hessFunction <- function(par) {
    # Use numeric Hessian approximation
    hessians <- hessApproxCpp(logLikFunction, par, delta = delta)
    print(hessians)
    return(hessians)
  }

  # Now run maxLik with both gradient and Hessian
  fit <- maxLik::maxLik(
    logLik = logLikFunction,
    grad = gradFunction,
    # hess = hessFunction,
    start = initial_params,
    method = method,
    printLevel = ifelse(verbose, 2, 0)
  )

  alpha_names <- paste0("alpha_", 1:(J-2))
  beta_names <- colnames(X)
  param_names <- c(alpha_names, beta_names)
  if (!is.null(Z_list)) {
    z_names <- unlist(lapply(seq_along(Z_list), function(i) {
      paste0(colnames(Z_list[[i]]), ":", categories[i+1])
    }))
    param_names <- c(param_names, z_names)
  }
  names(fit$estimate) <- param_names

  class(fit) <- c("orderedReg", class(fit))
  attr(fit, "formula") <- formula
  attr(fit, "partial_formula") <- partial_formula
  attr(fit, "data") <- data
  attr(fit, "family") <- family
  attr(fit, "categories") <- categories
  attr(fit, "X") <- X
  attr(fit, "y") <- y
  attr(fit, "weights") <- w
  if (!is.null(Z_list)) attr(fit, "Z_list") <- Z_list

  return(fit)
}
