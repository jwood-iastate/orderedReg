#' Estimate Ordinal and Partial Proportional Odds Models
#'
#' @name orderedReg
#' @param formula a formula specifying the outcome and variables for proportional odds
#' @param partial_formula a formula specifying variables for non-proportional odds (optional)
#' @param data the data frame
#' @param method the estimation method , called from \link[maxLik]{maxLik}. Options include: "NR", "CG", "BFGS", "BFGS-R", "NM", and "SANN".
#' @param family the family of the model ("logit" or "probit")
#' @return A list with model results and fit metrics
#'
#' @details
#' This function estimates ordered and generalized ordered (aka, partial proportional odds) regression models. This includes logit and probit ordinal regression models.
#'
#'
#' @import modelr maxLik stats
#' @export
orderedReg <- function(formula, partial_formula = NULL, data, method = "BHHH", family = "logit") {
  # Check family
  if (!family %in% c("logit", "probit")) stop("Family must be 'logit' or 'probit'")

  # Prepare the design matrix for the proportional odds model
  X <- as.matrix(modelr::model_matrix(data, formula))
  # If the intercept is present, it is typically the first column
  if ("(Intercept)" %in% colnames(X)) {
    # Remove the intercept (first column)
    X1 <- X[, -1]
  }
  else{
    X1 <- X
  }
  y <- stats::model.response(model.frame(formula, data))
  categories <- sort(unique(y))

  x_names_betas = colnames(X)[2:length(colnames(X))] # don't include the intercept
  x_names_intercepts <- c()

  for (i in 2:length(categories)){
    x_names_intercepts <- append(x_names_intercepts, paste0("Threshold ", categories[i-1],":", categories[i]))
  }

  x_names <- c(x_names_intercepts, x_names_betas)


  # Prepare the design matrix for the non-proportional odds model if specified
  if (!is.null(partial_formula)) {
    parterms <- terms(partial_formula)
    parpreds <- attributes(parterms)$term.labels
    partial_formula <- reformulate(parpreds, response=NULL, intercept=FALSE)

    Z <- as.matrix(modelr::model_matrix(data, partial_formula))

    # If the intercept is present, it is typically the first column
    if ("(Intercept)" %in% colnames(Z)) {
      # Remove the intercept (first column)
      Z <- Z[, -1]
    }

    z_names <- colnames(Z)

    for (i in 2:length(categories)){
      x_names <- c(x_names, paste0(z_names, ":", categories[i]))
    }
  }

  # Method of moments to determine initial values
  get_initial_values <- function(y, X, Z, family) {
    thresholds <- numeric(length(categories) - 1)
    beta <- numeric(ncol(X))

    for (i in seq_along(thresholds)) {
      current_category <- categories[i]
      mean_y <- mean(y <= current_category)
      if (family == "logit") {
        thresholds[i] <- stats::qlogis(mean_y)
      } else if (family == "probit") {
        thresholds[i] <- stats::qnorm(mean_y)
      }
    }

    beta <- rep(0,ncol(X1)) # Use 0 as starting values

    params <- c(thresholds, beta)
    if (!is.null(Z)) {
      Ncoef <- ncol(Z) * (length(categories)-1)
      gamma <- rep(0,Ncoef) # Use 0 as starting values
      params <- c(params, gamma)
    }
    return(params)
  }

  initial_params <- get_initial_values(y, X, Z, family)

  names(initial_params) <- x_names

  # Ensure ordered thresholds
  logLikFunction <- function(params) {
    k <- length(categories) - 1
    thresholds <- params[1:k]
    beta <- params[(k + 1):(k + ncol(X1))]

    if (!is.null(partial_formula)) {
      gamma <- matrix(params[(k + ncol(X1) + 1):length(params)], ncol = ncol(Z), byrow = TRUE)
    }

    logLik <- 0
    Prob <- rep(0,length(y))
    for (i in 1:length(categories)) {
      if (i==1){
        eta <- -thresholds[i] + X1 %*% beta

        if (!is.null(partial_formula)) {
          eta <- eta + Z %*% gamma[i, ]
        }

        if (family=="logit"){
          Pi <- 1/(1+exp(-eta))
        }
        else{
          Pi <- pnorm(-eta)
        }


        P <- (1-Pi)
        Prob <- ifelse(y==categories[i], P, Prob)
      }
      else {
        eta <- -thresholds[i-1] + X1 %*% beta

        if (!is.null(partial_formula)) {
          eta <- eta + Z %*% gamma[i-1, ]
        }

        if (family=="logit"){
          Pi <- 1/(1+exp(-eta))
        }
        else{
          Pi <- pnorm(-eta)
        }

        if (i<length(categories)){
          eta <- -thresholds[i] + X1 %*% beta

          if (!is.null(partial_formula)) {
            eta <- eta + Z %*% gamma[i, ]
          }

          if (family=="logit"){
            Pi2 <- 1/(1+exp(-eta))
          }
          else{
            Pi2 <- pnorm(-eta)
          }

          P <- Pi-Pi2
        }
        else{
          P <- Pi
        }

        Prob <- ifelse(y==categories[i], P, Prob)
      }
    }

    logLik <- log(max(Prob, .Machine$double.xmin))
    if (method=="BHHH") return(logLik) else return(sum(logLik))
  }

  # Use maxLik to estimate parameters
  result <- maxLik::maxLik(logLik = logLikFunction, start = initial_params, method = method)
  return(result)

}
