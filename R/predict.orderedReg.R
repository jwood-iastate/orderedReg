#' Predict method for orderedReg objects
#'
#' @param object A fitted orderedReg model
#' @param newdata New data for prediction (optional)
#' @param type Type of prediction: "class" or "prob"
#' @param ... Additional arguments
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
#'
#' # Predict probabilities for new data
#' predict(fitted_model, newdata = data, type = "prob")
#'
#' # Predict classes for new data
#' pred_classes <- predict(fitted_model, newdata = data, type = "class")
#'
#' # Generate table of predicted vs observed classes
#' table(pred_classes, data$y)
#'
#' @importFrom stats plogis pnorm model.matrix
#'
#' @export
predict.orderedReg <- function(object, newdata = NULL, type = "class", ...) {
  # Validate inputs
  if (missing(object)) {
    stop("Model object is required")
  }

  family <- attr(object, "family")
  N_thresholds <- attr(object, "N_thresholds")
  formula <- attr(object, "formula")
  Z_list <- if (!is.null(attr(object, "partial_formula"))) attr(object, "Z") else NULL

  # If no new data, use original data
  if (is.null(newdata)) {
    X <- attr(object, "X")
    if (!is.null(Z_list)) {
      Z_matrices <- Z_list
    } else {
      Z_matrices <- NULL
    }
  } else {
    # Prepare design matrix for new data
    X <- model.matrix(formula, data = newdata)
    if (!is.null(Z_list)) {
      Z_matrices <- lapply(Z_list, function(partial_formula) {
        model.matrix(partial_formula, data = newdata)
      })
    } else {
      Z_matrices <- NULL
    }
  }

  # Extract parameters
  thresholds <- object$estimate[1:N_thresholds]
  start_beta <- N_thresholds + 1
  end_beta <- start_beta + ncol(X) - 1
  beta <- object$estimate[start_beta:end_beta]

  # Extract gamma coefficients for each level
  gamma_list <- list()
  if (!is.null(Z_matrices)) {
    current_index <- end_beta + 1
    for (i in seq_along(Z_matrices)) {
      gamma_size <- ncol(Z_matrices[[i]])
      gamma_list[[i]] <- object$estimate[current_index:(current_index + gamma_size - 1)]
      current_index <- current_index + gamma_size
    }
  }

  # Compute linear predictors
  eta <- X %*% beta
  theta <- matrix(rep(eta, N_thresholds), nrow = nrow(X), ncol = N_thresholds)

  if (!is.null(Z_matrices)) {
    for (i in seq_along(Z_matrices)) {
      theta[, i] <- theta[, i] + Z_matrices[[i]] %*% gamma_list[[i]] + thresholds[i]
    }
  } else {
    for (i in 1:N_thresholds) {
      theta[, i] <- theta[, i] + thresholds[i]
    }
  }

  # Compute probabilities
  probs <- switch(family,
                  "logit" = plogis(theta, lower.tail = FALSE),
                  "probit" = pnorm(theta, lower.tail = FALSE),
                  {
                    warning("Invalid family argument: use 'logit' or 'probit'")
                    NULL
                  })

  # Compute actual class probabilities
  real_probs <- matrix(rep(1, nrow(X)), nrow = nrow(X), ncol = 1)
  real_probs <- cbind(real_probs, probs)

  for (i in 1:N_thresholds) {
    real_probs[, i] <- real_probs[, i] - real_probs[, i + 1]
  }

  # Return based on type
  switch(type,
         "class" = apply(real_probs, 1, which.max),
         "prob" = real_probs,
         {
           warning("Invalid type argument: use 'class' or 'prob'")
           NULL
         })
}
