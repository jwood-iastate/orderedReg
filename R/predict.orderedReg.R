#' Predict method for orderedReg objects
#'
#' @param object A fitted orderedReg model
#' @param newdata New data for prediction (optional)
#' @param type Type of prediction: "class" or "prob"
#' @param ... Additional arguments
#' @export
predict.orderedReg <- function(object, newdata = NULL, type = "class") {
  # Validate inputs
  if (missing(object)) {
    stop("Model object is required")
  }

  object$family <- attr(object, "family")

  # If no new data, use original data
  if (is.null(newdata)) {
    X <- attr(object, "X")
    Z <- if (!is.null(attr(object, "partial_formula"))) attr(object, "Z") else NULL
  } else {
    # Prepare design matrix for new data
    X <- model.matrix(attr(object, "formula"), data = newdata)
    Z <- if (!is.null(attr(object, "partial_formula"))) {
      model.matrix(attr(object, "partial_formula"), data = newdata)
    } else {
      NULL
    }
  }

  # Safely extract thresholds and other parameters
  N_thresholds <- attr(object, "N_thresholds")

  # Safely extract thresholds
  if (length(object$estimate) < N_thresholds) {
    stop("Insufficient estimates for the number of thresholds")
  }

  thresholds <- object$estimate[1:N_thresholds]

  start_beta <- N_thresholds + 1
  end_beta <- start_beta + ncol(X)-1
  beta <- object$estimate[start_beta:end_beta]

  # Extract beta coefficients
  if (!is.null(object$partial_formula)) {
    beta_partial <- object$estimate[(end_beta + 1):length(object$estimate)]
    }

  # Compute linear predictor
  eta <- X %*% beta

  # If partial formula exists, add Z matrix multiplied byt the coefficients for each level and save to different etas for each level
  # create theta as a matrix of N_thresholds = ncol and each columns starts with the vector of eta
  # then add the Z matrix multiplied by the coefficients for each level
  # then compute the probabilities
  theta <- matrix(rep(eta, N_thresholds), nrow = nrow(X), ncol = N_thresholds)

  if (!is.null(object$partial_formula)) {
    for (i in 1:(N_thresholds)) {
      if (i==1){
        gamma <- beta_partial[1:ncol(Z)]
      }else{
        gamma <- beta_partial[((i - 1) * ncol(Z)+1):(i * ncol(Z))]
      }
      theta[, i] <- eta + Z %*% gamma + thresholds[i]
    }
  } else {
    for (i in 1:(N_thresholds)) {
      theta[, i] <- eta + thresholds[i]
    }
  }

  # Compute Probabilities
  probs <- switch(object$family,
                  "logit" = plogis(theta, lower.tail = FALSE),
                  "probit" = pnorm(theta, lower.tail = FALSE),
                  {
                    warning("Invalid family argument: use 'logit' or 'probit'")
                    NULL
                  })

  # Compute Actual Class Probabilities
  real_probs <- matrix(rep(1, nrow(X)), nrow = nrow(X), ncol = 1)
  real_probs <- cbind(real_probs, probs)
  print(head(real_probs))

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

# S3 method registration
S3method(predict, orderedReg)

