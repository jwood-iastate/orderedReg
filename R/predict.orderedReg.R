#' Predict method for orderedReg objects
#'
#' @param object A fitted orderedReg model.
#' @param newdata Optional new data for which to predict.
#' @param type Type of prediction: "class" or "prob". Default is "class".
#' @param ... Additional arguments.
#'
#' @details
#' For probability predictions, the function returns a matrix of probabilities for each category.
#' For class predictions, the function returns the most likely category for each observation.
#'
#' @return Depending on `type`, either a vector of predicted classes or a matrix of predicted probabilities.
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
#' # Fit model
#' result <- orderedReg(formula = y ~ x1 + x2, data = data, family = "logit")
#'
#' # Predict classes on the training data
#' preds_class <- predict(result, newdata = data, type = "class")
#' table(preds_class, data$y)
#'
#' # Predict probabilities
#' preds_prob <- predict(result, newdata = data, type = "prob")
#' head(preds_prob)
#'
#' @importFrom stats model.matrix plogis pnorm
#' @export
predict.orderedReg <- function(object, newdata = NULL, type = "class", ...) {
  if (missing(object)) {
    stop("A fitted orderedReg model object is required.")
  }

  family <- attr(object, "family")
  formula <- attr(object, "formula")
  categories <- attr(object, "categories")
  N_thresholds <- length(categories)-1
  partial_formula <- attr(object, "partial_formula")
  X_orig <- attr(object, "X")

  # If no newdata provided, use original data
  if (is.null(newdata)) {
    X <- X_orig
    if (!is.null(partial_formula)) {
      # Use the original Z_list stored at fitting time
      Z_list <- attr(object, "Z_list")
    } else {
      Z_list <- NULL
    }
  } else {
    # Re-generate the design matrix for the new data
    # The fitted model used formula with no intercept (as thresholds act as intercepts),
    # so ensure consistent model.matrix call
    X <- model.matrix(update(formula, . ~ . - 1), data = newdata)

    # If partial formulas exist, reconstruct Z_list for new data
    if (!is.null(partial_formula)) {
      if (is.list(partial_formula)) {
        # Multiple formulas, one for each threshold
        if (length(partial_formula) != N_thresholds) {
          stop("partial_formula must have one formula per threshold (excluding the first level).")
        }
        Z_list <- vector("list", N_thresholds)
        for (i in seq_len(N_thresholds)) {
          # Remove intercept
          Z_mat <- model.matrix(update(partial_formula[[i]], . ~ . - 1), data = newdata)
          Z_list[[i]] <- Z_mat
        }
      } else {
        # Single partial formula applied to all thresholds
        Z_mat <- model.matrix(update(partial_formula, . ~ . - 1), data = newdata)
        Z_list <- replicate(N_thresholds, Z_mat, simplify = FALSE)
      }
    } else {
      Z_list <- NULL
    }
  }


  # Extract parameter estimates
  thresholds <- object$estimate[1:N_thresholds]
  start_beta <- N_thresholds + 1
  end_beta <- start_beta + ncol(X) - 1
  beta <- object$estimate[start_beta:end_beta]

  gamma_list <- list()
  if (!is.null(Z_list)) {
    current_index <- end_beta + 1
    for (i in seq_len(N_thresholds)) {
      gamma_size <- ncol(Z_list[[i]])
      gamma_list[[i]] <- object$estimate[current_index:(current_index + gamma_size - 1)]
      current_index <- current_index + gamma_size
    }
  }

  # Compute linear predictors for each threshold
  # The model uses CDF of (threshold_i - (Xb + Zg)), so define:
  # theta_i = threshold_i + (Xb + Zg)
  # and probability(Y<=category_i+1) = F(-theta_i)
  # Here, to ease calculations: P(Y<=j) = plogis(-theta_j) or pnorm(-theta_j)
  eta <- X %*% beta
  theta <- matrix(NA, nrow = nrow(X), ncol = N_thresholds)

  for (i in seq_len(N_thresholds)) {
    theta[, i] <- as.vector(eta) + thresholds[i]
    if (!is.null(Z_list)) {
      theta[, i] <- theta[, i] + as.vector(Z_list[[i]] %*% gamma_list[[i]])
    }
  }

  # Compute cumulative probabilities
  # P(Y <= category[i+1]) = F(-theta[,i]) where F is logistic or normal CDF
  if (family == "logit") {
    # F(-x) = plogis(-x) = 1 - plogis(x)
    # The code can either use plogis(-theta) directly or plogis(theta, lower.tail=FALSE)
    cum_probs <- plogis(-theta)
  } else if (family == "probit") {
    cum_probs <- pnorm(-theta)
  } else {
    stop("Invalid family: must be 'logit' or 'probit'.")
  }

  # Convert cumulative probs to category probabilities
  # P(Y=category_1) = cum_probs[,1]
  # P(Y=category_j) = cum_probs[,j] - cum_probs[,j-1] for j=2,...,N_thresholds
  # P(Y=category_{N_thresholds+1}) = 1 - cum_probs[,N_thresholds]
  real_probs <- cbind(cum_probs, 1) # Add a column for the upper bound
  real_probs <- cbind(real_probs[,1],
                      diff(real_probs),
                      1 - real_probs[,N_thresholds])

  # The above produces a (N x (N_thresholds+2)) matrix. We need only N_thresholds+1 columns.
  # Actually, let's do it step by step:
  real_probs <- matrix(nrow = nrow(X), ncol = N_thresholds + 1)
  real_probs[,1] <- cum_probs[,1]
  for (i in 2:N_thresholds) {
    real_probs[,i] <- cum_probs[,i] - cum_probs[,i-1]
  }
  real_probs[,N_thresholds+1] <- 1 - cum_probs[,N_thresholds]

  # Decide return type
  if (type == "class") {
    # Return the category index that has the highest probability
    # Convert numeric index to category levels if needed
    pred_indices <- max.col(real_probs, ties.method = "first")
    # categories is sorted, match prediction indices to actual categories
    return(categories[pred_indices])
  } else if (type == "prob") {
    return(real_probs)
  } else {
    stop("Invalid type argument: use 'class' or 'prob'.")
  }
}
