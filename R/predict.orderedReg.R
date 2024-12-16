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
#' # Fit model
#' result <- orderedReg(formula = y ~ x1 + x2, data = data, family = "logit")
#' summary(result)
#'
#' # Predict classes on the training data
#' preds_class <- predict(result, newdata = data, type = "class")
#' table(preds_class, data$y)
#'
#' # Predict probabilities
#' preds_prob <- predict(result, newdata = data, type = "prob")
#' head(preds_prob)
#'
#' @importFrom stats model.matrix plogis pnorm na.omit update
#' @export
predict.orderedReg <- function(object, newdata = NULL, type = "class", ...) {

  if (missing(object) || !inherits(object, "orderedReg")) {
    stop("A valid fitted orderedReg model object is required.")
  }

  type <- match.arg(type, choices = c("class", "prob"))

  if (!is.null(newdata)) {
    required_cols <- all.vars(attr(object, "formula"))
    missing_cols <- setdiff(required_cols, names(newdata))
    if (length(missing_cols) > 0) {
      stop(paste("Missing columns in newdata:", paste(missing_cols, collapse = ", ")))
    }
    if (sum(is.na(newdata)) > 0) {
      warning("Missing values in newdata. Removing rows with missing values.")
      newdata <- na.omit(newdata)
    }
  }

  family <- attr(object, "family")
  formula <- attr(object, "formula")
  categories <- attr(object, "categories")
  N_thresholds <- length(categories)-1
  partial_formula <- attr(object, "partial_formula")
  X_orig <- attr(object, "X")

  if (is.null(newdata)) {
    X <- X_orig
    if (!is.null(partial_formula)) {
      Z_list <- attr(object, "Z_list")
    } else {
      Z_list <- NULL
    }
  } else {
    X <- model.matrix(update(formula, . ~ . - 1), data = newdata)
    if (!is.null(partial_formula)) {
      if (is.list(partial_formula)) {
        if (length(partial_formula) != N_thresholds) {
          stop("partial_formula must have one formula per threshold (excluding the first level).")
        }
        Z_list <- vector("list", N_thresholds)
        for (i in seq_len(N_thresholds)) {
          Z_mat <- model.matrix(update(partial_formula[[i]], . ~ . - 1), data = newdata)
          Z_list[[i]] <- Z_mat
        }
      } else {
        Z_mat <- model.matrix(update(partial_formula, . ~ . - 1), data = newdata)
        Z_list <- replicate(N_thresholds, Z_mat, simplify = FALSE)
      }
    } else {
      Z_list <- NULL
    }
  }

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

  eta <- X %*% beta
  theta <- matrix(NA, nrow = nrow(X), ncol = N_thresholds)

  for (i in seq_len(N_thresholds)) {
    theta[, i] <- as.vector(eta) - thresholds[i]
    if (!is.null(Z_list)) {
      theta[, i] <- theta[, i] + as.vector(Z_list[[i]] %*% gamma_list[[i]])
    }
  }

  if (family == "logit") {
    cum_probs <- plogis(-theta)
  } else if (family == "probit") {
    cum_probs <- pnorm(-theta)
  } else {
    stop("Invalid family: must be 'logit' or 'probit'.")
  }

  real_probs <- matrix(nrow = nrow(X), ncol = N_thresholds + 1)
  real_probs[, 1] <- cum_probs[, 1]
  for (i in 2:N_thresholds) {
    real_probs[, i] <- cum_probs[, i] - cum_probs[, i-1]
  }
  real_probs[, N_thresholds + 1] <- 1 - cum_probs[, N_thresholds]

  epsilon <- .Machine$double.eps
  real_probs <- pmax(pmin(real_probs, 1 - epsilon), epsilon)

  if (type == "class") {
    pred_indices <- max.col(real_probs, ties.method = "first")
    return(categories[pred_indices])
  } else if (type == "prob") {
    return(real_probs)
  } else {
    stop("Invalid type argument: use 'class' or 'prob'.")
  }
}
