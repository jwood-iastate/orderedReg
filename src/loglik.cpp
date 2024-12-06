// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec logLikFunctionCpp(const arma::vec& params,
                            const arma::vec& categories,
                            const arma::mat& X1,
                            const arma::vec& y,
                            const std::string& family,
                            const Rcpp::List& Z_list = R_NilValue,
                            const Rcpp::Nullable<arma::vec>& weights_ = R_NilValue) {
  int k = categories.n_elem - 1; // Number of thresholds
  arma::vec thresholds = params.subvec(0, k - 1); // Extract thresholds
  arma::vec beta = params.subvec(k, k + X1.n_cols - 1); // Extract proportional odds coefficients

  // Parse Z_list into a vector of matrices for level-specific design matrices
  std::vector<arma::mat> Z_matrices;
  std::vector<arma::vec> gamma_vectors;
  int expected_params_size = k + X1.n_cols; // Start with thresholds and beta coefficients

  // Check if Z_list is not empty
  if (!Z_list.isNULL() && Z_list.size() > 0) {
    for (int i = 0; i < k; ++i) {
      arma::mat Z = as<arma::mat>(Z_list[i]);
      Z_matrices.push_back(Z);
      expected_params_size += Z.n_cols; // Add gamma size for each level
      arma::vec gamma = params.subvec(expected_params_size - Z.n_cols, expected_params_size - 1);
      gamma_vectors.push_back(gamma);
    }
  }

  // Check parameter size consistency
  if (params.n_elem != static_cast<arma::uword>(expected_params_size)) {
    stop("Mismatch in the size of 'params'. Expected size: " + std::to_string(expected_params_size) +
      ", but got: " + std::to_string(params.n_elem));
  }

  int n = y.n_elem; // Number of observations
  arma::vec Prob = arma::zeros(n); // Probabilities
  arma::vec eta(n); // Linear predictor

  for (unsigned int i = 0; i < categories.n_elem; ++i) {
    if (i == 0) {
      eta = -thresholds[i] + X1 * beta;
      if (!Z_matrices.empty()) {
        eta += Z_matrices[i] * gamma_vectors[i];
      }
      arma::vec Pi;
      if (family == "logit") {
        Pi = 1 / (1 + arma::exp(-eta));
      } else {
        Pi = arma::normcdf(-eta);
      }
      arma::vec P = 1 - Pi;
      arma::uvec idx = arma::find(y == categories[i]);
      Prob.elem(idx) = P.elem(idx);
    } else {
      eta = -thresholds[i - 1] + X1 * beta;
      if (!Z_matrices.empty()) {
        eta += Z_matrices[i - 1] * gamma_vectors[i - 1];
      }
      arma::vec Pi;
      if (family == "logit") {
        Pi = 1 / (1 + arma::exp(-eta));
      } else {
        Pi = arma::normcdf(-eta);
      }
      arma::vec P;
      if (i < categories.n_elem - 1) {
        arma::vec eta2 = -thresholds[i] + X1 * beta;
        if (!Z_matrices.empty()) {
          eta2 += Z_matrices[i] * gamma_vectors[i];
        }
        arma::vec Pi2;
        if (family == "logit") {
          Pi2 = 1 / (1 + arma::exp(-eta2));
        } else {
          Pi2 = arma::normcdf(-eta2);
        }
        P = Pi - Pi2;
      } else {
        P = Pi;
      }
      arma::uvec idx = arma::find(y == categories[i]);
      Prob.elem(idx) = P.elem(idx);
    }
  }

  // Ensure probabilities are above machine epsilon to avoid log(0)
  arma::vec eps = arma::vec(n).fill(std::numeric_limits<double>::min());
  arma::vec logProb = arma::log(arma::max(Prob, eps));

  // Handle weights
  if (weights_.isNotNull()) {
    arma::vec weights = as<arma::vec>(weights_);
    if (static_cast<int>(weights.n_elem) != n) { // Cast weights.n_elem for comparison
      stop("Length of weights must match the number of observations.");
    }
    logProb = logProb % weights; // Element-wise multiplication
  }

  return logProb;
}
