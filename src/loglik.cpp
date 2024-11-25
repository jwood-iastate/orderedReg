// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec logLikFunctionCpp(const arma::vec& params,
                            const arma::vec& categories,
                            const arma::mat& X1,
                            const arma::vec& y,
                            const std::string& family,
                            const Nullable<arma::mat>& Z_ = R_NilValue,
                            const bool partial = false,
                            const Nullable<arma::vec>& weights_ = R_NilValue) {
  int k = categories.n_elem - 1;
  arma::vec thresholds = params.subvec(0, k - 1);
  arma::vec beta = params.subvec(k, k + X1.n_cols - 1);

  arma::mat Z;
  arma::mat gamma;
  if (partial && Z_.isNotNull()) {
    Z = as<arma::mat>(Z_);
    int n_gamma_params = params.n_elem - (k + X1.n_cols);
    int n_gamma_rows = n_gamma_params / Z.n_cols;
    arma::vec gamma_params = params.subvec(k + X1.n_cols, params.n_elem - 1);
    gamma = arma::reshape(gamma_params, n_gamma_rows, Z.n_cols);
  }

  int n = y.n_elem;
  arma::vec Prob = arma::zeros(n);
  arma::vec eta(n);

  for (unsigned int i = 0; i < categories.n_elem; ++i) {
    if (i == 0) {
      eta = -thresholds[i] + X1 * beta;
      if (partial && Z_.isNotNull()) {
        eta += Z * gamma.row(i).t();
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
      if (partial && Z_.isNotNull()) {
        eta += Z * gamma.row(i - 1).t();
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
        if (partial && Z_.isNotNull()) {
          eta2 += Z * gamma.row(i).t();
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
  arma::vec weights;
  if (weights_.isNotNull()) {
    weights = as<arma::vec>(weights_);
    if (weights.n_elem != n) {
      stop("Length of weights must match the number of observations.");
    }
    logProb = logProb % weights; // Element-wise multiplication
  }

  return logProb;
}
