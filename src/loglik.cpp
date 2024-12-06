// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec logLikFunctionCpp(const arma::vec& params,
                            const arma::vec& categories,
                            const arma::mat& X1,
                            const arma::vec& y,
                            const std::string& family,
                            const List& Z_list,
                            const arma::vec& weights) {
  int k = categories.n_elem - 1;
  arma::vec thresholds = params.subvec(0, k - 1);
  arma::vec beta = params.subvec(k, k + X1.n_cols - 1);

  int current_index = k + X1.n_cols;
  std::vector<arma::mat> Z_matrices;
  std::vector<arma::vec> gamma_list;

  if (Z_list.size() > 0) {
    for (int i = 0; i < k; i++) {
      arma::mat Z = as<arma::mat>(Z_list[i]);
      Z_matrices.push_back(Z);
      int gamma_size = Z.n_cols;
      arma::vec gamma = params.subvec(current_index, current_index + gamma_size - 1);
      current_index += gamma_size;
      gamma_list.push_back(gamma);
    }
  }

  int n = y.n_elem;
  arma::vec Prob = arma::zeros(n);

  // Compute cumulative probabilities for each category
  // For category i:
  //   if i=1 (lowest category), P(Y=cat_1) = P(Y<=cat_1)
  //   if i=2,... P(Y=cat_i) = P(cat_{i-1}<Y<=cat_i)

  // Loop through categories
  for (int i = 0; i < (int)categories.n_elem; i++) {
    arma::vec eta(n, arma::fill::zeros);
    if (i == 0) {
      // P(Y=cat_1) = 1 - F(threshold_1 - Xb - Zg)
      eta = -thresholds[i] + X1 * beta;
      if (!Z_matrices.empty()) {
        eta += Z_matrices[i] * gamma_list[i];
      }
      arma::vec Pi;
      if (family == "logit") {
        Pi = 1.0 / (1.0 + arma::exp(-eta));
      } else {
        Pi = arma::normcdf(-eta);
      }
      // Probability of first category
      arma::vec P = 1.0 - Pi;
      arma::uvec idx = arma::find(y == categories[i]);
      Prob.elem(idx) = P.elem(idx);
    } else {
      // For category i (i>0), P(Y=cat_i) = F(threshold_i - Xb - Zg) - F(threshold_{i-1} - Xb - Zg) if i < last
      // For last category, P(Y=cat_last) = F(threshold_last - Xb - Zg)

      // upper cumulative probability
      eta = -thresholds[i - 1] + X1 * beta;
      if (!Z_matrices.empty()) {
        eta += Z_matrices[i - 1] * gamma_list[i - 1];
      }
      arma::vec Pi;
      if (family == "logit") {
        Pi = 1.0 / (1.0 + arma::exp(-eta));
      } else {
        Pi = arma::normcdf(-eta);
      }

      arma::vec P;
      if (i < (int)categories.n_elem - 1) {
        // Need lower cumulative probability for the next threshold
        arma::vec eta2 = -thresholds[i] + X1 * beta;
        if (!Z_matrices.empty()) {
          eta2 += Z_matrices[i] * gamma_list[i];
        }
        arma::vec Pi2;
        if (family == "logit") {
          Pi2 = 1.0 / (1.0 + arma::exp(-eta2));
        } else {
          Pi2 = arma::normcdf(-eta2);
        }
        P = Pi - Pi2;
      } else {
        // Last category
        P = Pi;
      }

      arma::uvec idx = arma::find(y == categories[i]);
      Prob.elem(idx) = P.elem(idx);
    }
  }

  // Avoid log(0)
  arma::vec eps_vec = arma::vec(n).fill(std::numeric_limits<double>::min());
  arma::vec logProb = arma::log(arma::max(Prob, eps_vec));
  logProb = logProb % weights; // apply weights

  return logProb;
}
