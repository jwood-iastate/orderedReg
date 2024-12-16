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

  for (int i = 1; i < k; i++) {
    // add absolute value of threshold to previous threshold for all but the first threshold
    thresholds[i] = thresholds[i - 1] + std::abs(thresholds[i]);
  }

  // Before main computation
  if (thresholds.has_nan() || beta.has_nan()) {
    return arma::vec(y.n_elem, arma::fill::zeros);
  }

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


// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// This function approximates gradient and Hessian of fun at x using central differences
// fun: an R function that takes a numeric vector and returns a scalar (double).
// x: numeric vector of parameters.
// delta: step size for finite differences.
// Returns a list with "gradient" and "Hessian".


Rcpp::List gradHessApproxCpp(Rcpp::Function fun, Rcpp::NumericVector x, double delta = 1e-4) {
  int nx = x.size();

  // Evaluate f(x)
  Rcpp::NumericVector fxVec = fun(x);
  if (fxVec.size() != 1) {
    Rcpp::stop("Function must return a scalar.");
  }
  double fx = fxVec[0];

  Rcpp::NumericVector grad(nx);
  Rcpp::NumericMatrix H(nx, nx);

  // Compute gradient and Hessian
  // Gradient and diagonal Hessian elements
  for (int j = 0; j < nx; j++) {
    Rcpp::NumericVector xadd = Rcpp::clone(x);
    Rcpp::NumericVector xsub = Rcpp::clone(x);
    xadd[j] = x[j] + delta;
    xsub[j] = x[j] - delta;

    double fadd = as<double>(fun(xadd));
    double fsub = as<double>(fun(xsub));

    // Diagonal Hessian element
    H(j,j) = (fadd - 2.0*fx + fsub)/(delta*delta);
    // Gradient element
    grad[j] = (fadd - fsub)/(2.0*delta);

    // Off-diagonal Hessian elements
    for (int i = 0; i < j; i++) {
      Rcpp::NumericVector xaa = Rcpp::clone(x);
      Rcpp::NumericVector xas = Rcpp::clone(x);
      Rcpp::NumericVector xsa = Rcpp::clone(x);
      Rcpp::NumericVector xss = Rcpp::clone(x);

      // i < j, we compute upper triangle
      xaa[i] = x[i] + delta; xaa[j] = x[j] + delta;
      xas[i] = x[i] + delta; xas[j] = x[j] - delta;
      xsa[i] = x[i] - delta; xsa[j] = x[j] + delta;
      xss[i] = x[i] - delta; xss[j] = x[j] - delta;

      double faa = as<double>(fun(xaa));
      double fas = as<double>(fun(xas));
      double fsa = as<double>(fun(xsa));
      double fss = as<double>(fun(xss));

      double Hij = (faa - fas - fsa + fss)/(4.0*delta*delta);
      H(i,j) = Hij;
      H(j,i) = Hij;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("Hessian") = H
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix hessApproxCpp(Rcpp::Function fun, Rcpp::NumericVector x, double delta=1e-4) {
  int nx = x.size();
  Rcpp::NumericVector fxVec = fun(x);
  if (fxVec.size()!=1) Rcpp::stop("Function must return a scalar.");
  double fx = fxVec[0];

  Rcpp::NumericMatrix H(nx, nx);

  // Diagonal and off-diagonal Hessian elements
  for (int j=0; j<nx; j++) {
    Rcpp::NumericVector xadd = Rcpp::clone(x);
    Rcpp::NumericVector xsub = Rcpp::clone(x);
    xadd[j] = x[j] + delta;
    xsub[j] = x[j] - delta;
    double fadd = as<double>(fun(xadd));
    double fsub = as<double>(fun(xsub));
    H(j,j) = (fadd - 2.0*fx + fsub)/(delta*delta);

    for (int i=0; i<j; i++) {
      Rcpp::NumericVector xaa = Rcpp::clone(x);
      Rcpp::NumericVector xas = Rcpp::clone(x);
      Rcpp::NumericVector xsa = Rcpp::clone(x);
      Rcpp::NumericVector xss = Rcpp::clone(x);

      xaa[i] = x[i]+delta; xaa[j] = x[j]+delta;
      xas[i] = x[i]+delta; xas[j] = x[j]-delta;
      xsa[i] = x[i]-delta; xsa[j] = x[j]+delta;
      xss[i] = x[i]-delta; xss[j] = x[j]-delta;

      double faa = as<double>(fun(xaa));
      double fas = as<double>(fun(xas));
      double fsa = as<double>(fun(xsa));
      double fss = as<double>(fun(xss));

      double Hij = (faa - fas - fsa + fss)/(4.0*delta*delta);
      H(i,j) = Hij;
      H(j,i) = Hij;
    }
  }

  return H;
}


 // [[Rcpp::export]]
 Rcpp::NumericMatrix gradApproxCpp(Rcpp::Function fun, Rcpp::NumericVector x, double delta=1e-4) {
   int nx = x.size();

   if (delta <= 0) {
     Rcpp::stop("Delta must be a positive number.");
   }

   // Evaluate f(x) to determine output size (n)
   Rcpp::NumericVector f0Vec = fun(x);
   int n = f0Vec.size();
   if (n == 0) {
     Rcpp::stop("Function must return a non-empty vector.");
   }

   // Create output matrix: n rows (observations), nx columns (parameters)
   Rcpp::NumericMatrix gradMat(n, nx);

   // For each parameter dimension, compute partial derivatives
   for (int i = 0; i < nx; i++) {
     Rcpp::NumericVector xadd = Rcpp::clone(x);
     Rcpp::NumericVector xsub = Rcpp::clone(x);
     xadd[i] = x[i] + delta;
     xsub[i] = x[i] - delta;

     // Evaluate function at perturbed parameters
     Rcpp::NumericVector fadd = fun(xadd);
     Rcpp::NumericVector fsub = fun(xsub);

     // Check that fadd and fsub have the same length as f0Vec
     if (fadd.size() != n || fsub.size() != n) {
       Rcpp::stop("Function must return vector of consistent length.");
     }

     // Compute central difference for each observation
     for (int obs = 0; obs < n; obs++) {
       gradMat(obs, i) = (fadd[obs] - fsub[obs]) / (2.0 * delta);
     }
   }

   return gradMat;
 }
