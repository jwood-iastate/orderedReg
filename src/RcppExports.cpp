// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logLikFunctionCpp
arma::vec logLikFunctionCpp(const arma::vec& params, const arma::vec& categories, const arma::mat& X1, const arma::vec& y, const std::string& family, const List& Z_list, const arma::vec& weights);
RcppExport SEXP _orderedReg_logLikFunctionCpp(SEXP paramsSEXP, SEXP categoriesSEXP, SEXP X1SEXP, SEXP ySEXP, SEXP familySEXP, SEXP Z_listSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type categories(categoriesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const List& >::type Z_list(Z_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikFunctionCpp(params, categories, X1, y, family, Z_list, weights));
    return rcpp_result_gen;
END_RCPP
}
// hessApproxCpp
Rcpp::NumericMatrix hessApproxCpp(Rcpp::Function fun, Rcpp::NumericVector x, double delta);
RcppExport SEXP _orderedReg_hessApproxCpp(SEXP funSEXP, SEXP xSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type fun(funSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(hessApproxCpp(fun, x, delta));
    return rcpp_result_gen;
END_RCPP
}
// gradApproxCpp
Rcpp::NumericMatrix gradApproxCpp(Rcpp::Function fun, Rcpp::NumericVector x, double delta);
RcppExport SEXP _orderedReg_gradApproxCpp(SEXP funSEXP, SEXP xSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type fun(funSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradApproxCpp(fun, x, delta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_orderedReg_logLikFunctionCpp", (DL_FUNC) &_orderedReg_logLikFunctionCpp, 7},
    {"_orderedReg_hessApproxCpp", (DL_FUNC) &_orderedReg_hessApproxCpp, 3},
    {"_orderedReg_gradApproxCpp", (DL_FUNC) &_orderedReg_gradApproxCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_orderedReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
