// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dmvnrm_arma_mc
arma::vec dmvnrm_arma_mc(arma::mat const& x, arma::rowvec const& mean, arma::mat const& sigma, bool const logd, int const cores);
RcppExport SEXP _fastmvtnorm_dmvnrm_arma_mc(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec const& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    Rcpp::traits::input_parameter< int const >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma_mc(x, mean, sigma, logd, cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastmvtnorm_dmvnrm_arma_mc", (DL_FUNC) &_fastmvtnorm_dmvnrm_arma_mc, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastmvtnorm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
