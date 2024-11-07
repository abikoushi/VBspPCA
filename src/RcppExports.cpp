// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ra1_norm
double ra1_norm(const double& a, const double& mean, const double& sd);
RcppExport SEXP _VBspPCA_ra1_norm(SEXP aSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(ra1_norm(a, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// ABsol
arma::mat ABsol(arma::mat A, arma::mat B);
RcppExport SEXP _VBspPCA_ABsol(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(ABsol(A, B));
    return rcpp_result_gen;
END_RCPP
}
// doVB_norm
List doVB_norm(arma::mat Z0, arma::mat W0, const arma::vec& y, const arma::uvec& rowi, const arma::uvec& coli, const int& Nr, const int& Nc, const int& L, const int& iter, const double& prior_prec, const double& a, const double& b);
RcppExport SEXP _VBspPCA_doVB_norm(SEXP Z0SEXP, SEXP W0SEXP, SEXP ySEXP, SEXP rowiSEXP, SEXP coliSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W0(W0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type rowi(rowiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type coli(coliSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nc(NcSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type prior_prec(prior_precSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_norm(Z0, W0, y, rowi, coli, Nr, Nc, L, iter, prior_prec, a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VBspPCA_ra1_norm", (DL_FUNC) &_VBspPCA_ra1_norm, 3},
    {"_VBspPCA_ABsol", (DL_FUNC) &_VBspPCA_ABsol, 2},
    {"_VBspPCA_doVB_norm", (DL_FUNC) &_VBspPCA_doVB_norm, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_VBspPCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
