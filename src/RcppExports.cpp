// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// read_mtx
List read_mtx(const std::string& readtxt, const arma::uvec& bag);
RcppExport SEXP _VBspPCA_read_mtx(SEXP readtxtSEXP, SEXP bagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type readtxt(readtxtSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type bag(bagSEXP);
    rcpp_result_gen = Rcpp::wrap(read_mtx(readtxt, bag));
    return rcpp_result_gen;
END_RCPP
}
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
// doVB_norm
List doVB_norm(const arma::vec& y, const arma::uvec& rowi, const arma::uvec& coli, const int& Nr, const int& Nc, const int& L, const int& iter, const double& prior_prec, const double& a, const double& b);
RcppExport SEXP _VBspPCA_doVB_norm(SEXP ySEXP, SEXP rowiSEXP, SEXP coliSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
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
    rcpp_result_gen = Rcpp::wrap(doVB_norm(y, rowi, coli, Nr, Nc, L, iter, prior_prec, a, b));
    return rcpp_result_gen;
END_RCPP
}
// doVB_norm_s_mtx
List doVB_norm_s_mtx(const std::string& file_path, const int& Nr, const int& Nc, const double& N1, const int& L, const int& ns, const int& iter, const int& subiter, const double& prior_prec, const double& a, const double& b, const double& delay, const double& forgetting);
RcppExport SEXP _VBspPCA_doVB_norm_s_mtx(SEXP file_pathSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP N1SEXP, SEXP LSEXP, SEXP nsSEXP, SEXP iterSEXP, SEXP subiterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP, SEXP delaySEXP, SEXP forgettingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nc(NcSEXP);
    Rcpp::traits::input_parameter< const double& >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type subiter(subiterSEXP);
    Rcpp::traits::input_parameter< const double& >::type prior_prec(prior_precSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< const double& >::type forgetting(forgettingSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_norm_s_mtx(file_path, Nr, Nc, N1, L, ns, iter, subiter, prior_prec, a, b, delay, forgetting));
    return rcpp_result_gen;
END_RCPP
}
// doVB_norm_woi
List doVB_norm_woi(const arma::vec& y, const arma::uvec& rowi, const arma::uvec& coli, const int& Nr, const int& Nc, const int& L, const int& iter, const double& prior_prec, const double& a, const double& b);
RcppExport SEXP _VBspPCA_doVB_norm_woi(SEXP ySEXP, SEXP rowiSEXP, SEXP coliSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
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
    rcpp_result_gen = Rcpp::wrap(doVB_norm_woi(y, rowi, coli, Nr, Nc, L, iter, prior_prec, a, b));
    return rcpp_result_gen;
END_RCPP
}
// doVB_norm_wo_s_mtx
List doVB_norm_wo_s_mtx(const std::string& file_path, const int& Nr, const int& Nc, const double& N1, const int& L, const int& ns, const int& iter, const int& subiter, const double& prior_prec, const double& a, const double& b, const double& delay, const double& forgetting);
RcppExport SEXP _VBspPCA_doVB_norm_wo_s_mtx(SEXP file_pathSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP N1SEXP, SEXP LSEXP, SEXP nsSEXP, SEXP iterSEXP, SEXP subiterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP, SEXP delaySEXP, SEXP forgettingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nc(NcSEXP);
    Rcpp::traits::input_parameter< const double& >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type subiter(subiterSEXP);
    Rcpp::traits::input_parameter< const double& >::type prior_prec(prior_precSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< const double& >::type forgetting(forgettingSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_norm_wo_s_mtx(file_path, Nr, Nc, N1, L, ns, iter, subiter, prior_prec, a, b, delay, forgetting));
    return rcpp_result_gen;
END_RCPP
}
// doVB_norm_woi_diag
List doVB_norm_woi_diag(const arma::vec& y, const arma::uvec& rowi, const arma::uvec& coli, const int& Nr, const int& Nc, const int& L, const int& iter, const double& prior_prec, const double& a, const double& b);
RcppExport SEXP _VBspPCA_doVB_norm_woi_diag(SEXP ySEXP, SEXP rowiSEXP, SEXP coliSEXP, SEXP NrSEXP, SEXP NcSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP prior_precSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
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
    rcpp_result_gen = Rcpp::wrap(doVB_norm_woi_diag(y, rowi, coli, Nr, Nc, L, iter, prior_prec, a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VBspPCA_read_mtx", (DL_FUNC) &_VBspPCA_read_mtx, 2},
    {"_VBspPCA_ra1_norm", (DL_FUNC) &_VBspPCA_ra1_norm, 3},
    {"_VBspPCA_doVB_norm", (DL_FUNC) &_VBspPCA_doVB_norm, 10},
    {"_VBspPCA_doVB_norm_s_mtx", (DL_FUNC) &_VBspPCA_doVB_norm_s_mtx, 13},
    {"_VBspPCA_doVB_norm_woi", (DL_FUNC) &_VBspPCA_doVB_norm_woi, 10},
    {"_VBspPCA_doVB_norm_wo_s_mtx", (DL_FUNC) &_VBspPCA_doVB_norm_wo_s_mtx, 13},
    {"_VBspPCA_doVB_norm_woi_diag", (DL_FUNC) &_VBspPCA_doVB_norm_woi_diag, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_VBspPCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
