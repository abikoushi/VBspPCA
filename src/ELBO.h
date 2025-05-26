#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

double calc_elbo(const double & R, const double & obs_prec,
                 const double & ahat, const double & bhat,
                 const double & a, const double & b);

double upsum(const arma::mat & f);
