#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void up_lambda(double & obs_prec, 
               double & bhat,
               const double & ahat,
               const double & R, const double & b);
