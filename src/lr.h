#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

double lr_default(const double & t,
                  const double & delay,
                  const double & forgetting);