#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

double KLnorm(double mq, double bq, double tau);

double KLnorm_mat(arma::mat Z, arma::mat W,
                  arma::mat prec_z,
                  arma::mat prec_w,
                  double prior_prec);