#include "KLnorm.h"

double KLnorm(double mq, double bq, double tau){
  double out = 0.5*(tau/bq - 1.0 - log(tau) + log(bq));
  return isnan(out)?0:out;
}

double KLnorm_mat(arma::mat Z, arma::mat W,
                  arma::mat prec_z,
                  arma::mat prec_w,
                  double prior_prec){
  double klv = 0;
  for(int l=0; l<Z.n_cols; l++){
    for(int i=0; i<Z.n_rows; i++){
      klv += KLnorm(Z(i,l), prec_z(i,l), prior_prec);
    }
    for(int i=0; i<W.n_rows; i++){
      klv += KLnorm(W(i,l), prec_w(i,l), prior_prec);
    }
  }
  return klv;
}