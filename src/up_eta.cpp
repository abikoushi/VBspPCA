#include "up_eta.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void up_eta_w_woi(arma::mat & num_w,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::mat & Z) {
  num_w.fill(0.0);
  for(arma::uword n = 0; n < y.n_rows; n++){
    num_w.row(coli(n)) += Z.row(rowi(n)) * y(n);
  }
}

void up_eta_z_woi(arma::mat & num_z,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::mat & W) {
  num_z.fill(0.0);
  for(arma::uword n = 0; n < y.n_rows; n++){
    num_z.row(rowi(n)) += W.row(coli(n)) * y(n);
  }
}

double residuals(const arma::vec y, const arma::vec y2, 
                 const arma::uvec & rowi, const arma::uvec & coli,
                 const arma::mat Z, const arma::mat W){
  double R = 0.0;
  for(arma::uword n = 0; n < y.n_rows; n++){
    double mn = dot(Z.row(rowi(n)), W.row(coli(n))); 
    R += -2.0*mn*y(n) + y2(n);
  }
  return R;
}