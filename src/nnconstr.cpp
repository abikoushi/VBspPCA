#include <RcppArmadillo.h>
#include <memory>
#include "truncnorm.h"
#include "KLnorm.h"
#include "nnconstr.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void set_constr(std::unique_ptr<nnconstr> & g,
                const std::string & constr_type){
  if(constr_type == "AN"){
    g.reset(new AN);
  }else if(constr_type == "NN"){
    g.reset(new NN);
  }else if(constr_type == "SN"){
    g.reset(new SN);
  }else{
    Rcpp::stop("This constr_type is not implemented\n");
  }
}

void up_eta_2D(arma::field<arma::mat> & eta,
               const arma::vec & xvlnk,
               const arma::vec & resid,
               const arma::umat X,
               const arma::uvec dims,
               const double & tau, const double & lambda,
               const int & k, const int & l){
  int n = dims(k);
  arma::vec num = arma::zeros<arma::vec>(n);
  num.rows(X.col(k)) += xvlnk % resid;
  eta(k).col(l) = num;
}