#include <RcppArmadillo.h>
#include "nnconstr.h"
#include "readline.h"
#include "lr.h"
#include "KLgamma.h"
#include "ELBO.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <memory>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;

//////
//SVB
//////

double up_etaH_2D_s(arma::field<arma::mat> & eta,
                    arma::mat & H,
                    arma::vec & xvl,
                    arma::vec & xvl2,
                    arma::vec & resid,
                    arma::field<arma::mat> & V,
                    arma::field<arma::mat> & V2,
                    const arma::umat & X, 
                    const double & tau, 
                    const double & lambda,
                    const arma::uvec dims,
                    const int & l, 
                    const double & NS){
  double klv = 0.0;
  int not_k = 1;
  for(int k = 0; k < 2; k++){
    arma::vec vkl = V(k).col(l);
    xvl /= vkl.rows(X.col(k));
    vkl = V2(k).col(l);
    xvl2 /= vkl.rows(X.col(k));
    up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, NS);
    H(k,l) = sum(V2(not_k).col(l));
    vkl = V(k).col(l);
    xvl %=  vkl.rows(X.col(k));
    vkl = V2(k).col(l);
    xvl2 %= vkl.rows(X.col(k));
    not_k = 0;
  }
  return klv;
}

double up_vpar_2D_s(arma::field<arma::mat> & eta,
                    arma::mat & H,
                    arma::field<arma::mat> & V,
                    arma::field<arma::mat> & V2,
                    arma::vec & f,
                    double & sumf2,
                    arma::mat & fm,
                    const arma::vec & y,
                    const arma::umat & X,
                    const arma::uvec dims,
                    const int & L,
                    const double & lambda, const double & tau,
                    const double & NS){
  double klv = 0;
  for(int l = 0; l < L; l++){
    arma::vec vl1 = V(0).col(l);
    arma::vec vl2 = V(1).col(l);
    arma::vec xvl = vl1.rows(X.col(0)) % vl2.rows(X.col(1));
    vl1 = V2(0).col(l);
    vl2 = V2(1).col(l);
    arma::vec xvl2 = vl1.rows(X.col(0)) % vl2.rows(X.col(1));
    f -= xvl; // sum of other than l
    sumf2 -= sum(xvl2);
    arma::vec resid = y - f;
    klv += up_etaH_2D_s(eta, H, xvl, xvl2, resid, V, V2,
                        X, tau, lambda, dims, l, NS);
    f += xvl; // sum of the all
    sumf2 += sum(xvl2);
    fm.col(l) = xvl;
  }
  return klv;
}

double up_bhat_s(double & bhat,
                 const double & ahat, 
                 const double & a, const double & b,
                 const arma::vec & y,
                 const arma::vec & f,
                 const double & sumf2,
                 const arma::mat & fm,
                 const double NS){
  double tmp = 0.0;
  tmp += NS*sum(y%y)*0.5;
  tmp -= NS*sum(y%f);
  tmp += sumf2*0.5;
  tmp +=  NS*upsum(fm);
  bhat = tmp + b;
  return - KLgamma(a, b, ahat, bhat);
}



// [[Rcpp::export]]
Rcpp::List doSVB_norm_woi_diag(arma::field<arma::mat> V,
                               double lambda,
                               const arma::vec & y,
                               const arma::umat & X,
                               const arma::uvec dims,
                               const int & L,
                               const std::string & constr_type,
                               const std::string & lr_type,
                               const arma::vec & lr_param,
                               const int & bsize,
                               const int & iter,
                               const double & tau,
                               const double & a,
                               const double & b,
                               const bool & display_progress){
  int N1 = y.n_rows;
  int K = 2;
  int N = prod(dims);
  arma::vec loglik = arma::zeros<arma::vec>(iter);
  double ahat = 0.5 * ((double) N) + a;
  double bhat = b;
  arma::field<arma::mat> V2(2);
  for(int k = 0; k < K; k++){
    arma::mat Vk = V(k);
    V2(k) = Vk%Vk;
  }
  arma::mat fm = V(0).rows(X.col(0)) % V(1).rows(X.col(1));
  arma::vec f = sum(fm, 1);
  double sumf2 = arma::accu( V2(0).rows(X.col(0)) % V2(1).rows(X.col(1)) );
  arma::field<arma::mat> eta = V;
  arma::mat H(K,L);
  H.fill(0.0);
  for(int k = 0; k < K; k++){
    eta(k).fill(0.0);
  }
  double NS = (double) N1 / (double) bsize;
  std::unique_ptr<nnconstr> constr;
  set_constr(constr, constr_type);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  Progress pb(iter, display_progress);
  for(int epoc = 0; epoc < iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    for(arma::uword step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      arma::vec val(bsize);
      arma::vec Sy = y.rows(bag);
      arma::umat SX(bag.n_rows, K);
      SX = X.rows(bag);
      arma::vec fs = f.rows(bag);
      arma::mat fms = fm.rows(bag);
      arma::field<arma::mat> eta_s = eta;
      arma::mat H_s = H;
      double bhat_s = bhat;
      loglik.row(epoc) += up_vpar_2D_s(eta_s, H_s, V, V2,
                 fs, sumf2, fms, Sy, SX, dims,
                 L, lambda, tau, NS);
      for(int k = 0; k < K; k++){
        eta(k) = rho2*eta(k) + rho*eta_s(k);
      }
      loglik.row(epoc) += up_bhat_s(bhat_s, ahat, a, b, Sy, fs, sumf2, fms, NS);
      H = rho2*H + rho*H_s;
      bhat = rho2*bhat + rho*bhat_s;
      //lambda = ahat/bhat;
      constr -> up_V_from_etaH_2D(V, V2, eta, H, tau, lambda, dims, L);
      fms = V(0).rows(SX.col(0)) % V(1).rows(SX.col(1));
      fs = sum(fms, 1);
      f.rows(bag) = fs;
      fm.rows(bag) = fms;
      sumf2 = arma::accu( V2(0).rows(X.col(0)) % V2(1).rows(X.col(1)) );
    }
    pb.increment();
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}
