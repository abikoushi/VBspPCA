#include <RcppArmadillo.h>
#include "nnconstr.h"
#include "readline.h"
#include "lr.h"
#include "KLgamma.h"
#include "ELBO.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <memory>
#include <iostream>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;

/////
//2D fun
////
double up_lambda_2d(double & lambda,
                    const double & ahat, 
                    const double & a, const double & b,
                    const arma::vec & Z,
                    const double & sumZ2,
                    const arma::vec & f,
                    const double & sumf2,
                    const arma::mat & fm){
  double tmp = 0.0;
  tmp += sumZ2*0.5;
  tmp -= sum(Z%f);
  tmp += sumf2*0.5;
  tmp += upsum(fm);
  double bhat = tmp + b;
  lambda = ahat/bhat;
  return (-lambda*tmp) + 0.5*std::log(lambda) - KLgamma(a, b, ahat, bhat);
}

double up_vpar_2D(arma::field<arma::mat> & eta,
                  arma::mat & H,
                  arma::field<arma::mat> & V,
                  arma::field<arma::mat> & V2,
                  arma::vec & f,
                  double & sumf2,
                  arma::mat & fm,
                  std::unique_ptr<nnconstr> & constr,
                  const arma::vec & y,
                  const arma::umat & X,
                  const arma::uvec dims,
                  const int & L,
                  const double & lambda, const double & tau){
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
    klv += constr -> up_V_eta_H_2D(eta, H, xvl, xvl2, resid, V, V2,
                                   X, tau, lambda, dims, l);
    f += xvl; // sum of the all
    sumf2 += sum(xvl2);
    fm.col(l) = xvl;
  }
  return klv;
}

//without intercept & diagonal posterior
// [[Rcpp::export]]
Rcpp::List doVB_norm_woi_diag(arma::field<arma::mat> V,
                              double lambda,
                              const arma::vec & y,
                              const arma::umat & X,
                              const arma::uvec dims,
                              const int & L,
                              const std::string & constr_type,
                              const int & iter,
                              const double & tau,
                              const double & a, const double & b,
                              const bool & display_progress){
  int N = y.n_rows;
  int K = 2;
  arma::vec loglik = arma::zeros<arma::vec>(iter);
  double ahat = 0.5 * N + a;
  arma::field<arma::mat> V2(2);
  for(int k = 0; k < K; k++){
    arma::mat Vk = V(k);
    V2(k) = Vk%Vk;
  }
  const double sumy2 = sum(y%y);
  arma::mat fm = V(0).rows(X.col(0)) % V(1).rows(X.col(1));
  arma::vec f = sum(fm, 1);
  double sumf2 = arma::accu( V2(0).rows(X.col(0)) % V2(1).rows(X.col(1)) );
  arma::field<arma::mat> eta = V;
  arma::mat H(K,L);
  H.fill(0.0);
  for(int k = 0; k < K; k++){
    eta(k).fill(0.0);
  }
  std::unique_ptr<nnconstr> constr;
  set_constr(constr, constr_type);
  Progress pb(iter, display_progress);
  for(int i = 0; i < iter; i++){
    double klv = up_vpar_2D(eta, H, V, V2,
                            f, sumf2, fm, constr, y, X, dims,
                            L, lambda, tau);
    loglik.row(i) = up_lambda_2d(lambda, ahat, a, b, y, sumy2, f, sumf2, fm);
    pb.increment();
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}

//////
//incremental VB
//////

double up_vpar_2D_mtx(arma::field<arma::mat> & eta,
                  arma::mat & H,
                  arma::field<arma::mat> & V,
                  arma::field<arma::mat> & V2,
                  arma::vec & f,
                  double & sumf2,
                  arma::mat & fm,
                  std::unique_ptr<nnconstr> & constr,
                  const std::string & readtxt,
                  const arma::uvec dims,
                  const int & L,
                  const double & lambda, const double & tau){
  double klv = 0;
  for(int l = 0; l < L; l++){
    arma::vec vl1 = V(0).col(l);
    arma::vec vl2 = V(1).col(l);
    
    arma::vec xvl = f;
    arma::vec xvl2 = f;
    
    arma::vec Z(f.n_rows);
    arma::umat X(f.n_rows, L);
    
    std::ifstream file(readtxt);
    std::string str;

    int n_header = 2;
    for(int i = 0; i < n_header; i++){
      std::getline(file, str);
    }
    int index = 0;
    while (std::getline(file, str)){
      std::stringstream ss(str);
      std::vector<std::string> svec;
      while( ss.good() ){
        std::string substr;
        std::getline(ss, substr, ' ');
        svec.push_back(substr);
      }
      int r_i = stoi(svec[0]) - 1;
      int c_i = stoi(svec[1]) - 1;
      X(index,0) = r_i;
      X(index,1) = c_i;
      double y = stod(svec[2]);
      Z(index) = y;
      arma::vec vl1 = V(0).col(l);
      arma::vec vl2 = V(1).col(l);
      xvl(index) = vl1(r_i) * vl2(c_i);
      vl1 = V2(0).col(l);
      vl2 = V2(1).col(l);
      xvl2(index) = vl1(r_i) * vl2(c_i);
      index++;
    }
    f -= xvl; // sum of other than l
    sumf2 -= sum(xvl2);
    arma::vec resid = Z - f;
    klv += constr -> up_V_eta_H_2D(eta, H,  xvl, xvl2, resid, V, V2, X, tau, lambda, dims, l);
    f += xvl; // sum of the all
    sumf2 += sum(xvl2);
    fm.col(l) = xvl;
  }
  return klv;
}


// [[Rcpp::export]]
Rcpp::List doVB_norm_woi_diag_mtx(arma::field<arma::mat> V,
                              double lambda,
                              const std::string & readtxt,
                              const arma::uvec dims,
                              int N1, int N,
                              const int & L,
                              const std::string & constr_type,
                              const int & iter,
                              const double & tau,
                              const double & a, const double & b,
                              const bool & display_progress){
  int K = 2;
  arma::vec loglik = arma::zeros<arma::vec>(iter);
  double ahat = 0.5 * N + a;
  arma::field<arma::mat> V2(2);
  for(int k = 0; k < K; k++){
    arma::mat Vk = V(k);
    V2(k) = Vk%Vk;
  }
  
  arma::vec Z(N1);
  
  double y;
  std::ifstream file(readtxt);
  std::string str;

  double sumy2 = 0;
  arma::mat fm(N1, L);
  double sumf2 = 0.0;

  int n_header = 2;
  for(int i = 0; i < n_header; i++){
    std::getline(file, str);
  }

  int index = 0;
  while (std::getline(file, str)){
    std::stringstream ss(str);
    std::vector<std::string> svec;
    while( ss.good() ){
      std::string substr;
      std::getline(ss, substr, ' ');
      svec.push_back(substr);
    }
    int r_i = stoi(svec[0]) - 1;
    int c_i = stoi(svec[1]) - 1;
    double y = stod(svec[2]);
    sumy2 += pow(y, 2);
    fm.row(index) = V(0).row(r_i) % V(1).row(c_i);
    sumf2 += arma::accu( V2(0).row(r_i) % V2(1).row(c_i) );
    Z(index) = y;
    index++;
  }
  
  arma::vec f = sum(fm, 1);
  
  arma::field<arma::mat> eta = V;
  arma::mat H(K,L);
  for(int k = 0; k < K; k++){
    eta(k).fill(0.0);
  }
  std::unique_ptr<nnconstr> constr;
  set_constr(constr, constr_type);
  Progress pb(iter, display_progress);
  for(int i = 0; i < iter; i++){
    double klv = up_vpar_2D_mtx(eta, H, V, V2,
                                f, sumf2, fm, constr, readtxt, dims,
                                L, lambda, tau);
    loglik.row(i) = up_lambda_2d(lambda, ahat, a, b, Z, sumy2, f, sumf2, fm);
    pb.increment();
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}