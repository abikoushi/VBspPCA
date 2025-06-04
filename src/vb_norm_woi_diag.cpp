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
                    const double & sumZ2,
                    const double & sumZf,
                    const double & upsumfm,
                    const double & sumf2){
  double tmp = sumZ2*0.5 - sumZf + sumf2*0.5 + upsumfm;
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
    loglik.row(i) = up_lambda_2d(lambda, ahat, a, b, sumy2, sum(y%f), upsum(fm), sumf2);
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

double sumdotproduct(const arma::mat V2_row, const arma::mat V2_col){
  double out = 0.0;
  for(arma::uword i = 0; i < V2_row.n_rows; i++){
    for(arma::uword j = 0; j < V2_col.n_rows; j++){
      out += dot(V2_row.row(i), V2_col.row(j));
    }
  }
  return out;
}

double sumupperproduct(const arma::mat V_row, const arma::mat V_col){
  double out = 0.0;
  arma::uword L = V_row.n_cols;
  for(arma::uword i = 0; i < V_row.n_rows; i++){
    for(arma::uword j = 0; j < V_col.n_rows; j++){
      for(arma::uword l = 0; l < (L-1); l++){
        for(arma::uword k = l+1; k < L; k++){
          out += V_row(i,l)*V_col(j,l) * V_row(i,k)*V_col(j,k);
        }
      }
    }
  }
  return out;
}

double up_lambda_2d_mtx(double & lambda,
                        const std::string & readtxt,
                        const arma::field<arma::mat> & V,
                        const arma::field<arma::mat> & V2,
                        const double sumZ2, 
                        const double & ahat, 
                        const double & a,
                        const double & b){
  double sumyf_new = 0.0;
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
    double y = stod(svec[2]);
    double fn = dot(V(0).row(r_i), V(1).row(c_i));
    sumyf_new += y*fn;
  }
  double tmp = sumZ2*0.5 - sumyf_new + 
    0.5*sumdotproduct(V2(0),V2(1)) + sumupperproduct(V(0), V(1));
  double bhat = tmp + b;
  lambda = ahat/bhat;
  return (-lambda*tmp) + 0.5*std::log(lambda) - KLgamma(a, b, ahat, bhat);
}

double up_vpar_2D_mtx(arma::field<arma::mat> & eta,
                  arma::mat & H,
                  arma::field<arma::mat> & V,
                  arma::field<arma::mat> & V2,
                  const std::unique_ptr<nnconstr> & constr,
                  const std::string & readtxt,
                  const arma::uvec dims,
                  const int & L,
                  const double & lambda, const double & tau){
  double klv = 0;
  for(int k=0; k<2; k++){
    eta(k).fill(0);    
  }
  for(int l = 0; l < L; l++){
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
      double y = stod(svec[2]);
      double fn = dot(V(0).row(r_i), V(1).row(c_i));
      arma::vec vl_r = V(0).col(l);
      arma::vec vl_c = V(1).col(l);
      double resid_n = y - (fn - vl_r(r_i)*vl_c(c_i));
      eta(0).col(l).row(r_i) += vl_c(c_i)*resid_n;
      eta(1).col(l).row(c_i) += vl_r(r_i)*resid_n;
      }
    //
      H(0,l) = sum(V2(1).col(l));
      H(1,l) = sum(V2(0).col(l));
      klv += constr -> up_V_from_etaH_2D(V, V2, eta, H, tau, lambda, dims, l);
      // for(int k = 0; k < 2; k++){
      //   arma::vec num = eta(k).col(l);
      //   arma::vec mu = num/(H(k,l) + tau/lambda);
      //   double B = lambda*H(k,l) + tau;
      //   double sigma2 = 1.0 / B;
      //   int n = dims(k);
      //   for(int i = 0; i < n; i++){
      //     klv += KLnorm(mu(i), B, tau);
      //     V(k).col(l).row(i) = mu(i);
      //     V2(k).col(l).row(i) = pow(mu(i), 2) + sigma2;
      //   }
      // }
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
  
  double y;
  std::ifstream file(readtxt);
  std::string str;

  double sumy2 = 0;

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
    double y = stod(svec[2]);
    sumy2 += pow(y, 2);
    index++;
  }
  
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
                                constr, readtxt, dims,
                                L, lambda, tau);
    loglik.row(i) = up_lambda_2d_mtx(lambda, readtxt, V, V2, sumy2, ahat, a, b);
    pb.increment();
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}