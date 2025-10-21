#include <RcppArmadillo.h>
#include "nnconstr.h"
#include "readline.h"
#include "lr.h"
#include "KLgamma.h"
#include "ELBO.h"
#include <memory>
#include <iostream>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;

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
        for(arma::uword k = (l+1); k < L; k++){
          out += V_row(i,l)*V_col(j,l) * V_row(i,k)*V_col(j,k);
        }
      }
    }
  }
  return out;
}

double up_lambda_2d_om(double & lambda,
                        const arma::vec y,
                        const arma::umat X,
                        const arma::field<arma::mat> & V,
                        const arma::field<arma::mat> & V2,
                        const double sumZ2, 
                        const double & ahat, 
                        const double & a,
                        const double & b){
  double sumyf_new = sum(y % sum(V(0).rows(X.col(0)) % V(1).rows(X.col(1)), 1));
  double tmp = sumZ2*0.5 - sumyf_new + 0.5*sumdotproduct(V2(0),V2(1)) + sumupperproduct(V(0), V(1));
  double bhat = tmp + b;
  lambda = ahat/bhat;
  return (-lambda*tmp) + 0.5*std::log(lambda) - KLgamma(a, b, ahat, bhat);
}

double up_vpar_2D_om(arma::field<arma::mat> & eta,
                      arma::mat & H,
                      arma::field<arma::mat> & V,
                      arma::field<arma::mat> & V2,
                      const std::unique_ptr<nnconstr> & constr,
                      const arma::vec y,
                      const arma::umat X,
                      const arma::uvec dims,
                      const int & L,
                      const double & lambda, const double & tau){
  double klv = 0;
  for(int k=0; k<2; k++){
    eta(k).fill(0);    
  }
  for(int l = 0; l < L; l++){
    arma::vec vl_r = V(0).col(l);
    arma::vec vl_c = V(1).col(l);
    arma::vec fn = sum(V(0).rows(X.col(0)) % V(1).rows(X.col(1)), 1);
    arma::vec resid_n = y - (fn - vl_r.rows(X.col(0)) % vl_c.rows(X.col(1)));
    arma::vec eta0_l = eta(0).col(l);
    arma::vec eta1_l = eta(1).col(l);
    eta0_l.rows(X.col(0)) += vl_c.rows(X.col(1))%resid_n;
    eta1_l.rows(X.col(1)) += vl_r.rows(X.col(0))%resid_n;
    eta(0).col(l) = eta0_l;
    eta(1).col(l) = eta1_l;
    //
    H(0,l) = sum(V2(1).col(l));
    H(1,l) = sum(V2(0).col(l));
    klv += constr -> up_V_from_etaH_2D(V, V2, eta, H, tau, lambda, dims, l);
  }
  return klv;
}

// [[Rcpp::export]]
Rcpp::List doVB_norm_woi_diag_om(arma::field<arma::mat> V,
                                  double lambda,
                                  const arma::vec y,
                                  const arma::umat X,
                                  const arma::uvec dims,
                                  const int & L,
                                  const std::string & constr_type,
                                  const int & maxit,
                                  const double & tau,
                                  const double & a, const double & b,
                                  const double & tol){
  int K = 2;
  double ahat = 0.5 * prod(dims) + a;
  arma::field<arma::mat> V2(2);
  for(int k = 0; k < K; k++){
    arma::mat Vk = V(k);
    V2(k) = Vk%Vk;
  }
  double sumy2 = sum(y%y);
  arma::field<arma::mat> eta = V;
  arma::mat H(K,L);
  for(int k = 0; k < K; k++){
    eta(k).fill(0.0);
  }
  std::unique_ptr<nnconstr> constr;
  set_constr(constr, constr_type);
  arma::vec loglik = arma::zeros<arma::vec>(maxit+1);
  for(int i = 1; i < maxit + 1; i++){
    double klv = up_vpar_2D_om(eta, H, V, V2, constr, y, X, dims, L, lambda, tau);
    loglik.row(i) = up_lambda_2d_om(lambda, y, X, V, V2, sumy2, ahat, a, b);
    if( abs(loglik(i)-loglik(i-1)) < tol ){
      loglik = loglik.rows(1,i);
      break;
    }
  }
  if(loglik.n_rows == maxit+1){
    loglik.shed_row(0);
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}

//////
//bin
//////

double up_lambda_2d_bin(double & lambda,
                        const int N1,
                        const std::string & readbin_x,
                        const std::string & readbin_y,
                        const arma::field<arma::mat> & V,
                        const arma::field<arma::mat> & V2,
                        const double sumZ2, 
                        const double & ahat, 
                        const double & a,
                        const double & b){
  double sumyf_new = 0.0;
  std::ifstream in_x(readbin_x, std::ios::binary);
  std::ifstream in_y(readbin_y, std::ios::binary);
  for (int i = 0; i < N1; ++i) {
    std::vector<unsigned int> buffer_x(2);
    in_x.read(reinterpret_cast<char*>(buffer_x.data()), 2 * sizeof(unsigned int));
    std::vector<double> buffer_y(1);
    in_y.read(reinterpret_cast<char*>(buffer_y.data()), sizeof(double));
    
    int r_i = buffer_x[0];
    int c_i = buffer_x[1];
    double y = buffer_y[0];
    double fn = dot(V(0).row(r_i), V(1).row(c_i));
    sumyf_new += y*fn;
  }
  
  in_x.close();
  in_y.close();

  double tmp = sumZ2*0.5 - sumyf_new + 
    0.5*sumdotproduct(V2(0),V2(1)) + sumupperproduct(V(0), V(1));
  double bhat = tmp + b;
  lambda = ahat/bhat;
  return (-lambda*tmp) + 0.5*std::log(lambda) - KLgamma(a, b, ahat, bhat);
}

double up_vpar_2D_bin(arma::field<arma::mat> & eta,
                      arma::mat & H,
                      arma::field<arma::mat> & V,
                      arma::field<arma::mat> & V2,
                      const std::unique_ptr<nnconstr> & constr,
                      const int N1,
                      const std::string & readbin_x,
                      const std::string & readbin_y,
                      const arma::uvec dims,
                      const int & L,
                      const double & lambda, const double & tau){
  double klv = 0;
  for(int k=0; k<2; k++){
    eta(k).fill(0);    
  }
  for(int l = 0; l < L; l++){
    std::ifstream in_x(readbin_x, std::ios::binary);
    std::ifstream in_y(readbin_y, std::ios::binary);
    for(int i = 0; i < N1; i++){
    std::vector<unsigned int> buffer_x(2);
    std::vector<double> buffer_y(1);    
    in_x.read(reinterpret_cast<char*>(buffer_x.data()), 2 * sizeof(unsigned int));
    in_y.read(reinterpret_cast<char*>(buffer_y.data()), sizeof(double));
    int r_i = buffer_x[0];
    int c_i = buffer_x[1];
    double y = buffer_y[0];
    double fn = dot(V(0).row(r_i), V(1).row(c_i));
    arma::vec vl_r = V(0).col(l);
    arma::vec vl_c = V(1).col(l);
    double resid_n = y - (fn - vl_r(r_i)*vl_c(c_i));
    eta(0).col(l).row(r_i) += vl_c(c_i)*resid_n;
    eta(1).col(l).row(c_i) += vl_r(r_i)*resid_n;
  }
    in_x.close();
    in_y.close();
    
    //
    H(0,l) = sum(V2(1).col(l));
    H(1,l) = sum(V2(0).col(l));
    klv += constr -> up_V_from_etaH_2D(V, V2, eta, H, tau, lambda, dims, l);
  }
  return klv;
}

// [[Rcpp::export]]
Rcpp::List doVB_norm_woi_diag_bin(arma::field<arma::mat> V,
                                  double lambda,
                                  const int N1,
                                  const std::string & readbin_x,
                                  const std::string & readbin_y,
                                  const arma::uvec dims,
                                  const int & L,
                                  const std::string & constr_type,
                                  const int & maxit,
                                  const double & tau,
                                  const double & a, const double & b,
                                  const double tol){
  int K = 2;
  int N = prod(dims);

  double ahat = 0.5 * N + a;
  arma::field<arma::mat> V2(2);
  for(int k = 0; k < K; k++){
    arma::mat Vk = V(k);
    V2(k) = Vk%Vk;
  }
  
  double sumy2 = 0;
  std::ifstream in_y(readbin_y, std::ios::binary);
  for (int i = 0; i < N1; ++i) {
    std::vector<double> buffer_y(1);
    in_y.read(reinterpret_cast<char*>(buffer_y.data()), sizeof(double));
    double y =  buffer_y[0];
    sumy2 += pow(y, 2);
  }
  in_y.close();

  arma::field<arma::mat> eta = V;
  arma::mat H(K,L);
  for(int k = 0; k < K; k++){
    eta(k).fill(0.0);
  }
  std::unique_ptr<nnconstr> constr;
  set_constr(constr, constr_type);
  arma::vec loglik = arma::zeros<arma::vec>(maxit+1);
  for(int i = 1; i < maxit+1; i++){
    double klv = up_vpar_2D_bin(eta, H, V, V2, 
                                constr, N1, readbin_x, readbin_y, dims,
                                L, lambda, tau);
    loglik.row(i) = up_lambda_2d_bin(lambda, N1, readbin_x, readbin_y, V, V2, sumy2, ahat, a, b);
    if( abs(loglik(i)-loglik(i-1)) < tol ){
      loglik = loglik.rows(1,i);
      break;
    }
  }
  if(loglik.n_rows == maxit+1){
    loglik.shed_row(0);
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}

//////
//mtx
//////

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
  }
  return klv;
}


// [[Rcpp::export]]
Rcpp::List doVB_norm_woi_diag_mtx(arma::field<arma::mat> V,
                                  double lambda,
                                  const std::string & readtxt,
                                  const arma::uvec dims,
                                  const int & L,
                                  const std::string & constr_type,
                                  const int & maxit,
                                  const double & tau,
                                  const double & a, const double & b,
                                  const double tol){
  int K = 2;
  int N = prod(dims);
  arma::vec loglik = arma::zeros<arma::vec>(maxit);
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

  for(int i = 1; i < maxit; i++){
    double klv = up_vpar_2D_mtx(eta, H, V, V2, 
                                constr, readtxt, dims,
                                L, lambda, tau);
    loglik.row(i) = up_lambda_2d_mtx(lambda, readtxt, V, V2, sumy2, ahat, a, b);
    if( abs(loglik(i)-loglik(i-1)) < tol ){
      loglik = loglik.rows(1,i);
      break;
    }
  }
  return Rcpp::List::create(Rcpp::Named("mean_row") = V(0),
                            Rcpp::Named("mean_col") = V(1),
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("H") = H,
                            Rcpp::Named("obs_prec") = lambda,
                            Rcpp::Named("logprob") = loglik);
}