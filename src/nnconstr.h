#include <RcppArmadillo.h>
#include <memory>
#include "truncnorm.h"
#include "KLnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void up_eta_2D(arma::field<arma::mat> & eta,
               const arma::vec & xvlnk,
               const arma::vec & resid,
               const arma::umat X,
               const arma::uvec dims,
               const double & tau, const double & lambda,
               const int & k, const int & l,
               const double & NS);

class nnconstr {
public:
  virtual double up_V_eta_H_2D(arma::field<arma::mat> & eta,
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
                                   const int & l) = 0;
  
  // virtual double up_V_eta_H_2D_mtx(arma::field<arma::mat> & eta,
  //                              arma::mat & H,
  //                              arma::vec & xvl,
  //                              arma::vec & xvl2,
  //                              arma::vec & resid,
  //                              arma::field<arma::mat> & V,
  //                              arma::field<arma::mat> & V2,
  //                              const std::string & readtxt,
  //                              const double & tau, 
  //                              const double & lambda,
  //                              const arma::uvec dims,
  //                              const int & l) = 0;
  
  virtual void up_V_from_etaH_2D(arma::field<arma::mat> & V,
                                   arma::field<arma::mat> & V2,
                                   const arma::field<arma::mat> & eta,
                                   const arma::mat & H,
                                   const double & tau, 
                                   const double & lambda,
                                   const arma::uvec dims,
                                   const int & L) = 0;
  virtual ~nnconstr(){}
};

//nonnegative
class NN : public nnconstr{
  double up_V_eta_H_2D(arma::field<arma::mat> & eta,
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
                           const int & l){
    int not_k = 1;
    double klv = 0.0;
    for(int k = 0; k < 2; k++){
      arma::vec vkl = V(k).col(l);
      xvl /= vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 /= vkl.rows(X.col(k));
      up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, 1.0);
      H(k,l) = sum(V2(not_k).col(l));
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      double den = H(k,l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma = 1.0 / sqrt(B);
      for(int i = 0; i < n; i++){
        klv += KLtruncnorm(mu(i), B, tau);
        V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
        V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
      }
      vkl = V(k).col(l);
      xvl %=  vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 %= vkl.rows(X.col(k));
      not_k = 0;
    }
    return klv;
  }
  
  // double up_V_eta_H_2D_mtx(arma::field<arma::mat> & eta,
  //                                  arma::mat & H,
  //                                  arma::vec & xvl,
  //                                  arma::vec & xvl2,
  //                                  arma::vec & resid,
  //                                  arma::field<arma::mat> & V,
  //                                  arma::field<arma::mat> & V2,
  //                                  const std::string & readtxt,
  //                                  const double & tau, 
  //                                  const double & lambda,
  //                                  const arma::uvec dims,
  //                                  const int & l){
  //   int not_k = 1;
  //   double klv = 0.0;
  //   for(int k = 0; k < 2; k++){
  //     arma::vec vkl = V(k).col(l);
  //     xvl /= vkl.rows(X.col(k));
  //     vkl = V2(k).col(l);
  //     xvl2 /= vkl.rows(X.col(k));
  //     up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, 1.0);
  //     H(k,l) = sum(V2(not_k).col(l));
  //     int n = dims(k);
  //     arma::vec num = eta(k).col(l);
  //     double den = H(k,l);
  //     arma::vec mu = num/(H(k,l) + tau/lambda);
  //     double B = lambda*H(k,l) + tau;
  //     double sigma = 1.0 / sqrt(B);
  //     for(int i = 0; i < n; i++){
  //       klv += KLtruncnorm(mu(i), B, tau);
  //       V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
  //       V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
  //     }
  //     vkl = V(k).col(l);
  //     xvl %=  vkl.rows(X.col(k));
  //     vkl = V2(k).col(l);
  //     xvl2 %= vkl.rows(X.col(k));
  //     not_k = 0;
  //   }
  //   return klv;
  // }
  
  void up_V_from_etaH_2D(arma::field<arma::mat> & V,
                           arma::field<arma::mat> & V2,
                           const arma::field<arma::mat> & eta,
                           const arma::mat & H,
                           const double & tau, 
                           const double & lambda,
                           const arma::uvec dims,
                           const int & L){
  double klv = 0;
  for(int l = 0; l < L; l++){
    for(int k = 0; k < 2; k++){
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma = 1.0 / sqrt(B);
      for(int i = 0; i < n; i++){
        klv += KLtruncnorm(mu(i), B, tau);
        V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
        V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
      }
    }
  }
  }
};

//arrow negative
class AN : public nnconstr{
  double up_V_eta_H_2D(arma::field<arma::mat> & eta,
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
                           const int & l){
    int not_k = 1;
    double klv = 0.0;
    for(int k = 0; k < 2; k++){
      arma::vec vkl = V(k).col(l);
      xvl /= vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 /= vkl.rows(X.col(k));
      up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, 1.0);
      H(k,l) = sum(V2(not_k).col(l));
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma2 = 1.0 / B;
      for(int i = 0; i < n; i++){
        klv += KLnorm(mu(i), B, tau);
        V(k).col(l).row(i) = mu(i);
        V2(k).col(l).row(i) = pow(mu(i),2) + sigma2;
      }
      vkl = V(k).col(l);
      xvl %=  vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 %= vkl.rows(X.col(k));
      not_k = 0;
    }
    return klv;
  }
  void up_V_from_etaH_2D(arma::field<arma::mat> & V,
                                   arma::field<arma::mat> & V2,
                                   const arma::field<arma::mat> & eta,
                                   const arma::mat & H,
                                   const double & tau, 
                                   const double & lambda,
                                   const arma::uvec dims,
                                   const int & L){
    double klv = 0.0;
    for(int l = 0; l < L; l++){
    for(int k = 0; k < 2; k++){
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma2 = 1.0 / B;
      for(int i = 0; i < n; i++){
        klv += KLnorm(mu(i), B, tau);
        V(k).col(l).row(i) = mu(i);
        V2(k).col(l).row(i) = pow(mu(i),2) + sigma2;
      }
    }
    }
  }
};

//semi-nonnegative
class SN : public nnconstr{
  double up_V_eta_H_2D(arma::field<arma::mat> & eta,
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
                           const int & l){
    int not_k = 1;
    double klv = 0.0;
    for(int k = 0; k < 1; k++){
      arma::vec vkl = V(k).col(l);
      xvl /= vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 /= vkl.rows(X.col(k));
      up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, 1.0);
      H(k,l) = sum(V2(not_k).col(l));
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      double den = H(k,l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma2 = 1.0 / B;
      for(int i = 0; i < n; i++){
        klv += KLnorm(mu(i), B, tau);
        V(k).col(l).row(i) = mu(i);
        V2(k).col(l).row(i) = pow(mu(i),2) + sigma2;
      }
      vkl = V(k).col(l);
      xvl %=  vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 %= vkl.rows(X.col(k));
      not_k = 0;
    }
    for(int k = 1; k < 2; k++){
      arma::vec vkl = V(k).col(l);
      xvl /= vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 /= vkl.rows(X.col(k));
      up_eta_2D(eta, xvl, resid, X, dims, tau, lambda, k, l, 1.0);
      H(k,l) = sum(V2(not_k).col(l));
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      double den = H(k,l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma = 1.0 / sqrt(B);
      for(int i = 0; i < n; i++){
        klv += KLtruncnorm(mu(i), B, tau);
        V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
        V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
      }
      vkl = V(k).col(l);
      xvl %=  vkl.rows(X.col(k));
      vkl = V2(k).col(l);
      xvl2 %= vkl.rows(X.col(k));
    }
    return klv;
  }
  void up_V_from_etaH_2D(arma::field<arma::mat> & V,
                                   arma::field<arma::mat> & V2,
                                   const arma::field<arma::mat> & eta,
                                   const arma::mat & H,
                                   const double & tau, 
                                   const double & lambda,
                                   const arma::uvec dims,
                                   const int & L){
    double klv = 0.0;
    for(int l = 0; l < L; l++){
    for(int k = 0; k < 1; k++){
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      double den = H(k,l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma2 = 1.0 / B;
      for(int i = 0; i < n; i++){
        klv += KLnorm(mu(i), B, tau);
        V(k).col(l).row(i) = mu(i);
        V2(k).col(l).row(i) = pow(mu(i),2) + sigma2;
      }
    }
    for(int k = 1; k < 2; k++){
      arma::vec vkl = V(k).col(l);
      int n = dims(k);
      arma::vec num = eta(k).col(l);
      double den = H(k,l);
      arma::vec mu = num/(H(k,l) + tau/lambda);
      double B = lambda*H(k,l) + tau;
      double sigma = 1.0 / sqrt(B);
      for(int i = 0; i < n; i++){
        klv += KLtruncnorm(mu(i), B, tau);
        V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
        V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
      }
    }
    }
  }
};

void set_constr(std::unique_ptr<nnconstr> & g,
                const std::string & constr_type);
