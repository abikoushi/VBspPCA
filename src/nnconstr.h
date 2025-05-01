#include <RcppArmadillo.h>
#include <memory>
#include "truncnorm.h"
#include "KLnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class nnconstr {
public:
  virtual double up_V_from_etaH_2D(const arma::field<arma::mat> & eta,
                                   const arma::mat & H,
                                   const double & tau, 
                                   const double & lambda,
                                   arma::field<arma::mat> & V,
                                   arma::field<arma::mat> & V2,
                                   const arma::uvec dims,
                                   const int & k, const int & l) = 0;
  virtual ~nnconstr(){}
};

//non negative
class NN : public nnconstr{
  double up_V_from_etaH_2D(const arma::field<arma::mat> & eta,
                           const arma::mat & H,
                           const double & tau, 
                           const double & lambda,
                           arma::field<arma::mat> & V,
                           arma::field<arma::mat> & V2,
                           const arma::uvec dims,
                           const int & k, const int & l){
    int n = dims(k);
    arma::vec num = eta(k).col(l);
    double den = H(k,l);
    arma::vec mu = num/(H(k,l) + tau/lambda);
    double B = lambda*H(k,l) + tau;
    double klv = 0.0;
    double sigma = 1.0 / sqrt(B);
    for(int i = 0; i < n; i++){
      klv += KLtruncnorm(mu(i), B, tau);
      V(k).col(l).row(i) = truncmoment1(mu(i), sigma);
      V2(k).col(l).row(i) = truncmoment2(mu(i), sigma);
  }
    return klv;
  }
};

//Arrow negative
class AN : public nnconstr{
  double up_V_from_etaH_2D(const arma::field<arma::mat> & eta,
                           const arma::mat & H,
                           const double & tau, 
                           const double & lambda,
                           arma::field<arma::mat> & V,
                           arma::field<arma::mat> & V2,
                           const arma::uvec dims,
                           const int & k, const int & l){
    int n = dims(k);
    arma::vec num = eta(k).col(l);
    double den = H(k,l);
    arma::vec mu = num/(H(k,l) + tau/lambda);
    double B = lambda*H(k,l) + tau;
    double klv = 0.0;
    double sigma2 = 1.0 / B;
    for(int i = 0; i < n; i++){
      klv += KLnorm(mu(i), B, tau);
      V(k).col(l).row(i) = mu(i);
      V2(k).col(l).row(i) = pow(mu(i),2) + sigma2;
    }
    return klv;
  }
};

void set_constr(std::unique_ptr<nnconstr> & g,
                const std::string & constr_type);
