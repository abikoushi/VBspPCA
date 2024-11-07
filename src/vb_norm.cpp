#include "RcppArmadillo.h"
#include "KLgamma.h"
#include "KLnorm.h"
#include "truncnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat ABsol(arma::mat A, arma::mat B){
  return solve(A, B, arma::solve_opts::likely_sympd);
}


//eta :1st order
double up_eta(arma::mat & num_z,
              arma::mat & num_w,
              const arma::vec & y,
              const arma::uvec & rowi,
              const arma::uvec & coli,
              arma::mat & Z,
              arma::mat & W) {
  double lp = 0.0;
  num_z.fill(0.0);
  num_w.fill(0.0);
  for(int n=0; n<y.n_rows; n++){
    lp -= y(n)*arma::dot(Z.row(rowi(n)), W.row(coli(n)));
    num_w.row(coli(n)) += Z.row(rowi(n))*y(n);
    num_z.row(rowi(n)) += W.row(coli(n))*y(n);
  }
  return lp;
}

//up_H :2nd order
double up_H(arma::mat & Z,
            arma::mat & W,
            arma::mat & prec_z,
            arma::mat & prec_w,
            arma::mat & num_z,
            arma::mat & num_w,
            const double & obs_prec,
            const double & prior_prec) {
  int L = W.n_cols;
  /*
  arma::rowvec WW = sum(pow(W,2), 0);
  prec_z = obs_prec*(WW) + prior_prec;
  Z = obs_prec*(num_z.each_row()/prec_z);
  arma::rowvec ZZ = sum(pow(Z,2), 0);
  prec_w = obs_prec*(ZZ) + prior_prec;
  W = obs_prec*(num_w.each_row()/prec_w);
  */
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L)); 
  arma::mat ZZ = Z.t() * Z;
  prec_w = obs_prec*ZZ + prior;
  W = obs_prec*arma::trans(solve(prec_w, num_w.t(), arma::solve_opts::likely_sympd));
  
  arma::mat WW = W.t() * W;
  prec_z = obs_prec*WW + prior;
  Z = obs_prec*arma::trans(solve(prec_z, num_z.t(), arma::solve_opts::likely_sympd));
  //Z.row(0).print();
  //Rprintf("\n");
  return -0.5*accu(ZZ%WW);
}

double up_theta(arma::mat & Z, 
                arma::mat & W,
         arma::mat & prec_z,
         arma::mat & prec_w,
         double & R,
         const arma::vec & y,
         const arma::uvec & rowi,
         const arma::uvec & coli,
         const double & obs_prec,
         const double & prior_prec){
  R = 0;
  arma::mat num_z(Z.n_rows, Z.n_cols);
  arma::mat num_w(W.n_rows, W.n_cols);
  R += up_eta(num_z, num_w, y, rowi, coli, Z, W);
  R += up_H(Z, W, prec_z, prec_w, num_z, num_w, obs_prec, prior_prec);
  return obs_prec*R;
}

void up_lambda(double & lambda,
               const double & ahat, 
               const double & bhat){
  double tmp = 0;
  lambda = bhat/ahat;
}

// [[Rcpp::export]]
List doVB_norm(arma::mat Z0,
               arma::mat W0,
               const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr,
               const int & Nc,
               const int & L,
               const int & iter,
               const double & prior_prec,
               const double & a,
               const double & b){
  //arma::mat Z = arma::randn<arma::mat>(Nr, L);
  //arma::mat W = arma::randn<arma::mat>(Nc, L);
  arma::mat Z = Z0;
  arma::mat W = W0;
  arma::mat prec_z;
  arma::mat prec_w;
  double obs_prec = a/b;
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double R = 0;
  const double sumy2 = 0.5*accu(pow(y,2));
  for (int i=0; i<iter; i++) {
    double lp = up_theta(Z, W, prec_z, prec_w, R, y, rowi, coli, obs_prec, prior_prec);
    //up_lambda(obs_prec, ahat, sumy2+R+b);
    logprob(i) = lp*obs_prec-obs_prec*sumy2;
  }
  return List::create(Named("mean_z") = Z,
                      Named("mean_w") = W,
                      Named("prec_z") = prec_z,
                      Named("prec_w") = prec_w,
                      Named("obs_prec") = obs_prec,
                      Named("logprob") = logprob);
}
