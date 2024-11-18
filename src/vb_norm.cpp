#include "RcppArmadillo.h"
#include "KLgamma.h"
#include "KLnorm.h"
#include "truncnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void up_eta_w(arma::mat & num_w,
              const arma::vec & y,
              const arma::uvec & rowi,
              const arma::uvec & coli,
              const arma::mat & Z) {

  num_w.fill(0.0);
  for(int n=0; n<y.n_rows; n++){
    num_w.row(coli(n)) += Z.row(rowi(n))*y(n);
  }
}

double up_eta_z(arma::mat & num_z,
              const arma::vec & y,
              const arma::uvec & rowi,
              const arma::uvec & coli,
              const arma::mat & Z,
              const arma::mat & W) {
  num_z.fill(0.0);
  double R = 0;
  for(int n=0; n<y.n_rows; n++){
    R += dot(Z.row(rowi(n)),W.row(coli(n)))*y(n);
    num_z.row(rowi(n)) += W.row(coli(n))*y(n);
  }
  return 2.0*R;
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
  arma::mat num_z(Z.n_rows, Z.n_cols);
  arma::mat num_w(W.n_rows, W.n_cols);
  int L = W.n_cols;
  
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L)); 
  arma::mat ZZ = Z.t() * Z;
  prec_w = obs_prec*ZZ + prior;
  up_eta_w(num_w, y, rowi, coli, Z);
  W = obs_prec*arma::trans(solve(prec_w, num_w.t(), arma::solve_opts::likely_sympd));
  
  arma::mat WW = W.t() * W;
  prec_z = obs_prec*WW + prior;
  R = up_eta_z(num_z,y,rowi,coli, Z, W);
  Z = obs_prec*arma::trans(solve(prec_z, num_z.t(), arma::solve_opts::likely_sympd));
  R += accu(WW%ZZ); //!!!
  return obs_prec*R;
}

// [[Rcpp::export]]
List doVB_norm(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr,
               const int & Nc,
               const int & L,
               const int & iter,
               const double & prior_prec,
               const double & a,
               const double & b){
  arma::mat Z = arma::randn<arma::mat>(Nr, L);
  arma::mat W = arma::randn<arma::mat>(Nc, L);
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
    logprob(i) = -0.5*obs_prec*(sumy2+R);
  }
  return List::create(Named("mean_z") = Z,
                      Named("mean_w") = W,
                      Named("prec_z") = prec_z,
                      Named("prec_w") = prec_w,
                      Named("obs_prec") = obs_prec,
                      Named("logprob") = logprob);
}


//////
//mini-batche version
/////

double up_theta_s(arma::mat & Z, 
                arma::mat & W,
                arma::mat & prec_z,
                arma::mat & prec_w,
                double & R,
                const arma::vec & y,
                const arma::uvec & rowi,
                const arma::uvec & coli,
                const arma::uvec & uid_r,
                const arma::uvec & uid_c,
                const double & obs_prec,
                const double & prior_prec,
                const double & NS){
  arma::mat num_z(Z.n_rows, Z.n_cols);
  arma::mat num_w(W.n_rows, W.n_cols);
  int L = W.n_cols;
  
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L)); 
  arma::mat ZZ = Z.t() * Z;
  prec_w = obs_prec*ZZ + prior;
  up_eta_w(num_w, y, rowi, coli, Z);
  W.rows(uid_c) = obs_prec*arma::trans(solve(NS*prec_w, num_w.rows(uid_c).t(), arma::solve_opts::likely_sympd));
  
  arma::mat WW = W.t() * W;
  prec_z = obs_prec*WW + prior;
  R = up_eta_z(num_z, y, rowi, coli, Z, W);
  Z.rows(uid_r) = obs_prec*arma::trans(solve(NS*prec_z, num_z.rows(uid_r).t(), arma::solve_opts::likely_sympd));
  R += accu(WW%ZZ); //!!!
  return obs_prec*R;
}

// [[Rcpp::export]]
List doVB_norm_s(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr,
               const int & Nc,
               const int & L,
               const int & iter,
               const double & prior_prec,
               const double & a,
               const double & b,
               const double & N1, 
               arma::mat Z, arma::mat W,
               arma::mat prec_z, arma::mat prec_w){
  const arma::uvec uid_r = unique(rowi);
  const arma::uvec uid_c = unique(coli);
  double obs_prec = a/b;
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double R = 0;
  double S = y.n_rows;
  double NS = S/N1;
  const double sumy2 = 0.5*accu(pow(y,2));
  for (int i=0; i<iter; i++) {
    double lp = up_theta_s(Z, W, prec_z, prec_w, R, y, rowi, coli, uid_r, uid_c, obs_prec, prior_prec, NS);
    //up_lambda(obs_prec, ahat, sumy2+R+b);
    logprob(i) = -0.5*obs_prec*(sumy2+R);
  }
  return List::create(Named("mean_z") = Z,
                      Named("mean_w") = W,
                      Named("prec_z") = prec_z,
                      Named("prec_w") = prec_w,
                      Named("obs_prec") = obs_prec,
                      Named("logprob") = logprob);
}
