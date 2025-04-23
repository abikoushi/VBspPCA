#include "RcppArmadillo.h"
#include "KLgamma.h"
#include "KLnorm.h"
#include "readline.h"
#include "lr.h"
#include "ELBO.h"
#include "lambda.h"
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


void up_theta_woi(arma::mat & Z,
              arma::mat & W,
              arma::mat & ZZ,
              arma::mat & WW,
              arma::mat & cov_z,
              arma::mat & cov_w,
              double & R,
              const arma::vec & y,
              const arma::vec & y2,
              const arma::uvec & rowi,
              const arma::uvec & coli,
              const double & obs_prec,
              const double & prior_prec,
              const double & N){
  arma::mat num_z(Z.n_rows, Z.n_cols);
  arma::mat num_w(W.n_rows, W.n_cols);
  int L = W.n_cols;
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  //up W
  cov_w = inv(obs_prec*ZZ + prior);
  up_eta_w_woi(num_w, y, rowi, coli, Z);
  W = obs_prec*num_w*cov_w;
  WW = W.t() * W + cov_w;
  //up Z
  cov_z = inv(obs_prec*WW + prior);
  up_eta_z_woi(num_z, y, rowi, coli, W);
  Z = obs_prec*num_z*cov_z;
  ZZ = Z.t() * Z + cov_z;
  R = arma::trace(WW*ZZ);
  R += residuals(y, y2, rowi, coli, Z, W);
}

// [[Rcpp::export]]
List doVB_norm_woi(const arma::vec & y,
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
  arma::mat cov_z = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat cov_w = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat ZZ = Z.t() * Z + cov_z;
  arma::mat WW = W.t() * W + cov_w;
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double bhat = b;
  double obs_prec = a/b;
  double R = 0;
  const arma::vec y2 = pow(y,2);
  for (int i=0; i<iter; i++) {
    up_theta_woi(Z, W, ZZ, WW, cov_z, cov_w, R, y, y2, rowi, coli, obs_prec, prior_prec, N);
    up_lambda(obs_prec, bhat, ahat, R, b);
    logprob(i) = calc_elbo(R, obs_prec, ahat, bhat, a, b);
  }
  return List::create(Named("mean_row") = Z,
                      Named("mean_col") = W,
                      Named("cov_row") = cov_z,
                      Named("cov_col") = cov_w,
                      Named("obs_prec") = obs_prec,
                      Named("prec_shape") = ahat,
                      Named("prec_rate") = bhat,
                      Named("logprob") = logprob);
}

//////
//mini-batches version
/////
void up_U(arma::mat U, 
          const arma::vec & y,
          const arma::mat & Z, 
          const arma::mat & W,
          const arma::uvec & rowi,
          const arma::uvec & coli){
  arma::mat Rl = Z.rows(rowi) % W.rows(coli);
  arma::vec R = sum(Rl, 1);
  for(arma::uword l = 0; l < Z.n_cols; l++){
    U.col(l) = y - (R - Rl.col(l));    
  }
}

void up_eta_woi(arma::mat & Z, arma::mat & W,
                arma::mat & cov_z,
                arma::mat & cov_w,
                const arma::mat & U,
                const arma::uvec & rowi, const arma::uvec & coli,
                const arma::uvec & uid_r,
                const arma::uvec & uid_c,
                const double & obs_prec,
                const double & prior_prec,
                const double & NS) {
  arma::mat num_z = Z;
  arma::mat num_w = W;
  num_z.fill(0.0);
  num_w.fill(0.0);
  for(arma::uword n=0; n < U.n_rows; n++){
    for(arma::uword l=0; l < Z.n_cols; l++){
      num_z(rowi(n), l) += W(coli(n), l) * U(n, l);
      num_w(coli(n), l) += Z(rowi(n), l) * U(n, l);
    }
  }
  W.rows(uid_c) = NS*obs_prec*num_w.rows(uid_c)*cov_w;
  Z.rows(uid_r) = NS*obs_prec*num_z.rows(uid_r)*cov_z;
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(Z.n_cols));
  arma::mat WW = W.t() * W + cov_w;
  arma::mat ZZ = Z.t() * Z + cov_z;
  cov_z = inv(obs_prec*WW + prior);
  cov_w = inv(obs_prec*ZZ + prior);
}


// [[Rcpp::export]]
List doVB_norm_wo_s_mtx(const std::string & file_path,
                        const int & Nr,
                        const int & Nc,
                        const double & N1,
                        const int & L,
                        const int & ns,
                        const int & iter,
                        const double & prior_prec,
                        const double & a,
                        const double & b,
                        const arma::vec & lr_param,
                        const std::string & lr_type){
  arma::mat Z = arma::randn<arma::mat>(Nr, L);
  arma::mat W = arma::randn<arma::mat>(Nc, L);
  arma::mat cov_z = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat cov_w = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat ZZ = Z.t() * Z;
  arma::mat WW = W.t() * W;
  double obs_prec = a/b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double bhat = b;
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  for(int epoc = 0; epoc < iter; epoc++){
    arma::uvec row_i(ns);
    arma::uvec col_i(ns);
    arma::vec val(ns);
    arma::umat bags = randpick_c(N1, ns);
    double NS = (double) N1 / (double) ns; 
    for(arma::uword step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      readmtx(row_i, col_i, val, file_path, bag);
      arma::uvec uid_r = unique(row_i);
      arma::uvec uid_c = unique(col_i);
      arma::mat Zs = Z;
      arma::mat Ws = W;
      arma::mat cov_zs = cov_z;
      arma::mat cov_ws = cov_w;
      double bhat_s = bhat;
      arma::mat U = arma::zeros<arma::mat>(bag.n_rows, Z.n_cols);
      up_U(U, val, Zs, Ws, row_i, col_i);
      up_eta_woi(Zs, Ws, cov_zs, cov_ws, U, row_i, col_i, 
                 uid_r, uid_c, obs_prec, prior_prec, NS);
      double rho = g -> lr_t(epoc, lr_param);
      double rho2 = 1.0 - rho;
      Z.rows(uid_r) = rho2*Z.rows(uid_r) + rho*Zs.rows(uid_r);
      W.rows(uid_c) = rho2*W.rows(uid_c) + rho*Ws.rows(uid_c);
      cov_z = rho2*cov_z + rho*cov_zs;
      cov_w = rho2*cov_w + rho*cov_ws;
      bhat = rho2*bhat + rho*bhat_s;
    }
    //pb.increment();
  }
  return List::create(Named("mean_row") = Z,
                      Named("mean_col") = W,
                      Named("cov_row") = cov_z,
                      Named("cov_col") = cov_w,
                      Named("obs_prec") = obs_prec,
                      Named("prec_shape") = ahat,
                      Named("prec_rate") = bhat,
                      Named("logprob") = lp);
}


///
//diagonal covariance
///
