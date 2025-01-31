#include "RcppArmadillo.h"
#include "KLgamma.h"
#include "KLnorm.h"
#include "readline.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double lr_default(const double & t,
                  const double & delay,
                  const double & forgetting){
  return pow(t+delay, -forgetting);
}

double calc_elbo(const double & R, const double & obs_prec,
                 const double & ahat, const double & bhat,
                 const double & a, const double & b){
  return -0.5*obs_prec*R + 0.5*(R::digamma(ahat)-log(bhat)) - kl2gamma(ahat, bhat, a, b);
}

void up_lambda(double & obs_prec, double & bhat, const double & ahat, const double & R, const double & b){
  bhat = 0.5*R+b;
  obs_prec = ahat/bhat;
}

void up_eta_w(arma::mat & num_w,
              const arma::vec & y,
              const arma::uvec & rowi,
              const arma::uvec & coli,
              const arma::mat & Z,
              const arma::vec & B) {
  num_w.fill(0.0);
  for(int n=0; n<y.n_rows; n++){
    num_w.row(coli(n)) += Z.row(rowi(n))*(y(n) - B(rowi(n)));
  }
}

void up_eta_z(arma::mat & num_z,
                const arma::vec & y,
                const arma::uvec & rowi,
                const arma::uvec & coli,
                const arma::mat & W,
                const arma::vec & B) {
  num_z.fill(0.0);
  for(int n=0; n<y.n_rows; n++){
    num_z.row(rowi(n)) += W.row(coli(n))*(y(n) - B(rowi(n)));
  }
}

double up_eta_B(arma::vec & num_B,
                const arma::vec & y,
                const arma::vec & y2,
                const arma::uvec & rowi,
                const arma::uvec & coli,
                const arma::mat & Z,
                const arma::mat & W,
                const arma::vec & B) {
  num_B.fill(0.0);
  double R = 0.0;
  for(int n=0; n<y.n_rows; n++){
    double mn = dot(Z.row(rowi(n)), W.row(coli(n))); 
    R += -2.0*mn*(y(n) - B(rowi(n))) +2.0*mn*B(rowi(n))+ y2(n);
    num_B.row(rowi(n)) += y(n) - mn;
  }
  return R;
}

void up_theta(arma::mat & Z,
                arma::mat & W,
                arma::vec & B,
                arma::mat & cov_z,
                arma::mat & cov_w,
                double & cov_B,
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
  arma::vec num_B(B.n_rows);
  int L = W.n_cols;
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L)); 
  arma::mat ZZ = Z.t() * Z + cov_z;
  cov_w = inv(obs_prec*ZZ + prior);
  up_eta_w(num_w, y, rowi, coli, Z, B);
  W = obs_prec*num_w*cov_w;
  arma::mat WW = W.t() * W + cov_w;
  cov_z = inv(obs_prec*WW + prior);
  up_eta_z(num_z, y, rowi, coli, W, B);
  Z = obs_prec*num_z*cov_z;
  R = up_eta_B(num_B, y, y2, rowi, coli, Z, W, B);
  cov_B = 1.0/(N*obs_prec+prior_prec);
  B = (num_B*obs_prec)*cov_B;
  R += sum(arma::diagvec(WW*ZZ)) + sum(B%B + cov_B); //!!!
  //return obs_prec*R;
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
  arma::vec B = arma::randn<arma::vec>(Nr);
  arma::mat cov_z = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat cov_w = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  double cov_B = prior_prec;
  double obs_prec = a/b;
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double bhat = b;
  double R = 0;
  const arma::vec y2 = pow(y,2);
  for (int i=0; i<iter; i++) {
    up_theta(Z, W, B, cov_z, cov_w, cov_B, R, y, y2, rowi, coli, obs_prec, prior_prec, N);
    up_lambda(obs_prec, bhat, ahat, R, b);
    logprob(i) = calc_elbo(R, obs_prec, ahat, bhat, a, b);
  }
  return List::create(Named("mean_row") = Z,
                      Named("mean_col") = W,
                      Named("mean_bias") = B,
                      Named("cov_row") = cov_z,
                      Named("cov_col") = cov_w,
                      Named("cov_bias") = cov_B,
                      Named("obs_prec") = obs_prec,
                      Named("prec_shape") = ahat,
                      Named("prec_rate") = bhat,
                      Named("logprob") = logprob);
}

//////
//mini-batche version
/////
void up_theta_s(arma::mat & Z, 
                  arma::mat & W,
                  arma::vec & B,
                  arma::mat & cov_z,
                  arma::mat & cov_w,
                  double & cov_B,
                  double & R,
                  const arma::vec & y,
                  const arma::vec & y2,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::uvec & uid_r,
                  const arma::uvec & uid_c,
                  const double & obs_prec,
                  const double & prior_prec,
                  const double & NS,
                  const double & N){
  arma::mat num_z(Z.n_rows, Z.n_cols);
  arma::mat num_w(W.n_rows, W.n_cols);
  arma::vec num_B(W.n_rows);
  int L = W.n_cols;
  
  const arma::mat prior = arma::diagmat(prior_prec*arma::ones<arma::vec>(L)); 
  arma::mat ZZ = Z.t() * Z + cov_z;
  arma::mat WW = W.t() * W + cov_w;
  cov_w = NS*inv((obs_prec*WW + prior));
  up_eta_w(num_w, y, rowi, coli, Z, B);
  //W.rows(uid_c) = obs_prec*arma::trans(solve(NS*prec_w, num_w.rows(uid_c).t(), arma::solve_opts::likely_sympd));
  W.rows(uid_c) = obs_prec*num_w.rows(uid_c)*cov_w;
  
  cov_z = NS*inv((obs_prec*WW + prior));
  up_eta_z(num_z, y, rowi, coli, W, B);
  Z.rows(uid_r) = obs_prec*num_z.rows(uid_r)*cov_z;
  //Z.rows(uid_r) = obs_prec*arma::trans(solve(NS*prec_z, num_z.rows(uid_r).t(), arma::solve_opts::likely_sympd));
  
  R = up_eta_B(num_B, y, y2, rowi, coli, Z, W, B);
  cov_B = NS/(N*obs_prec+prior_prec);
  B.rows(uid_r) = num_B.rows(uid_r)*obs_prec*cov_B;
  R += sum(diagvec(WW*ZZ)) + sum(B%B + cov_B); //!!!!
  //return obs_prec*R;
}

double doVB_norm_s_sub(const arma::vec & y,
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
                     arma::mat & Z,
                     arma::mat & W,
                     arma::vec & B,
                     arma::mat & cov_z, 
                     arma::mat & cov_w,
                     double & cov_B,
                     double & obs_prec,
                     double & ahat,
                     double & bhat){
  const arma::uvec uid_r = unique(rowi);
  const arma::uvec uid_c = unique(coli);
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double R = 0;
  double S = y.n_rows;
  double NS = S/N1;
  const arma::vec y2 = pow(y,2);
  double lp;
  for (int i=0; i<iter; i++) {
    up_theta_s(Z, W,  B, cov_z, cov_w, cov_B, R, y, y2, rowi, coli, uid_r, uid_c, obs_prec, prior_prec, NS, N);
    up_lambda(obs_prec, bhat, ahat, R, b);
    lp = calc_elbo(R, obs_prec, ahat, bhat, a, b);
  }
  return lp;
}

void rankindex(arma::uvec & x, const arma::uvec & uid){
  for(int i=0; i<uid.n_rows; i++){
    x.rows(find(uid(i) == x)).fill(i);
  }
}

// [[Rcpp::export]]
List doVB_norm_s_mtx(const std::string & file_path,
                       const int & Nr,
                       const int & Nc,
                       const double & N1,
                       const int & L,
                       const int & ns,
                       const int & iter,
                       const int & subiter,
                       const double & prior_prec,
                       const double & a,
                       const double & b,
                       const double & delay,
                       const double & forgetting){
  arma::mat Z = arma::randn<arma::mat>(Nr, L);
  arma::mat W = arma::randn<arma::mat>(Nc, L);
  arma::vec B = arma::randn<arma::vec>(Nr);
  arma::mat cov_z = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  arma::mat cov_w = arma::diagmat(prior_prec*arma::ones<arma::vec>(L));
  double cov_B = prior_prec;
  double obs_prec = a/b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double bhat = b;
  for(int epoc=0; epoc<iter; epoc++){
    arma::uvec row_i(ns);
    arma::uvec col_i(ns);
    arma::vec val(ns);
    arma::umat bags = randpick_c(N1, ns);
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      readmtx(row_i, col_i, val, file_path, bag);
      
      arma::uvec uid_r = unique(row_i);
      arma::uvec uid_c = unique(col_i);
      
      arma::mat Zs = Z.rows(uid_r);
      arma::mat Ws = W.rows(uid_c);
      arma::mat cov_zs = cov_z;
      arma::mat cov_ws = cov_w;
      arma::vec Bs = B.rows(uid_r);
      double cov_Bs = cov_B;
      double bhat_s = bhat;
      rankindex(row_i, uid_r);
      rankindex(col_i, uid_c);
      lp(epoc) += doVB_norm_s_sub(val, row_i, col_i,
         Nr, Nc, N1, L, subiter, 
         prior_prec, a, b, 
         Zs, Ws, Bs, cov_zs, cov_ws, cov_Bs, obs_prec, ahat, bhat_s);
      double rho = lr_default(epoc, delay, forgetting);
      double rho2 = 1-rho;
      Z.rows(uid_r) = rho2*Z.rows(uid_r) + rho*Zs;
      W.rows(uid_c) = rho2*W.rows(uid_c)+ rho*Ws;
      B.rows(uid_r) = rho2*B.rows(uid_r) + rho*Bs;
      cov_z = rho2*cov_z + rho*cov_zs;
      cov_w = rho2*cov_w + rho*cov_ws;
      cov_B = rho2*cov_B + rho*cov_Bs;
      bhat = rho2*bhat_s + rho*bhat_s;
    }
    //pb.increment();
  }
  return List::create(Named("mean_row") = Z,
                    Named("mean_col") = W,
                    Named("mean_bias") = B,
                    Named("cov_row") = cov_z,
                    Named("cov_col") = cov_w,
                    Named("cov_bias") = cov_B,
                    Named("obs_prec") = obs_prec,
                    Named("prec_shape") = ahat,
                    Named("prec_rate") = bhat,
                    Named("logprob") = lp);
}





/*
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
                 arma::mat Z,
                 arma::mat W,
                 arma::vec B,
                 arma::mat cov_z, 
                 arma::mat cov_w,
                 double cov_B,
                 double obs_prec){
  const arma::uvec uid_r = unique(rowi);
  const arma::uvec uid_c = unique(coli);
  arma::vec logprob = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc; //int to double
  double ahat = 0.5*N+a;
  double R = 0;
  double S = y.n_rows;
  double NS = S/N1;
  double bhat = b;
  const arma::vec y2 = pow(y,2);
  for (int i=0; i<iter; i++) {
    up_theta_s(Z, W,  B, cov_z, cov_w, cov_B, R, y, y2, rowi, coli, uid_r, uid_c, obs_prec, prior_prec, NS, N);
    up_lambda(obs_prec, bhat, ahat, R, b);
    logprob(i) = calc_elbo(R, obs_prec, ahat, bhat, a, b);
  }
  return List::create(Named("mean_row") = Z,
                      Named("mean_col") = W,
                      Named("mean_bias") = B,
                      Named("cov_row") = cov_z,
                      Named("cov_col") = cov_w,
                      Named("cov_bias") = cov_B,
                      Named("obs_prec") = obs_prec,
                      Named("prec_shape") = ahat,
                      Named("prec_rate") = bhat,
                      Named("logprob") = logprob);
}
*/




