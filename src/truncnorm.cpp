#include <RcppArmadillo.h>
#include "truncnorm.h"

// static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(const double & a) {
  double ainv = 1.0 / a;
  double x;
  double rho;
  do {
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0, 1) > rho);
  return x;
}

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(const double & a) {
  double x = -DBL_MAX;
  while (x < a) {
    x = R::rnorm(0, 1);
  }
  return x;
}

/* Exponential rejection sampling (a, b) */
double ers_ab(const double & a, const double & b) {
  double ainv = 1.0 / a;
  double x;
  double rho;
  do {
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0, 1) > rho || x > b);
  return x;
}

double hnrs_ab(double a, double b) {
  double x = a - 1.0;
  while (x < a || x > b) {
    x = R::rnorm(0, 1);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling (a, b) */
double urs_ab(const double & a, const double & b) {
  const double phi_a = R::dnorm4(a, 0.0, 1.0, false);
  double x = 0.0;
  /* Upper bound of normal density on [a, a+1] */
  const double ub = (a < 0 && b > 0) ? M_1_SQRT_2PI : phi_a;
  do {
    x = R::runif(a, b);
  } while (R::runif(0, 1) * ub > R::dnorm4(x, 0, 1, 0));
  return x;
}

double rhalf_norm(const double & mean, const double & sd) {
  const double alpha = (- mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha);
  } else {
    return mean + sd * ers_a_inf(alpha);
  }
}

// [[Rcpp::export]]
double ra1_norm(const double & a, const double & mean, const double & sd) {
  const double alpha = (a - mean) / sd;
  const double beta = (a + 1 - mean) / sd;
  // const double Phi_a = R::pnorm(lower, 0, 1, true, false);
  // const double Phi_b = R::pnorm(upper, 0, 1, true, false);
  // double u = R::runif(Phi_a, Phi_b);
  // return mean+sd*R::qnorm(u, 0, 1, true, false);
  const double phi_a = R::dnorm(alpha, 0, 1, false);
  const double phi_b = R::dnorm(beta, 0, 1, false);
  if (alpha > 0) {      /* 3 */
  if (phi_a / phi_b <= t2) { /* 3 (a) */
  return mean + sd * urs_ab(alpha, beta);
  } else {
    if (alpha < t3) { /* 3 (b) */
  return mean + sd * hnrs_ab(alpha, beta);
    } else { /* 3 (c) */
  return mean + sd * ers_ab(alpha, beta);
    }
  }
  } else {                     /* 3s */
  if (phi_b / phi_a <= t2) { /* 3s (a) */
  return mean - sd * urs_ab(-beta, -alpha);
  } else {
    if (beta > -t3) { /* 3s (b) */
  return mean - sd * hnrs_ab(-beta, -alpha);
    } else { /* 3s (c) */
  return mean - sd * ers_ab(-beta, -alpha);
    }
  }
  }
}

//truncated x>0
double moment1(const double & mu, const double & sigma){
  return mu + sigma*exp(R::dnorm(-mu/sigma, 0, 1, 1)-R::pnorm(-mu/sigma, 0, 1, 0, 1));
}

double moment2(const double & mu, const double & sigma) {
  return pow(mu,2) + pow(sigma,2)+mu*sigma*exp(R::dnorm(-mu/sigma, 0, 1, 1)-R::pnorm(-mu/sigma,0, 1, 0, 1));
}

/*
 arma::vec EV(const arma::vec & mu, const double & sigma){
 arma::vec V = mu;
 for(int i=0;i<mu.n_rows;i++){
 V(i) = moment1(mu(i), sigma);
 }
 return V;
 }
 
 arma::vec EV2(const arma::vec & mu, const double & sigma){
 arma::vec V2 = mu;
 for(int i=0;i<mu.n_rows;i++){
 V2(i) = moment2(mu(i), sigma);
 }
 return V2;
 }
 */
