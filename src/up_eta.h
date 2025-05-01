#include "RcppArmadillo.h"

void up_eta_w_woi(arma::mat & num_w,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::mat & Z);

void up_eta_z_woi(arma::mat & num_z,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::mat & W);

double residuals(const arma::vec y, const arma::vec y2, 
                 const arma::uvec & rowi, const arma::uvec & coli,
                 const arma::mat Z, const arma::mat W);
