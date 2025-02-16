#include "ELBO.h"
#include "KLgamma.h"

double calc_elbo(const double & R, const double & obs_prec,
                 const double & ahat, const double & bhat,
                 const double & a, const double & b){
  return -0.5*obs_prec*R + 0.5*(R::digamma(ahat)-log(bhat)) - kl2gamma(ahat, bhat, a, b);
}
