
void up_lambda(double & obs_prec, 
               double & bhat, const double & ahat,
               const double & R, const double & b){
  bhat = 0.5*R + b;
  obs_prec = ahat/bhat;
}