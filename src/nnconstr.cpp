#include <RcppArmadillo.h>
#include <memory>
#include "truncnorm.h"
#include "KLnorm.h"
#include "nnconstr.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void set_constr(std::unique_ptr<nnconstr> & g,
                const std::string & constr_type){
  if(constr_type == "AN"){
    g.reset(new AN);
  }else if(constr_type == "NN"){
    g.reset(new NN);
  }else{
    Rcpp::stop("This constr_type is not implemented\n");
  }
}
