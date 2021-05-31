#ifndef Model1_Par_pxvb_hpp
#define Model1_Par_pxvb_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;

List Para_PXVb1(arma::mat R, arma::umat block_inf, uword nblocks,
                      arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                      arma::vec& mu,  double& beta0, double& sgga2, int coreNum,
                      const int& constr, const double& epsStopLogLik, const int& maxIter);

#endif /* Model1_Par_pxvb_hpp */