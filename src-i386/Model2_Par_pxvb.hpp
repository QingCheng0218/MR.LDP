#ifndef Model2_Par_pxvb_hpp
#define Model2_Par_pxvb_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;


List Para_PXVb2(arma::mat R, arma::umat block_inf, uword nblocks,
                arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                arma::vec& mu, arma::vec& muA, double& beta0, double& sgga2, double& sgal2, int coreNum,
                const int& constr, const double& epsStopLogLik, const int& maxIter);

#endif /* Model2_Par_pxvb_hpp */