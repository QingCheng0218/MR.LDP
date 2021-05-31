#ifndef Par_pxvb_hpp
#define Par_pxvb_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;

List PXvbfunM1_par(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                   const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1,const  arma::vec& se2,
                   arma::vec& mu, double& sgga2, double& beta0, 
                   const int& constr, const double& epsStopLogLik, const int& maxIter, int coreNum);
List PXvbfunM2_par(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                     const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1, const arma::vec& se2,
                     arma::vec& mu, arma::vec& muA, double& sgga2, double& sgal2, double& beta0, 
                     const int& constr, const double& epsStopLogLik, const int& maxIter, int coreNum);
#endif /* Par_pxvb_hpp */