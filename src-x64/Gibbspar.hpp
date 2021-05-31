#ifndef Gibbspar_hpp
#define Gibbspar_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;

Rcpp::List gibbsres1(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                     arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                     arma::vec &gamma,  double &beta0, double &sgga2, 
                     int coreNum, int IterMax);
Rcpp::List gibbsres2(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                     arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                     arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                     int coreNum, int IterMax);
  
void Varres(field<mat> &F4H2a, arma::mat R, arma::umat block_inf, uword nblocks, 
            arma::vec bh2, arma::vec se2, arma::ivec N2);
arma::vec Herit_iMax(arma::mat R, arma::umat block_inf, uword nblocks, 
                     arma::vec bh2, arma::vec se2, arma::ivec N2, arma::mat ALres);
double heritability(field<mat> F4H2a, arma::umat block_inf, uword nblocks,  arma::vec  alpha);

#endif /* Gibbspar_hpp */