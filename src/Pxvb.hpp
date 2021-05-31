#ifndef Pxvb_hpp
#define Pxvb_hpp

#include "RcppArmadillo.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

List vbfunM1(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
	arma::vec mu, double& sgga2, double beta0, arma::mat R,
	const int&constr, const double &epsStopLogLik, const int& maxIter);

List PXvbfunM1(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
	arma::vec mu, double& sgga2, double beta0, arma::mat R,
	const int&constr, const double &epsStopLogLik, const int& maxIter);

List vbfunM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
	arma::vec mu, arma::vec muA, double& sgga2, double& sgal2, double beta0,
	arma::mat R, const int&constr, const double &epsStopLogLik, const int& maxIter);

List PXvbfunM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
	arma::vec mu, arma::vec muA, double& sgga2, double& sgal2, double beta0, arma::mat R,
	const int&constr, const double &epsStopLogLik, const int& maxIter);

List rcpparma_SCLMM2_VB(vec &hatmur, vec &hatsr, vec &hatmu2r, vec &hats2r, mat &R,
                        int max_iter, double tol,  bool fix_alphag);

List VBMgM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
            arma::vec mu, double& sgga2, double sgal2, double beta0, arma::mat R,
            const int&constr, const double &epsStopLogLik, const int& maxIter);

List MgMgib2(arma::vec gammah, arma::vec Gammah, arma::vec sg2, arma::vec sG2,
             arma::vec mu, double& sgga2, double sgal2, double beta0,
             arma::mat R, int IterMax, int agm, double bgm);

#endif /* Pxvb_hpp */