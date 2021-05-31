#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>
#include "time.h"
#include <iostream>
#include <thread>
#include <vector>
#include "pdsoft.hpp"
#include "Pxvb.hpp"
#include "Pxvb_block.hpp"
#include "Gibbspar.hpp"
#include "ReadGeneFile.hpp"
#include "CalCorr.hpp"
#include "data_loader.hpp"
#include "Model1_Par_pxvb.hpp"
#include "Model2_Par_pxvb.hpp"
#include "Par_pxvb.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
Rcpp::List MRLDP_SimVb(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                 arma::vec gamma, arma::vec alpha, double beta0, double& sgga2, double& sgal2, arma::mat R,
                 const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){

  double tstat, Iteration, diff;
  vec loglik;

  if(model==1){
    List VbResult1 = vbfunM1(bh1, bh2, se1, se2, gamma, sgga2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = VbResult1["beta0"];
    tstat = VbResult1["tstat"];
    sgga2 = VbResult1["sgga2"];
    diff = VbResult1["diff"];
    Iteration = VbResult1["Iteration"];
    loglik = as<vec>(VbResult1["loglik"]);
  }else if(model==2){
    List VbResult2 = vbfunM2(bh1, bh2, se1, se2, gamma, alpha, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = VbResult2["beta0"];
    tstat = VbResult2["tstat"];
    sgga2 = VbResult2["sgga2"];
    sgal2 = VbResult2["sgal2"];
    diff = VbResult2["diff"];
    Iteration = VbResult2["Iteration"];
    loglik = as<vec>(VbResult2["loglik"]);
  }

  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}

// [[Rcpp::export]]
Rcpp::List MRLDP_SimPXvb(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                 arma::vec gamma, arma::vec alpha, double beta0, double& sgga2, double& sgal2, arma::mat R,
                 const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){

  double tstat, Iteration, diff;
  vec loglik;
  if(model==1){
    List PXvbResult1 = PXvbfunM1(bh1, bh2, se1, se2, gamma, sgga2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = PXvbResult1["beta0"];
    tstat = PXvbResult1["tstat"];
    sgga2 = PXvbResult1["sgga2"];
    diff = PXvbResult1["diff"];
    Iteration = PXvbResult1["Iteration"];
    loglik = as<vec>(PXvbResult1["loglik"]);
  }else if(model==2){
    List PXvbResult2 = PXvbfunM2(bh1, bh2, se1, se2, gamma, alpha, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = PXvbResult2["beta0"];
    tstat = PXvbResult2["tstat"];
    sgga2 = PXvbResult2["sgga2"];
    sgal2 = PXvbResult2["sgal2"];
    diff = PXvbResult2["diff"];
    Iteration = PXvbResult2["Iteration"];
    loglik = as<vec>(PXvbResult2["loglik"]);
  }
  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_SimVBMg(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
            arma::vec gamma, double sgga2, double sgal2, double beta0, arma::mat R,
            const int&constr, const double &epsStopLogLik, const int& maxIter){
  double tstat, Iteration, diff;
  vec loglik;
  List VBMgresult = VBMgM2(bh1, bh2, se1, se2, gamma, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
  beta0 = VBMgresult["beta0"];
  tstat = VBMgresult["tstat"];
  sgga2 = VBMgresult["sgga2"];
  diff = VBMgresult["diff"];
  Iteration = VBMgresult["Iteration"];
  loglik = as<vec>(VBMgresult["loglik"]);

  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_SimMggibbs(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                      arma::vec gamma, double& sgga2, double sgal2, double beta0, arma::mat R,
                      int maxIter, int agm, double bgm){
  vec BETAres, SGGMres;
  List Mggibbsresult = MgMgib2(bh1, bh2, se1, se2, gamma, sgga2, sgal2, beta0, R,
                            maxIter, agm, bgm);
  BETAres = as<vec>(Mggibbsresult["BETAres"]);
  SGGMres = as<vec>(Mggibbsresult["SGGMres"]);
  
  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_SimGibbspar(arma::mat R, arma::umat block_inf, uword nblocks,
                       arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                       arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                       int coreNum, int IterMax, arma::ivec N2, int model){
  
  field<mat> F4Rblock(nblocks, 1);
  for(int i=0; i< (int)(nblocks); i++){
    F4Rblock(i, 0) = R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1));
  }
  
  vec BETAres, SGGMres, SGALres;
  mat GAres, ALres;
  vec pve_a = zeros(IterMax, 1);
  
  if(model==1){
    // simulation result for model 1.
    List GibbsResult = gibbsres1(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                                 gamma,  beta0, sgga2,  coreNum, IterMax);
    
    BETAres = as<vec>(GibbsResult["BETAres"]);
    SGGMres = as<vec>(GibbsResult["SGGMres"]);
    GAres = as<mat>(GibbsResult["GAres"]);
  }else if(model==2){
    // simulation result for model 2.
    List GibbsResult = gibbsres2(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                                 gamma, alpha, beta0, sgga2, sgal2, coreNum, IterMax);
    
    BETAres = as<vec>(GibbsResult["BETAres"]);
    SGGMres = as<vec>(GibbsResult["SGGMres"]);
    SGALres = as<vec>(GibbsResult["SGALres"]);
    GAres = as<mat>(GibbsResult["GAres"]);
    ALres = as<mat>(GibbsResult["ALres"]);
    
    // ------------------------------------------------------------------------------------
    // caculate the PVE(proportion of phenotypic variation explained by direct effect(alpha))
    
    pve_a = Herit_iMax(R, block_inf, nblocks, bh2, se2, N2, ALres);
  }
  
  
  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres,
    Rcpp::Named("GAres") = GAres,
    Rcpp::Named("SGALres") = SGALres,
    Rcpp::Named("ALres") = ALres,
    Rcpp::Named("pve_a") = pve_a
  );
  
  return output;
}

//--------------------------------------------------------------------------------//

// [[Rcpp::export]]
Rcpp::List MRLDP_RealVb(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                  std::string stringname3,
                  arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                  arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                  int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];

  arma::mat R = zeros(p, p);
  for(int i=0; i< (int)(nblocks); i++){
    R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
  }

  double tstat, Iteration, diff;
  vec loglik;
  if(model==1){
    List VbResult1 = vbfunM1(bh1, bh2, se1, se2, gamma, sgga2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = VbResult1["beta0"];
    tstat = VbResult1["tstat"];
    sgga2 = VbResult1["sgga2"];
    diff = VbResult1["diff"];
    Iteration = VbResult1["Iteration"];
    loglik = as<vec>(VbResult1["loglik"]);
  }else if(model==2){
    List VbResult2 = vbfunM2(bh1, bh2, se1, se2, gamma, alpha, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = VbResult2["beta0"];
    tstat = VbResult2["tstat"];
    sgga2 = VbResult2["sgga2"];
    sgal2 = VbResult2["sgal2"];
    diff = VbResult2["diff"];
    Iteration = VbResult2["Iteration"];
    loglik = as<vec>(VbResult2["loglik"]);
  }

  List output = List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_RealPXvb(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                  std::string stringname3,
                  arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                  arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                  int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];

  arma::mat R = zeros(p, p);
  for(int i=0; i< (int)(nblocks); i++){
    R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
  }

  double tstat, Iteration, diff;
  vec loglik;
  if(model==1){
    List PXvbResult1 = PXvbfunM1(bh1, bh2, se1, se2, gamma, sgga2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = PXvbResult1["beta0"];
    tstat = PXvbResult1["tstat"];
    sgga2 = PXvbResult1["sgga2"];
    diff = PXvbResult1["diff"];
    Iteration = PXvbResult1["Iteration"];
    loglik = as<vec>(PXvbResult1["loglik"]);
  }else if(model==2){
    List PXvbResult2 = PXvbfunM2(bh1, bh2, se1, se2, gamma, alpha, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
    beta0 = PXvbResult2["beta0"];
    tstat = PXvbResult2["tstat"];
    sgga2 = PXvbResult2["sgga2"];
    sgal2 = PXvbResult2["sgal2"];
    diff = PXvbResult2["diff"];
    Iteration = PXvbResult2["Iteration"];
    loglik = as<vec>(PXvbResult2["loglik"]);
  }

  List output = List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_RealPXvb_block(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                          std::string stringname3,
                          arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                          arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                          int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  
  double tstat, Iteration, diff;
  vec loglik;
  if(model==1){
    List PXvbResult1_block = PXvbfunM1_block(F4Rblock, block_inf,  nblocks,
                                             bh1, bh2, se1, se2, gamma, sgga2, beta0, 
                                             constr, epsStopLogLik, maxIter);
    beta0 = PXvbResult1_block["beta0"];
    tstat = PXvbResult1_block["tstat"];
    sgga2 = PXvbResult1_block["sgga2"];
    diff = PXvbResult1_block["diff"];
    Iteration = PXvbResult1_block["Iteration"];
    loglik = as<vec>(PXvbResult1_block["loglik"]);
  }else if(model==2){
    List PXvbResult2_block = PXvbfunM2_block(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                                             gamma, alpha, sgga2, sgal2, beta0, constr, epsStopLogLik, maxIter); 
                                             
    beta0 = PXvbResult2_block["beta0"];
    tstat = PXvbResult2_block["tstat"];
    sgga2 = PXvbResult2_block["sgga2"];
    sgal2 = PXvbResult2_block["sgal2"];
    diff = PXvbResult2_block["diff"];
    Iteration = PXvbResult2_block["Iteration"];
    loglik = as<vec>(PXvbResult2_block["loglik"]);
  }
  
  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}





// [[Rcpp::export]]
Rcpp::List MRLDP_RealPXvb_par(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                              std::string stringname3, arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                              arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                              int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  
  double tstat, Iteration, diff;
  vec loglik;
  
  if(model==1){
    List ParapxvbResult1 = PXvbfunM1_par(F4Rblock, block_inf, nblocks,
                                         bh1, bh2, se1, se2, gamma, sgga2, beta0, 
                                         constr, epsStopLogLik, maxIter, coreNum);
    beta0 = ParapxvbResult1["beta0"];
    tstat = ParapxvbResult1["tstat"];
    sgga2 = ParapxvbResult1["sgga2"];
    diff = ParapxvbResult1["diff"];
    Iteration = ParapxvbResult1["Iteration"];
    loglik = as<vec>(ParapxvbResult1["loglik"]);
  }else if(model==2){
    List PXvbResult2_block = PXvbfunM2_par(F4Rblock,  block_inf, nblocks, bh1, bh2, se1, se2,
                                           gamma, alpha, sgga2, sgal2, beta0,
                                           constr, epsStopLogLik, maxIter,coreNum);
    
    beta0 = PXvbResult2_block["beta0"];
    tstat = PXvbResult2_block["tstat"];
    sgga2 = PXvbResult2_block["sgga2"];
    sgal2 = PXvbResult2_block["sgal2"];
    diff = PXvbResult2_block["diff"];
    Iteration = PXvbResult2_block["Iteration"];
    loglik = as<vec>(PXvbResult2_block["loglik"]);
  }
  
  
  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_RealParPXvb(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                    std::string stringname3,
                    arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                    arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                    int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  
  arma::mat R = zeros(p, p);
  for(int i=0; i< (int)(nblocks); i++){
    R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
  }
  
  double tstat, Iteration, diff;
  vec loglik;
  if(model==1){
    List ParapxvbResult1 = Para_PXVb1(R, block_inf, nblocks, bh1, bh2, se1, se2, gamma,
                                           beta0, sgga2, coreNum, constr, epsStopLogLik, maxIter);
    beta0 = ParapxvbResult1["beta0"];
    tstat = ParapxvbResult1["tstat"];
    sgga2 = ParapxvbResult1["sgga2"];
    diff = ParapxvbResult1["diff"];
    Iteration = ParapxvbResult1["Iteration"];
    loglik = as<vec>(ParapxvbResult1["loglik"]);
  }else if(model==2){
    List ParapxvbResult2 = Para_PXVb2(R, block_inf, nblocks, bh1, bh2, se1, se2, gamma, alpha,
                                      beta0, sgga2, sgal2, coreNum, constr, epsStopLogLik, maxIter);
    beta0 = ParapxvbResult2["beta0"];
    tstat = ParapxvbResult2["tstat"];
    sgga2 = ParapxvbResult2["sgga2"];
    sgal2 = ParapxvbResult2["sgal2"];
    diff = ParapxvbResult2["diff"];
    Iteration = ParapxvbResult2["Iteration"];
    loglik = as<vec>(ParapxvbResult2["loglik"]);
  }
  
  List output = List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_RealVBMg(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel,  std::string block_file,
                    std::string stringname3,
                    arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                    arma::vec &gamma, double &beta0, double &sgga2, double sgal2,
                    int coreNum, double lam, const int&constr, const double &epsStopLogLik, const int& maxIter, const int& model){
  int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  
  arma::mat R = zeros(p, p);
  for(int i=0; i< (int)(nblocks); i++){
    R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
  }
  
  double tstat, Iteration, diff;
  vec loglik;
  List VBMgresult = VBMgM2(bh1, bh2, se1, se2, gamma, sgga2, sgal2, beta0, R, constr, epsStopLogLik, maxIter);
  beta0 = VBMgresult["beta0"];
  tstat = VBMgresult["tstat"];
  sgga2 = VBMgresult["sgga2"];
  sgal2 = VBMgresult["sgal2"];
  diff = VBMgresult["diff"];
  Iteration = VBMgresult["Iteration"];
  loglik = as<vec>(VBMgresult["loglik"]);
  
  List output = List::create(
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("sgga2") = sgga2,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("Iteration") = Iteration
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List MRLDP_RealMggibbs(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                    std::string stringname3,
                    arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                    arma::vec &gamma, double &beta0, double &sgga2, double sgal2,
                    int coreNum, double lam, int maxIter, int agm, double bgm){
  int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  
  arma::mat R = zeros(p, p);
  for(int i=0; i< (int)(nblocks); i++){
    R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
  }
  
  vec BETAres, SGGMres;
  List Mggibbsresult = MgMgib2(bh1, bh2, se1, se2, gamma, sgga2, sgal2, beta0, R,
                               maxIter, agm, bgm);
  BETAres = as<vec>(Mggibbsresult["BETAres"]);
  SGGMres = as<vec>(Mggibbsresult["SGGMres"]);
  
  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres
  );
  return output;
}



// [[Rcpp::export]]
Rcpp::List MRLDP_RealGibbspar(arma::ivec bp, arma::ivec chr, arma::uvec avbIndex, arma::uvec &idx4panel, std::string block_file,
                        std::string stringname3,
                        arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                        arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                        int coreNum, double lam, int IterMax, arma::ivec N2, int model){

  // int p = bp.size();
  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];

  vec BETAres;
  vec SGGMres;
  vec SGALres;
  mat GAres;
  mat ALres;
  vec pve_a = zeros(IterMax, 1);
  if(model==1){
    List GibbsResult1 = gibbsres1(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                                  gamma, beta0, sgga2, coreNum, IterMax);

    BETAres = as<vec>(GibbsResult1["BETAres"]);
    SGGMres = as<vec>(GibbsResult1["SGGMres"]);
    GAres = as<vec>(GibbsResult1["GAres"]);

  }else if(model==2){
    List GibbsResult2 = gibbsres2(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                                  gamma, alpha, beta0, sgga2, sgal2, coreNum, IterMax);

    BETAres = as<vec>(GibbsResult2["BETAres"]);
    SGGMres = as<vec>(GibbsResult2["SGGMres"]);
    SGALres = as<vec>(GibbsResult2["SGALres"]);
    GAres = as<mat>(GibbsResult2["GAres"]);
    ALres = as<mat>(GibbsResult2["ALres"]);

    // ------------------------------------------------------------------------------------
    // caculate the PVE(proportion of phenotypic variation explained by direct effect(alpha))
    // arma::mat R = zeros(p, p);
    // for(int i=0; i< (int)(nblocks); i++){
    //   R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1)) = F4Rblock(i, 0);
    // }
    // 
    // pve_a = Herit_iMax(R, block_inf, nblocks, bh2, se2, N2, ALres);
  }


  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres,
    Rcpp::Named("GAres") = GAres,
    Rcpp::Named("SGALres") = SGALres,
    Rcpp::Named("ALres") = ALres,
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("nblocks") = nblocks
    // Rcpp::Named("pve_a") = pve_a
  );
  return output;
}





