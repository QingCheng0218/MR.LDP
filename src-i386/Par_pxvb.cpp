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
// #include "CalCorr_match_exposure.hpp"
// #include "data_loader_match_exposure.hpp"
#include "Model1_Par_pxvb.hpp"
#include "Model2_Par_pxvb.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]


class paraBlock_GamAlp2{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  
  vec betasig;
  vec betamean;
  vec Mu2;
  vec MuA2;
  vec xi_nu;
  vec xi_de;
  double beta0;
  double sgga2;
  double sgal2;
  double xi;
  int constr;
  
  
  field<vec> F4se1, F4se2, F4mu, F4muA, F4v2, F4v2A;
  field<vec> F4sg2, F4sG2, F4GinvsG2, F4ginvsg2;
  field<vec> F4RinsGmu, F4Rinsgmu, F4RinsGmuA;
  field<vec> F4RinsginsG;
  field<mat> F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg;
  
  paraBlock_GamAlp2(const uword &nblocks, field<vec> &F4se1, field<vec> &F4se2, field<vec> &F4ginvsg2, field<vec> &F4GinvsG2,  
                    field<vec> &F4v2, field<vec> &F4v2A,
                   field<mat> &F4Rins, field<mat> &F4Rins2, const field<mat> &F4Rblock, field<mat> &F4insGRinsG, field<mat> &F4insgRinsg,
                   field<vec> &F4mu, field<vec> &F4muA, field<vec> &F4Rinsgmu, field<vec> &F4RinsGmu, field<vec> &F4RinsGmuA,
                   arma::vec &betasig, arma::vec &betamean, arma::vec &Mu2, arma::vec &MuA2, arma::vec &xi_nu, arma::vec &xi_de,
                   const int& constr, double &xi, double &beta0, double &sgga2, double &sgal2){

    this -> nblocks = nblocks;
    this -> F4se1 = F4se1;
    this -> F4se2 = F4se2;
    this -> F4ginvsg2 = F4ginvsg2;
    this -> F4GinvsG2 = F4GinvsG2;
    this -> F4v2 = F4v2;
    this -> F4v2A = F4v2A;
    this -> F4Rins = F4Rins;
    this -> F4Rins2 = F4Rins2;
    this -> F4Rblock = F4Rblock;
    this -> F4insGRinsG = F4insGRinsG;
    this -> F4insgRinsg = F4insgRinsg;
    this -> F4mu = F4mu;
    this -> F4muA = F4muA;
    this -> F4Rinsgmu = F4Rinsgmu;
    this -> F4RinsGmu = F4RinsGmu;
    this -> F4RinsGmuA = F4RinsGmuA;
    this -> betasig = betasig;
    this -> betamean = betamean;
    this -> Mu2 = Mu2;
    this -> MuA2 = MuA2;
    this -> constr = constr;
    this -> xi = xi;
    this -> xi_nu = xi_nu;
    this -> xi_de = xi_de;
    this -> beta0 = beta0;
    this -> sgga2 = sgga2;
    this -> sgal2 = sgal2;
  }
  
  int  next_GamAlp2();
  void loop_by_block_gibbs_GamAlp2(int i);
  void update_by_thread_GamAlp2(int thread_id);
  
};


//
void paraBlock_GamAlp2::loop_by_block_gibbs_GamAlp2(int i){
  double xi2 = xi * xi;
  double beta02 = pow(beta0, 2);

  vec se1 = F4se1(i, 0);
  vec se2 = F4se2(i, 0);
  mat Rins = F4Rins(i, 0);
  mat Rins2 = F4Rins2(i, 0);
  mat Rblock = F4Rblock(i, 0);
  vec ginvsg2 = F4ginvsg2(i, 0);
  vec GinvsG2 = F4GinvsG2(i, 0);
  mat insgRinsg = F4insgRinsg(i,  0);
  mat insGRinsG = F4insGRinsG(i,  0);
  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);

  vec invse1 = 1. / se1;
  vec invse2 = 1. / se2;

  vec mu = F4mu(i, 0);
  vec muA = F4muA(i, 0);
  vec Rinsgmu = F4Rinsgmu(i, 0);
  vec RinsGmu = F4RinsGmu(i, 0);
  vec RinsGmuA = F4RinsGmuA(i, 0);
  vec Rdiag = diagvec(F4Rblock(i, 0));

  int p_block = F4Rblock(i, 0).n_rows;
  vec v2 = 1. / (beta02*diaginsGRinsG + xi2*diaginsgRinsg + 1. / sgga2);
  F4v2(i, 0) = v2;

  for (int j = 0; j < p_block; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2, RinSmuAjj;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];


    RinSmujj = Rinsgmu[j] - Rdiag[j]*mu[j] * invse1[j];
    RinSmujj2 = RinsGmu[j] - Rdiag[j]*mu[j] * invse2[j];
    RinSmuAjj = RinsGmuA[j];
  
    mu[j] =(beta0*GinvsG2[j] - beta02 * invse2[j]*RinSmujj2 - beta0 * invse2[j]*RinSmuAjj + xi*ginvsg2[j] - xi2 * invse1[j]*RinSmujj)*v2[j];
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];

  }

  vec v2A = 1. / (invse2%invse2 + 1. / sgal2);
  F4v2A(i, 0) = v2A;

  for (int k = 0; k < p_block; k = k + 1){
    vec tmp3;
    double RinSmuAkk;
    tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
    RinSmuAkk = RinsGmuA[k] - Rdiag[k]*muA[k] * invse2[k];
    muA[k] = (GinvsG2[k] - beta0 * invse2[k]*RinsGmu[k] - invse2[k]*RinSmuAkk)*v2A[k];
    RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
  }
  // -----------------------------------------------------------------------
  F4mu(i, 0) = mu;
  F4muA(i, 0) = muA;
  F4Rinsgmu(i, 0) = Rinsgmu;
  F4RinsGmu(i, 0) = RinsGmu;
  F4RinsGmuA(i, 0) = RinsGmuA;
  // -----------------------------------------------------------------------
  // M step
  if (constr != 1){
    betasig[i] = as_scalar(mu.t()*insGRinsG*mu + v2.t()*diaginsGRinsG);
    betamean[i] = as_scalar(GinvsG2.t()*mu - muA.t()*insGRinsG*mu);
  }
  // -----------------------------------------------------------------------
  // for sgga2, sgal2 and xi iteration
  // cout << "-----------------------------check erro 3.4 " << endl;
  Mu2[i] = sum(mu%mu) + sum(v2);
  MuA2[i] = sum(muA%muA) + sum(v2A);
  xi_nu[i] = as_scalar(ginvsg2.t()*mu);
  xi_de[i] = as_scalar(mu.t()*insgRinsg*mu + v2.t()*diaginsgRinsg);


  se1.reset();
  se2.reset();
  invse1.reset();
  invse2.reset();

  mu.reset();
  muA.reset();
  Rinsgmu.reset();
  RinsGmu.reset();
  RinsGmuA.reset();

  Rins.reset();
  Rins2.reset();
  Rblock.reset();
  Rdiag.reset();

  ginvsg2.reset();
  GinvsG2.reset();
  insgRinsg.reset();
  insGRinsG.reset();
  diaginsGRinsG.reset();
  diaginsgRinsg.reset();

}


std::mutex _mtx5;
int paraBlock_GamAlp2::next_GamAlp2(){
  std::lock_guard<std::mutex> lockGuard(_mtx5);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}




void paraBlock_GamAlp2::update_by_thread_GamAlp2(int thread_id){
  while(true){
    int idx = next_GamAlp2();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_GamAlp2(idx);
  }
}

// [[Rcpp::export]]
List PXvbfunM2_par(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                     const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1, const arma::vec& se2,
                     arma::vec& mu, arma::vec& muA, double& sgga2, double& sgal2, double& beta0, 
                     const int& constr, const double& epsStopLogLik, const int& maxIter, int coreNum){
  
  int p = bh1.n_elem;
  ivec NB = zeros<ivec>(nblocks, 1);
  
  field<vec> F4se1(nblocks, 1), F4se2(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1), F4v2(nblocks, 1), F4v2A(nblocks, 1);
  field<vec> F4sg2(nblocks, 1), F4sG2(nblocks, 1), F4GinvsG2(nblocks, 1), F4ginvsg2(nblocks, 1), F4diaginsGRinsG(nblocks, 1);
  field<vec> F4diaginsgRinsg(nblocks, 1), F4RinsGmu(nblocks, 1), F4Rinsgmu(nblocks, 1), F4RinsGmuA(nblocks, 1);
  field<mat> F4insGRinsG(nblocks, 1), F4insgRinsg(nblocks, 1), F4Rins(nblocks, 1), F4Rins2(nblocks, 1);
 
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    vec se1_block = se1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec se2_block = se2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    vec mu_block = mu.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec muA_block = muA.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh1_blcok = bh1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh2_blcok = bh2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    mat R_block = F4Rblock(nn, 0);
    NB[nn] = se1_block.n_elem;
    F4mu(nn, 0) = mu_block;
    F4muA(nn, 0) = muA_block;
    F4se1(nn, 0) = se1_block;
    F4se2(nn, 0) = se2_block;
    
    F4sg2(nn, 0) = sg2_block;
    F4sG2(nn, 0) = sG2_block;
    
    F4GinvsG2(nn, 0) = bh2_blcok / sG2_block;
    F4ginvsg2(nn, 0) = bh1_blcok / sg2_block;
    // cout << "!!!!" << F4GinvsG2(nn, 0) << endl;
    F4insGRinsG(nn, 0) = diagmat(1. / se2_block)*R_block*diagmat(1. / se2_block);
    F4insgRinsg(nn, 0) = diagmat(1. / se1_block)*R_block*diagmat(1. / se1_block);
    
    F4diaginsGRinsG(nn, 0) = diagvec(F4insGRinsG(nn, 0));
    F4diaginsgRinsg(nn, 0) = diagvec(F4insgRinsg(nn, 0));
    
    F4Rins(nn, 0) = R_block*diagmat(1 / se1_block);
    F4Rins2(nn, 0) = R_block*diagmat(1 / se2_block);
    
    F4RinsGmu(nn, 0) = R_block*diagmat(1 / se2_block)*mu_block;
    F4Rinsgmu(nn, 0) = R_block*diagmat(1 / se1_block)*mu_block;
    F4RinsGmuA(nn, 0) = R_block*diagmat(1. / se2_block)*muA_block;
    
    
  }
  
  double xi = 1;
  double xi2 = xi*xi;
  double beta02 = beta0*beta0;
  
  vec v2,v2A;
  
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
 
  for(int iter = 2; iter <= maxIter; iter ++){
    // E step
    vec betamean = zeros(nblocks, 1); // mean for each block.
    vec betasig = zeros(nblocks, 1); // sigma for each block.
    vec Mu2 = zeros(nblocks, 1); // sgga2 for each block.
    vec MuA2 = zeros(nblocks, 1); // sgal2 for each block.
    vec xi_nu = zeros(nblocks, 1); // xi for each block.
    vec xi_de = zeros(nblocks, 1); // xi for each block.
    // parallel for gamma
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    
    
    paraBlock_GamAlp2 parobj_alp2(nblocks, F4se1, F4se2, F4ginvsg2, F4GinvsG2, F4v2, F4v2A, 
                                    F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg,
                                    F4mu, F4muA, F4Rinsgmu, F4RinsGmu, F4RinsGmuA,
                                    betasig, betamean, Mu2, MuA2, xi_nu, xi_de,
                                    constr, xi, beta0, sgga2, sgal2);

    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_GamAlp2::update_by_thread_GamAlp2, &parobj_alp2, i_thread);
    }
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    
    betasig = parobj_alp2.betasig;
    betamean = parobj_alp2.betamean;
    Mu2 = parobj_alp2.Mu2;
    MuA2 = parobj_alp2.MuA2;
    xi_nu = parobj_alp2.xi_nu;
    xi_de = parobj_alp2.xi_de;
    F4v2 = parobj_alp2.F4v2;
    F4v2A = parobj_alp2.F4v2A;
    F4RinsGmu = parobj_alp2.F4RinsGmu;
    F4Rinsgmu = parobj_alp2.F4Rinsgmu;
    F4mu = parobj_alp2.F4mu;
    F4muA = parobj_alp2.F4muA;
    F4RinsGmuA = parobj_alp2.F4RinsGmuA;

    // -----------------------------------------------------------------------
    // M step 
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      double sig2b;
      sig2b = 1. / sum(betasig);
      beta0 = sum(betamean)*sig2b;
    }
    

    // update sgga2
    sgga2 = sum(Mu2) / p;
    // update sgal2
    sgal2 = sum(MuA2) / p;
    // update xi.
    xi = sum(xi_nu)/ sum(xi_de);
    // xi = 1;

    // Reduction step.
   
    xi2 = xi*xi;
    mu = mu*xi;
    beta0 = beta0 / xi;
    
    beta02 = beta0*beta0;
    sgga2 = sgga2*xi2;


    // lower bound to check convergence.
    vec low_b = zeros(nblocks, 1);
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){

      F4v2(ll, 0) = F4v2(ll, 0)*xi2;
      F4mu(ll, 0) = F4mu(ll, 0)*xi;
      F4RinsGmu(ll, 0) = F4RinsGmu(ll, 0)*xi;
      F4Rinsgmu(ll, 0) = F4Rinsgmu(ll, 0)*xi;

      // cout << "-----------------------------check erro 6.1**** " << endl;
      double term1, term2;
      int p_b = NB[ll];
      term1 = as_scalar((beta0*F4GinvsG2(ll, 0) + F4ginvsg2(ll, 0)).t()*F4mu(ll, 0)) -
        0.5*as_scalar(F4mu(ll, 0).t()*(beta02*F4insGRinsG(ll, 0) + F4insgRinsg(ll, 0) + 1./sgga2*diagmat(ones(p_b,1)))*F4mu(ll, 0)) -
        0.5*as_scalar(F4v2(ll, 0).t()*(beta02*F4diaginsGRinsG(ll, 0)  + F4diaginsgRinsg(ll, 0) + 1/sgga2)) + 0.5*sum(log(F4v2(ll, 0)));


      term2 = as_scalar(- beta0*F4muA(ll, 0).t()*F4insGRinsG(ll, 0).t()*F4mu(ll, 0) + F4GinvsG2(ll, 0).t()*F4muA(ll, 0) -
        0.5*F4muA(ll, 0).t()*(F4insGRinsG(ll, 0) + 1/sgal2 * diagmat(ones(p_b)))*F4muA(ll, 0)) -
        0.5* as_scalar(F4v2A(ll, 0).t()*(F4diaginsGRinsG(ll, 0) + 1/sgal2)) + 0.5*sum(log(F4v2A(ll, 0)));
      low_b[ll] = term1 + term2;
    }

    
    xi = 1;
    xi2 = 1;
    
    double low = sum(low_b) - 0.5*p*log(sgga2) - 0.5*p*log(sgal2);
    

    loglik(iter-1) = low;

    if(loglik(iter -1) - loglik(iter -2) < -1e-7){
      perror("The likelihood failed to increase!");
    }
    //
    Iteration = iter;
    if(iter > 2){
      if(abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik){
        break;
      }
    }
    
  }
  
  vec loglik_out;
  int to = Iteration - 1;
  loglik_out = loglik.subvec(0, to);
  double diff = loglik_out(to) - loglik_out(to - 1);

  // for loglikelihood ratio test.
  vec tstat_b = zeros(nblocks, 1);
  double invsgga2 = 1./sgga2;
  double invsgal2 = 1./sgal2;
  for(int bb = 0; bb < (int)(nblocks); bb = bb + 1){
    int pb = NB[bb];
    mat MG = beta02*F4insGRinsG(bb, 0) + xi2*F4insgRinsg(bb, 0) + invsgga2*diagmat(ones(pb, 1));
    mat MA = F4insGRinsG(bb, 0) + invsgal2*diagmat(ones(pb, 1));
    mat SigG = inv(MG);
    mat SigA = inv(MA);

    double t1, t2;
    t1 = as_scalar((beta0* F4GinvsG2(bb, 0) + xi*F4ginvsg2(bb, 0)).t()*F4mu(bb, 0) -
      0.5*(F4mu(bb, 0).t()*MG*F4mu(bb, 0)))  + sum(log(diagvec(chol(SigG))));

    t2 = as_scalar(-beta0*F4mu(bb, 0).t()*F4insGRinsG(bb, 0)*F4muA(bb, 0)  + F4GinvsG2(bb, 0).t()*F4muA(bb, 0) -
      0.5 * F4muA(bb, 0).t()*MA*F4muA(bb, 0))  + sum(log(diagvec(chol(SigA))));

    tstat_b[bb] = t1 + t2;
  }
  // 
  double tstat = sum(tstat_b) - 0.5*p*(1 + log(sgga2)) - 0.5*p*(1 + log(sgal2)) + 0.5*p;

   
  List output = List::create(
    Rcpp::Named("Iteration") = Iteration,
    Rcpp::Named("loglik") = loglik_out,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("sgal2") = sgal2,
    Rcpp::Named("sgga2") = sgga2
  );
  return output;
}

class paraBlock_Gam1{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  
  vec betasig;
  vec betamean;
  vec Mu2;
  vec MuA2;
  vec xi_nu;
  vec xi_de;
  double beta0;
  double sgga2;
  double xi;
  int constr;
  
  
  field<vec> F4se1, F4se2, F4mu,  F4v2;
  field<vec> F4GinvsG2, F4ginvsg2;
  field<vec> F4RinsGmu, F4Rinsgmu;
  field<vec> F4RinsginsG;
  field<mat> F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg;
  
  paraBlock_Gam1(const uword &nblocks, field<vec> &F4se1, field<vec> &F4se2, field<vec> &F4ginvsg2, field<vec> &F4GinvsG2,  
                 field<vec> &F4v2,
                 field<mat> &F4Rins, field<mat> &F4Rins2, const field<mat> &F4Rblock, field<mat> &F4insGRinsG, field<mat> &F4insgRinsg,
                 field<vec> &F4mu,  field<vec> &F4Rinsgmu, field<vec> &F4RinsGmu,
                 arma::vec &betasig, arma::vec &betamean, arma::vec &Mu2,  arma::vec &xi_nu, arma::vec &xi_de,
                 const int& constr, double &xi, double &beta0, double &sgga2){
    
    this -> nblocks = nblocks;
    this -> F4se1 = F4se1;
    this -> F4se2 = F4se2;
    this -> F4ginvsg2 = F4ginvsg2;
    this -> F4GinvsG2 = F4GinvsG2;
    this -> F4v2 = F4v2;
    this -> F4Rins = F4Rins;
    this -> F4Rins2 = F4Rins2;
    this -> F4Rblock = F4Rblock;
    this -> F4insGRinsG = F4insGRinsG;
    this -> F4insgRinsg = F4insgRinsg;
    this -> F4mu = F4mu;
    this -> F4Rinsgmu = F4Rinsgmu;
    this -> F4RinsGmu = F4RinsGmu;
    this -> betasig = betasig;
    this -> betamean = betamean;
    this -> Mu2 = Mu2;
    this -> constr = constr;
    this -> xi = xi;
    this -> xi_nu = xi_nu;
    this -> xi_de = xi_de;
    this -> beta0 = beta0;
    this -> sgga2 = sgga2;
  }
  
  int  next_Gam1();
  void loop_by_block_gibbs_Gam1(int i);
  void update_by_thread_Gam1(int thread_id);
  
};

void paraBlock_Gam1::loop_by_block_gibbs_Gam1(int i){
  double xi2 = xi * xi;
  double beta02 = pow(beta0, 2);
  
  vec se1 = F4se1(i, 0);
  vec se2 = F4se2(i, 0);
  mat Rins = F4Rins(i, 0);
  mat Rins2 = F4Rins2(i, 0);
  mat Rblock = F4Rblock(i, 0);
  vec ginvsg2 = F4ginvsg2(i, 0);
  vec GinvsG2 = F4GinvsG2(i, 0);
  mat insgRinsg = F4insgRinsg(i,  0);
  mat insGRinsG = F4insGRinsG(i,  0);
  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);
  
  vec invse1 = 1. / se1;
  vec invse2 = 1. / se2;
  
  vec mu = F4mu(i, 0);
  vec Rinsgmu = F4Rinsgmu(i, 0);
  vec RinsGmu = F4RinsGmu(i, 0);
  vec Rdiag = diagvec(F4Rblock(i, 0));
  
  int p_block = F4Rblock(i, 0).n_rows;
  vec v2 = 1. / (beta02*diaginsGRinsG + xi2*diaginsgRinsg + 1. / sgga2);
  F4v2(i, 0) = v2;
  for (int j = 0; j < p_block; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];
    
    RinSmujj = Rinsgmu[j] - Rdiag[j]*mu[j] * invse1[j];
    RinSmujj2 = RinsGmu[j] - Rdiag[j]*mu[j] * invse2[j];
    
    mu[j] = (beta0*GinvsG2[j] - beta02 * invse2[j]*RinSmujj2  + xi*ginvsg2[j] - xi2 * invse1[j]*RinSmujj)*v2[j];
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];
  }
  
  // -----------------------------------------------------------------------
  F4mu(i, 0) = mu;
  F4Rinsgmu(i, 0) = Rinsgmu;
  F4RinsGmu(i, 0) = RinsGmu;
  // -----------------------------------------------------------------------
  // M step
  if (constr != 1){
    betasig[i] = as_scalar(mu.t()*insGRinsG*mu + v2.t()*diaginsGRinsG);
    betamean[i] = as_scalar(GinvsG2.t()*mu);
  }
  // -----------------------------------------------------------------------
  // for sgga2, sgal2 and xi iteration
  Mu2[i] = sum(mu%mu) + sum(v2);
  xi_nu[i] = as_scalar(ginvsg2.t()*mu);
  xi_de[i] = as_scalar(mu.t()*insgRinsg*mu + v2.t()*diaginsgRinsg);
  
  
  se1.reset();
  se2.reset();
  invse1.reset();
  invse2.reset();
  
  mu.reset();
  Rinsgmu.reset();
  RinsGmu.reset();
  
  Rins.reset();
  Rins2.reset();
  Rblock.reset();
  Rdiag.reset();
  
  ginvsg2.reset();
  GinvsG2.reset();
  insgRinsg.reset();
  insGRinsG.reset();
  diaginsGRinsG.reset();
  diaginsgRinsg.reset();
  
}

std::mutex _mtx6;
int paraBlock_Gam1::next_Gam1(){
  std::lock_guard<std::mutex> lockGuard(_mtx6);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_Gam1::update_by_thread_Gam1(int thread_id){
  while(true){
    int idx = next_Gam1();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_Gam1(idx);
  }
}

// [[Rcpp::export]]
List PXvbfunM1_par(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                   const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1,const  arma::vec& se2,
                   arma::vec& mu, double& sgga2, double& beta0, 
                   const int& constr, const double& epsStopLogLik, const int& maxIter, int coreNum){
  
  int p = bh1.n_elem;
  ivec NB = zeros<ivec>(nblocks, 1);
  
  field<vec> F4se1(nblocks, 1), F4se2(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1);
  field<vec> F4sg2(nblocks, 1), F4sG2(nblocks, 1), F4GinvsG2(nblocks, 1), F4ginvsg2(nblocks, 1), F4diaginsGRinsG(nblocks, 1);
  field<vec> F4diaginsgRinsg(nblocks, 1), F4RinsGmu(nblocks, 1), F4Rinsgmu(nblocks, 1);
  field<mat> F4insGRinsG(nblocks, 1), F4insgRinsg(nblocks, 1), F4Rins(nblocks, 1), F4Rins2(nblocks, 1);
  
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    vec se1_block = se1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec se2_block = se2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    vec mu_block = mu.subvec(block_inf(nn, 0), block_inf(nn, 1));
    // vec muA_block = muA.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh1_blcok = bh1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh2_blcok = bh2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    mat R_block = F4Rblock(nn, 0);
    NB[nn] = se1_block.n_elem;
    F4mu(nn, 0) = mu_block;
    // F4muA(nn, 0) = muA_block;
    F4se1(nn, 0) = se1_block;
    F4se2(nn, 0) = se2_block;
    
    F4sg2(nn, 0) = sg2_block;
    F4sG2(nn, 0) = sG2_block;
    
    F4GinvsG2(nn, 0) = bh2_blcok / sG2_block;
    F4ginvsg2(nn, 0) = bh1_blcok / sg2_block;
    
    F4insGRinsG(nn, 0) = diagmat(1. / se2_block)*R_block*diagmat(1. / se2_block);
    F4insgRinsg(nn, 0) = diagmat(1. / se1_block)*R_block*diagmat(1. / se1_block);
    
    F4diaginsGRinsG(nn, 0) = diagvec(F4insGRinsG(nn, 0));
    F4diaginsgRinsg(nn, 0) = diagvec(F4insgRinsg(nn, 0));
    
    F4Rins(nn, 0) = R_block*diagmat(1 / se1_block);
    F4Rins2(nn, 0) = R_block*diagmat(1 / se2_block);
    
    F4RinsGmu(nn, 0) = R_block*diagmat(1 / se2_block)*mu_block;
    F4Rinsgmu(nn, 0) = R_block*diagmat(1 / se1_block)*mu_block;
    // F4RinsGmuA(nn, 0) = R_block*diagmat(1. / se2_block)*muA_block;
    
    
  }
  
  double xi = 1;
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  vec v2;
  double xi2;
  double beta02 = beta0*beta0;
  field<vec> F4v2(nblocks, 1);
  for(int iter = 2; iter <= maxIter; iter++){
    // -----------------------------------------------------------------------
    // E step
    vec betamean = zeros(nblocks, 1); // mean for each block.
    vec betasig = zeros(nblocks, 1); // sigma for each block.
    vec Mu2 = zeros(nblocks, 1); // sgga2 for each block.
    vec xi_nu = zeros(nblocks, 1); // xi for each block.
    vec xi_de = zeros(nblocks, 1); // xi for each block.
    // parallel for gamma
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    
    
    paraBlock_Gam1 parobj1(nblocks, F4se1, F4se2, F4ginvsg2, F4GinvsG2,  
                           F4v2, F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg,
                           F4mu,F4Rinsgmu, F4RinsGmu,
                           betasig, betamean, Mu2,  xi_nu, xi_de,
                           constr, xi, beta0, sgga2);
    
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_Gam1::update_by_thread_Gam1, &parobj1, i_thread);
    }
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    
    betasig = parobj1.betasig;
    betamean = parobj1.betamean;
    Mu2 = parobj1.Mu2;
    xi_nu = parobj1.xi_nu;
    xi_de = parobj1.xi_de;
    F4v2 = parobj1.F4v2;
    F4RinsGmu = parobj1.F4RinsGmu;
    F4Rinsgmu = parobj1.F4Rinsgmu;
    F4mu = parobj1.F4mu;
    
    // -----------------------------------------------------------------------
    
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      double sig2b;
      sig2b = 1. / sum(betasig);
      beta0 = sum(betamean)*sig2b;
    }
    
    // update xi.
    xi = sum(xi_nu)/ sum(xi_de);
    //update sgga2
    sgga2 = sum(Mu2) / p;
    // ------------------------------------------------------------------
    // Reduction step.
    xi2 = xi*xi;
    beta0 = beta0 / xi;
    beta02 = beta0*beta0;
    sgga2 = sgga2*xi2;
    
    // lower bound to check convergence.
    vec low_b = zeros(nblocks, 1);
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){
      F4v2(ll, 0) = F4v2(ll, 0)*xi2;
      F4mu(ll, 0) = F4mu(ll, 0)*xi;
      F4RinsGmu(ll, 0) = F4RinsGmu(ll, 0)*xi;
      F4Rinsgmu(ll, 0) = F4Rinsgmu(ll, 0)*xi;
      
      double term1;
      int p_b = NB[ll];
      term1 = as_scalar((beta0*F4GinvsG2(ll, 0) + F4ginvsg2(ll, 0)).t()*F4mu(ll, 0)) -
        0.5*as_scalar(F4mu(ll, 0).t()*(beta02*F4insGRinsG(ll, 0) + F4insgRinsg(ll, 0) + 1./sgga2*diagmat(ones(p_b,1)))*F4mu(ll, 0)) -
        0.5*as_scalar(F4v2(ll, 0).t()*(beta02*F4diaginsGRinsG(ll, 0)  + F4diaginsgRinsg(ll, 0) + 1/sgga2)) + 0.5*sum(log(F4v2(ll, 0)));
      
      low_b[ll] = term1;
    }
    
    double low = sum(low_b) - 0.5*p*log(sgga2);
    // ------------------------------------------------------------------
    // Reduction step.
    xi = 1;
    xi2 = 1;
    
    
    
    loglik(iter-1) = low;
    if(loglik(iter -1) - loglik(iter -2) < -1e-7){
      perror("The likelihood failed to increase!");
      cout << loglik(iter -1) - loglik(iter -2) << endl;
    }
    
    Iteration = iter;
    if(iter > 2){
      if(abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik){
        break;
      }
    }
    
  }
  
  vec loglik_out;
  int to = Iteration - 1;
  loglik_out = loglik.subvec(0, to);
  double diff = loglik_out(to) - loglik_out(to - 1);
  
  // for loglikelihood ratio test.
  vec tstat_b = zeros(nblocks, 1);
  for(int bb = 0; bb < (int)(nblocks); bb = bb + 1){
    int pb = NB[bb];
    mat M = beta02*F4insGRinsG(bb, 0) + xi2*F4insgRinsg(bb, 0) + 1./sgga2*diagmat(ones(pb, 1));
    mat SigG = inv(M);
    
    double t1, t2;
    t1 = as_scalar((beta0* F4GinvsG2(bb, 0) + xi*F4ginvsg2(bb, 0)).t()*F4mu(bb, 0) - 
      0.5*(F4mu(bb, 0).t()*M*F4mu(bb, 0)));
    t2 = sum(log(diagvec(chol(SigG))));
    tstat_b[bb] = t1 + t2;
  }
  
  double tstat = sum(tstat_b) - 0.5*p*(1 + log(sgga2));
  
  
  List output = List::create(
    Rcpp::Named("Iteration") = Iteration,
    Rcpp::Named("loglik") = loglik_out,
    Rcpp::Named("diff") = diff,
    Rcpp::Named("tstat") = tstat,
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("sgga2") = sgga2
  );
  
  return output;
}




