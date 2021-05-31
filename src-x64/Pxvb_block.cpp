#include "RcppArmadillo.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]



// [[Rcpp::export]]
List PXvbfunM1_block(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                     const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1,const  arma::vec& se2,
                     arma::vec& mu, double& sgga2, double& beta0, 
                     const int& constr, const double& epsStopLogLik, const int& maxIter){
 
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
  double xi2, beta02;
  field<vec> F4v2(nblocks, 1);
  for(int iter = 2; iter <= maxIter; iter++){
    
    for(int nb1 = 0; nb1 < (int)(nblocks); nb1 = nb1 + 1){
      
      v2 = 1. / (pow(beta0, 2)*F4diaginsGRinsG(nb1, 0) + pow(xi,2)*F4diaginsgRinsg(nb1, 0) + 1. / sgga2);
      F4v2(nb1, 0) = v2;
      int p_block = NB[nb1];
      
      for (int j = 0; j < (int)p_block; j= j+1){
        vec tmp1, tmp2;
        double RinSmujj, RinSmujj2;
        tmp1 = F4Rinsgmu(nb1, 0) - F4Rins(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        tmp2 = F4RinsGmu(nb1, 0) - F4Rins2(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        
        RinSmujj = F4Rinsgmu(nb1, 0)[j] - F4Rblock(nb1, 0)(j, j)*F4mu(nb1, 0)[j] / F4se1(nb1, 0)[j];
        RinSmujj2 = F4RinsGmu(nb1, 0)[j] - F4Rblock(nb1, 0)(j, j)*F4mu(nb1, 0)[j] / F4se2(nb1, 0)[j];
        
        
        F4mu(nb1, 0)[j] =(beta0*F4GinvsG2(nb1, 0)[j] - pow(beta0, 2) / F4se2(nb1, 0)[j]*RinSmujj2 + 
          xi*F4ginvsg2(nb1, 0)[j] - pow(xi, 2) / F4se1(nb1, 0)[j]*RinSmujj)*v2[j];
        F4Rinsgmu(nb1, 0) = tmp1 + F4Rins(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        F4RinsGmu(nb1, 0) = tmp2 + F4Rins2(nb1, 0).col(j)*F4mu(nb1, 0)[j];
      }
    }
    
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      vec betanu = zeros(nblocks, 1); // numerator for each block.
      vec betade = zeros(nblocks, 1); // denominator for each block.
      for (int jj = 0; jj < (int)(nblocks); jj = jj + 1){
        betade[jj] = as_scalar(F4mu(jj, 0).t()*F4insGRinsG(jj, 0)*F4mu(jj, 0) + F4v2(jj, 0).t()*F4diaginsGRinsG(jj, 0));
        betanu[jj] = as_scalar(F4GinvsG2(jj, 0).t()*F4mu(jj, 0));
      }
      beta0 = sum(betanu) / sum(betade);
    }
    
    vec sgga2_b = zeros(nblocks, 1); 
    vec xi_nu = zeros(nblocks, 1), xi_de = zeros(nblocks, 1); 
    
    for(int kk = 0; kk < (int)(nblocks); kk = kk + 1){
      sgga2_b[kk] = sum(F4mu(kk, 0)%F4mu(kk, 0)) + sum(F4v2(kk, 0));
      xi_nu[kk] = as_scalar(F4ginvsg2(kk, 0).t()*F4mu(kk, 0));
      xi_de[kk] = as_scalar(F4mu(kk, 0).t()*F4insgRinsg(kk, 0)*F4mu(kk, 0) + F4v2(kk, 0).t()*F4diaginsgRinsg(kk, 0));
    }
    //update sgga2
    sgga2 = sum(sgga2_b) / p;
    // update xi.
    xi = sum(xi_nu)/sum(xi_de);
    
    // Reduction step.
    xi2 = xi*xi;
    beta0 = beta0 / xi;
    beta02 = beta0*beta0;
    sgga2 = sgga2*xi2;
    
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){
      F4v2(ll, 0) = F4v2(ll, 0)*xi2;
      F4mu(ll, 0) = F4mu(ll, 0)*xi;
      F4RinsGmu(ll, 0) = F4RinsGmu(ll, 0)*xi;
      F4Rinsgmu(ll, 0) = F4Rinsgmu(ll, 0)*xi;
    }
    
    xi = 1;
    xi2 = 1;

    
    //lower bound to check convergence.
    vec low_b = zeros(nblocks, 1);
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){
      double term1;
      int p_b = NB[ll];
      term1 = as_scalar((beta0*F4GinvsG2(ll, 0) + xi*F4ginvsg2(ll, 0)).t()*F4mu(ll, 0)) -
        0.5*as_scalar(F4mu(ll, 0).t()*(beta02*F4insGRinsG(ll, 0) + xi2*F4insgRinsg(ll, 0) + 1./sgga2*diagmat(ones(p_b,1)))*F4mu(ll, 0)) -
        0.5*as_scalar(F4v2(ll, 0).t()*(beta02*F4diaginsGRinsG(ll, 0)  + xi2*F4diaginsgRinsg(ll, 0) + 1/sgga2)) + 0.5*sum(log(F4v2(ll, 0)));
      
      low_b[ll] = term1;
    }
    
    double low = sum(low_b) - 0.5*p*log(sgga2);
    
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

// [[Rcpp::export]]
List PXvbfunM2_block(const arma::field<mat>& F4Rblock, const arma::umat& block_inf, const uword& nblocks,
                     const arma::vec& bh1, const arma::vec& bh2, const arma::vec& se1, const arma::vec& se2,
                     arma::vec& mu, arma::vec& muA, double& sgga2, double& sgal2, double& beta0, 
                     const int& constr, const double& epsStopLogLik, const int& maxIter){
  
  int p = bh1.n_elem;
  ivec NB = zeros<ivec>(nblocks, 1);
  
  field<vec> F4se1(nblocks, 1), F4se2(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1);
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

  vec v2,v2A;
  
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  double xi = 1;
  double xi2, beta02;
  field<vec> F4v2(nblocks, 1), F4v2A(nblocks, 1);
  for(int iter = 2; iter <= maxIter; iter ++){
    
    for(int nb1 = 0; nb1 < (int)(nblocks); nb1 = nb1 + 1){
      
      v2 = 1. / (pow(beta0, 2)*F4diaginsGRinsG(nb1, 0) + pow(xi,2)*F4diaginsgRinsg(nb1, 0) + 1. / sgga2);
      F4v2(nb1, 0) = v2;
      int p_block = NB[nb1];
      for (int j = 0; j < (int)p_block; j= j+1){
        vec tmp1, tmp2;
        double RinSmujj, RinSmujj2, RinSmuAjj;
        tmp1 = F4Rinsgmu(nb1, 0) - F4Rins(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        tmp2 = F4RinsGmu(nb1, 0) - F4Rins2(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        
        RinSmujj = F4Rinsgmu(nb1, 0)[j] - F4Rblock(nb1, 0)(j, j)*F4mu(nb1, 0)[j] / F4se1(nb1, 0)[j];
        RinSmujj2 = F4RinsGmu(nb1, 0)[j] - F4Rblock(nb1, 0)(j, j)*F4mu(nb1, 0)[j] / F4se2(nb1, 0)[j];
        RinSmuAjj = F4RinsGmuA(nb1, 0)[j];
        
        F4mu(nb1, 0)[j] =(beta0*F4GinvsG2(nb1, 0)[j] - pow(beta0, 2) / F4se2(nb1, 0)[j]*RinSmujj2 - 
          beta0 / F4se2(nb1, 0)[j]*RinSmuAjj + xi*F4ginvsg2(nb1, 0)[j] - pow(xi, 2) / F4se1(nb1, 0)[j]*RinSmujj)*v2[j];
        F4Rinsgmu(nb1, 0) = tmp1 + F4Rins(nb1, 0).col(j)*F4mu(nb1, 0)[j];
        F4RinsGmu(nb1, 0) = tmp2 + F4Rins2(nb1, 0).col(j)*F4mu(nb1, 0)[j];
      }
      
      v2A = 1. / (1. / (F4sG2(nb1, 0))+1. / sgal2);
      F4v2A(nb1, 0) = v2A;
      int pa_block = NB[nb1];
      for (int k = 0; k < (int)pa_block; k = k + 1){
        vec tmp3;
        double RinSmuAkk;
        tmp3 = F4RinsGmuA(nb1, 0) - F4Rins2(nb1, 0).col(k)*F4muA(nb1, 0)[k];
        RinSmuAkk = F4RinsGmuA(nb1, 0)[k] - F4Rblock(nb1, 0)(k, k)*F4muA(nb1, 0)[k] / F4se2(nb1, 0)[k];
        F4muA(nb1, 0)[k] = (F4GinvsG2(nb1, 0)[k] - beta0*(1. / F4se2(nb1, 0)[k])*F4RinsGmu(nb1, 0)[k] - 1. / F4se2(nb1, 0)[k]*RinSmuAkk)*v2A[k];
        F4RinsGmuA(nb1, 0) = tmp3 + F4Rins2(nb1, 0).col(k)*F4muA(nb1, 0)[k];
      }
      
    }
   
    
    // M step
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      vec betanu = zeros(nblocks, 1); // numerator for each block.
      vec betade = zeros(nblocks, 1); // denominator for each block.
      for (int jj = 0; jj < (int)(nblocks); jj = jj + 1){
        betade[jj] = as_scalar(F4mu(jj, 0).t()*F4insGRinsG(jj, 0)*F4mu(jj, 0) + F4v2(jj, 0).t()*F4diaginsGRinsG(jj, 0));
        betanu[jj] = as_scalar(F4GinvsG2(jj, 0).t()*F4mu(jj, 0) - F4muA(jj, 0).t()*F4insGRinsG(jj, 0)*F4mu(jj, 0));
      }
      beta0 = sum(betanu) / sum(betade);
    }

    vec sgga2_b = zeros(nblocks, 1); 
    vec sgal2_b = zeros(nblocks, 1);
    vec xi_nu = zeros(nblocks, 1), xi_de = zeros(nblocks, 1); 
    
    for(int kk = 0; kk < (int)(nblocks); kk = kk + 1){
      sgga2_b[kk] = sum(F4mu(kk, 0)%F4mu(kk, 0)) + sum(F4v2(kk, 0));
      sgal2_b[kk] = sum(F4muA(kk, 0)%F4muA(kk, 0)) + sum(F4v2A(kk, 0));
      xi_nu[kk] = as_scalar(F4ginvsg2(kk, 0).t()*F4mu(kk, 0));
      xi_de[kk] = as_scalar(F4mu(kk, 0).t()*F4insgRinsg(kk, 0)*F4mu(kk, 0) + F4v2(kk, 0).t()*F4diaginsgRinsg(kk, 0));
    }
    //update sgga2
    sgga2 = sum(sgga2_b) / p;
    // //update sgal2
    sgal2 = sum(sgal2_b) / p;
    
    // update xi.
    xi = sum(xi_nu)/sum(xi_de);
    // xi = 1;
    
    // Reduction step.
    xi2 = xi*xi;
    beta0 = beta0 / xi;
    beta02 = beta0*beta0;
    sgga2 = sgga2*xi2;
    
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){
      F4v2(ll, 0) = F4v2(ll, 0)*xi2;
      F4mu(ll, 0) = F4mu(ll, 0)*xi;
      F4RinsGmu(ll, 0) = F4RinsGmu(ll, 0)*xi;
      F4Rinsgmu(ll, 0) = F4Rinsgmu(ll, 0)*xi;
    }

    xi = 1;
    xi2 = 1;
    
    //lower bound to check convergence.
    vec low_b = zeros(nblocks, 1);
    for(int ll = 0; ll < (int)(nblocks); ll = ll + 1){
      double term1, term2;
      int p_b = NB[ll];
      term1 = as_scalar((beta0*F4GinvsG2(ll, 0) + xi*F4ginvsg2(ll, 0)).t()*F4mu(ll, 0)) -
        0.5*as_scalar(F4mu(ll, 0).t()*(beta02*F4insGRinsG(ll, 0) + xi2*F4insgRinsg(ll, 0) + 1./sgga2*diagmat(ones(p_b,1)))*F4mu(ll, 0)) -
        0.5*as_scalar(F4v2(ll, 0).t()*(beta02*F4diaginsGRinsG(ll, 0)  + xi2*F4diaginsgRinsg(ll, 0) + 1/sgga2)) + 0.5*sum(log(F4v2(ll, 0)));
      
      term2 = as_scalar(- beta0*F4muA(ll, 0).t()*F4insGRinsG(ll, 0).t()*F4mu(ll, 0) + F4GinvsG2(ll, 0).t()*F4muA(ll, 0) - 
        0.5*F4muA(ll, 0).t()*(F4insGRinsG(ll, 0) + 1/sgal2 * diagmat(ones(p_b)))*F4muA(ll, 0)) -
        0.5* as_scalar(F4v2A(ll, 0).t()*(F4diaginsGRinsG(ll, 0) + 1/sgal2)) + 0.5*sum(log(F4v2A(ll, 0)));
      
      low_b[ll] = term1 + term2;
    }
    
    double low = sum(low_b) - 0.5*p*log(sgga2) - 0.5*p*log(sgal2);
    
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
    mat MG = beta02*F4insGRinsG(bb, 0) + xi2*F4insgRinsg(bb, 0) + 1./sgga2*diagmat(ones(pb, 1));
    mat MA = F4insGRinsG(bb, 0) + 1./sgal2*diagmat(ones(pb, 1));
    mat SigG = inv(MG);
    mat SigA = inv(MA);
    
    double t1, t2;
    t1 = as_scalar((beta0* F4GinvsG2(bb, 0) + xi*F4ginvsg2(bb, 0)).t()*F4mu(bb, 0) - 
      0.5*(F4mu(bb, 0).t()*MG*F4mu(bb, 0)))  + sum(log(diagvec(chol(SigG))));
    
    t2 = as_scalar(-beta0*F4mu(bb, 0).t()*F4insGRinsG(bb, 0)*F4muA(bb, 0)  + F4GinvsG2(bb, 0).t()*F4muA(bb, 0) -
      0.5 * F4muA(bb, 0).t()*MA*F4muA(bb, 0))  + sum(log(diagvec(chol(SigA))));
    
    tstat_b[bb] = t1 + t2;
  }
 
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


