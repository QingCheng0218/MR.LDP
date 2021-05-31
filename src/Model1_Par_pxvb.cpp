#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]

class paraBlock_VBgamma1{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  vec GinvsG2;
  vec ginvsg2, v2, RinsGmu, Rinsgmu, se1, se2, mu;
  mat insGRinsG, insgRinsg, Rins, Rins2, R;
  double sgga2 =  1;
  double beta0 = 0;
  
  paraBlock_VBgamma1(uword nblocks, umat block_inf, arma::vec &mu,  double beta0, double sgga2, 
                     arma::vec& RinsGmu, arma::vec& Rinsgmu, arma::vec& GinvsG2, arma::vec &ginvsg2, arma::vec &v2,
                     arma::vec& se1, arma::vec& se2, arma::mat& insGRinsG, arma::mat& insgRinsg, arma::mat& R, arma::mat& Rins, arma::mat& Rins2){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> mu = mu;
    this -> RinsGmu = RinsGmu;
    this -> Rinsgmu = Rinsgmu;
    this -> GinvsG2 = GinvsG2;
    this -> ginvsg2 = ginvsg2;
    this -> v2 = v2;
    this -> insGRinsG = insGRinsG;
    this -> insgRinsg = insgRinsg;
    this -> R = R;
    this -> Rins = Rins;
    this -> Rins2 = Rins2;
    this -> beta0 = beta0;
    this -> sgga2 = sgga2;
    this -> se1 = se1;
    this -> se2 = se2;
  }
  
  int  next_gamma1();
  void loop_by_block_PXVb_gamma1(int i);
  void update_by_thread_VBgamma1(int thread_id);
  
  
};

void PXVb1_update_gamma(arma::vec& GinvsG2, arma::vec& ginvsg2, arma::vec& v2, 
                        arma::mat& R, arma::mat& Rins, arma::mat& Rins2, arma::vec& RinsGmu, arma::vec& Rinsgmu,
                        arma::vec& se1, arma::vec& se2, arma::vec& mu, double beta0, double sgga2){
  int p = Rins.n_rows;
  double xi = 1;
  for (int j = 0; j < p; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];
    
    RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
    RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
    
    mu[j] = (beta0* GinvsG2[j] - pow(beta0, 2) / se2[j]*RinSmujj2  + xi*ginvsg2[j] - 1. / se1[j]*RinSmujj)*v2[j];
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];
  }
  
}


void paraBlock_VBgamma1::loop_by_block_PXVb_gamma1(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  
  
  vec GinvsG2_b, ginvsg2_b,  RinsGmu_b, Rinsgmu_b, se1_b, se2_b, v2_b, mu_b;
  mat Rins_b, Rins2_b, R_b;
  
  mu_b = mu.rows(index);
  GinvsG2_b = GinvsG2.rows(index);
  ginvsg2_b = ginvsg2.rows(index);
  RinsGmu_b = RinsGmu.rows(index);
  Rinsgmu_b = Rinsgmu.rows(index);
  se1_b = se1.rows(index);
  se2_b = se2.rows(index);
  v2_b = v2.rows(index);
  
  R_b = R.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  Rins_b = Rins.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  Rins2_b = Rins2.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  
  PXVb1_update_gamma(GinvsG2_b, ginvsg2_b, v2_b, R_b, Rins_b, Rins2_b, RinsGmu_b, Rinsgmu_b, se1_b, se2_b, mu_b, beta0, sgga2);
  
  mu.rows(index) = mu_b;
  RinsGmu.rows(index) = RinsGmu_b;
  Rinsgmu.rows(index) = Rinsgmu_b;
}

std::mutex _mtx11;
int paraBlock_VBgamma1::next_gamma1(){
  std::lock_guard<std::mutex> lockGuard(_mtx11);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_VBgamma1::update_by_thread_VBgamma1(int thread_id){
  while(true){
    int idx = next_gamma1();
    if(idx == -1){
      break;
    }
    loop_by_block_PXVb_gamma1(idx);
  }
}

// [[Rcpp::export]]
Rcpp::List Para_PXVb1(arma::mat R, arma::umat block_inf, uword nblocks,
                      arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                      arma::vec& mu,  double& beta0, double& sgga2, int coreNum,
                      const int& constr, const double& epsStopLogLik, const int& maxIter){
  
  int p = bh1.size();
  double xi, beta02, xi2;
  
  vec sg2 = pow(se1, 2);
  vec sG2 = pow(se2, 2);
  
  vec GinvsG2 = bh2 / sG2;
  vec ginvsg2 = bh1 / sg2;
  
  mat insGRinsG = diagmat(1. / se2)*R*diagmat(1. / se2);
  mat insgRinsg = diagmat(1. / se1)*R*diagmat(1. / se1);
  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);
  
  mat Rins = R*diagmat(1 / se1);
  mat Rins2 = R*diagmat(1 / se2);
  vec RinsGmu = R*diagmat(1 / se2)*mu;
  vec Rinsgmu = R*diagmat(1 / se1)*mu;
  vec v2;
  
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  
  for(int iter = 2; iter <= maxIter; iter ++){
    v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);
    // E step parallel.
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    paraBlock_VBgamma1 parobj_gamma1(nblocks, block_inf, mu,  beta0, sgga2,
                                     RinsGmu, Rinsgmu, GinvsG2, ginvsg2, v2,
                                     se1, se2, insGRinsG, insgRinsg, R, Rins, Rins2);
    
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_VBgamma1::update_by_thread_VBgamma1, &parobj_gamma1, i_thread);
    }
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    mu = parobj_gamma1.mu;
    RinsGmu = parobj_gamma1.RinsGmu;
    Rinsgmu = parobj_gamma1.Rinsgmu;
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      double sig2b;
      sig2b = 1. / as_scalar(mu.t()*insGRinsG*mu + v2.t()*diaginsGRinsG);
      beta0 = as_scalar(GinvsG2.t()*mu)*sig2b;
    }
    //update sgga2
    sgga2 = as_scalar(sum(mu%mu + v2)) / p;
    // update xi
    xi = as_scalar(ginvsg2.t()*mu)/ as_scalar(mu.t()*insgRinsg*mu + v2.t()*diaginsgRinsg);
    // Reduction step.
    xi2 = xi*xi;
    mu = mu*xi;
    beta0 = beta0 / xi;
    beta02 = beta0*beta0;
    sgga2 = sgga2*xi2;
    v2 = v2*xi2;
    RinsGmu = RinsGmu*xi;
    Rinsgmu = Rinsgmu*xi;
    
    xi = 1;
    xi2 = 1;
    //lower bound to check convergence.
    double term1, term2, low;
    term1 = as_scalar((beta0*GinvsG2 + xi*ginvsg2).t()*mu) -
      0.5*as_scalar(mu.t()*(beta02*insGRinsG + xi2*insgRinsg + 1./sgga2*diagmat(ones(p,1)))*mu) -
      0.5*as_scalar(v2.t()*(beta02*diaginsGRinsG  + xi2*diaginsgRinsg + 1/sgga2))- 0.5*p*log(sgga2);
    
    term2 = 0.5*sum(log(v2));
    low = term1 + term2;
    loglik(iter-1) = low;
    if(loglik(iter -1) - loglik(iter -2) < -1e-7){
      perror("The likelihood failed to increase!");
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
  // 
  // for loglikelihood ratio test.
  mat M = pow(beta0, 2)*insGRinsG + xi2*insgRinsg + 1./sgga2*diagmat(ones(p, 1));
  mat SigG = inv(M);
  
  // 
  double lb1, tstat;
  lb1 = as_scalar((xi*ginvsg2 + beta0*GinvsG2).t()*mu - 0.5*(mu.t()*M*mu)) - 0.5*p*(1 +log(sgga2));
  
  tstat = lb1 + sum(log(diagvec(chol(SigG))));
  
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

