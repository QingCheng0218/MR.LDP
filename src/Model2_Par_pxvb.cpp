#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>



using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

class paraBlock_VBgamma2{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  uword nblocks;
  vec GinvsG2, ginvsg2, v2, RinsGmu, Rinsgmu, RinsGmuA, se1, se2,mu;
  mat insGRinsG, insgRinsg, Rins, Rins2, R;
  
  double sgga2 =  1;
  double beta0 = 0;
  
  paraBlock_VBgamma2(uword nblocks, umat block_inf, vec& mu,  double beta0, double sgga2, 
                     arma::vec& RinsGmu, arma::vec& Rinsgmu, arma::vec& RinsGmuA, arma::vec& GinvsG2, arma::vec& ginvsg2, arma::vec& v2,
                     arma::vec& se1, arma::vec& se2, arma::mat& insGRinsG, arma::mat& insgRinsg, arma::mat& R, arma::mat& Rins, arma::mat& Rins2){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> mu = mu;
    this -> RinsGmu = RinsGmu;
    this -> Rinsgmu = Rinsgmu;
    this -> RinsGmuA = RinsGmuA;
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
  
  int  next_gamma2();
  void loop_by_block_PXVb_gamma2(int i);
  void update_by_thread_VBgamma2(int thread_id);
  
};

class paraBlock_VBalpha{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  uword nblocks;
  vec GinvsG2, v2A, RinsGmu, RinsGmuA, se1, se2, muA;
  mat Rins2, R;
  
  double sgga2 =  1;
  double beta0 = 0;
  
  paraBlock_VBalpha(uword nblocks, umat block_inf, arma::vec& muA,  double beta0, 
                    arma::vec& RinsGmu, arma::vec& RinsGmuA, arma::vec& GinvsG2, arma::vec& v2A, arma::vec& se2, arma::mat& R,arma::mat& Rins2){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> muA = muA;
    this -> RinsGmu = RinsGmu;
    this -> RinsGmuA = RinsGmuA;
    this -> GinvsG2 = GinvsG2;
    this -> v2A = v2A;
    this -> R = R;
    this -> Rins2 = Rins2;
    this -> beta0 = beta0;
    this -> se2 = se2;
  }
  
  int  next_alpha();
  void loop_by_block_PXVb_alpha(int i);
  void update_by_thread_VBalpha(int thread_id);
  
};
void PXVb2_update_gamma(arma::vec& GinvsG2, arma::vec& ginvsg2, arma::vec& v2, 
                        arma::mat& R, arma::mat& Rins, arma::mat& Rins2, arma::vec& RinsGmu, arma::vec& Rinsgmu, arma::vec& RinsGmuA,
                        arma::vec& se1, arma::vec& se2, arma::vec& mu, double beta0, double sgga2){
  int p = Rins.n_rows;
  double xi = 1;
  for (int j = 0; j < p; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2, RinSmuAjj;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];
    
    RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
    RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
    RinSmuAjj = RinsGmuA[j];
    
    mu[j] =(beta0*GinvsG2[j] - pow(beta0, 2) / se2[j]*RinSmujj2 - 
      beta0 / se2[j]*RinSmuAjj + xi*ginvsg2[j] - pow(xi, 2) / se1[j]*RinSmujj)*v2[j];
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];
  }
  
}

void PXVb2_update_alpha(arma::vec& GinvsG2, arma::vec& v2A, 
                        arma::mat& R, arma::mat& Rins2, arma::vec& RinsGmu, arma::vec& RinsGmuA,
                        arma::vec& se2, arma::vec& muA, double beta0){
  int p = Rins2.n_rows;
  // double xi = 1;
  for (int k = 0; k < p; k = k + 1){
    vec tmp3;
    double RinSmuAkk;
    tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
    RinSmuAkk = RinsGmuA[k] - R(k, k)*muA[k] / se2[k];
    muA[k] = (GinvsG2[k] - beta0*(1. / se2[k])*RinsGmu[k] - 1. / se2[k]*RinSmuAkk)*v2A[k];
    RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
  }
  
}

void paraBlock_VBgamma2::loop_by_block_PXVb_gamma2(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  
  
  vec GinvsG2_b, ginvsg2_b,  RinsGmu_b, Rinsgmu_b, RinsGmuA_b, se1_b, se2_b, v2_b, mu_b;
  mat Rins_b, Rins2_b, R_b;
  
  mu_b = mu.rows(index);
  
  GinvsG2_b = GinvsG2.rows(index);
  ginvsg2_b = ginvsg2.rows(index);
  RinsGmu_b = RinsGmu.rows(index);
  Rinsgmu_b = Rinsgmu.rows(index);
  RinsGmuA_b = RinsGmuA.rows(index);
  
  se1_b = se1.rows(index);
  se2_b = se2.rows(index);
  v2_b = v2.rows(index);
  
  Rins_b = Rins.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  Rins2_b = Rins2.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  R_b = R.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  
  PXVb2_update_gamma(GinvsG2_b, ginvsg2_b, v2_b, R_b, Rins_b, Rins2_b, RinsGmu_b, Rinsgmu_b, RinsGmuA_b,
                     se1_b, se2_b, mu_b, beta0, sgga2);
  mu.rows(index) = mu_b;
  RinsGmu.rows(index) = RinsGmu_b;
  Rinsgmu.rows(index) = Rinsgmu_b;
}

void paraBlock_VBalpha::loop_by_block_PXVb_alpha(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  
  vec GinvsG2_b, RinsGmu_b, RinsGmuA_b, se2_b, v2A_b, muA_b;
  mat Rins2_b, R_b;
  
  muA_b = muA.rows(index);
  GinvsG2_b = GinvsG2.rows(index);
  RinsGmu_b = RinsGmu.rows(index);
  RinsGmuA_b = RinsGmuA.rows(index);
  se2_b = se2.rows(index);
  v2A_b = v2A.rows(index);
  
  Rins2_b = Rins2.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  R_b = R.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1));
  
  PXVb2_update_alpha(GinvsG2_b, v2A_b, R_b, Rins2_b, RinsGmu_b, RinsGmuA_b, se2_b, muA_b, beta0);
  muA.rows(index) = muA_b;
  RinsGmuA.rows(index) = RinsGmuA_b;
}

std::mutex _mtx21;
int paraBlock_VBgamma2::next_gamma2(){
  std::lock_guard<std::mutex> lockGuard(_mtx21);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_VBgamma2::update_by_thread_VBgamma2(int thread_id){
  while(true){
    int idx = next_gamma2();
    if(idx == -1){
      break;
    }
    loop_by_block_PXVb_gamma2(idx);
  }
}

std::mutex _mtx22;
int paraBlock_VBalpha::next_alpha(){
  std::lock_guard<std::mutex> lockGuard(_mtx22);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_VBalpha::update_by_thread_VBalpha(int thread_id){
  while(true){
    int idx = next_alpha();
    if(idx == -1){
      break;
    }
    loop_by_block_PXVb_alpha(idx);
  }
}

// [[Rcpp::export]]
Rcpp::List Para_PXVb2(arma::mat R, arma::umat block_inf, uword nblocks,
                      arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
                      arma::vec& mu, arma::vec& muA, double& beta0, double& sgga2, double& sgal2, int coreNum,
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
  vec RinsGmuA = R* diagmat(1. / se2)*muA;
  
  vec v2,v2A;
  
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  for(int iter = 2; iter <= maxIter; iter ++){
    // E step
    v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + pow(xi,2)*diaginsgRinsg + 1. / sgga2);
    // -----------------------------------------------------------------------
    
    // parallel for gamma
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    paraBlock_VBgamma2 parobj_gamma2(nblocks, block_inf, mu, beta0, sgga2, 
                                     RinsGmu, Rinsgmu, RinsGmuA, GinvsG2, ginvsg2, v2,
                                     se1, se2, insGRinsG, insgRinsg, R, Rins, Rins2);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_VBgamma2::update_by_thread_VBgamma2, &parobj_gamma2, i_thread);
    }
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    mu = parobj_gamma2.mu;
    RinsGmu = parobj_gamma2.RinsGmu;
    Rinsgmu = parobj_gamma2.Rinsgmu;
    v2A = 1. / (1. / sG2 + 1. / sgal2);
    // cout << mu.t()<< endl;
    // -----------------------------------------------------------------------
    // parallel for alpha
    paraBlock_VBalpha parobj_alpha(nblocks, block_inf, muA,  beta0, RinsGmu, RinsGmuA, GinvsG2, v2A, se2, R, Rins2);
    
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_VBalpha::update_by_thread_VBalpha, &parobj_alpha, i_thread);
    }
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    muA = parobj_alpha.muA;
    RinsGmuA = parobj_alpha.RinsGmuA;
    
    // M step 
    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      double sig2b;
      sig2b = 1. / as_scalar(mu.t()*insGRinsG*mu + v2.t()*diaginsGRinsG);
      beta0 = as_scalar(GinvsG2.t()*mu - muA.t()*insGRinsG*mu)*sig2b;
    }
    // update sgga2
    sgga2 = as_scalar(sum(mu%mu + v2)) / p;
    // update sgal2
    sgal2 = as_scalar(sum(muA%muA + v2A)) / p;
    // update xi.
    xi = as_scalar(ginvsg2.t()*mu)/ as_scalar(mu.t()*insgRinsg*mu + v2.t()*diaginsgRinsg);
    // xi = 1;
    
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
      0.5*as_scalar(v2.t()*(beta02*diaginsGRinsG  + xi2*diaginsgRinsg + 1/sgga2))- 0.5*p*log(sgga2) + 0.5*sum(log(v2));
    
    term2 = as_scalar(- beta0*muA.t()*insGRinsG.t()*mu + GinvsG2.t()*muA - 0.5*muA.t()*(insGRinsG + 1/sgal2 * diagmat(ones(p)))*muA) -
      0.5* as_scalar(v2A.t()*(diaginsGRinsG + 1/sgal2)) - 0.5*p*log(sgal2) + 0.5*sum(log(v2A));
    
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
  
  // for loglikelihood ratio test.
  mat MG = beta02*insGRinsG + xi2*insgRinsg + 1./sgga2*diagmat(ones(p, 1));
  mat MA = insGRinsG + 1./sgal2*diagmat(ones(p, 1));
  mat SigG = inv(MG);
  mat SigA = inv(MA);
  double t1, t2, tstat;
  t1 = as_scalar((beta0* GinvsG2 + xi*ginvsg2).t()*mu - 0.5*(mu.t()*MG*mu)) - 
    0.5*p*(1 + log(sgga2))  + sum(log(diagvec(chol(SigG))));
  
  t2 = as_scalar(-beta0*mu.t()*insGRinsG*muA  + GinvsG2.t()*muA - 
    0.5 * muA.t()*MA*muA)- 0.5*p*(1 + log(sgal2))  + sum(log(diagvec(chol(SigA)))) + 0.5*p;
  
  tstat = t1 + t2;
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

