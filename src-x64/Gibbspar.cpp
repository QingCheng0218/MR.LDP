#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

class paraBlock_gamma1{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  field<mat> F4Rblock;
  
  vec bh1;
  vec bh2;
  vec s12;
  vec s22;
  
  vec gamma;
  
  double sgga2 =  1;
  double agm = 0.001, bgm = 0.001;
  double beta0 = 0;
  
  vec GinvsG2;
  mat insGRinsG;
  paraBlock_gamma1(uword nblocks, umat block_inf, vec &gamma, field<mat> F4Rblock, vec bh1, vec bh2, vec s12, vec s22,
                   double beta0, double sgga2, double agm, double bgm, arma::vec &GinvsG2, arma::mat &insGRinsG){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> gamma = gamma;
    this -> F4Rblock = F4Rblock;
    this -> bh1 = bh1;
    this -> bh2 = bh2;
    this -> s12 = s12;
    this -> s22 = s22;
    this -> beta0 = beta0;
    this -> sgga2 = sgga2;
    this -> agm = agm;
    this -> bgm = bgm;
    this -> GinvsG2 = GinvsG2;
    this -> insGRinsG = insGRinsG;
  }
  
  int  next_gamma1();
  void loop_by_block_gibbs_gamma1(int i);
  void update_by_thread_Gamma1(int thread_id);
  
  
};

class paraBlock_gamma{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  field<mat> F4Rblock;
  
  vec bh1;
  vec bh2;
  vec s12;
  vec s22;
  
  vec gamma;
  vec alpha;
  double sgga2 =  1;
  double sgal2 =  1;
  double agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001;
  double beta0 = 0;
  
  vec GinvsG2;
  mat insGRinsG;
  paraBlock_gamma(uword nblocks, umat block_inf, vec &gamma, vec &alpha, field<mat> F4Rblock, vec bh1, vec bh2, vec s12, vec s22,
                  double beta0, double sgga2, double sgal2, double agm, double bgm, arma::vec &GinvsG2, arma::mat &insGRinsG){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> gamma = gamma;
    this -> alpha = alpha;
    this -> F4Rblock = F4Rblock;
    this -> bh1 = bh1;
    this -> bh2 = bh2;
    this -> s12 = s12;
    this -> s22 = s22;
    this -> beta0 = beta0;
    this -> sgga2 = sgga2;
    this -> sgal2 = sgal2;
    this -> agm = agm;
    this -> bgm = bgm;
    this -> GinvsG2 = GinvsG2;
    this -> insGRinsG = insGRinsG;
  }
  
  int  next_gamma();
  void loop_by_block_gibbs_gamma(int i);
  void update_by_thread_Gamma(int thread_id);
  
  
  
};

class paraBlock_alpha{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  field<mat> F4Rblock;
  
  vec bh1;
  vec bh2;
  vec s12;
  vec s22;
  
  vec gamma;
  vec alpha;
  
  double sgal2;
  double aal, bal;
  double beta0;
  
  vec GinvsG2;
  
  paraBlock_alpha(uword nblocks, umat block_inf, vec &gamma, vec &alpha, field<mat> F4Rblock, vec bh1, vec bh2, vec s12, vec s22,
                  double beta0,  double sgal2, const double aal, const double bal, arma::vec GinvsG2){
    this -> nblocks = nblocks;
    this -> block_inf = block_inf;
    this -> gamma = gamma;
    this -> alpha = alpha;
    this -> F4Rblock = F4Rblock;
    this -> bh1 = bh1;
    this -> bh2 = bh2;
    this -> s12 = s12;
    this -> s22 = s22;
    this -> beta0 = beta0;
    this -> sgal2 = sgal2;
    this -> aal = aal;
    this -> bal = bal;
    this -> GinvsG2 = GinvsG2;
  }
  int  next_alpha();
  void loop_by_block_gibbs_alpha(int i);
  void update_by_thread_Alpha(int thread_id);
  
};
// 

std::mutex _mtx0;
int paraBlock_gamma1::next_gamma1(){
  std::lock_guard<std::mutex> lockGuard(_mtx0);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

std::mutex _mtx2;
int paraBlock_gamma::next_gamma(){
  std::lock_guard<std::mutex> lockGuard(_mtx2);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

std::mutex _mtx3;
int paraBlock_alpha::next_alpha(){
  std::lock_guard<std::mutex> lockGuard(_mtx3);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

//-----------------------------------------------------------------------------
void gibbs1_update_gamma(arma::vec gammah, arma::vec Gammah, arma::vec sg2, arma::vec sG2, arma::mat R,
                         arma::vec &mu, double &beta0, double &sgga2,double agm, 
                         double bgm, arma::vec &GinvsG2, arma::mat &insGRinsG){
  int p = R.n_rows;
  
  GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;
  
  insGRinsG = diagmat(1. / sqrt(sG2))*R*diagmat(1. / sqrt(sG2));
  mat insgRinsg = diagmat(1. / sqrt(sg2))*R*diagmat(1. / sqrt(sg2));
  
  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);
  
  mat Rins = R*diagmat(1 / sqrt(sg2));
  mat Rins2 = R*diagmat(1 / sqrt(sG2));
  vec RinsGmu = R*diagmat(1 / sqrt(sG2))*mu;
  vec Rinsgmu = R*diagmat(1 / sqrt(sg2))*mu;
  
  vec v2,v2A = zeros(p, 1);
  
  v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);
  for (int j = 0; j < p; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2,  tmpga;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];
    
    RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / sqrt(sg2[j]);
    RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / sqrt(sG2[j]);
    // RinSmuAjj = RinsGmuA[j];
    
    tmpga = (beta0*GinvsG2[j] - pow(beta0, 2) / sqrt(sG2[j])*RinSmujj2  + ginvsg2[j] - 1. / sqrt(sg2[j])*RinSmujj)*v2[j];
    mu[j] = tmpga + randn()*sqrt(v2[j]);
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];
  }
  
}

// [[Rcpp::export]]
void gibbs2_update_gamma(arma::vec gammah, arma::vec Gammah, arma::vec sg2, arma::vec sG2, arma::mat R,
                         arma::vec &mu, arma::vec &muA, double &beta0, double &sgga2, double & sgal2,
                         double agm, double bgm, arma::vec &GinvsG2, arma::mat &insGRinsG)
{
  int p = R.n_rows;
  
  GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;
  
  insGRinsG = diagmat(1. / sqrt(sG2))*R*diagmat(1. / sqrt(sG2));
  mat insgRinsg = diagmat(1. / sqrt(sg2))*R*diagmat(1. / sqrt(sg2));
  
  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);
  
  mat Rins = R*diagmat(1 / sqrt(sg2));
  mat Rins2 = R*diagmat(1 / sqrt(sG2));
  vec RinsGmu = R*diagmat(1 / sqrt(sG2))*mu;
  vec Rinsgmu = R*diagmat(1 / sqrt(sg2))*mu;
  vec RinsGmuA = R* diagmat(1. / sqrt(sG2))*muA;
  
  // cout<<"GinvsG2:"<< GinvsG2.t()<<endl;
  
  //E step
  vec v2 = zeros(p, 1);
  
  v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);
  for (int j = 0; j < p; j= j+1){
    vec tmp1, tmp2;
    double RinSmujj, RinSmujj2, RinSmuAjj, tmpga;
    tmp1 = Rinsgmu - Rins.col(j)*mu[j];
    tmp2 = RinsGmu - Rins2.col(j)*mu[j];
    
    RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / sqrt(sg2[j]);
    RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / sqrt(sG2[j]);
    RinSmuAjj = RinsGmuA[j];
    
    tmpga = (beta0*GinvsG2[j] - pow(beta0, 2) / sqrt(sG2[j])*RinSmujj2 - beta0 / sqrt(sG2[j])*RinSmuAjj + ginvsg2[j] - 1. / sqrt(sg2[j])*RinSmujj)*v2[j];
    mu[j] = tmpga + randn()*sqrt(v2[j]);
    Rinsgmu = tmp1 + Rins.col(j)*mu[j];
    RinsGmu = tmp2 + Rins2.col(j)*mu[j];
  }
  // cout<<"mu:"<< mu.t()<<endl;
  // return output;
}

//-----------------------------------------------------------------------------
// [[Rcpp::export]]
void gibbs2_update_alpha(arma::vec gammah, arma::vec Gammah, arma::vec sg2, arma::vec sG2, arma::mat R,
                         arma::vec mu, arma::vec &muA, double beta0,  double sgal2, arma::vec GinvsG2)
{
  
  int p = R.n_rows;
  
  mat Rins2 = R*diagmat(1 / sqrt(sG2));
  vec RinsGmu = R*diagmat(1 / sqrt(sG2))*mu;
  vec Rinsgmu = R*diagmat(1 / sqrt(sg2))*mu;
  
  vec RinsGmuA = R* diagmat(1. / sqrt(sG2))*muA;
  vec v2A = zeros(p, 1);
  
  v2A = 1. / (1. / (sG2)+1. / sgal2);
  for (int k = 0; k < p; k = k + 1){
    vec tmp3;
    double tmpal, RinSmuAkk;
    tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
    RinSmuAkk = RinsGmuA[k] - R(k, k)*muA[k] / sqrt(sG2[k]);
    tmpal = (GinvsG2[k] - beta0*(1. / sqrt(sG2[k]))*RinsGmu[k] - 1. / sqrt(sG2[k])*RinSmuAkk)*v2A[k];
    muA[k] = tmpal + randn()*sqrt(v2A[k]);
    RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
  }
  
}
//-----------------------------------------------------------------------------

void paraBlock_gamma1::loop_by_block_gibbs_gamma1(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  vec bh1_b, bh2_b, s12_b, s22_b;
  mat R_b;
  
  bh1_b = bh1.rows(index);
  s12_b = pow(s12.rows(index), 2);
  
  bh2_b = bh2.rows(index);
  s22_b = pow(s22.rows(index), 2);
  
  R_b = F4Rblock(i, 0);
  vec mu = gamma.rows(index);
  
  int p1 = R_b.n_rows;
  vec GinvsG2_b = zeros(p1, 1);
  mat insGRinsG_b = zeros(p1, p1);
  
  gibbs1_update_gamma(bh1_b, bh2_b, s12_b, s22_b, R_b, mu,  beta0, sgga2, agm, bgm, GinvsG2_b, insGRinsG_b);
  
  gamma.rows(index) = mu;
  
  GinvsG2(index) = GinvsG2_b;
  
  insGRinsG.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1)) = insGRinsG_b;
  
}
void paraBlock_gamma::loop_by_block_gibbs_gamma(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  vec bh1_b, bh2_b, s12_b, s22_b;
  mat R_b;
  
  bh1_b = bh1.rows(index);
  s12_b = pow(s12.rows(index), 2);
  
  bh2_b = bh2.rows(index);
  s22_b = pow(s22.rows(index), 2);
  
  R_b = F4Rblock(i, 0);
  vec mu = gamma.rows(index);
  vec muA = alpha.rows(index);
  
  int p1 = R_b.n_rows;
  vec GinvsG2_b = zeros(p1, 1);
  mat insGRinsG_b = zeros(p1, p1);
  
  gibbs2_update_gamma(bh1_b, bh2_b, s12_b, s22_b, R_b, mu, muA, beta0, sgga2, sgal2, agm, bgm, GinvsG2_b, insGRinsG_b);
  
  gamma.rows(index) = mu;
  
  GinvsG2(index) = GinvsG2_b;
  
  insGRinsG.submat(index(0), index(0), index(index.size() - 1), index(index.size() - 1)) = insGRinsG_b;
  
}
void paraBlock_alpha::loop_by_block_gibbs_alpha(int i){
  
  uword tmp1 = block_inf(i, 1) - block_inf(i, 0) + 1;
  uvec index =linspace<uvec>(block_inf(i, 0), block_inf(i, 1), tmp1);
  
  vec bh1_b, bh2_b, s12_b, s22_b;
  mat R_b;
  
  bh1_b = bh1.rows(index);
  s12_b = pow(s12.rows(index), 2);
  
  bh2_b = bh2.rows(index);
  s22_b = pow(s22.rows(index), 2);
  
  vec GinvsG2_b = GinvsG2.rows(index);
  
  R_b = F4Rblock(i, 0);
  vec mu = gamma.rows(index);
  vec muA = alpha.rows(index);
  
  gibbs2_update_alpha(bh1_b, bh2_b, s12_b, s22_b, R_b, mu, muA, beta0, sgal2, GinvsG2_b);
  
  alpha.rows(index) = muA;
  
  
}

void paraBlock_gamma1::update_by_thread_Gamma1(int thread_id){
  while(true){
    int idx = next_gamma1();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_gamma1(idx);
  }
}

void paraBlock_gamma::update_by_thread_Gamma(int thread_id){
  while(true){
    int idx = next_gamma();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_gamma(idx);
  }
}
void paraBlock_alpha::update_by_thread_Alpha(int thread_id){
  while(true){
    int idx = next_alpha();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_alpha(idx);
  }
}

// [[Rcpp::export]]
void Para_Gibbs1(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                 arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                 arma::vec &gamma,  double &beta0, double &sgga2, int coreNum){
  
  uword P = bh1.size();
  vec GinvsG2 = zeros(P, 1);
  mat insGRinsG = zeros(P, P);
  
  double agm = 0.001, bgm = 0.001;
  
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  
  paraBlock_gamma1 parobj_gamma1(nblocks, block_inf, gamma,  F4Rblock,
                                 bh1, bh2, s12, s22, beta0, sgga2, agm, bgm, GinvsG2, insGRinsG);
  
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_gamma1::update_by_thread_Gamma1, &parobj_gamma1, i_thread);
  }
  
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  gamma = parobj_gamma1.gamma;
  GinvsG2 = parobj_gamma1.GinvsG2;
  insGRinsG = parobj_gamma1.insGRinsG;
  vec mu = gamma;
  // -----------------------------------------------------------------------
  //update beta0
  double sig2b, mub;
  sig2b = 1. / as_scalar(mu.t()*insGRinsG*mu);
  mub = as_scalar(GinvsG2.t()*mu)*sig2b;
  beta0 = mub + randn()*sqrt(sig2b);
  //update sgga2
  double tagm, tbgm;
  tagm = agm + P / 2;
  tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
  sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
  
}
// [[Rcpp::export]]
void Para_Gibbs2(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                int coreNum){
  
  uword P = bh1.size();
  vec GinvsG2 = zeros(P, 1);
  mat insGRinsG = zeros(P, P);
  
  double agm = 0.001, bgm = 0.001;
  double aal = 0.001, bal = 0.001;
  
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  
  paraBlock_gamma parobj_gamma(nblocks, block_inf, gamma, alpha, F4Rblock,
                               bh1, bh2, s12, s22, beta0, sgga2,sgal2, agm, bgm, GinvsG2, insGRinsG);
  
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_gamma::update_by_thread_Gamma, &parobj_gamma, i_thread);
  }
  
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  gamma = parobj_gamma.gamma;
  GinvsG2 = parobj_gamma.GinvsG2;
  insGRinsG = parobj_gamma.insGRinsG;
  
  // -----------------------------------------------------------------------
  // parallel for alpha
  paraBlock_alpha parobj_alpha(nblocks, block_inf, gamma, alpha, F4Rblock,
                               bh1, bh2, s12, s22, beta0, sgal2, aal, bal, GinvsG2);
  
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_alpha::update_by_thread_Alpha, &parobj_alpha, i_thread);
  }
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  alpha = parobj_alpha.alpha;
  vec mu = gamma;
  vec muA = alpha;
  // -----------------------------------------------------------------------
  //update beta0
  double sig2b, mub;
  sig2b = 1. / as_scalar(mu.t()*insGRinsG*mu);
  mub = as_scalar(GinvsG2.t()*mu - muA.t()*insGRinsG*mu)*sig2b;
  beta0 = mub + randn()*sqrt(sig2b);
  //update sgga2
  double tagm, tbgm, taal, tbal;
  tagm = agm + P / 2;
  tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
  sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
  //update sgal2
  taal = aal + P / 2;
  tbal = as_scalar(muA.t()*muA) / 2 + bal;
  sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
  
}

// [[Rcpp::export]]
Rcpp::List gibbsres1(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                     arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                     arma::vec &gamma,  double &beta0, double &sgga2, 
                     int coreNum, int IterMax){
  
  vec BETAres = zeros(IterMax, 1);
  vec SGGMres = zeros(IterMax, 1);
  vec SGALres = zeros(IterMax, 1);
  int p = bh1.size();
  mat GAres = zeros(p, IterMax);
  mat ALres = zeros(p, IterMax);
  
  for (int iter = 0; iter < IterMax; iter++){
    Para_Gibbs1(F4Rblock, block_inf, nblocks, bh1, bh2, s12, s22, gamma, beta0, sgga2, coreNum);
    
    BETAres[iter] = beta0;
    SGGMres[iter] = sgga2;
    GAres.col(iter) = gamma;
    
  };
  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres,
    Rcpp::Named("GAres") = GAres
  );
  return output;
}

// [[Rcpp::export]]
Rcpp::List gibbsres2(arma::field<mat> F4Rblock, arma::umat block_inf, uword nblocks,
                    arma::vec bh1, arma::vec bh2, arma::vec s12, arma::vec s22,
                    arma::vec &gamma, arma::vec &alpha, double &beta0, double &sgga2, double &sgal2,
                    int coreNum, int IterMax){
  
  vec BETAres = zeros(IterMax, 1);
  vec SGGMres = zeros(IterMax, 1);
  vec SGALres = zeros(IterMax, 1);
  int p = bh1.size();
  mat GAres = zeros(p, IterMax);
  mat ALres = zeros(p, IterMax);
  
  for (int iter = 0; iter < IterMax; iter++){
    Para_Gibbs2(F4Rblock, block_inf, nblocks, bh1, bh2, s12, s22,
               gamma, alpha, beta0, sgga2, sgal2, coreNum);
    
    BETAres[iter] = beta0;
    SGGMres[iter] = sgga2;
    SGALres[iter] = sgal2;
    GAres.col(iter) = gamma;
    ALres.col(iter) = alpha;
  };
  List output = List::create(
    Rcpp::Named("BETAres") = BETAres,
    Rcpp::Named("SGGMres") = SGGMres,
    Rcpp::Named("GAres") = GAres,
    Rcpp::Named("SGALres") = SGALres,
    Rcpp::Named("ALres") = ALres
  );
  return output;
}

// [[Rcpp::export]]
void Varres(field<mat> &F4H2a, arma::mat R, arma::umat block_inf, uword nblocks, 
            arma::vec bh2, arma::vec se2, arma::ivec N2){
  // bh2 and se2 are summary statistic for outcome, se2 means standard error 2
  // field<mat> F4H2a(nblocks, 1);
  for(int ii = 0; ii < (int)(nblocks); ii++){
    int idex1 = block_inf(ii, 0);
    int idex2 = block_inf(ii, 1);
    int nb = idex2 - idex1 + 1;
    arma::mat hha = zeros(nb, nb);
    arma::ivec bidex = linspace<ivec>(idex1, idex2, nb); 
    for(int j = 0; j < nb; j++){
      double tmp1 = N2[bidex[j]]*se2[bidex[j]]*se2[bidex[j]] + bh2[bidex[j]]*bh2[bidex[j]];
      for(int k = 0; k < nb; k++){
        double tmp2 = N2[bidex[k]]*se2[bidex[k]]*se2[bidex[k]]  + bh2[bidex[k]]*bh2[bidex[k]];
        hha(j, k) = R(bidex[j], bidex[k])/sqrt(tmp1 * tmp2);
      }
    }
    F4H2a(ii, 0) = hha;
  }
  
  
}

// [[Rcpp::export]]
arma::vec Herit_iMax(arma::mat R, arma::umat block_inf, uword nblocks, 
                     arma::vec bh2, arma::vec se2, arma::ivec N2, arma::mat ALres){
  // caculate heritability for all iterations.
  field<mat> F4H2a(nblocks, 1);
  Varres(F4H2a, R, block_inf, nblocks, bh2, se2, N2);
  int iMax = ALres.n_cols;
  vec H2a = zeros(iMax, 1);
  for(int itr = 0; itr < iMax; itr ++){
    vec alpha = ALres.col(itr);
    double h2a = 0;
    for(int ii = 0; ii < (int)(nblocks); ii++){
      int idex1 = block_inf(ii, 0);
      int idex2 = block_inf(ii, 1);
      int nb = idex2 - idex1 + 1;
      // arma::mat hha = zeros(nb, nb);
      arma::uvec bidex = linspace<uvec>(idex1, idex2, nb);
      
      
      vec alpha_b = alpha.rows(bidex);
      mat varres_b = F4H2a(ii, 0);
      h2a += as_scalar(alpha_b.t()*varres_b*alpha_b);
    }
    H2a[itr] = h2a;
  }
  
  return H2a;
}

// [[Rcpp::export]]
double heritability(field<mat> F4H2a, arma::umat block_inf, uword nblocks,  arma::vec  alpha){
  // caculate heritability for 1 iteration.
  double h2a = 0;
  for(int ii = 0; ii < (int)(nblocks); ii++){
    int idex1 = block_inf(ii, 0);
    int idex2 = block_inf(ii, 1);
    int nb = idex2 - idex1 + 1;
    // arma::mat hha = zeros(nb, nb);
    arma::uvec bidex = linspace<uvec>(idex1, idex2, nb);
    
    
    vec alpha_b = alpha.rows(bidex);
    mat varres_b = F4H2a(ii, 0);
    h2a += as_scalar(alpha_b.t()*varres_b*alpha_b);
  }
  
  return h2a;
}
