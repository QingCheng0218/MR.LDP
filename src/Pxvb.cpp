#include "RcppArmadillo.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List vbfunM1(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
             arma::vec mu, double& sgga2, double beta0,arma::mat R,  
             const int&constr, const double &epsStopLogLik, const int& maxIter){
  
  int p = bh1.n_elem;
  
  
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
  
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  vec v2 = zeros(p, 1);
  for(int iter = 2; iter <= maxIter; iter++){
    v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);
    // cout << "v2" << v2.subvec(0, 4).t() << endl;
    for (int j = 0; j < p; j= j+1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];
      
      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
      
      mu[j] = (beta0*GinvsG2[j] - beta0*beta0*RinSmujj2 / se2[j]  + ginvsg2[j] - RinSmujj / se1[j])*v2[j];
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
      
    }
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
    sgga2 = (sum(mu%mu) + sum(v2)) / p;
    
    double term1, term2, term3, term4, term5, low;
    
    term1 = as_scalar((beta0*GinvsG2 + ginvsg2).t()*mu);
    term2 = - 0.5*as_scalar(mu.t()*(pow(beta0, 2)*insGRinsG + insgRinsg + 1./sgga2*diagmat(ones(p,1)))*mu);
    term3 = - 0.5*as_scalar(v2.t()*(pow(beta0, 2)*diaginsGRinsG  + diaginsgRinsg + 1/sgga2));
    term4 = - 0.5*p*log(sgga2);
    term5 = 0.5*sum(log(v2));
    low = term1 + term2 + term3 + term4 + term5;
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
  // 
  // for loglikelihood ratio test.
  mat M = pow(beta0, 2)*insGRinsG + insgRinsg + 1./sgga2*diagmat(ones(p, 1));
  mat SigG = inv(M);
  
  
  double lb1, tstat;
  lb1 = as_scalar((ginvsg2 + beta0*GinvsG2).t()*mu - 0.5*(mu.t()*M*mu)) - 0.5*p*(1 +log(sgga2));
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

// [[Rcpp::export]]
List PXvbfunM1(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
               arma::vec mu, double& sgga2, double beta0, arma::mat R,
               const int&constr, const double &epsStopLogLik, const int& maxIter){
  int p = bh1.n_elem;
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
  
  double xi = 1;
  vec loglik(maxIter);
  // initialization of likelihood.
  loglik(0) = NAN;
  
  int Iteration = 1;
  vec v2 = zeros(p, 1);
  double xi2, beta02;
  for(int iter = 2; iter <= maxIter; iter++){
    v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);
    for (int j = 0; j < p; j= j+1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];
      
      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
      
      mu[j] = (beta0*GinvsG2[j] - pow(beta0, 2) / se2[j]*RinSmujj2  + xi*ginvsg2[j] - 1. / se1[j]*RinSmujj)*v2[j];
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
    }
    
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
    // sgga2 = as_scalar(sum(mu%mu + v2)) / p;
    sgga2 = (sum(mu%mu) + sum(v2)) / p;
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


// [[Rcpp::export]]
List vbfunM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
             arma::vec mu, arma::vec muA, double& sgga2, double& sgal2, double beta0,
             arma::mat R, const int&constr, const double &epsStopLogLik, const int& maxIter){

  int p = bh1.n_elem;

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
  double beta02 = beta0*beta0;
  for(int iter = 2; iter <= maxIter; iter ++){


    v2 = 1. / (beta02*diaginsGRinsG + diaginsgRinsg + 1. / sgga2);

    for (int j = 0; j < p; j= j + 1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2, RinSmuAjj;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];

      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
      RinSmuAjj = RinsGmuA[j];

      mu[j] =(beta0*GinvsG2[j] - beta02*RinSmujj2/ se2[j] - beta0*RinSmuAjj / se2[j] + ginvsg2[j] - RinSmujj / se1[j])*v2[j];
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
    }



    v2A = 1. / (1. / (sG2)+1. / sgal2);
    for (int k = 0; k < p; k = k + 1){
      vec tmp3;
      double RinSmuAkk;
      tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
      RinSmuAkk = RinsGmuA[k] - R(k, k)*muA[k] / se2[k];
      muA[k] = (GinvsG2[k] - beta0*(1. / se2[k])*RinsGmu[k] - 1. / se2[k]*RinSmuAkk)*v2A[k];
      RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
    }

    //update beta0
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
    beta02 = beta0*beta0;
    //update sgga2
    sgga2 = (sum(mu%mu) + sum(v2)) / p;
    // //update sgal2
    sgal2 = (sum(muA%muA) + sum(v2A)) / p;

    //lower bound to check convergence.
    double term1, term2, term3, term4, term5, term6, low;
    term1 = as_scalar((beta0*GinvsG2 + ginvsg2).t()*mu + GinvsG2.t()*muA);
    term2 = - 0.5*(as_scalar(mu.t()*(beta02*insGRinsG + insgRinsg)*mu) + sum(mu%mu) / sgga2);
    term3 = - 0.5*(as_scalar(v2.t()*(beta02*diaginsGRinsG  + diaginsgRinsg)) + sum(v2)/sgga2);
    term4 = - 0.5*p*log(sgga2) + 0.5*sum(log(v2)) - 0.5*p*log(sgal2) + 0.5*sum(log(v2A));
    term5 = - as_scalar(beta0*muA.t()*insGRinsG.t()*mu + 0.5*muA.t()*insGRinsG*muA) - 0.5*sum(muA%muA)/sgal2;
    term6 = - 0.5* (as_scalar(v2A.t()*diaginsGRinsG) + sum(v2A)/sgal2);

    low = term1 + term2 + term3 + term4 + term5 + term6;
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
  mat MG = pow(beta0, 2)*insGRinsG + insgRinsg + 1./sgga2*diagmat(ones(p, 1));
  mat MA = insGRinsG + 1./sgal2*diagmat(ones(p, 1));
  mat SigG = inv(MG);
  mat SigA = inv(MA);

  double t1, t2, tstat;
  t1 = as_scalar((beta0* GinvsG2 + ginvsg2).t()*mu - 0.5*(mu.t()*MG*mu)) -
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

// [[Rcpp::export]]
List PXvbfunM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
               arma::vec mu, arma::vec muA, double& sgga2, double& sgal2, double beta0, arma::mat R,
               const int&constr, const double &epsStopLogLik, const int& maxIter){

  int p = bh1.n_elem;

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
  double xi = 1;
  double xi2, beta02;
  for(int iter = 2; iter <= maxIter; iter ++){

    v2 = 1. / (pow(beta0, 2)*diaginsGRinsG + pow(xi,2)*diaginsgRinsg + 1. / sgga2);
    for (int j = 0; j < p; j= j+1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2, RinSmuAjj;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];

      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - R(j, j)*mu[j] / se2[j];
      RinSmuAjj = RinsGmuA[j];

      mu[j] =(beta0*GinvsG2[j] - pow(beta0, 2) / se2[j]*RinSmujj2 - beta0 / se2[j]*RinSmuAjj + xi*ginvsg2[j] - pow(xi, 2) / se1[j]*RinSmujj)*v2[j];
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
    }

    v2A = 1. / (1. / (sG2)+1. / sgal2);
    for (int k = 0; k < p; k = k + 1){
      vec tmp3;
      double RinSmuAkk;
      tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
      RinSmuAkk = RinsGmuA[k] - R(k, k)*muA[k] / se2[k];
      muA[k] = (GinvsG2[k] - beta0*(1. / se2[k])*RinsGmu[k] - 1. / se2[k]*RinSmuAkk)*v2A[k];
      RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
    }

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
    //update sgga2
    sgga2 = (sum(mu%mu) + sum(v2)) / p;
    // //update sgal2
    sgal2 = (sum(muA%muA) + sum(v2A)) / p;

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

// [[Rcpp::export]]
List VBMgM2(arma::vec bh1, arma::vec bh2, arma::vec se1, arma::vec se2,
            arma::vec mu, double& sgga2, double sgal2, double beta0, arma::mat R,
            const int&constr, const double &epsStopLogLik, const int& maxIter){

  int p = bh1.n_elem;

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


  mat sGRsG = diagmat(sqrt(sG2))*R*diagmat(sqrt(sG2));
  mat sGRinsG = diagmat(sqrt(sG2))*R*diagmat(1. / sqrt(sG2));
  mat A = sGRsG + sGRinsG*sGRinsG.t()*sgal2;
  mat inA = inv(A);
  mat B = sGRinsG.t()*inA*sGRinsG;
  mat Btemp = inA*sGRinsG;
  vec B1 = inA*sGRinsG*bh2;
  vec v2 = zeros(p, 1);
  vec Bmu = B*mu;
  double xi = 1;
  vec loglik(maxIter);
  //Initialiazation of likelihood.
  loglik(0) = NAN;
  int Iteration = 1;
  double xi2, beta02;

  for(int iter = 2; iter <= maxIter; iter++){

    v2 = 1. / (pow(beta0, 2)*diagvec(B) + pow(xi,2)*diaginsgRinsg + 1. / sgga2);
    for (int j = 0; j < p; j= j+1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = Bmu - B.col(j)*mu[j];

      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / sqrt(sg2[j]);
      RinSmujj2 = Bmu[j] - mu[j]*B(j, j);

      mu[j] =(beta0*B1[j] - pow(beta0, 2)*RinSmujj2 +
        xi*ginvsg2[j] - pow(xi, 2) / sqrt(sg2[j])*RinSmujj)*v2[j];
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      Bmu = tmp2 + B.col(j)*mu[j];
    }

    //update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      double sig2b;
      sig2b = 1. / as_scalar(mu.t()*B*mu + v2.t()*diagvec(B));
      beta0 = as_scalar(B1.t()*mu)*sig2b;
    }

    //update sgga2
    sgga2 = as_scalar(sum(mu%mu + v2)) / p;
    //update xi.
    // xi = as_scalar(ginvsg2.t()*mu)/ as_scalar(mu.t()*insgRinsg*mu + v2.t()*diaginsgRinsg);
    xi = 1;
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
    double term1, low;
    term1 = as_scalar((beta0*B1 + xi*ginvsg2).t()*mu) -
      0.5*as_scalar(mu.t()*(beta02*B + xi2*insgRinsg + 1./sgga2*diagmat(ones(p,1)))*mu) -
      0.5*as_scalar(v2.t()*(beta02*diagvec(B)  + xi2*diaginsgRinsg + 1/sgga2))- 0.5*p*log(sgga2) + 0.5*sum(log(v2));

    low = term1;
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

  mat MG = beta02*B + xi2*insgRinsg + 1./sgga2*diagmat(ones(p, 1));
  mat SigG = inv(MG);

  double t1,  t2, tstat;
  t1 = as_scalar((beta0* B1 + xi*ginvsg2).t()*mu - 0.5*(mu.t()*MG*mu)) - 0.5*p*(1 + log(sgga2)) + sum(log(diagvec(chol(SigG))));
  t2 = sum(log(diagvec(chol(A))));
  tstat = t1 - t2;

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
List MgMgib2(arma::vec gammah, arma::vec Gammah, arma::vec sg2, arma::vec sG2,
             arma::vec mu, double& sgga2, double sgal2, double beta0,
             arma::mat R, int IterMax, int agm, double bgm){

  int p = gammah.n_elem;

  vec Betares = zeros(IterMax, 1);
  vec Sgga2res = zeros(IterMax, 1);

  vec GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;

  mat insGRinsG = diagmat(1. / sqrt(sG2))*R*diagmat(1. / sqrt(sG2));
  mat insgRinsg = diagmat(1. / sqrt(sg2))*R*diagmat(1. / sqrt(sg2));

  vec diaginsGRinsG = diagvec(insGRinsG);
  vec diaginsgRinsg = diagvec(insgRinsg);

  mat Rins = R*diagmat(1 / sqrt(sg2));
  mat Rins2 = R*diagmat(1 / sqrt(sG2));
  vec RinsGmu = R*diagmat(1 / sqrt(sG2))*mu;
  vec Rinsgmu = R*diagmat(1 / sqrt(sg2))*mu;

  mat sGRsG = diagmat(sqrt(sG2))*R*diagmat(sqrt(sG2));
  mat sGRinsG = diagmat(sqrt(sG2))*R*diagmat(1. / sqrt(sG2));
  mat A = sGRsG + sGRinsG*sGRinsG.t()*sgal2;
  mat inA = inv(A);
  mat B = sGRinsG.t()*inA*sGRinsG;
  mat Btemp = inA*sGRinsG;
  vec B1 = inA*sGRinsG*Gammah;
  vec v2 = zeros(p, 1);
  vec Bmu = B*mu;

  double tmpga;

  for(int iter = 1; iter <= IterMax; iter++){

    v2 = 1. / (pow(beta0, 2)*diagvec(B) + diaginsgRinsg + 1. / sgga2);
    for (int j = 0; j < p; j= j+1){
      vec tmp1, tmp2;
      double RinSmujj, RinSmujj2;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = Bmu - B.col(j)*mu[j];

      RinSmujj = Rinsgmu[j] - R(j, j)*mu[j] / sqrt(sg2[j]);
      RinSmujj2 = Bmu[j] - mu[j]*B(j, j);

      tmpga =(beta0*B1[j] - pow(beta0, 2)*RinSmujj2 + ginvsg2[j] - 1. / sqrt(sg2[j])*RinSmujj)*v2[j];
      mu[j] = tmpga + randn()*sqrt(v2[j]);
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      Bmu = tmp2 + B.col(j)*mu[j];
    }

    //update beta0
    double sig2b, mub;
    sig2b = 1. / as_scalar(mu.t()*B*mu);
    mub = as_scalar(B1.t()*mu)*sig2b;
    beta0 = mub + randn()*sqrt(sig2b);

    Betares[iter] = beta0;
    //update sgga2
    double tagm, tbgm;
    tagm = agm + p / 2;
    tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    Sgga2res[iter] = sgga2;

  }

  List output = List::create(
    Rcpp::Named("BETAres") = Betares,
    Rcpp::Named("SGGMres") = Sgga2res
  );
  return output;
}


