#ifndef MRLDPdefaultInitialValue_hpp
#define MRLDPdefaultInitialValue_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>

using namespace Rcpp;
using namespace arma;
using namespace std;

class Options_MRLD{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_MRLD(int p){
    this -> gamma = zeros(p, 1);
    this -> sgga2 = 0.01;
    this -> beta0 = 0;
    this -> epsStopLogLik = 1e-7;
    this -> maxIter = 10000; 
  }
  
  Options_MRLD(vec gamma, double sgga2, double beta0, double epsStopLogLik, int maxIter){
    this -> gamma = gamma;
    this -> sgga2 = sgga2;
    this -> beta0 = beta0;
    this -> epsStopLogLik = epsStopLogLik;
    this -> maxIter = maxIter;
  }
  vec gamma;
  double sgga2;
  double beta0;
  double epsStopLogLik;
  int maxIter;

};
struct ObjMRLD{
  int Iteration;
  vec loglik_out;
  double diff;
  double tstat;
  double beta0;
  double sgga2;
};


class Options_MRLDP{
public:
  // Constructor definition

  Options_MRLDP(int p){
    this -> gamma = zeros(p, 1);
    this -> alpha = zeros(p, 1);
    this -> sgga2 = 0.01;
    this -> sgal2 = 0.01;
    this -> beta0 = 0;
    this -> epsStopLogLik = 1e-7;
    this -> maxIter = 10000; 
  }
  
  Options_MRLDP(vec gamma, vec alpha, double sgga2, double sgal2, double beta0, double epsStopLogLik, int maxIter){
    this -> gamma = gamma;
    this -> alpha = alpha;
    this -> sgga2 = sgga2;
    this -> sgal2 = sgal2;
    this -> beta0 = beta0;
    this -> epsStopLogLik = epsStopLogLik;
    this -> maxIter = maxIter;
  }
  vec gamma;
  vec alpha;
  double sgga2;
  double sgal2;
  double beta0;
  double epsStopLogLik;
  int maxIter;
  
};
struct ObjMRLDP{
  int Iteration;
  vec loglik_out;
  double diff;
  double tstat;
  double beta0;
  double sgga2;
  double sgal2;
};

#endif /* MRLDPdefaultInitialValue_hpp */