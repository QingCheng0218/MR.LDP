#ifndef CalCorr_hpp
#define CalCorr_hpp

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
#include "ReadGeneFile.hpp"



using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost;

List Cal_blockR(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                std::string stringname3, int coreNum, double lam);

List Cal_block_Rmatrix(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                       std::string stringname3, int coreNum, double lam);

#endif /* CalCorr_hpp */