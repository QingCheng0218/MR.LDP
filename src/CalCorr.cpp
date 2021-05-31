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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]


class DataBlock{
public:
  int chrom;
  int block_no;
  int start;
  int end;
  DataBlock(int chrom, int block_no);
  DataBlock();
  DataBlock(const DataBlock& block);
  fmat corr;
  double heritability_beta;
  double heritability;
};

DataBlock::DataBlock(int chrom, int block_no){
  this->chrom = chrom;
  this->block_no = block_no;
  this->heritability_beta = 0;
}

DataBlock::DataBlock(){

}

DataBlock::DataBlock(const DataBlock& block){
  this->chrom = block.chrom;
  this->block_no = block.block_no;
  this->start = block.start;
  this->end = block.end;
}


// [[Rcpp::export]]
vector<umat> load_block_file(string block_file){
  ivec chrom_vec;

  vector<umat> blocks_mat;
  vector<int *> blocks;
  std::ifstream ifs(block_file.c_str());

  std::string line;
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  std::getline(ifs, line);

  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    i++;
  }

  ifs.close();
  ifs.open(block_file.c_str());
  std::getline(ifs, line);
  chrom_vec.resize(i);
  umat block_mat(i, 2);
  i = 0;
  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::algorithm::split(fields, line, boost::is_any_of("\t"));
    chromsome = atoi(fields[0].substr(3).c_str());
    chrom_vec[i] = chromsome - 1;
    block_mat(i, 0) = atof(fields[1].c_str());
    block_mat(i, 1) = atof(fields[2].c_str());
    i++;
  }

  for (int k = 0; k < 22; k++){
    blocks_mat.push_back(block_mat.rows(find(chrom_vec == k)));
  }
  ifs.close();
  return blocks_mat;
}

void cal_blocks(vector<DataBlock> &datablocks, arma::ivec bp, arma::ivec chr, std::string block_file){
  vector<arma::fmat> corr_blocks;
  vector<umat>  blocks = load_block_file(block_file.c_str());

  Col<int> chrom_list;
  chrom_list = arma::unique(chr);
  ivec chr_index(chrom_list.size());

  double sum2 = 0;
  for(int i = 1; i<=(int)(chrom_list.size()); i++){
    int chr_i = chrom_list[i - 1];
    chr_index[i - 1] = sum(chr==chr_i);
    // chr_index[i - 1] = sum(chr==i);
    sum2 += chr_index[i - 1];
  }

  ivec tmpsum = cumsum(chr_index);
  uword n_chrom = chr_index.size();
  vector<int> block_vector;
  int last_block_no = -1;

  for (uword i = 0; i < n_chrom; i++) {
    ivec u_i = get_interval(tmpsum, i);
    int chri_idx = chrom_list[i] - 1;
    umat block_i = blocks[chri_idx];
    // umat block_i = blocks[i];
    uword n_block = block_i.n_rows;
    int last_chrom;

    for (int kk = 0; kk < (int)(u_i.size()); kk++){
      bool in = false;
      int block_no = -1;
      int k = u_i[kk];
      for (int j = 0; j < (int)(n_block); j++){
        in = (bp[k] <= (int)(block_i(j, 1))) && (bp[k] >= (int)(block_i(j, 0)));
        if (in){
          block_no = j;
          break;
        }
      }
      if (block_no != -1){
        if (block_no < last_block_no){
          // int last_chrom = datablocks[datablocks.size() - 1].chrom;
          if (last_chrom == (int)(i)){
            cout << "Error Order" << endl;
          }
        }
        if (block_no != last_block_no && last_block_no != -1){
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[0];
          block.end = block_vector[block_vector.size() - 1];
          datablocks.push_back(block);
          block_vector.clear();
        }
        block_vector.push_back(k);
        last_block_no = block_no;
      }
      last_chrom = i;
    }
  }

  if (block_vector.size() > 0){
    DataBlock block((int)(n_chrom - 1), last_block_no);
    block.start = block_vector[0];
    block.end = block_vector[block_vector.size() - 1];
    datablocks.push_back(block);
  }

  // Col<int> chrom_vec;
  // chrom_vec.resize(datablocks.size());
  // chrom_vec.zeros();
  // for (int i = 0; i < (int)(datablocks.size()); i++){
  //   int chrom_idx = datablocks[i].chrom;
  //   uvec idx_u = find(chrom_list == chrom_idx + 1);
  //   uword idx = idx_u[0];
  //   chrom_vec[i] = (int)idx;
  // }

}
// [[Rcpp::export]]
List test_blocks(arma::ivec bp, arma::ivec chr, std::string block_file){
  vector<DataBlock> datablocks;
  vector<arma::fmat> corr_blocks;
  vector<umat>  blocks = load_block_file(block_file.c_str());

  Col<int> chrom_list;

  chrom_list = arma::unique(chr);
  ivec chr_index(chrom_list.size());

  double sum2 = 0;
  for(int i = 1; i<=(int)(chrom_list.size()); i++){
    int chr_i = chrom_list[i - 1];
    chr_index[i - 1] = sum(chr==chr_i);
    // chr_index[i - 1] = sum(chr==i);
    sum2 += chr_index[i - 1];
  }
  cout<< "chr_index:"<<chr_index.t()<<endl;
  ivec tmpsum = cumsum(chr_index);
  uword n_chrom = chr_index.size();
  vector<int> block_vector;
  int last_block_no = -1;
  // cout << "check error 0 "<< endl;
  cout << "n_chrom: "<< n_chrom << endl;
  cout <<"tmpsum:"<< tmpsum<<endl;


  for (uword i = 0; i < n_chrom; i++) {
    ivec u_i = get_interval(tmpsum, i);
    // cout << "check error" <<endl;
    // cout <<"u_i:"<<u_i.t()<<endl;
    int chri_idx = chrom_list[i] - 1;
    umat block_i = blocks[chri_idx];
    uword n_block = block_i.n_rows;
    int last_chrom;

    for (int kk = 0; kk < (int)(u_i.size()); kk++){
      bool in = false;
      int block_no = -1;
      int k = u_i[kk];
      for (int j = 0; j < (int)(n_block); j++){
        in = (bp[k] <= (int)(block_i(j, 1))) && (bp[k] >= (int)(block_i(j, 0)));
        if (in){
          block_no = j;
          break;
        }
      }
      if (block_no != -1){
        if (block_no < last_block_no){
          if (last_chrom == chri_idx){
            cout << "Error Order" << endl;
          }
        }
        if (block_no != last_block_no && last_block_no != -1){
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[0];
          block.end = block_vector[block_vector.size() - 1];
          datablocks.push_back(block);
          block_vector.clear();
        }
        block_vector.push_back(k);
        last_block_no = block_no;
      }
      last_chrom = i;

    }
  }

  // cout << "check error 1"<< endl;
  if (block_vector.size() > 0){
    DataBlock block((int)(n_chrom - 1), last_block_no);
    block.start = block_vector[0];
    block.end = block_vector[block_vector.size() - 1];
    datablocks.push_back(block);
  }


  Col<int> chrom_vec;
  chrom_vec.resize(datablocks.size());
  chrom_vec.zeros();


  // for (int i = 0; i < (int)(datablocks.size()); i++){
  //   int chrom_idx = datablocks[i].chrom;
  //   uvec idx_u = find(chrom_list == chrom_idx + 1);
  //   // cout <<"i:" << i << "idx_u: " << idx_u<<endl;
  //   cout<<"chrom_idx :"<<chrom_idx <<endl;
  //   uword idx = idx_u[0];
  //   chrom_vec[i] = (int)idx;
  // }


  // -------------------------------------------------------
  // cout << "check error2 "<< endl;
  uword nblocks = datablocks.size();
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }
  List output = List::create(
    Rcpp::Named("block_inf") = block_inf
  );
  return output;

}


class paraBlock_CorR{

public:
  int current_idx=0;
  uword Ngene_active;
  int n_thread = 1;
  vector<DataBlock> datablocks;
  arma::Mat<unsigned>* X;
  // Chroms* chroms;
  Col<int> chrom_vec;

  uword nblocks;
  field<mat> F4Rblock;
  arma::ivec bp;
  arma::ivec chr;
  arma::uvec avbIndex;
  std::string block_file;
  std::string stringname3;
  double lam;
  arma::uvec Nb;


  paraBlock_CorR(uword &nblocks, vector<DataBlock> &datablocks, arma::Mat<unsigned>* X, field<mat>& F4Rblock, arma::uvec &Nb, double lam){
    this -> nblocks = nblocks;
    this -> datablocks = datablocks;
    this -> X = X;
    this -> F4Rblock = F4Rblock;
    this -> Nb = Nb;
    this -> lam = lam;
  }



  arma::fmat calCorr();
  arma::mat cal_blockcor(Mat<unsigned>& X);

  int  next_gibbs();
  void update_by_thread_gibbs(int thread_id);
  mat loop_by_block_CarCorr(int i);
  void update_by_thread_CarCorr(int thread_id);


};

arma::mat paraBlock_CorR::cal_blockcor(Mat<unsigned>& X){
  uword size1 = X.n_rows;
  uword size2 = X.n_cols;
  vec meanX(size2);
  vec sqrtsum(size2);
  mat Xnorm = conv_to<mat>::from(X);
  for (int i = 0; i < (int)(size2); i++) { //calculate the mean of the vector and sqrt sum
    meanX[i] = sum(Xnorm.col(i))*1.0 / size1;
    Xnorm.col(i) -= meanX[i];
    vec v_i = Xnorm.col(i);
    mat pd = v_i.t() * v_i;
    sqrtsum[i] = sqrt(pd.at(0));
    if (sqrtsum[i] > 0){
      Xnorm.col(i) /= sqrtsum[i];
    }
  }
  arma::mat corr(size2, size2);
  arma::mat eyeI(size2, size2);
  eyeI.eye();
  mat cor_ij(1, 1);
  corr.eye();
  for (int i = 0; i < (int)(size2); i++){
    for (int j = i + 1; j < (int)(size2); j++) {
      cor_ij = ((Xnorm.col(i)).t() * (Xnorm.col(j)));
      double value = cor_ij.at(0);
      corr(i, j) = value;
      corr(j, i) = value;
    }
  }
  //
   // corr *= 0.9;
   // corr += 0.1*eyeI;

  return corr;

}
std::mutex _mtx1;
int paraBlock_CorR::next_gibbs(){
  std::lock_guard<std::mutex> lockGuard(_mtx1);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

mat paraBlock_CorR::loop_by_block_CarCorr(int i){
  Mat<unsigned> subX = X->cols(datablocks[i].start, datablocks[i].end);
  arma::mat corr0 = cal_blockcor(subX);
  arma::mat corr;
  int p1 = corr0.n_rows;
  arma::mat LAM = zeros(p1, p1);
  LAM.fill(lam);
  
  if(lam > 0.5)
  {
    corr = corr0;
    arma::mat eyeI(p1, p1);
    eyeI.eye();
    corr = corr0;
    corr *= lam;
    corr += (1 - lam)*eyeI;
  }
  else
  {
    pdsoftObj out = pdsoft(corr0, LAM);
    corr = out.theta;
  }

  F4Rblock(i, 0) = corr;
  Nb(i) = datablocks[i].end - datablocks[i].start + 1;

  return(corr);
}

void paraBlock_CorR::update_by_thread_CarCorr(int thread_id){
  while(true){
    int idx = next_gibbs();
    if(idx == -1){
      break;
    }
    loop_by_block_CarCorr(idx);
  }
}

// [[Rcpp::export]]
List Cal_blockR(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                std::string stringname3, int coreNum, double lam){

  vector<DataBlock> datablocks;

  cal_blocks(datablocks, bp, chr, block_file);

  string famfile = stringname3;
  famfile +=".fam";
  string bimfile = stringname3;
  bimfile += ".bim";
  int N = getLineNum(famfile);
  int P = getLineNum(bimfile);
  long long Nlong = (long long)N;
  long long Plong = (long long)P;
  unsigned* X0 = new unsigned[Nlong * Plong];
  readPlink(stringname3, N, P, X0);
  arma::Mat<unsigned>* Xdata =  new arma::Mat<unsigned>(X0, N, P, false, false);
  //
  uword M = Xdata -> n_cols;
  arma::Mat<unsigned>* X = Xdata;
  arma::Mat<unsigned> subX0;
  if(avbIndex.size() != 0 && avbIndex.size() < M){
    subX0 = Xdata->cols(avbIndex);
    X = &subX0;
  }
  
  // revise the genotype data for matching with exposure.
  if(idx4panel.n_elem > 0){
    // for(int i = 0; i < idx4panel.n_elem; i++){
    //   X->cols(i) = 2 - X->cols(idx4panel.elem(i));
    // }
    X->cols(idx4panel) = 2 - X->cols(idx4panel);
  }

  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);

  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }

  field<mat> F4Rblock(nblocks, 1);
  // -----------------------------------------------------------------------
  // parallel for correlation matrix
  // set parallele structure object
  paraBlock_CorR parobj(nblocks, datablocks, X, F4Rblock, Nb, lam);
  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_CorR::update_by_thread_CarCorr, &parobj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  F4Rblock = parobj.F4Rblock;
  // ---------------------------------------------------------------------------
  // resize the F4Block
  // mat R;
  // R.zeros(avbIndex.size(), avbIndex.size());
  // 
  // for(int i = 0; i< (int)(nblocks); i++){
  //   R.submat(datablocks[i].start, datablocks[i].start, datablocks[i].end, datablocks[i].end) = F4Rblock(i, 0);
  // }
  // 

  List output = List::create(
    // Rcpp::Named("genotype") = subX0,
    // Rcpp::Named("R") = R,
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("F4Rblock") = F4Rblock,
    Rcpp::Named("Nb") = parobj.Nb,
    Rcpp::Named("nblocks") = nblocks
  );
  return output;
}


// [[Rcpp::export]]
List Cal_block_Rmatrix(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                std::string stringname3, int coreNum, double lam){
  
  vector<DataBlock> datablocks;
  
  cal_blocks(datablocks, bp, chr, block_file);
  
  string famfile = stringname3;
  famfile +=".fam";
  string bimfile = stringname3;
  bimfile += ".bim";
  int N = getLineNum(famfile);
  int P = getLineNum(bimfile);
  long long Nlong = (long long)N;
  long long Plong = (long long)P;
  unsigned* X0 = new unsigned[Nlong * Plong];
  readPlink(stringname3, N, P, X0);
  arma::Mat<unsigned>* Xdata =  new arma::Mat<unsigned>(X0, N, P, false, false);
  //
  uword M = Xdata -> n_cols;
  arma::Mat<unsigned>* X = Xdata;
  arma::Mat<unsigned> subX0;
  if(avbIndex.size() != 0 && avbIndex.size() < M){
    subX0 = Xdata->cols(avbIndex);
    X = &subX0;
  }
  
  // revise the genotype data for matching with exposure.
  if(idx4panel.n_elem > 0){

    X->cols(idx4panel) = 2 - X->cols(idx4panel);
  }
  
  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);
  
  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }
  
  field<mat> F4Rblock(nblocks, 1);
  // -----------------------------------------------------------------------
  // parallel for correlation matrix
  // set parallele structure object
  paraBlock_CorR parobj(nblocks, datablocks, X, F4Rblock, Nb, lam);
  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_CorR::update_by_thread_CarCorr, &parobj, i_thread);
  }
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  F4Rblock = parobj.F4Rblock;
  // ---------------------------------------------------------------------------
  // resize the F4Block
  mat R;
  R.zeros(avbIndex.size(), avbIndex.size());

  for(int i = 0; i< (int)(nblocks); i++){
    R.submat(datablocks[i].start, datablocks[i].start, datablocks[i].end, datablocks[i].end) = F4Rblock(i, 0);
  }

  
  List output = List::create(
    // Rcpp::Named("genotype") = subX0,
    Rcpp::Named("R") = R
    // Rcpp::Named("block_inf") = block_inf,
    // Rcpp::Named("F4Rblock") = F4Rblock,
    // Rcpp::Named("Nb") = parobj.Nb,
    // Rcpp::Named("nblocks") = nblocks
  );
  return output;
}
