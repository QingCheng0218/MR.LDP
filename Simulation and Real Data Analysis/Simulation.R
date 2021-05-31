library("MR.LDP");
library("mvtnorm");
library("ggplot2");
library("PDSCE");
library("mr.raps")
library("MendelianRandomization")
# ----------------------------------------------------------------
rho <- 0.4;
h2z <- 0.05;
L <- 10
b0 <- 0.1; # b0 = 0.1 for point estimate; b0 = 0 for power.
lam <- 0.055;
n3 <- 500;
Alrate <- 1;
Garate <- 1;
h2y <- 0.1;
M = 50;
p = M*L;
n1 <- 20000;
n2 <- 20000;
m <- p;
bw <- ifelse(m>p,p,m);

filename <- paste0("r", 10*rho,  "hz", h2z*100, "hy", h2y*10,
                   "L", L, "b", b0, "lam", lam*100,"n3",
                   n3,"alrat", Alrate*10 , ".Rdata");

if(b0!=0){nrep = 100;}else{nrep = 500}
coreNum = 20;
epsStopLogLik = 1e-7;
maxIter = 10000;
# ----------------------------------------------------------------
set.seed(2019);
# # -------------------------------------------------------------------
# Catergorize x into genoypte data 0, 1, 2
genRawGeno <- function(maf, L, M, rho, n){
  SIGMA = matrix(nrow=M,ncol=M)
  for (i in 1:M){
    for (j in 1:M){
      SIGMA[i,j] = rho^(abs(i-j));
    }
  }
  
  nsnp = L*M;
  X = NULL;
  for ( l in 1:L ){
    #set.seed(1000);
    index = (M*(l-1)+1): (M*l);
    AAprob = maf[index]^2.;
    Aaprob = 2*maf[index]*(1-maf[index]);
    quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),M,2);
    Xt = rmvnorm(n, mean=rep(0,M), sigma=SIGMA, method="chol")
    Xt2 = matrix(0,n,M);
    for (j in 1:M){
      cutoff = qnorm(quanti[j,]);
      Xt2[Xt[,j] < cutoff[1],j] = 0;
      Xt2[Xt[,j] >= cutoff[1] & Xt[,j] < cutoff[2],j] = 1;  ## attention
      Xt2[Xt[,j] >= cutoff[2],j] = 2;
    }
    X <- cbind(X,Xt2);
  }
  return(X)
}
p = M*L;
maf = runif(p,0.05,0.5);
x = genRawGeno(maf, L, M, rho, n1 + n2 + n3);
# -------------------------------------------------------------------
x1 = x[1:n1,];
x2 = x[(n1+1):(n1+n2),];
x12 = x[1:(n1+n2),];
x3 = x[(n1+n2+1):(n1+n2+n3),];

R0 <- cor(x3);
R = pdsoft(R0, lam)$theta;
diag(R) <- rep(1, p);
mask = kronecker(diag(L),matrix(1,M,M));
R = R*mask;
# ------------------------------------------------------------
S.time <- H2g.res <- H2a.res <- rep(0, nrep);
irep = 1;
bhat3 = matrix(nrow=nrep, ncol=3);
se3 = matrix(nrow=nrep, ncol=3);
bhatNp3 = seNp3 = matrix(nrow=nrep, ncol=3);

bres = bres2 = matrix(0, ncol = 2, nrow = nrep);
tstat = tstat2 = matrix(0, ncol = 2, nrow = nrep);
Iterres = Iterres2 = matrix(0, ncol = 2, nrow = nrep);
Diff = Diff2 = matrix(0, ncol = 2, nrow = nrep);
H2a.res <- rep(0, nrep);
H2g.res <- rep(0, nrep);

if(rho==0.4){
  id4ld <- seq(1, p, 5);
  m = p;
  index = id4ld;
}else if(rho==0.8){
  id4ld <- seq(1, p, 20);
  m = p;
  index = id4ld;
}else if(rho==0){
  index = 1:m;
}

R2.pa = 0.05
bhatM1 = seM1 = nldM1 = nusedM1 = matrix(0, nrow = nrep, ncol = length(R2.pa));
bhatM2 = seM2 = nldM2 = nusedM2 = matrix(0, nrow = nrep, ncol = length(R2.pa));
#-------------------------------------------------------------------------#

for(irep in 1:nrep){
  
  set.seed(irep);
  Stm <- proc.time();
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  
  # ------------------------------------------------------------------------
  sigma2g = 0.005;
  gano = floor(p * Garate);
  gamma.nz = rnorm(gano)*sqrt(sigma2g);
  indxGA = sample(1:p,gano);
  gamma = numeric(p);
  gamma[indxGA] = gamma.nz;
  
  Su = matrix(c(1,0.8,0.8,1),nrow=2)
  bu = rmvnorm(q,mean=rep(0,2), sigma = Su,method="chol")
  by = bu[,1]; bz = bu[,2];
  uby = u%*%by; ubz = u%*%bz;
  uby = uby/sqrt(as.numeric(var(uby)/0.6));
  ubz = ubz/sqrt(as.numeric(var(ubz)/0.2));
  
  x12g = x12%*%gamma;
  
  if(b0!=0){
    h2ga = (h2y *( 1 + b0^2))/(b0^2 * (1 - h2y));
    gamma0 = gamma/sqrt(as.numeric(var(x12g)/h2ga));
    x12g = x12%*%gamma0;
  }
  
  
  yall = x12g + uby + rnorm(n1+n2)*as.numeric(sqrt(1-var(uby)));
  # ------------------------------------------------------------------------
  # The direct effects on Z
  h2yb = var(b0*yall);
  h2al = (h2z + h2z*h2yb)/(1 - h2z)
  
  sigma2a <- 0.005;
  if(h2z==0){
    alpha0 = rep(0, m);
    x12a = x12%*%alpha0;
  }else{
    alno = floor(p*Alrate);
    alpha.nz <- rnorm(alno)*sqrt(sigma2a);
    # sparse setting for pleiotropy
    indxAL = sample(1:p,alno);
    alpha = numeric(p);
    alpha[indxAL] = alpha.nz;
    
    x12a = x12%*%alpha;
    alpha0 = alpha/sqrt(as.numeric(var(x12a)/(h2al)));
    x12a = x12%*%alpha0;
  }
  # ------------------------------------------------------------------------
  #fm0 = lmm_pxem(y,matrix(rep(1,n1),ncol=1),x1, 100)
  resz = ubz + rnorm(n1+n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b0*yall  + x12a +  resz;
  H2a.res[irep] <- var(x12a)/var(zall);
  H2g.res[irep] <- var(b0*x12g)/var(zall);
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  # create summary statistics
  gammahall = numeric(p); Gammahall = numeric(p);
  sg2all = numeric(p); sG2all = numeric(p);
  pval = numeric(p);
  for (i in 1:p){
    fm = lm(y~1+x1[,i]);
    gammahall[i] = summary(fm)$coefficients[2,1];
    sg2all[i] = summary(fm)$coefficients[2,2]^2;
    pval[i] = summary(fm)$coefficients[2,4];
    
    fm = lm(z~1+x2[,i]);
    Gammahall[i] = summary(fm)$coefficients[2,1];
    sG2all[i] = summary(fm)$coefficients[2,2]^2;
  }
  
  
  sig.indx = index#[pval<0.05/p]
  sg2 = sg2all[sig.indx];
  sG2 = sG2all[sig.indx];
  gammah = gammahall[sig.indx];
  Gammah = Gammahall[sig.indx];
  # plot(gammahall,Gammahall)
  # summary(lm(Gammahall~gammahall))
  # -------------------------------------------------------------------------------------
  se1 = sqrt(sg2); se2 = sqrt(sG2);
  se1all = sqrt(sg2all); se2all = sqrt(sG2all);
  # # # #-------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------
  # initial value.
  beta0 <- 0;
  sgga2 <- 1;
  agm <- bgm <- aal <- bal <- 0.001;
  gamma <- rep(0.01, p);
  alpha <- rep(0.01, p);
  sgga2 <- 1;
  sgal2 <- 1;
  nblocks = L;
  block_inf <- cbind(seq(0, p-M, M), seq(M-1, p, M));
  
  # -------------------------------------------------------------------------------------
  
  pxresult = MRLDP_SimPXvb(gammahall, Gammahall, se1all, se2all, gamma, alpha,  beta0, sgga2, sgal2, R,
                           0, epsStopLogLik, maxIter, model = 1);
  tstat[irep, 1] = pxresult$tstat;
  bres[irep, 1] = pxresult$beta0;
  Iterres[irep, 1] = pxresult$Iteration;
  Diff[irep, 1] = pxresult$diff;
  
  if(b0==0){
    pxresult1 = MRLDP_SimPXvb(gammahall, Gammahall, se1all, se2all, gamma, alpha,  beta0, sgga2, sgal2, R,
                              1, epsStopLogLik, maxIter, model = 1);
    tstat[irep, 2] = pxresult1$tstat;
    bres[irep, 2] = pxresult1$beta0;
    Iterres[irep, 2] = pxresult1$Iteration;
    Diff[irep, 2] = pxresult1$diff;
    
  }
  
  pxresult = MRLDP_SimPXvb(gammahall, Gammahall, se1all, se2all, gamma, alpha,  beta0, sgga2, sgal2, R,
                           0, epsStopLogLik, maxIter, model = 2);
  tstat2[irep, 1] = pxresult$tstat;
  bres2[irep, 1] = pxresult$beta0;
  Iterres2[irep, 1] = pxresult$Iteration;
  Diff2[irep, 1] = pxresult$diff;
  
  if(b0==0){
    pxresult1 = MRLDP_SimPXvb(gammahall, Gammahall, se1all, se2all, gamma, alpha,  beta0, sgga2, sgal2, R,
                              1, epsStopLogLik, maxIter, model = 2);
    tstat2[irep, 2] = pxresult1$tstat;
    bres2[irep, 2] = pxresult1$beta0;
    Iterres2[irep, 2] = pxresult1$Iteration;
    Diff2[irep, 2] = pxresult1$diff;
    
  }
  
  
  # # The RAPs method.
  if(h2z==0){
    raps = mr.raps(gammah, Gammah, se1, se2);
  }else if(h2z!=0){
    raps = mr.raps(gammah, Gammah, se1, se2, over.dispersion = TRUE);
  }
  
  bhat3[irep,1] = raps$beta.hat;
  se3[irep,1] = raps$beta.se;
  # #-------------------------------------------------------------------------------------
  # # The inverse-variance weighted method.
  ivm.res <- mr_ivw(mr_input(bx = gammah, bxse = se1, by = Gammah, byse = se2));
  bhat3[irep, 2] <- ivm.res$Estimate;
  se3[irep, 2] <- ivm.res$StdError;
  # # -------------------------------------------------------------------------------------
  # # # The MR-Egger method.
  mregger.res <- mr_egger(mr_input(bx = gammah, bxse = se1, by = Gammah, byse = se2));
  bhat3[irep, 3] <- mregger.res$Estimate;
  se3[irep, 3] <- mregger.res$StdError.Est;
  # -------------------------------------------------------------------------------------
  
  # # The RAPs method for robust.
  if(h2z==0){
    raps = mr.raps(gammahall, Gammahall, se1all, se2all);
  }else if(h2z!=0){
    raps = mr.raps(gammahall, Gammahall, se1all, se2all, over.dispersion = TRUE);
  }
  
  bhatNp3[irep,1] = raps$beta.hat;
  seNp3[irep,1] = raps$beta.se;
  # #-------------------------------------------------------------------------------------
  # # The inverse-variance weighted method.
  ivm.res <- mr_ivw(mr_input(bx = gammahall, bxse = se1all, by = Gammahall, byse = se2all));
  bhatNp3[irep, 2] <- ivm.res$Estimate;
  seNp3[irep, 2] <- ivm.res$StdError;
  # # -------------------------------------------------------------------------------------
  # # # The MR-Egger method.
  mregger.res <- mr_egger(mr_input(bx = gammahall, bxse = se1all, by = Gammahall, byse = se2all));
  bhatNp3[irep, 3] <- mregger.res$Estimate;
  seNp3[irep, 3] <- mregger.res$StdError.Est;
  # -------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------
  bh1_in = gammah;
  bh2_in = Gammah;
  s12_in = sqrt(sg2);
  s22_in = sqrt(sG2);
  # GSMR method.
  R_in = R
  n_ref = n3
  
  colnames(R_in) = 1:p;
  rownames(R_in) = 1:p;
  snp_coeff_id = 1:p;
  
  
  gwas_thresh = 10   # GWAS threshold to select SNPs as the instruments for the GSMR analysis
  single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
  multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
  nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
  ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
  gsmr_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development
  
  for(rr in 1:length(R2.pa)){
    ld_r2_thresh = R2.pa[rr];
    if(h2z==0){
      heidi_outlier_flag = F;
      gsmrM1 = try(gsmr(bh1_in, s12_in, pval_x, bh2_in, s22_in, pval_y,
                        R_in, snp_coeff_id, n_ref, heidi_outlier_flag,
                        gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh,
                        nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta));
      if(class(gsmrM1) != 'try-error'){
        bhatM1[irep, rr] = gsmrM1$bxy;
        seM1[irep, rr] = gsmrM1$bxy_se;
        nldM1[irep, rr] = p - length(gsmrM1$linkage_snps);
        nusedM1[irep, rr] = length(gsmrM1$used_index);
        
        filename <- paste0("r", 10*rho,  "hz", h2z*100, "hy", h2y*10, "L", L, "b", b0, "lam",
                           lam*100,"n3", n3,"ld",10*ld_r2_thresh,"alrat", Alrate*10, "M1.Rdata");
        save(bhatM1, seM1, nldM1, nusedM1, file = filename);
      }
      
      
      heidi_outlier_flag = T;
      gsmrM2 = try(gsmr(bh1_in, s12_in, pval_x, bh2_in, s22_in, pval_y,
                        R_in, snp_coeff_id, n_ref, heidi_outlier_flag,
                        gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh,
                        nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta));
      if(class(gsmrM2) != 'try-error'){
        bhatM2[irep, rr] = gsmrM2$bxy;
        seM2[irep, rr] = gsmrM2$bxy_se;
        nldM2[irep, rr] = p - length(gsmrM2$linkage_snps);
        nusedM2[irep, rr] = length(gsmrM2$used_index);
        filename <- paste0("r", 10*rho,  "hz", h2z*100, "hy", h2y*10, "L", L, "b", b0, "lam",
                           lam*100,"n3", n3,"ld",10*ld_r2_thresh,"alrat", Alrate*10, "M2.Rdata");
        
        save(bhatM2, seM2, nldM2, nusedM2, file = filename);
      }
    }else{
      heidi_outlier_flag = T;
      gsmrM2 = try(gsmr(bh1_in, s12_in, pval_x, bh2_in, s22_in, pval_y,
                        R_in, snp_coeff_id, n_ref, heidi_outlier_flag,
                        gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh,
                        nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta));
      if(class(gsmrM2)!="try-error"){
        bhatM2[irep, rr] = gsmrM2$bxy;
        seM2[irep, rr] = gsmrM2$bxy_se;
        nldM2[irep, rr] = p - length(gsmrM2$linkage_snps);
        nusedM2[irep, rr] = length(gsmrM2$used_index);
        filename <- paste0("r", 10*rho,  "hz", h2z*100, "hy", h2y*10, "L", L, "b", b0, "lam",
                           lam*100,"n3", n3,"ld",10*ld_r2_thresh,"alrat", Alrate*10,"M2.Rdata");
        
        save(bhatM2, seM2, nldM2, nusedM2, file = filename);
      }
      
    }
  }
  
  # -------------------------------------------------------------------------------------
  Ftm <- proc.time();
  S.time[irep] <- Ftm[3] - Stm[3];
  #print(Sys.time())
  if(irep%%5==0){
    print(Sys.time())
    cat(irep,"/",nrep," time: ",S.time[irep],"\n",sep="")
    save(bres, bres2, tstat, tstat2, Iterres, Iterres2, Diff, Diff2,
         bhat3, se3, bhatNp3, seNp3, H2a.res, H2g.res, S.time, file = filename);
  }
  
  # -------------------------------------------------------------------------------------
  #
}
