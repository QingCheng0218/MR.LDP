rm(list = ls());
library(MR.LDP);
library("MendelianRandomization");
library(mr.raps);
library(gsmr)
# For CAD-CAD study.
filescreen= "heart attack_myocardial infarction.txt";
fileexposure = "c4d.txt";
fileoutcome = "cardiogram.txt";
# -------------------------------------------------
# For Height-Height study.
# filescreen = "UKbiobank_height_used.txt";
# fileexposure = "GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N_used.txt";
# fileoutcome = "GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N_used.txt";
# -------------------------------------------------
stringname3 = "all_chr_1000G";
block_file = "fourier_ls-all.bed";
nrefname3 = paste0(stringname3,".fam");

QC.pa = c(1, 0);
LAM.pa = c(0.1, 0.15);

para.all <- expand.grid(QC.pa, LAM.pa)
for(i in 1:nrow(para.all)){
  assign(paste0("res",i), para.all[i, ])
}

para.iter <- res1;
QCindex <-  0;
lam <- 0.1;

filename = paste0("cad_cad","Lam", lam*100, stringname3, "QC", QCindex, ".Rdata");

# ---------------------------------------------------------------------
# imputed data
pva_cutoff = 2;
scrres = matchscreen(filescreen, fileexposure, fileoutcome, stringname3,  pva_cutoff)

bh1 = scrres$bh1;
bh2 = scrres$bh2;
s12 = scrres$s12;
s22 = scrres$s22;
chr = scrres$chr;
bp = scrres$bp;
rsname = scrres$rsname;
avbIndex = scrres$idxin;
idx4panel = scrres$idx4panel;

screendata = read.table(filescreen, header = T);
expdata = read.table(fileexposure, header = T);
outdata = read.table(fileoutcome, header = T);

# ------------------------------------------------------------------------------
# Select instrumental SNPs using screening data under different thresholds
idx = match(rsname, screendata$SNP);
Xdata0 = screendata[idx, ];


PVA = c(5e-06, 1e-05, 5e-05, 1e-04, 5e-04, 1e-03, 5e-03, 1e-02)
pres = rep(1, length(PVA));
bres = rep(1, length(PVA));
Tstat = rep(1, length(PVA));
bhat = se = tstatres = matrix(0, ncol = 3, nrow = length(PVA));
MHCpres = Xpres = Ypres = rep(1, length(PVA));


LDR2 = c(0.001, 0.05);
Bres = Sres = array(0, dim = c(length(LDR2), 5, length(PVA)));
Nused = Nldres = Pres = matrix(0, nrow = length(LDR2), ncol = length(PVA));
pres = rep(0, length(PVA));

for(irep in 1:length(PVA)){
  
  Stm <- proc.time();
  
  pvalue = PVA[irep];
  idx1 = which(Xdata0$pvalue < pvalue);
  Xdata = Xdata0[idx1, ];
  idx0 = match(Xdata$SNP, rsname);
  # ------------------------------------------------------------------------------
  bh1_in = bh1[idx0];
  bh2_in = bh2[idx0];
  s12_in = s12[idx0];
  s22_in = s22[idx0];
  bp_in = bp[idx0];
  chr_in = chr[idx0];
  avbIndex_in = avbIndex[idx0];
  rsname_in = rsname[idx0];
  idx4panel_in = match(avbIndex[intersect((idx4panel + 1), idx0)], avbIndex_in) -1;
  # ------------------------------------------------------------------------------
  if(QCindex){
    QCresult = summaryQC(mhcstart, mhcend, bh1_in, bh2_in, s12_in, s22_in, bp_in,
                         chr_in, rsname_in, avbIndex_in, idx4panel_in,Inf, Inf)
    bh1new = QCresult$bh1new;
    bh2new = QCresult$bh2new;
    s12new = QCresult$s12new;
    s22new = QCresult$s22new;
    bpnew = QCresult$bpnew;
    chrnew = QCresult$chrnew;
    avbIndexnew = QCresult$avbIndexnew;
    idx4panelnew = QCresult$idx4panelnew;
    rsnamenew = QCresult$rsnamenew;
    pmhc = QCresult$pmhc;
    px = QCresult$px;
    py = QCresult$py;
  }else{
    bh1new = bh1_in;
    bh2new = bh2_in;
    s12new = s12_in;
    s22new = s22_in;
    bpnew = bp_in;
    chrnew = chr_in;
    rsnamenew = rsname_in;
    idx4panelnew = idx4panel_in;
    avbIndexnew = avbIndex_in;
    pmhc = 0;
    px = 0;
    py = 0;
  }
  
  p = length(avbIndexnew);
  
  MHCpres[irep] = pmhc;
  Xpres[irep] = px;
  Ypres[irep] = py;
  # ------------------------------------------------------------------------------
  
  gamma = rep(0.01, p);
  alpha = rep(0.01, p);
  
  coreNum = 24;
  # fit = lm(bh2_in~bh1_in)
  # limfit = fit$coefficients[2]
  
  sgga2 =  0.01;
  sgal2 =  0.01;
  beta0 = 1;
  IterMax = 10000
  N2 = rep(1000, p);
  coreNum = 24;
  
  
  pres[irep] = p;
  epsStopLogLik = 1e-6;
  
  if(p >2){
    result0 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                                   gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 0, epsStopLogLik, IterMax, model = 1)
    
    result1 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                                   gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 1, epsStopLogLik, IterMax, model = 1)
    
    
    Tstat = 2*(result0$tstat - result1$tstat)
    bhat[irep, 1] = result0$beta0;
    se[irep, 1] = bhat[irep, 1]/sqrt(Tstat);
    tstatres[irep, 1] = Tstat;
    
    #-------------------------------------------------------------------------------------
    result0 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                                   gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 0, epsStopLogLik, IterMax, model = 2)
    
    result1 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                                   gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 1, epsStopLogLik, IterMax, model = 2)
    
    
    Tstat = 2*(result0$tstat - result1$tstat)
    bhat[irep, 2] = result0$beta0;
    se[irep, 2] = bhat[irep, 1]/sqrt(Tstat);
    tstatres[irep, 2] = Tstat;
    
    #-------------------------------------------------------------------------------------
    # The linear regression method.
    fit = summary(lm(bh2new~bh1new))
    bhat[irep, 3] <- fit$coefficients[2, 1];
    se[irep, 3] <- fit$coefficients[2, 2];
    #-------------------------------------------------------------------------------------
    
  }
  
  
  # ------------------------------------------------------------------------------
  # obtain the pvalue of exposure data and outcome data
  id_x = match(rsnamenew, expdata$SNP);
  exp_in = expdata[id_x, ];
  pval_x = exp_in$pvalue;
  id_y = match(rsnamenew, outdata$SNP);
  out_in = outdata[id_y, ];
  pval_y = out_in$pvalue;
  # ------------------------------------------------------------------------------
  
  
  Rblockres = Cal_block_Rmatrix(bpnew, chrnew, avbIndexnew-1,idx4panel, block_file, stringname3, coreNum, lam);
  
  Rnew = Rblockres$R
  n_ref = getLineNum(nrefname3);

  colnames(Rnew) = rsnamenew;
  rownames(Rnew) = rsnamenew;
  snp_coeff_id = rsnamenew;
  
  
  gwas_thresh = 1   # GWAS threshold to select SNPs as the instruments for the GSMR analysis
  single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
  multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
  nsnps_thresh = 1   # the minimum number of instruments required for the GSMR analysis
  heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
  
  ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
  gsmr_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
  
  bhat5 = se5 = matrix(0, ncol = 5, nrow = length(LDR2))
  # bhat = se = nbres = rep(0, length(LDR2));
  nld = nused = rep(0, length(LDR2));
  
  for(rr in 1:length(LDR2)){
    
    ld_r2_thresh = LDR2[rr]   # LD r2 threshold to remove SNPs in high LD
    
    gsmr_results = gsmr(bh1new, s12new, pval_x, bh2new, s22new, pval_y,
                        Rnew, snp_coeff_id, n_ref, heidi_outlier_flag, 
                        gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, 
                        nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta) 
    
    
    bhat5[rr, 1] = gsmr_results$bxy;
    se5[rr, 1] = gsmr_results$bxy_se;
    nld[rr] = p - length(gsmr_results$linkage_snps);
    nused[rr] = length(gsmr_results$used_index);
    
    id00 = match(gsmr_results$linkage_snps, rsnamenew);
    idx = 1:p;
    id4ld = idx[-id00];
    
    bh14ld = bh1_in[id4ld];
    bh24ld = bh2_in[id4ld];
    s124ld = s12_in[id4ld];
    s224ld = s22_in[id4ld];
    
    #-------------------------------------------------------------------------------------
    # The linear regression method.
    fit = summary(lm(bh24ld~bh14ld))
    bhat5[rr, 2] <- fit$coefficients[2, 1];
    se5[rr, 2] <- fit$coefficients[2, 2];
    #-------------------------------------------------------------------------------------
    # The Raps method.
    raps = mr.raps(bh14ld, bh24ld, s124ld, s224ld);
    bhat5[rr, 3] <- raps$beta.hat;
    se5[rr, 3] <- raps$beta.se;
    #-------------------------------------------------------------------------------------
    # The inverse-variance weighted method.
    ivm.res <- mr_ivw(mr_input(bx = bh14ld, bxse = s124ld, by = bh24ld, byse = s224ld));
    bhat5[rr, 4] <- ivm.res$Estimate;
    se5[rr, 4] <- ivm.res$StdError;
    # -------------------------------------------------------------------------------------
    # # The MR-Egger method.
    mregger.res <- mr_egger(mr_input(bx = bh14ld, bxse = s124ld, by = bh24ld, byse = s224ld));
    bhat5[rr, 5] <- mregger.res$Estimate;
    se5[rr, 5] <- mregger.res$StdError.Est;
  }
  
  colnames(bhat5) = colnames(se5) = c("GSMR", "LS","Raps","IVW","MREgger");
  Bres[, , irep] = bhat5;
  Sres[, , irep] = se5;
  Nused[, irep] = nused;
  Nldres[, irep] = nld;
  
  
  Ftm <- proc.time();
  S.time <- Ftm[3] - Stm[3];
  print(Sys.time())
  cat( irep,"/",length(PVA)," time: ",S.time,"\n",sep="");

  save(PVA, pres, MHCpres, Xpres, Ypres, bhat, se, tstatres,
       Bres, Sres, pres, Nused, Nldres,file = filename);
  
}
