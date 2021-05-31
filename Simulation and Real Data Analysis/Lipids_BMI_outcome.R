rm(list = ls());
library(MR.LDP)
library(mr.raps);
library(gsmr);
library("MendelianRandomization");

# -------------------------------------------------
# BMI-outcome analysis
filescreen= "All_BMI_BBJ2017.txt";
fileexposure = "SNP_gwas_mc_merge_nogc.tbl.uniq_matchpa.txt";
fileoutcome = "cardiogram.txt";
# -------------------------------------------------
# # Lipids-outcome analysis
# filescreen= "HDL2010.txt"; # "LDL2010.txt"; "TC2010.txt";
# fileexposure = "HDL2013.txt"; # "LDL2013.txt"; "TC2013.txt";
# fileoutcome = "cardiogram.txt";
# -------------------------------------------------
stringname3 = "F:/share/Nutstore/Model3/overlap/RealData/UK10K_hm3_EUR_QC_group1";
block_file = "fourier_ls-all.bed";
nrefname3 = paste0(stringname3,".fam");
filename = paste0("BmiYLam", lam*100, name3, ".Rdata");
# ---------------------------------------------------------------------
pva_cutoff <- 1*10^(-4);
lam =  0.1;
QCindex = 0;
# ------------------parameters for gsmr-----------------#
R2.pa = c(0.001, 0.05)
gwas_thresh = 1  # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis

ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
# ------------------parameters for gsmr-----------------#
RESbeta = RESse = matrix(0, ncol = 5, nrow = length(R2.pa));
Nldres = Ngsmr = rep(0, length(R2.pa));
scrres = matchscreen(filescreen, fileexposure,  fileoutcome, stringname3,  pva_cutoff)
bh1 = as.numeric(scrres$bh1);
bh2 = as.numeric(scrres$bh2);
s12 = as.numeric(scrres$s12);
s22 = as.numeric(scrres$s22);
chr = as.numeric(scrres$chr);
bp = scrres$bp;
rsname = scrres$rsname
avbIndex = scrres$idxin;
p = length(avbIndex)
idx4panel = scrres$idx4panel;

# ------------------------------------------------------------------------------

if(QCindex){
  QCresult = summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp,
                       chr, rsname, avbIndex, idx4panel, Inf, Inf)
  bh1new = QCresult$bh1new;
  bh2new = QCresult$bh2new;
  s12new = QCresult$s12new;
  s22new = QCresult$s22new;
  bpnew = QCresult$bpnew;
  chrnew = QCresult$chrnew;
  avbIndexnew = QCresult$avbIndexnew;
  rsnamenew = QCresult$rsnamenew;
  idx4panelnew = QCresult$idx4panelnew;
  
  pmhc = QCresult$pmhc;
  px = QCresult$px;
  py = QCresult$py;
  p = length(avbIndexnew);
}else{
  bh1new = bh1;
  bh2new = bh2;
  s12new = s12;
  s22new = s22;
  bpnew = bp;
  chrnew = chr;
  rsnamenew = rsname;
  idx4panelnew = idx4panel;
  avbIndexnew = avbIndex;
  pmhc = 0;
  px = 0;
  py = 0;
}



# ------------------------------------------------------------------------------


gamma = rep(0.01, p);
alpha = rep(0.01, p);


sgga2 =  0.01;
sgal2 =  0.01;
beta0 = 1;
model = 2
IterMax = 10000
coreNum = 24;
epsStopLogLik = 1e-6;



result0 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                               gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 0, epsStopLogLik, IterMax, model = 2)

result1 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, bh1new, bh2new, s12new, s22new,
                               gamma, alpha, beta0, sgga2, sgal2, coreNum, lam, 1, epsStopLogLik, IterMax, model = 2)


Iter0 <- result0$Iteration;
Iter1 <- result1$Iteration;
bhat1 = result0$beta0;
Tstat <- 2*(result0$tstat - result1$tstat);
se1 = abs(bhat1/sqrt(Tstat))

#-------------------------------------------------------------------------------------
expdata = read.table(fileexposure, header = T);
outdata = read.table(fileoutcome, header = T);

# ------------------------------------------------------------------------------
# obtain the pvalue of exposure data and outcome data
id_x = match(rsnamenew, expdata$SNP);
exp_in = expdata[id_x, ];
pval_x = exp_in$pvalue;
id_y = match(rsnamenew, outdata$SNP);
out_in = outdata[id_y, ];
pval_y = out_in$pvalue;
# ------------------------------------------------------------------------------
Rblockres = Cal_block_Rmatrix(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file, stringname3, coreNum, lam);

Rnew = Rblockres$R
n_ref = getLineNum(nrefname3);

colnames(Rnew) = rsnamenew;
rownames(Rnew) = rsnamenew;
snp_coeff_id = rsnamenew;


for(rr in 1:length(R2.pa)){
  ld_r2_thresh = R2.pa[rr];
  gsmr_results = try(gsmr(bh1new, s12new, pval_x, bh2new, s22new, pval_y,
                          Rnew, snp_coeff_id, n_ref, heidi_outlier_flag, 
                          gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, 
                          nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta));
  if(class(gsmr_results) != 'try-error'){
    bhat2 = gsmr_results$bxy;
    se2 = gsmr_results$bxy_se;
    Nldres[rr] = p - length(gsmr_results$linkage_snps);
    Ngsmr[rr] = length(gsmr_results$used_index);
    
    id00 = match(gsmr_results$linkage_snps, rsname);
    idx = 1:p;
    id4ld = idx[-id00];
    
    bh14ld = bh1[id4ld];
    bh24ld = bh2[id4ld];
    s124ld = s12[id4ld];
    s224ld = s22[id4ld];
    #-------------------------------------------------------------------------------------
    # The Raps method.
    raps = mr.raps(bh14ld, bh24ld, s124ld, s224ld, over.dispersion = TRUE);
    bhat3 <- raps$beta.hat;
    se3 <- raps$beta.se;
    #-------------------------------------------------------------------------------------
    # The inverse-variance weighted method.
    ivm.res <- mr_ivw(mr_input(bx = bh14ld, bxse = s124ld, by = bh24ld, byse = s224ld));
    bhat4 <- ivm.res$Estimate;
    se4 <- ivm.res$StdError;
    # -------------------------------------------------------------------------------------
    # # The MR-Egger method.
    mregger.res <- mr_egger(mr_input(bx = bh14ld, bxse = s124ld, by = bh24ld, byse = s224ld));
    bhat5 <- mregger.res$Estimate;
    se5<- mregger.res$StdError.Est;
    brestmp = c(bhat1, bhat2, bhat3, bhat4, bhat5);
    setmp = c(se1, se2, se3, se4, se5);
    RESbeta[rr, ] = brestmp;
    RESse[rr, ] = setmp;
    
  }
  
}


save(RESbeta, RESse, Nldres, Ngsmr, 
     p,  fileYname, fileX, n_ref, file = filename);

