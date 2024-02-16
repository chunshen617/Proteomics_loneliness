#Two-sample MR, using loneliness to proteins as an example
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(foreach)
library(doParallel)

gwas_path <- '/home1/prot_GCTA_result/'

#GWAS summary statistics of loneliness (calculated by PLINK2)
LO_gwas <- fread('LO_GWAS_plink2.csv')
#significant threshold
LO_gwas_sig <- LO_gwas[LO_gwas$P<1e-6,]
#revise SNP name
LO_snpN <- strsplit(LO_gwas_sig$ID,':')
LO_gwas_sig$SNP <- unlist(lapply(LO_snpN,function(x) x[1]))
LO_gwas_sig$Phenotype <- 'LO'
LO_gwas_sig$BETA <- log(LO_gwas_sig$OR)

LO_gwas_sig$A2 <- ifelse(LO_gwas_sig$A1==LO_gwas_sig$REF,'ALT','REF')
LO_gwas_sig$noneffect[LO_gwas_sig$A2=='ALT'] <- LO_gwas_sig$ALT[LO_gwas_sig$A2=='ALT']
LO_gwas_sig$noneffect[LO_gwas_sig$A2=='REF'] <- LO_gwas_sig$REF[LO_gwas_sig$A2=='REF']

#format
LO_exp_data <- format_data(LO_gwas_sig, type = "exposure",
                           phenotype_col = "Phenotype",
                           snp_col = "SNP",
                           beta_col = "BETA",
                           se_col = "LOG(OR)_SE",
                           effect_allele_col = "A1",
                           other_allele = 'noneffect',
                           eaf_col = "A1_FREQ",
                           pval_col = "P",
                           chr_col = "CHR",
                           pos_col = "POS")

#clump
LO_snp_clump <- ld_clump_local(data.frame(rsid = LO_exp_data$SNP,pval = LO_exp_data$pval.exposure),
                             clump_kb=10000,
                             clump_r2=0.001,
                             clump_p = 1,
                             bfile = '/home1/shenchun/DATA/Genome1000_EUR/EUR',
                             plink_bin = '/home1/shenchun/software/plink_1.9/plink')

LO_exp_clump <- LO_exp_data[which(LO_exp_data$SNP %in% LO_snp_clump$rsid),]

#Exposure
#significant proteins associated with loneliness
Result_LO_m2 <- fread('Result_LO_protein_assoc_M2.csv',data.table = F)
bonf_LO_m2 <- Result_LO_m2[Result_LO_m2$pval<0.025/nrow(Result_LO_m2),]
protName <- bonf_LO_m2$V1

multiResultClass <- function(prot_mr_res=NULL,prot_het_res=NULL,prot_ple_res=NULL){
  me <- list(
    prot_mr_res = prot_mr_res,
    prot_het_res = prot_het_res,
    prot_ple_res = prot_ple_res
  )

  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

registerDoParallel(20)

result_mr_all <- foreach(i=1:length(protName)) %dopar% {

  prot <- protName[i]

  #read protein gwas results
  protFile <- paste0(prot,'_assoc.fastGWA')
  prot_gwas <- fread(paste0(gwas_path,protFile))

  #revise SNP name
  snpN <- strsplit(prot_gwas$SNP,':')
  snpN2 <- unlist(lapply(snpN,function(x) x[1]))
  #extract outcome data based on instruments
  prot_gwas1 <- prot_gwas[which(snpN2 %in% LO_exp_clump$SNP),]
  prot_gwas1$SNP2 <- snpN2[which(snpN2 %in% LO_exp_clump$SNP)]
  prot_gwas1$Phenotype <- prot

  out_data1 <- format_data(prot_gwas1, type = "outcome",
                           phenotype_col = "Phenotype",
                           snp_col = "SNP2",
                           beta_col = "BETA",
                           se_col = "SE",
                           eaf_col = "AF1",
                           effect_allele_col = "A1",
                           other_allele = 'A2',
                           pval_col = "P",
                           chr_col = "CHR",
                           pos_col = "POS")

  #harmonise data
  dat1 <- harmonise_data(exposure_dat = LO_exp_clump,
                         outcome_dat = out_data1)

  #perform MR
  result <- multiResultClass()
  result$prot_mr_res <- mr(dat1,method_list = c("mr_wald_ratio","mr_ivw","mr_weighted_median","mr_egger_regression"))
  #Sensitivity analyses
  #Heterogeneity statistics
  result$prot_het_res <- mr_heterogeneity(dat1)
  #Horizontal pleiotropy
  result$prot_ple_res <- mr_pleiotropy_test(dat1)

  return(result)
}

save(result_mr_all, file = 'Result_MR_loTOprotein.RData')
