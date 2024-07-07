library(data.table)
library(ieugwasr)
library(cause)

#GWAS summary statistics of loneliness (calculated by PLINK2)
LO_gwas <- fread('LO_GWAS_plink2.csv')
#revise SNP name
LO_snpN <- strsplit(LO_gwas$ID,':')
LO_gwas$SNP <- unlist(lapply(LO_snpN,function(x) x[1]))

LO_gwas$beta <- log(LO_gwas$OR)

LO_gwas$A2 <- ifelse(LO_gwas$A1==LO_gwas$REF,'ALT','REF')
LO_gwas$noneffect[LO_gwas$A2=='ALT'] <- LO_gwas$ALT[LO_gwas$A2=='ALT']
LO_gwas$noneffect[LO_gwas$A2=='REF'] <- LO_gwas$REF[LO_gwas$A2=='REF']

#protein gwas
gwas_path <- '/home1/prot_GCTA_result/'

Result_LO_m2 <- fread('/Result_LO_protein_assoc_M2.csv',data.table = F)
bonf_LO_m2 <- Result_LO_m2[Result_LO_m2$pval<0.025/nrow(Result_LO_m2),]
protName <- bonf_LO_m2$V1

#perform CAUSE
Result_cause_LOtoProt <- list()

set.seed(666)
for (i in 1:length(protName)){

    print(i)

    prot <- protName[i]
    #read protein gwas results
    protFile <- paste0(prot,'_assoc.fastGWA')
    prot_gwas <- fread(paste0(gwas_path,protFile))
    #revise SNP name
    prot_snpN <- strsplit(prot_gwas$SNP,':')
    prot_gwas$SNP2 <- unlist(lapply(prot_snpN,function(x) x[1]))

    #format
    gwas <- gwas_merge(LO_gwas, prot_gwas, snp_name_cols = c('SNP','SNP2'),
                       beta_hat_cols = c('beta','BETA'),
                       se_cols = c('LOG(OR)_SE','SE'),
                       A1_cols = c('A1','A1'),
                       A2_cols = c('noneffect','A2'))
    
    gwas <- merge(gwas,LO_gwas[,c('SNP','P')],by.x='snp',by.y='SNP',sort=F)
    colnames(gwas)[ncol(gwas)] <- 'p1'
    gwas <- merge(gwas,prot_gwas[,c('SNP2','P')],by.x='snp',by.y='SNP2',sort=F)
    colnames(gwas)[ncol(gwas)] <- 'p2'
    
    gwas2 <- gwas[which(gwas$p1<=1),]
    
    #Calculate nuisance parameters
    varlist <- with(gwas2, sample(snp, size=1000000, replace=FALSE))
    params <- est_cause_params(gwas2, varlist)
    
    #LD Pruning
    gwas_clump <- ld_clump_local(data.frame(rsid = gwas2$snp,pval = gwas2$p1),
                                   clump_kb=10000,
                                   clump_r2=0.01,
                                   clump_p = 1e-3,
                                   bfile = '/home1/shenchun/DATA/Genome1000_EUR/EUR',
                                   plink_bin = '/home1/shenchun/software/plink_1.9/plink')
    
    top_vars <- gwas_clump$rsid
    
    #Fit CAUSE
    res <- cause(X=gwas2, variants = top_vars, param_ests = params)
    
    #summary(res, ci_size=0.95)
    #plot(res)
    Result_cause_LOtoProt[[i]] <- res

}
names(Result_cause_LOtoProt) <- protName

save(Result_cause_LOtoProt,
     file = '/home1/shenchun/Projects/SI_protein/NHB_revision_v3/S6_causal_analysis/CAUSE/Result_cause_LOtoProt_2.RData')

#plot(Result_cause_LOtoProt[["GDF15"]])


sig_toprot <- fread('/home1/shenchun/Projects/SI_protein/NHB_revision_v3/S6_causal_analysis/MR_results/Result_IVW_sig_toProt.csv')
protName <- sig_toprot$outcome

#perform CAUSE
Result_cause_LOtoProt <- list()

set.seed(666)
for (i in 1:length(protName)){
  
  print(i)
  
  prot <- protName[i]
  #read protein gwas results
  protFile <- paste0(prot,'_assoc.fastGWA')
  prot_gwas <- fread(paste0(gwas_path,protFile))
  #revise SNP name
  prot_snpN <- strsplit(prot_gwas$SNP,':')
  prot_gwas$SNP2 <- unlist(lapply(prot_snpN,function(x) x[1]))
  
  #format
  gwas <- gwas_merge(LO_gwas, prot_gwas, snp_name_cols = c('SNP','SNP2'),
                     beta_hat_cols = c('beta','BETA'),
                     se_cols = c('LOG(OR)_SE','SE'),
                     A1_cols = c('A1','A1'),
                     A2_cols = c('noneffect','A2'))
  
  gwas <- merge(gwas,LO_gwas[,c('SNP','P')],by.x='snp',by.y='SNP',sort=F)
  colnames(gwas)[ncol(gwas)] <- 'p1'
  gwas <- merge(gwas,prot_gwas[,c('SNP2','P')],by.x='snp',by.y='SNP2',sort=F)
  colnames(gwas)[ncol(gwas)] <- 'p2'
  
  gwas2 <- gwas[which(gwas$p1<=1),]
  
  #Calculate nuisance parameters
  varlist <- with(gwas2, sample(snp, size=1000000, replace=FALSE))
  params <- est_cause_params(gwas2, varlist)
  
  #LD Pruning
  gwas_clump <- ld_clump_local(data.frame(rsid = gwas2$snp,pval = gwas2$p1),
                               clump_kb=10000,
                               clump_r2=0.01,
                               clump_p = 1e-3,
                               bfile = '/home1/shenchun/DATA/Genome1000_EUR/EUR',
                               plink_bin = '/home1/shenchun/software/plink_1.9/plink')
  
  top_vars <- gwas_clump$rsid
  
  #Fit CAUSE
  res <- cause(X=gwas2, variants = top_vars, param_ests = params)
  
  Result_cause_LOtoProt[[i]] <- res
  
}
names(Result_cause_LOtoProt) <- protName

save(Result_cause_LOtoProt, file = 'Result_cause_LOtoProt_MRsig.RData')

