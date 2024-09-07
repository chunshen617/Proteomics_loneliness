# colocalization analysis
# Chun Shen, 2024

# Load R packages
library(coloc)
library(data.table)
library(ieugwasr)

## 500 Kb on both side of a significant GWAS SNP
GWAS_WINDOW <- 500000

## LO gwas
LO_GWAS <- fread('LO_GWAS_plink2.csv')
colnames(LO_GWAS)[1] <- 'CHR'
LO_GWAS$beta <- log(LO_GWAS$OR)

##====== LO significant snps
LO_GWAS_sig <- LO_GWAS[which(LO_GWAS$P < 5e-8), ]
snpN <- strsplit(LO_GWAS_sig$ID,':')
LO_GWAS_sig$SNP <- unlist(lapply(snpN,function(x) x[1]))
LO_sigsnp <- ld_clump_local(data.frame(rsid = LO_GWAS_sig$SNP,pval = LO_GWAS_sig$P),
                               clump_kb=10000,
                               clump_r2=0.001,
                               clump_p = 1,
                               bfile = '/Genome1000_EUR/EUR',
                               plink_bin = '/software/plink_1.9/plink')
LO_sigsnp_clump <- LO_GWAS_sig[LO_GWAS_sig$SNP %in% LO_sigsnp$rsid,]

LO_sigsnp_clump$start <- LO_sigsnp_clump$POS-GWAS_WINDOW
LO_sigsnp_clump$end <- LO_sigsnp_clump$POS+GWAS_WINDOW

#Proteins
#significant proteins in IVW
sig_toprot <- fread('Result_IVW_sig_toProt.csv')
protName <- sig_toprot$outcome

LO_coloc_result <- vector(mode='list', length=length(protName))
names(LO_coloc_result) <- protName

for (jj in 1:length(protName)){
  print(jj)
  
  pQTLData <- fread(paste0('/home1/prot_GCTA_result/',protName[jj],'_assoc.fastGWA'))
  prot_sig <- pQTLData[which(pQTLData$P < 5e-8), ]
  prot_snpN <- strsplit(prot_sig$SNP,':')
  prot_sig$SNP <- unlist(lapply(prot_snpN,function(x) x[1]))
  prot_sigsnp <- ld_clump_local(data.frame(rsid = prot_sig$SNP,pval = prot_sig$P),
                              clump_kb=10000,
                              clump_r2=0.001,
                              clump_p = 1,
                              bfile = '/Genome1000_EUR/EUR',
                              plink_bin = '/software/plink_1.9/plink')
  prot_sigsnp_clump <- prot_sig[prot_sig$SNP %in% prot_sigsnp$rsid,]
  prot_sigsnp_clump$start <- prot_sigsnp_clump$POS-GWAS_WINDOW
  prot_sigsnp_clump$end <- prot_sigsnp_clump$POS+GWAS_WINDOW
  
  #all snps
  GWASData_sig_clump <- rbind(LO_sigsnp_clump[,c('CHR','start','end','SNP')],prot_sigsnp_clump[,c('CHR','start','end','SNP')])
  
  for (i in 1:nrow(GWASData_sig_clump)){
    currchr <- GWASData_sig_clump$CHR[i]
    
    GWASD_currchr <- LO_GWAS[which(LO_GWAS$CHR == currchr),]
    pQTLData_currchr <- pQTLData[which(pQTLData$CHR == currchr), ]
    
    ## gwas data for the current loci
    currloci_GWASdata <- GWASD_currchr[which((GWASD_currchr[, c("POS")] >= GWASData_sig_clump$start[i]) & 
                                               (GWASD_currchr[, c("POS")] <= GWASData_sig_clump$end[i])), ]
    ## pQTLs for the current loci
    currloci_pqtldata <- pQTLData_currchr[which((pQTLData_currchr[, c("POS")] >= GWASData_sig_clump$start[i]) & 
                                                  (pQTLData_currchr[, c("POS")] <= GWASData_sig_clump$end[i])), ]
    currloci_merge <- merge(currloci_GWASdata,currloci_pqtldata,by.x='ID',by.y='SNP',suffixes=c("_gwas","_eqtl"))
    
    #make sure that the betas refer to the same allele!
    idx <- which(currloci_merge$A1_gwas!=currloci_merge$A1_eqtl)
    currloci_merge$beta_eqtl_revise <- currloci_merge$BETA
    currloci_merge$beta_eqtl_revise[idx] <- currloci_merge$BETA[idx]*-1
    
    result <- coloc.abf(dataset1=list(snp = currloci_merge$ID,pvalues=currloci_merge$P_gwas, type="cc", beta=currloci_merge$beta, varbeta=currloci_merge$`LOG(OR)_SE`^2, position=currloci_merge$POS_gwas), 
                        dataset2=list(snp = currloci_merge$ID,pvalues=currloci_merge$P_eqtl, type="quant", beta=currloci_merge$beta_eqtl_revise, varbeta = currloci_merge$SE^2, N = currloci_merge$N, position=currloci_merge$POS_eqtl), MAF=currloci_merge$A1_FREQ)
    #sensitivity(result,rule="H4 > 0.5")
    
    LO_coloc_result[[jj]][[i]] <- result
  }
  names(LO_coloc_result[[jj]]) <- GWASData_sig_clump$SNP
  
}

#clean results
LO_prot_coloc <- data.frame()

for (i in 1:length(protName)){
  
  rr <- LO_coloc_result[[protName[i]]]
  ss <- lapply(rr,function(x) x[['summary']])
  tb <- as.data.frame(do.call(rbind,ss))
  tb$protName <- protName[i]
  tb$snpName <- rownames(tb)
  rownames(tb) <- NULL
  
  LO_prot_coloc <- rbind(LO_prot_coloc,tb)
}

write.csv(LO_prot_coloc,file = 'Result_loTOprot_coloc.csv', row.names = F)
