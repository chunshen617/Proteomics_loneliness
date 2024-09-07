# GSMR analysis, using loneliness to proteins as an example
# Chun Shen, 2024

# Load R packages
library(gsmr2)
library(data.table)

#LO gwas
LO_gwas <- fread('LO_GWAS_plink2.csv')
#significant threshold
LO_gwas_sig <- LO_gwas[LO_gwas$P<1e-6,]
#revise SNP name
LO_snpN <- strsplit(LO_gwas_sig$ID,':')
LO_gwas_sig$SNP <- unlist(lapply(LO_snpN,function(x) x[1]))

LO_gwas_sig$beta <- log(LO_gwas_sig$OR)

LO_gwas_sig$A2 <- ifelse(LO_gwas_sig$A1==LO_gwas_sig$REF,'ALT','REF')
LO_gwas_sig$noneffect[LO_gwas_sig$A2=='ALT'] <- LO_gwas_sig$ALT[LO_gwas_sig$A2=='ALT']
LO_gwas_sig$noneffect[LO_gwas_sig$A2=='REF'] <- LO_gwas_sig$REF[LO_gwas_sig$A2=='REF']

#clumping
LO_snp_clump <- ld_clump_local(data.frame(rsid = LO_gwas_sig$SNP,pval = LO_gwas_sig$P),
                               clump_kb=10000,
                               clump_r2=0.001,
                               clump_p = 1,
                               bfile = '/Genome1000_EUR/EUR',
                               plink_bin = '/software/plink_1.9/plink')

LO_exp_clump <- LO_gwas_sig[which(LO_gwas_sig$SNP %in% LO_snp_clump$rsid),]

# protein GWAS folder
gwas_path <- '/home1/prot_GCTA_result/'

#significant proteins associated with loneliness
Result_LO_m2 <- fread('Result_LO_protein_assoc_M2.csv',data.table = F)
bonf_LO_m2 <- Result_LO_m2[Result_LO_m2$pval<0.025/nrow(Result_LO_m2),]
protName <- bonf_LO_m2$V1

# individual genetic data 
snp_coeff_id <- scan('/home1/GSMR/gsmr_allSNPs.xmat.gz', what='', nlines=1)
snp_coeff <- read.table('/home1/GSMR/gsmr_allSNPs.xmat.gz', header=F, skip=2)

# Estimate LD correlation matrix using R 
snp_id <- Reduce(intersect,	list(LO_exp_clump$ID, snp_coeff_id))
snp_order <- match(snp_id, snp_coeff_id)
snp_coeff_id <- snp_coeff_id[snp_order]
snp_coeff <- snp_coeff[,snp_order]

ldrho <- cor(snp_coeff,use='p')
colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id

Result_gsmr_LOtoProt <- data.frame(matrix(nrow=length(protName),ncol=7))
colnames(Result_gsmr_LOtoProt) <- c('bxy','bxy_se','pval','Nlinkage_snps','Npleio_snps','exposure','outcome')
for (i in 1:length(protName)){

    print(i)

    prot <- protName[i]
    #read protein gwas results
    protFile <- paste0(prot,'_assoc.fastGWA')
    prot_gwas <- fread(paste0(gwas_path,protFile))
  
    prot_gwas_out <- prot_gwas[which(prot_gwas$SNP %in% LO_exp_clump$ID),]
    identical(LO_exp_clump$ID,prot_gwas_out$SNP)

    #gsmr data
    gsmr_data <- data.frame(SNP = LO_exp_clump$ID,
                            a1 = LO_exp_clump$A1,
                            a2 = LO_exp_clump$noneffect,
                            a1_freq = LO_exp_clump$A1_FREQ,
                            bzx = LO_exp_clump$beta,
                            bzx_se = LO_exp_clump$`LOG(OR)_SE`,
                            bzx_pval = LO_exp_clump$P,
                            bzx_n = LO_exp_clump$OBS_CT,
                            bzy_se = prot_gwas_out$SE,
                            bzy_pval = prot_gwas_out$P,
                            bzy_n = prot_gwas_out$N)
    gsmr_data$bzy <- ifelse(prot_gwas_out$A1==gsmr_data$a1,prot_gwas_out$BETA,-1*prot_gwas_out$BETA)

    gsmr_data <- gsmr_data[match(snp_id, gsmr_data$SNP),]

    # GSMR
    bzx	=	gsmr_data$bzx				 #	SNP	effects	on	the	risk	factor 
    bzx_se	=	gsmr_data$bzx_se				 #	standard	errors	of	bzx 
    bzx_pval	=	gsmr_data$bzx_pval			 #	p-values	for	bzx 
    bzy	=	gsmr_data$bzy				 #	SNP	effects	on	the	disease 
    bzy_se	=	gsmr_data$bzy_se				 #	standard	errors	of	bzy 
    bzy_pval	=	gsmr_data$bzy_pval				 #	p-values	for	bzy 
    n_ref	=	10000				#	Sample	size	of	the	reference	sample 
    gwas_thresh	=	1e-6				#	GWAS	threshold	to	select	SNPs	as	the	instruments	for	the	GSMR	analysis 
    single_snp_heidi_thresh	=	 0.05				#	p-value	threshold	for	single-SNP-based	HEIDI-outlier	analysis 
    multi_snp_heidi_thresh	=	 0.05				#	p-value	threshold	for	multi-SNP-based	HEIDI-outlier	analysis 
    nsnps_thresh	=	 1			#	the	minimum	number	of	instruments	required	for	the	GSMR	analysis 
    heidi_outlier_flag	=	 T				#	flag	for	HEIDI-outlier	analysis 
    ld_r2_thresh	=	 0.05				#	LD	r2	threshold	to	remove	SNPs	in	high	LD 
    ld_fdr_thresh	=	 0.05			#	FDR	threshold	to	remove	the	chance	correlations	between	the	SNP	instruments 
    gsmr_beta	=	0					#	0	-	the	original	HEIDI-outlier	method;	1	-	the	new	HEIDI-outlier	method	that	is	currently	under	development	 
    gsmr_results	=	gsmr(bzx,	bzx_se,	bzx_pval,	bzy,	bzy_se,	bzy_pval,	ldrho,	snp_coeff_id,	n_ref,	heidi_outlier_flag,	gwas_thresh,	single_snp_heidi_thresh,	multi_snp_heidi_thresh,	nsnps_thresh,	ld_r2_thresh,	ld_fdr_thresh,	gsmr_beta)				 #	GSMR	analysis	 

    Result_gsmr_LOtoProt[i,1] <- gsmr_results[["bxy"]]
    Result_gsmr_LOtoProt[i,2] <- gsmr_results[["bxy_se"]]
    Result_gsmr_LOtoProt[i,3] <- gsmr_results[["bxy_pval"]]
    Result_gsmr_LOtoProt[i,4] <- length(gsmr_results$linkage_snps)
    Result_gsmr_LOtoProt[i,5] <- length(gsmr_results$pleio_snps)
    Result_gsmr_LOtoProt[i,6] <- 'LO'
    Result_gsmr_LOtoProt[i,7] <- prot
}

write.csv(Result_gsmr_LOtoProt,file='Result_GSMR_LOtoProt.csv',row.names = F)
