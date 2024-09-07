#NMR data preprocessing
library(ukbnmr)
library(data.table)

dat <- fread('blood_raw_prot.csv')
nmr <- extract_biomarkers(dat)
biomarker_qc_flags <- extract_biomarker_qc_flags(dat)
sample_qc_flags <- extract_sample_qc_flags(dat)

processed <- remove_technical_variation(dat) 
processed.nmr <- processed[['biomarkers']]

write.csv(processed.nmr,file='nmr_processed.csv',row.names = F)

#############
#Association to blood biomarkers 
#load protein data
data_prot <- fread('protein_UKB.csv',data.table = F)
#MR significant proteins
sig_toprot <- fread('Result_IVW_sig_toProt.csv')
prot_use <- sig_toprot$outcome
#extract protein data
data_prot_use <- data_prot[,c('ID',prot_used,'Batch','age','sex','site','timeGap',
                               'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                               paste0('PC',1:20))]

#load blood biomarkers data
dat_bld <- fread('blood_processed.csv')
dat_prot_bld <- merge(data_prot_use,dat_bld,by.x='ID',by.y='eid',sort=F)
var_bld <- colnames(dat_prot_bld)[-1]

result <- vector(mode = 'list',length = length(prot_use))
names(result) <- prot_use

for (i in 1:length(prot_use)){
  print(i)
  
  rr <- data.frame(matrix(nrow=length(var_bld),ncol=3))
  colnames(rr) <- c('tval','pval','N')
  rownames(rr) <- var_bld
  
  for (ii in 1:length(var_bld)){
    dat <- na.omit(subset(dat_prot_bld,select = c('ID',prot_use[i],var_bld[ii],'sex','age','site','Batch','eth2','edu_4c','inc_2c','smokeNow','alc_2c','bmi','timeGap',
                                             paste0('PC',1:20))))
    dat[,prot_use[i]] <- RankNorm(dat[,prot_use[i]])
    
    fit <- lm(as.formula(paste0(var_bld[ii],'~',prot_use[i],'+sex+age+site+Batch+eth2+edu_4c+inc_2c+smokeNow+alc_2c+bmi+timeGap+',
                                paste0("PC", 1:20, collapse = "+"))),
              data=dat)
    sum <- summary(fit)
    
    rr[ii,1] <- sum$coefficients[2,3]
    rr[ii,2] <- sum$coefficients[2,4]
    rr[ii,3] <- nrow(dat)
  }
  result[[i]] <- rr
}

save(result,file = 'phewas_bld.RData')
