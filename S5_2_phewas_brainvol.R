library(data.table)
library(RNOmni)

##Load data
data_prot <- fread('protein_UKB.csv',data.table = F)
#MR significant proteins
sig_toprot <- fread('Result_IVW_sig_toProt.csv')
prot_use <- sig_toprot$outcome
data_prot_used <- data_prot[,'ID',prot_used,'Batch','age','sex','site','timeGap',
                               'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                               paste0('PC',1:20)]

#Brain volume data
aparc_vol <- fread('UKB_2_0_aparc_GrayVol.csv',data.table = F)
aparc_vol_merge <- merge(data_prot_used,aparc_vol,by.x = 'eid',by.y = 'SubID')
aparc_vol_merge$age_dif <- aparc_vol_merge$age_image-aparc_vol_merge$age

Result_aparcVol <- vector(mode='list', length=length(prot_used))
names(Result_aparcVol) <- prot_used

for (i in 1:length(prot_used)){
  print(i)
  
  result <- data.frame(matrix(nrow = length(aparcName),ncol = 3))
  rownames(result) <- aparcName
  colnames(result) <- c('t_M1','p_M1','n')
  
  for (ii in 1:length(aparcName)){
    dat <- na.omit(subset(aparc_vol_merge,select = c('eid',prot_used[i],aparcName[ii],'sex','age','site2','ICV,'Batch','eth2','edu_4c','inc_2c','smokeNow','alc_2c','bmi','timeGap',
                                             paste0('PC',1:20))))
    dat[,prot_used[i]] <- RankNorm(dat[,prot_used[i]])
    
    fit <- lm(as.formula(paste0(aparcName[ii],'~',prot_used[i],'+sex+age+site2+ICV+Batch+eth2+edu_4c+inc_2c+smokeNow+alc_2c+bmi+timeGap+',
                                paste0("PC", 1:20, collapse = "+"))),
              data=dat)
    sum <- summary(fit)
    
    result[ii,1] <- sum$coefficients[2,3]
    result[ii,2] <- sum$coefficients[2,4]
    result[ii,3] <- nrow(dat2)
  }
  
  Result_aparcVol[[i]] <- result
  
}

save(Result_aparcVol, file = 'Result_prot_aparc_vol.RData')
