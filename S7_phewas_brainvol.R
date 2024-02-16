library(data.table)
library(RNOmni)

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
#MR significant proteins
mr_toProt <- fread('Result_MR_toProt_ME.csv')
fdrsig_mr_ivw <- mr_toProt[mr_toProt$pfdr_com<0.05,] #9 protein and protein modules
prot_use <- fdrsig_mr_ivw$outcome[-9]
data_prot_used <- data_prot[,c('ID',prot_used,'sex','age','Batch','ethnicity','edu_4c','smokeNow','alc_2c','bmi_3c','inc_2c','timeGap')]

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
    dat <- na.omit(subset(aparc_vol_merge,select = c('eid',prot_used[i],aparcName[ii],'sex','age','age_dif','site2','Batch','ICV','ethnicity','edu_4c','inc_2c','smokeNow','alc_2c','bmi_3c','timeGap')))
    dat[,prot_used[i]] <- RankNorm(dat[,prot_used[i]])
    
    fit <- lm(as.formula(paste0(aparcName[ii],'~',prot_used[i],'+sex+age+age_dif+site2+Batch+ICV+ethnicity+edu_4c+inc_2c+smokeNow+alc_2c+bmi_3c+timeGap')),
              data=dat)
    sum <- summary(fit)
    
    result[ii,1] <- sum$coefficients[2,3]
    result[ii,2] <- sum$coefficients[2,4]
    result[ii,3] <- nrow(dat2)
  }
  
  Result_aparcVol[[i]] <- result
  
}

save(Result_aparcVol, file = 'Result_prot_aparc_vol.RData')
