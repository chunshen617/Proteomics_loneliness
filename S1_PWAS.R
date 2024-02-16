library(data.table)
library(RNOmni)
library(foreach)
library(doParallel)

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)

##Model1:controlling for batch, time gap between blood collection and protein measurement, age, sex and site
model1logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+site')
  fit <- glm(as.formula(f),data=data,family=binomial)
  or <- exp(coef(fit))
  s <- summary(fit)
  result <- data.frame(N=nrow(data),
                       OR=or[2],
                       beta=s[["coefficients"]][2,1],
                       se=s[["coefficients"]][2,2],
                       z=s[["coefficients"]][2,3],
                       pval=s[["coefficients"]][2,4])
  return(result)
}
#Model2:controlling for batch, time gap between blood collection and protein measurement, age, sex and site,
#ethnicity, qualification, smoking, alcohol intake, TDI, and BMI
model2logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+ethnicity+edu_4c+smokeNow+alc_2c+bmi_3c+inc_2c+site')
  fit <- glm(as.formula(f),data=data,family=binomial)
  or <- exp(coef(fit))
  s <- summary(fit)
  result <- data.frame(N=nrow(data),
                       OR=or[2],
                       beta=s[["coefficients"]][2,1],
                       se=s[["coefficients"]][2,2],
                       z=s[["coefficients"]][2,3],
                       pval=s[["coefficients"]][2,4])
  return(result)
}

##setup parallel backend to use many processors
registerDoParallel(20)

##perform analysis
Result_SI_m1 <- foreach(i=1:length(protName), .combine='rbind') %dopar% {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site'))
  dat_t <- subset(timeGap,select = c('eid',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by.x='ID',by.y='eid')
  dat2 <- na.omit(dat2)
  
  r <- model1logistic(dat2,protName2[i],'SI_2c')
  r
}
Result_LO_m1 <- foreach(i=1:length(protName), .combine='rbind') %dopar% {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site'))
  dat_t <- subset(timeGap,select = c('eid',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by.x='ID',by.y='eid')
  dat2 <- na.omit(dat2)
  
  r <- model1logistic(dat2,protName2[i],'LO_2c')
  r
}

Result_SI_m2 <- foreach(i=1:length(protName), .combine='rbind') %dopar% {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',
                                     'ethnicity','edu_4c','smokeNow','alc_2c','bmi_3c','inc_2c'))
  dat_t <- subset(timeGap,select = c('eid',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by.x='ID',by.y='eid')
  dat2 <- na.omit(dat2)
  
  r <- model2logistic(dat2,protName2[i],'SI_2c')
  r
}
Result_LO_m2 <- foreach(i=1:length(protName), .combine='rbind') %dopar% {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',
                                     'ethnicity','edu_4c','smokeNow','alc_2c','bmi_3c','inc_2c'))
  dat_t <- subset(timeGap,select = c('eid',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by.x='ID',by.y='eid')
  dat2 <- na.omit(dat2)
  
  r <- model2logistic(dat2,protName2[i],'LO_2c')
  r
}

##save results
write.csv(Result_SI_m1, file = 'Result_SI_protein_assoc_M1.csv')
write.csv(Result_LO_m1, file = 'Result_LO_protein_assoc_M1.csv')
write.csv(Result_SI_m2, file = 'Result_SI_protein_assoc_M2_noIRN.csv')
write.csv(Result_LO_m2, file = 'Result_LO_protein_assoc_M2_noIRN.csv')
