library(parallel)
library(data.table)
library(RNOmni)

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)

##################
# Simple model: controlling for age，sex，site，technical，first 20 genetic PCs
model1logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+site+',
              paste0("PC", 1:20, collapse = "+"))
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

# Full model:  controlling for age, sex, site, batch, time gap, 20 PC, ethnicity, education level, 
#household income, smoking, alcohol consumption, BMI
model2logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))
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

##################
# SI, simple model
SI.model1 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model1logistic(dat2,protName2[i],'SI_2c')
}, mc.cores = 10)

SI_model1 <- do.call(rbind,SI.model1)
write.csv(SI_model1,file = 'Result_SI_protein_assoc_M1.csv')

##################
# LO, simple model
LO.model1 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model1logistic(dat2,protName2[i],'LO_2c')
}, mc.cores = 10)

LO_model1 <- do.call(rbind,LO.model1)
write.csv(LO_model1,file = 'Result_LO_protein_assoc_M1.csv')

##################
# SI, full model
SI.model2 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',
                                     'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                     paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model2logistic(dat2,protName2[i],'SI_2c')
}, mc.cores = 10)

SI_model2 <- do.call(rbind,SI.model2)
write.csv(SI_model2,file = 'Result_SI_protein_assoc_M2.csv')

##################
# LO, full model
LO.model2 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',
                                     'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                     paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model2logistic(dat2,protName2[i],'LO_2c')
}, mc.cores = 8)

LO_model2 <- do.call(rbind,LO.model2)
write.csv(LO_model2,file = 'Result_LO_protein_assoc_M2.csv')
