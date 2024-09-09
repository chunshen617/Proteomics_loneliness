# Protein-wide association analysis using logistic regression
# Chun Shen, 2024

# Load R packages
library(parallel)
library(data.table)
library(RNOmni)

###Load data
# time gap: N participants * M proteins
timeGap <- fread('olink_timeGap.csv')
colnames(timeGap)[1] <- 'ID'
#protein data and other covariates
data_prot <- fread('protein_UKB.csv',data.table = F)
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)
colnames(data_prot)[2:2921] <- protName2
#covariates as factor
data_prot$SI_2c <- as.factor(data_prot$SI_2c) # Binary
data_prot$LO_2c <- as.factor(data_prot$LO_2c) # Binary
data_prot$sex <- as.factor(data_prot$sex) # Binary
data_prot$site <- as.factor(data_prot$site) # 22 levels
data_prot$Batch <- as.factor(data_prot$Batch) # 8 levels
data_prot$eth2 <- as.factor(data_prot$eth2) # 5 levels
data_prot$edu_4c <- as.factor(data_prot$edu_4c) # 4 levels
data_prot$inc_2c <- as.factor(data_prot$inc_2c) # Binary
data_prot$smokeNow <- as.factor(data_prot$smokeNow) # Binary
data_prot$alc_2c <- as.factor(data_prot$alc_2c) # Binary

##################
# Simple model, controlling for age，sex，site，batch, time gap between blood collection and protein measurement，20PC
model1logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+site+',
              paste0("PC", 1:20, collapse = "+"))
  fit <- glm(as.formula(f),data=data,family=binomial)
  or <- exp(coef(fit))
  s <- summary(fit)
  result <- data.frame(protname=protName,
                       N=nrow(data),
                       OR=or[2],
                       beta=s[["coefficients"]][2,1],
                       se=s[["coefficients"]][2,2],
                       z=s[["coefficients"]][2,3],
                       pval=s[["coefficients"]][2,4])
  return(result)
}

# Full model, controlling for age, sex, site, batch, time gap, 20 PC, ethnicity, education level, 
#household income, smoking, alcohol consumption, and BMI
model2logistic <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])
  f <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))
  fit <- glm(as.formula(f),data=data,family=binomial)
  or <- exp(coef(fit))
  s <- summary(fit)
  result <- data.frame(protname=protName,
                       N=nrow(data),
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
write.csv(SI_model1, file = 'Result_SI_protein_assoc_M1.csv')

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
write.csv(LO_model1, file = 'Result_LO_protein_assoc_M1.csv')

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
write.csv(SI_model2, file = 'Result_SI_protein_assoc_M2.csv')

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
}, mc.cores = 10)

LO_model2 <- do.call(rbind,LO.model2)
write.csv(LO_model2, file = 'Result_LO_protein_assoc_M2.csv')
