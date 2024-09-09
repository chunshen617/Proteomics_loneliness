# Sensitivity analysis: social isolation and loneliness raw scores were used and modeled with ordered logistic regression
# Chun Shen, 2024

# Load R packages
library(parallel)
library(data.table)
library(RNOmni)
library(MASS)
library(lmtest)

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
data_prot$SI <- as.factor(data_prot$SI) # 4 levels
data_prot$LO <- as.factor(data_prot$LO) # 3 levels
data_prot$sex <- as.factor(data_prot$sex) # Binary
data_prot$site <- as.factor(data_prot$site) # 22 levels
data_prot$Batch <- as.factor(data_prot$Batch) # 8 levels
data_prot$eth2 <- as.factor(data_prot$eth2) # 5 levels
data_prot$edu_4c <- as.factor(data_prot$edu_4c) # 4 levels
data_prot$inc_2c <- as.factor(data_prot$inc_2c) # Binary
data_prot$smokeNow <- as.factor(data_prot$smokeNow) # Binary
data_prot$alc_2c <- as.factor(data_prot$alc_2c) # Binary

# Full model, controlling for age, sex, site, batch, time gap, 20 PC, ethnicity, education level, 
#household income, smoking, alcohol consumption, and BMI
model2order <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])

  f1 <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))
  # baseline model
  f2 <- paste0(phenotype,'~age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))

  fit1 <- polr(as.formula(f1),data=data,Hess=TRUE)
  fit2 <- polr(as.formula(f2),data=data,Hess=TRUE)

  # convergence
  convf1 <- fit1$convergence
  convf2 <- fit2$convergence

  ctable <- coef(summary(fit1))
  ct <- coeftest(fit1)

  # compare model to baseline model
  lres <- anova(fit1,fit2)
  modelP <- pchisq(lres[2,"LR stat."], df=lres[2,"   Df"], lower.tail=FALSE)

  result <- data.frame(protN=protName,
                       N=nrow(data),
                       modelP=modelP,
                       beta=ctable[protName,'Value'],
                       se=ctable[protName,'Std. Error'],
                       t=ctable[protName,'t value'],
                       pval=ct[1,4])
  return(result)
}

##################
# SI, full model
SI.model2 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI','LO','Batch','age','sex','site',
                                     'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                     paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model2order(dat2,protName2[i],'SI')
}, mc.cores = 10)

SI_model2 <- do.call(rbind,SI.model2)
write.csv(SI_model2, file = 'Result_SI_protein_order_M2.csv',row.names = F)

# LO, full model
LO.model2 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SI','LO','Batch','age','sex','site',
                                     'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                     paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model2order(dat2,protName2[i],'LO')
}, mc.cores = 10)

LO_model2 <- do.call(rbind,LO.model2)
write.csv(LO_model2, file = 'Result_LO_protein_order_M2.csv',row.names = F)
