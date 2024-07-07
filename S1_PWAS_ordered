library(parallel)
library(data.table)
library(RNOmni)
library(MASS)
library(lmtest)

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)

# Full model:  age, sex, site, batch, time gap, 20 PC, ethnicity, education level, 
#household income, smoking, alcohol consumption, BMI
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
write.csv(SI_model2,file = 'Result_SI_protein_order_M2.csv',row.names = F)

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
write.csv(LO_model2,file = 'Result_LO_protein_order_M2.csv',row.names = F)
