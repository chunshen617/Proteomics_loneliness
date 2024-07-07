library(parallel)
library(data.table)
library(RNOmni)
library(nnet)

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)

#recode group
data_prot$SL_group <- NA
data_prot$SL_group[data_prot$SI_2c==0 & data_prot$LO_2c==0] <- 1
data_prot$SL_group[data_prot$SI_2c==1 & data_prot$LO_2c==0] <- 2
data_prot$SL_group[data_prot$SI_2c==0 & data_prot$LO_2c==1] <- 3
data_prot$SL_group[data_prot$SI_2c==1 & data_prot$LO_2c==1] <- 4
data_prot$SL_group <- as.factor(data_prot$SL_group)

# Full model:  age, sex, site, batch, time gap, 20 PC, ethnicity, education level, 
#household income, smoking, alcohol consumption, BMI
model2multinom <- function(data,protName,phenotype){
  data[,protName] <- RankNorm(data[,protName])

  f1 <- paste0(phenotype,'~',protName,'+age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))
  # baseline model
  f2 <- paste0(phenotype,'~age+sex+Batch+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+site+',
              paste0("PC", 1:20, collapse = "+"))

  fit1 <- multinom(as.formula(f1),data=data,maxit=1000)
  fit2 <- multinom(as.formula(f2),data=data,maxit=1000)

  # convergence
  convf1 <- fit1$convergence
  convf2 <- fit2$convergence

  # compare model to baseline model
  lres <- anova(fit1,fit2)
  modelP <- pchisq(lres[2,"LR stat."], df=lres[2,"   Df"], lower.tail=FALSE)

  sumx <- summary(fit1)

  z <- sumx$coefficients/sumx$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1))*2
  or_prot <- exp(sumx$coefficients[,protName])
  
  result <- data.frame(protN=protName,
                       N=nrow(data),
                       or_S1L0=or_prot[1],
                       z_S1L0=z[1,2],
                       p_S1L0=p[1,2],
                       or_S0L1=or_prot[2],
                       z_S0L1=z[2,2],
                       p_S0L1=p[2,2],
                       or_S1L1=or_prot[3],
                       z_S1L1=z[3,2],
                       p_S1L1=p[3,2])
  return(result)
}

##################
SLgroup.model2 <- mclapply(1:length(protName2), function(i) {
  dat <- subset(data_prot,select = c('ID',protName2[i],'SL_group','Batch','age','sex','site',
                                     'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                     paste0('PC',1:20)))
  dat_t <- subset(timeGap,select = c('ID',protName[i]))
  colnames(dat_t)[2] <- 'timeGap'
  dat2 <- merge(dat,dat_t,by='ID')
  dat2 <- na.omit(dat2)
  
  model2multinom(dat2,protName2[i],'SL_group')
}, mc.cores = 10)

SLgroup_model2 <- do.call(rbind,SLgroup.model2)
write.csv(SLgroup_model2,file = 'Result_SLgroup_protein_assoc_M2_multinom.csv')

