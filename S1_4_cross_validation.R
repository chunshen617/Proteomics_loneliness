# Sensitivity analysis: cross-validation by randomly splitting socially isolated, lonely, and control participants seperately
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

##################
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

####################
# Number of randomization
Nt <- 100
data_prot_SI1 <- data_prot[data_prot$SI_2c==1,] 
data_prot_SI0 <- data_prot[data_prot$SI_2c==0,] 

# Social isolation
SI.result <- list()
for (n in 1:Nt){

    resam_1 <- sample(1:nrow(data_prot_SI1))
    resam_0 <- sample(1:nrow(data_prot_SI0))

    dat_resam <- list()
    dat_resam[[1]] <- rbind(data_prot_SI1[which(resam_1<=nrow(data_prot_SI1)/2),],data_prot_SI0[which(resam_0<=nrow(data_prot_SI0)/2),])
    dat_resam[[2]] <- rbind(data_prot_SI1[which(resam_1>nrow(data_prot_SI1)/2),],data_prot_SI0[which(resam_0>nrow(data_prot_SI0)/2),])
    
    SI.result.resam <- list()
    
    for (s in 1:2){
      
      dat_resam2 <- dat_resam[[s]]
      
      SI.result.resam[[s]] <- mclapply(1:length(protName2), function(i) {
        dat <- subset(dat_resam2,select = c('ID',protName2[i],'SI_2c','LO_2c','Batch','age','sex','site',
                                            'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                            paste0('PC',1:20)))
        dat_t <- subset(timeGap,select = c('ID',protName[i]))
        colnames(dat_t)[2] <- 'timeGap'
        dat2 <- merge(dat,dat_t,by='ID')
        dat2 <- na.omit(dat2)
        
        #cov as factor
        dat2$SI_2c <- as.factor(dat2$SI_2c)
        dat2$sex <- as.factor(dat2$sex)
        dat2$site <- as.factor(dat2$site)
        dat2$Batch <- as.factor(dat2$Batch)
        dat2$eth2 <- as.factor(dat2$eth2)
        dat2$edu_4c <- as.factor(dat2$edu_4c)
        dat2$inc_2c <- as.factor(dat2$inc_2c)
        dat2$smokeNow <- as.factor(dat2$smokeNow)
        dat2$alc_2c <- as.factor(dat2$alc_2c)
        
        if (nrow(dat2)!=0){
          model2logistic(dat2,protName2[i],'SI_2c')
        }
        
      }, mc.cores = 10)
      
    }
    
    SI.result[[n]] <- SI.result.resam
    
}
save(SI.result, file = "SI_PWAS_randomsplit_100t.RData")

# Loneliness
data_prot_LO1 <- data_prot[data_prot$LO_2c==1,]
data_prot_LO0 <- data_prot[data_prot$LO_2c==0,] 

LO.result <- list()
for (n in 1:Nt){
  
  resam_1 <- sample(1:nrow(data_prot_LO1))
  resam_0 <- sample(1:nrow(data_prot_LO0))
  
  dat_resam <- list()
  dat_resam[[1]] <- rbind(data_prot_LO1[which(resam_1<=nrow(data_prot_LO1)/2),],data_prot_LO0[which(resam_0<=nrow(data_prot_LO0)/2),])
  dat_resam[[2]] <- rbind(data_prot_LO1[which(resam_1>nrow(data_prot_LO1)/2),],data_prot_LO0[which(resam_0>nrow(data_prot_LO0)/2),])
  
  LO.result.resam <- list()
  
  for (s in 1:2){
    
    dat_resam2 <- dat_resam[[s]]
    
    LO.result.resam[[s]] <- mclapply(1:length(protName2), function(i) {
      dat <- subset(dat_resam2,select = c('ID',protName2[i],'LO_2c','LO_2c','Batch','age','sex','site',
                                          'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                          paste0('PC',1:20)))
      dat_t <- subset(timeGap,select = c('ID',protName[i]))
      colnames(dat_t)[2] <- 'timeGap'
      dat2 <- merge(dat,dat_t,by='ID')
      dat2 <- na.omit(dat2)
      
      #cov as factor
      dat2$LO_2c <- as.factor(dat2$LO_2c)
      dat2$sex <- as.factor(dat2$sex)
      dat2$site <- as.factor(dat2$site)
      dat2$Batch <- as.factor(dat2$Batch)
      dat2$eth2 <- as.factor(dat2$eth2)
      dat2$edu_4c <- as.factor(dat2$edu_4c)
      dat2$inc_2c <- as.factor(dat2$inc_2c)
      dat2$smokeNow <- as.factor(dat2$smokeNow)
      dat2$alc_2c <- as.factor(dat2$alc_2c)
      
      if (nrow(dat2)!=0){
        model2logistic(dat2,protName2[i],'LO_2c')
      }
      
    }, mc.cores = 10)
    
  }
  
  LO.result[[n]] <- LO.result.resam
  
}
save(LO.result, file = "LO_PWAS_randomsplit_100t2.RData")

