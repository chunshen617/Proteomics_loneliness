# Cox-based mediation analysis, using CVD as an example
# Chun Shen, 2024

# Load R packages
library('data.table')
library('lubridate')                                        
library('RNOmni')
library('survival')
library('foreach')
library('doParallel')

###Load data
#protein data and other covariates
data_prot <- fread('protein_UKB.csv',data.table = F)
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)
colnames(data_prot)[2:2921] <- protName2

# MR sig proteins
sig_toprot <- fread('Result_IVW_sig_toProt.csv')
prot_used <- sig_toprot$outcome

data_prot_use <- data_prot[,c('ID',prot_used,'SI_2c','LO_2c','sex','eth2','age','site',
                              'edu_4c','smokeNow','alc_2c','bmi','inc_2c')]

#Health outcome data from UKB first-occurence
phenoData <- fread('COXphenoData_raw.csv')

timedif <- function(eventD,attendD,deathD,endD){
  result <- vector()
  #event exist
  result[is.na(eventD)==FALSE] <- 
    time_length(difftime(eventD[is.na(eventD)==FALSE],
                         attendD[is.na(eventD)==FALSE]),'years')
  #no event & isn't dead
  result[is.na(eventD) & is.na(deathD)] <- 
    time_length(difftime(endD,
                         attendD[is.na(eventD) & is.na(deathD)]),'years')
  #no event but dead
  result[is.na(eventD) & is.na(deathD)==FALSE] <- 
    time_length(difftime(deathD[is.na(eventD) & is.na(deathD)==FALSE],
                         attendD[is.na(eventD) & is.na(deathD)==FALSE]),'years')
  return(result)
}

#Disease: CVD
Data_merge <- merge(data_prot_use,phenoData[,c('eid','131296-0.0','131298-0.0','131300-0.0','131302-0.0','131304-0.0',
                                             '131306-0.0','40000-0.0','53-0.0','20002-0.0')],
                    by.x = 'ID', by.y='eid')
colnames(Data_merge)[196:ncol(Data_merge)] <- c('I20Date','I21Date','I22Date','I23Date','I24Date','I25Date','DeathDate','AttendDate','SelfReport')

#Outcome: all-cause CVD diagnosis
Data_merge$CVD1Date <- apply(Data_merge[,c('I20Date','I21Date','I22Date','I23Date','I24Date','I25Date')],1,function(x) min(x,na.rm = T))
#self-reported CVD at baseline
CVD_selfD <- unique(c(grep('1066',Data_merge$SelfReport),grep('1075',Data_merge$SelfReport),
                      grep('1076',Data_merge$SelfReport),grep('1077',Data_merge$SelfReport),
                      grep('1078',Data_merge$SelfReport),grep('1485',Data_merge$SelfReport)))
Data_merge2 <- Data_merge[-CVD_selfD,]
#prevalent cases
Data_merge3 <- Data_merge2[-which(time_length(difftime(Data_merge2$CVD1Date,Data_merge2$AttendDate),'years')<3),]
#outcome
Data_merge3$outcome <- ifelse(is.na(Data_merge3$CVD1Date)==FALSE,1,0)
Data_merge3$fultime <- timedif(Data_merge3$CVD1Date,Data_merge3$AttendDate,Data_merge3$DeathDate,'2022-11-30')

##############
registerDoParallel(10)

Result_CVD_cmed2 <- foreach(i=1:length(prot_used), .combine='rbind') %dopar% {
  # Covariate
  mydata <- na.omit(subset(Data_merge3,select=c('ID',prot_used[i],'outcome','fultime','LO_2c','age','sex','site',
                                           'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c')))
  mydata$mediator <- RankNorm(mydata[,2])

  # category
  mydata$LO_2c <- as.factor(mydata$LO_2c)
  mydata$sex <- as.factor(mydata$sex)
  mydata$eth2 <- as.factor(mydata$eth2)
  mydata$edu_4c <- as.factor(mydata$edu_4c)
  mydata$inc_2c <- as.factor(mydata$inc_2c)
  mydata$smokeNow <- as.factor(mydata$smokeNow)
  mydata$alc_2c <- as.factor(mydata$alc_2c)
  mydata$site <- as.factor(mydata$site)
  
  ##############
  f1 <- 'Surv(fultime,outcome)~LO_2c+age+sex+eth2+edu_4c+smokeNow+inc_2c+alc_2c+bmi+site'
  f2 <- 'Surv(fultime,outcome)~LO_2c+mediator+age+sex+eth2+edu_4c+smokeNow+inc_2c+alc_2c+bmi+site'

  M1 <- coxph(as.formula(f1),mydata)
  M2 <- coxph(as.formula(f2),mydata)

  s_M1 <- summary(M1)
  s_M2 <- summary(M2)
  
  result <- data.frame(protname = prot_used[i],
                       disease='CVD',
                       n=nrow(mydata),
                       OR1=s_M1[["coefficients"]][1,2],
                       OR2=s_M2[["coefficients"]][1,2],
                       pval1=s_M1[["coefficients"]][1,5],
                       pval2=s_M2[["coefficients"]][1,5],
                       med=((s_M1[["coefficients"]][1,2]-s_M2[["coefficients"]][1,2])/(s_M1[["coefficients"]][1,2]-1))*100)

  result
}

write.csv(Result_CVD_cmed2, file='Result_PERM_CVD_M2_v2.csv', row.names = F)
