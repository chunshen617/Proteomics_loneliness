#Mediation analysis caculated by PERM, using CVD as an example
library('data.table')
library('lubridate')                                        
library('RNOmni')
library('survival')
library('foreach')
library('doParallel')

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]

#MR significant proteins
mr_toProt <- fread('Result_MR_toProt_ME.csv')
fdrsig_mr_ivw <- mr_toProt[mr_toProt$pfdr_com<0.05,] #9 protein and protein modules
prot_use <- fdrsig_mr_ivw$outcome
#extract protein data
data_prot_use <- data_prot[,c('ID',prot_use,'SI_2c','LO_2c','age','sex','site','ethnicity','edu_4c','smokeNow','inc_2c','alc_2c','bmi_3c')]

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

#######Association between protein and disease
#Perform cox model
registerDoParallel(20)

Result_CVD_med1 <- foreach(i=1:length(prot_use), .combine='rbind') %dopar% {
  
  #simple model: adjust for age, sex
  dat_simple <- na.omit(subset(Data_merge3,select=c(prot_use[i],'outcome','fultime','LO_2c','age','sex')))
  dat_simple$exposure <- RankNorm(dat_simple[,1])
  
  M11 <- coxph(Surv(fultime,outcome)~LO_2c+age+sex,dat_simple)
  M12 <- coxph(Surv(fultime,outcome)~LO_2c+exposure+age+sex,dat_simple)
  
  s_M11 <- summary(M11)
  s_M12 <- summary(M12)
  
  result <- data.frame(OR1=s_M11[["coefficients"]][1,2],
                       OR2=s_M12[["coefficients"]][1,2],
                       pval1=s_M11[["coefficients"]][1,5],
                       pval2=s_M12[["coefficients"]][1,5],
                       med=((s_M11[["coefficients"]][1,2]-s_M12[["coefficients"]][1,2])/(s_M11[["coefficients"]][1,2]-1))*100)
  
  result
}

save(Result_CVD_med1, file='Result_prot_CVD_LOmed_M1.RData')
