#Prospective association between protein and diseases and mortality, using CVD as an example
library('data.table')
library('lubridate')                                        
library('RNOmni')
library('survival')
library('foreach')
library('doParallel')

##Load data
data_prot <- fread('olink_ins0_cov.csv',data.table = F) #plasma protein data with covariates
protName <- colnames(data_prot)[2:2921]

##Model 2 result
Result_SI_m2 <- fread('Result_SI_protein_assoc_M2.csv',data.table = F)
Result_LO_m2 <- fread('Result_LO_protein_assoc_M2.csv',data.table = F)
#Bonferroni correction
bonf_SI_m2 <- Result_SI_m2[Result_SI_m2$pval<0.05/(nrow(Result_SI_m2)+nrow(Result_LO_m2)),]#171
bonf_LO_m2 <- Result_LO_m2[Result_LO_m2$pval<0.05/(nrow(Result_SI_m2)+nrow(Result_LO_m2)),]#43

prot_use <- unique(c(bonf_SI_m2$V1,bonf_LO_m2$V1))#183
data_prot_use <- data_prot[,c('ID',prot_used,'age','sex','timeGap',
                               'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                               paste0('PC',1:20)]

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

Result_CVD_M2 <- foreach(i=1:length(prot_use), .combine='rbind') %dopar% {
  dat_full <- na.omit(subset(Data_merge3,select=c('ID',prot_use[i],'outcome','fultime',
                                                  'age','sex','ethnicity','edu_4c','smokeNow','alc_2c','bmi_3c','inc_2c')))
  dat_full$exposure <- RankNorm(dat_full[,2])

  f <- paste0('Surv(fultime,outcome)~exposure+age+sex+timeGap+eth2+edu_4c+smokeNow+alc_2c+bmi+inc_2c+',
              paste0("PC", 1:20, collapse = "+"))
  M2 <- coxph(as.formula(f),dat_full)
  s_M2 <- summary(M2)
  c_M2 <- cox.zph(M2)
  
  result <- data.frame(n=s_M2[["n"]],
                       nevent=s_M2[["nevent"]],
                       OR=s_M2[["coefficients"]][1,2],
                       lower=s_M2[["conf.int"]][1,3],
                       upper=s_M2[["conf.int"]][1,4],
                       pval=s_M2[["coefficients"]][1,5],
                       zph=c_M2[["table"]][10,3])
  result
}

rownames(Result_CVD_M2) <- prot_use
write.csv(Result_CVD_M2,file='Result_prot_CVD_cox_M2.csv')
