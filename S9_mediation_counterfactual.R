#Mediation analysis caculated by counterfactual-based mediation analysis, using CVD as an example
library('data.table')
library('lubridate')                                        
library('RNOmni')
library('survival')
library('foreach')
library('doParallel')
library('timereg')

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

##############
registerDoParallel(20)

Result_CVD_cmed1 <- foreach(i=1:length(prot_use), .combine='rbind') %dopar% {
  #simple model: adjust for age, sex, site
  dat_simple <- na.omit(subset(Data_merge3,select=c('ID',prot_use[i],'outcome','fultime','LO_2c','age','sex')))
  dat_simple$exposure <- RankNorm(dat_simple[,2])
  
  ##############
  #Aalen model
  #Step 1: protein~LO_2c+covs
  dat_simple$LOtemp <- dat_simple$LO_2c
  fitM_1 <- lm(exposure ~ factor(LOtemp)+age+sex, data=dat_simple) 
  fitMvar_1 <- summary(fitM_1)$sigma^2
  
  #Step 2: construct new variables id and LO_2c where the last corresponds to the value of the exposure relative to the indirect path.
  dat_simple1 <- dat_simple
  dat_simple2 <- dat_simple
  dat_simple1$LOstar <- dat_simple$LO_2c
  dat_simple2$LOstar <- 1-dat_simple$LO_2c
  newDat_simple <- rbind(dat_simple1, dat_simple2)
  
  #Step 3: compute weights
  newDat_simple$LOtemp <- newDat_simple$LO_2c
  tempDir_1 <- dnorm(newDat_simple$exposure, 
                     mean=predict(fitM_1, type = "response", newdata=newDat_simple), 
                     sd=sqrt(fitMvar_1) )
  newDat_simple$LOtemp <- newDat_simple$LOstar
  tempIndir_1 <- dnorm(newDat_simple$exposure, 
                       mean=predict(fitM_1, type = "response", newdata=newDat_simple), 
                       sd=sqrt(fitMvar_1) )
  newDat_simple$weight <- tempIndir_1/tempDir_1
  
  #Step 4: the MSM model for direct and indirect effects can be estimated
  #Aalens Model
  fitYaalen_1 <- aalen(Surv(fultime, outcome)~const(factor(LO_2c))+const(factor(LOstar))+const(age)+const(sex), 
                       data=newDat_simple, weights=newDat_simple$weight, clusters=newDat_simple$eid) 
  #clean results
  sumAalen_1 <- as.data.frame(capture.output(summary(fitYaalen_1)))
  sumAalen_2 <- sumAalen_1[17:20,]
  sumAalen_3 <- strsplit(sumAalen_2, " ")
  filtered_data <- lapply(sumAalen_3, function(x) x[x != ""])
  max_cols <- max(sapply(filtered_data, length))
  padded_data <- lapply(filtered_data, function(x) c(x, rep(NA, max_cols - length(x))))
  result_Aalen_1 <- as.data.frame(do.call(rbind, padded_data))
  
  #direct and indirect effect
  result <- data.frame(dirEff=result_Aalen_1[1,2],dirSE=result_Aalen_1[1,3],dirRSE=result_Aalen_1[1,4],dirZ=result_Aalen_1[1,5],dirP=result_Aalen_1[1,6],
                       indirEff=result_Aalen_1[2,2],indirSE=result_Aalen_1[2,3],indirRSE=result_Aalen_1[2,4],indirZ=result_Aalen_1[2,5],indirP=result_Aalen_1[2,6],)
  
  result
}

save(Result_CVD_cmed1,file='Result_causalMediation_CVD.RData')
