# Counterfactual-based mediation analysis, using CVD as an example
# Chun Shen, 2024

# Load R packages
library('data.table')
library('lubridate')                                        
library('RNOmni')
library('survival')
library('foreach')
library('doParallel')
library('timereg')

###Load data
# time gap: N participants * M proteins
timeGap <- fread('olink_timeGap.csv')
colnames(timeGap)[1] <- 'ID'
#protein data and other covariates
data_prot <- fread('protein_UKB.csv',data.table = F)
protName <- colnames(data_prot)[2:2921]
protName2 <- chartr('-','_',protName)
colnames(data_prot)[2:2921] <- protName2

# MR sig proteins
sig_toprot <- fread('Result_IVW_sig_toProt.csv')
prot_used <- sig_toprot$outcome

data_prot_use <- data_prot[,c('ID',prot_used,'SI_2c','LO_2c','sex','eth2','age','site','Batch',
                              'edu_4c','smokeNow','alc_2c','bmi','inc_2c',paste0("PC", 1:20))]

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
  dat <- subset(Data_merge3,select=c('ID',prot_used[i],'outcome','fultime','LO_2c','age','sex','site','Batch',
                                           'eth2','edu_4c','smokeNow','alc_2c','bmi','inc_2c',
                                           paste0("PC", 1:20)))
  dat_t <- subset(timeGap,select = c('ID',prot_used[i]))
  colnames(dat_t)[2] <- 'timeGap'
  mydata <- na.omit(merge(dat,dat_t,by='ID'))
  mydata$mediator <- RankNorm(mydata[,2])

  # category
  mydata$sex <- as.factor(mydata$sex)
  mydata$eth2 <- as.factor(mydata$eth2)
  mydata$edu_4c <- as.factor(mydata$edu_4c)
  mydata$inc_2c <- as.factor(mydata$inc_2c)
  mydata$smokeNow <- as.factor(mydata$smokeNow)
  mydata$alc_2c <- as.factor(mydata$alc_2c)
  mydata$site <- as.factor(mydata$site)
  mydata$Batch <- as.factor(mydata$Batch)
  
  ##############
  #Step 1: Mediator (protein) ~ LO+covs
  mydata$LOtemp <- mydata$LO_2c
  f1 <- paste0('mediator ~ LOtemp+age+sex+eth2+edu_4c+smokeNow+inc_2c+alc_2c+bmi+timeGap+site+Batch+',paste0("PC", 1:20, collapse = "+"))
  fitM <- lm(as.formula((f1)), data=mydata) 
  fitMvariance <- summary(fitM)$sigma^2
  
  #Step 2: construct new variables id and LO_2c where the last corresponds to the value of the exposure relative to the indirect path.
  mydata1 <- mydata
  mydata2 <- mydata
  mydata1$LOstar <- mydata$LO_2c
  mydata2$LOstar <- 1-mydata$LO_2c
  newdata <- rbind(mydata1,mydata2)
  
  #Step 3: compute weights
  newdata$LOtemp <- newdata$LO_2c
  tempDir <- dnorm(newdata$mediator,
                  mean=predict(fitM,type='response',newdata=newdata),
                  sd=sqrt(fitMvariance))
  newdata$LOtemp <- newdata$LOstar
  tempIndir <- dnorm(newdata$mediator,
                    mean=predict(fitM,type='response',newdata=newdata),
                    sd=sqrt(fitMvariance))
  newdata$weightM <- tempIndir/tempDir

  #Step 4: Aalens model
  f2 <- 'Surv(fultime, outcome) ~ const(factor(LO_2c))+const(factor(LOstar))+const(age)+const(sex)+const(eth2)+const(edu_4c)+const(smokeNow)+const(inc_2c)+const(alc_2c)+const(bmi)+const(timeGap)+const(site)'
  fitYaalen <- aalen(as.formula(f2),data=newdata,weights=newdata$weightM,clusters=newdata$ID)

  sumAalen_1 <- as.data.frame(capture.output(summary(fitYaalen)))
  sumAalen_2 <- sumAalen_1[17:52,]
  sumAalen_3 <- strsplit(sumAalen_2, " ")
  filtered_data <- lapply(sumAalen_3, function(x) x[x != ""])
  max_cols <- max(sapply(filtered_data, length))
  padded_data <- lapply(filtered_data, function(x) c(x, rep(NA, max_cols - length(x))))
  result_Aalen <- as.data.frame(do.call(rbind, padded_data))
  
  #direct and indirect effect
  result <- data.frame(protname = prot_used[i],disease = 'CVD',
                      dirEff=result_Aalen[1,2],dirRSE=result_Aalen[1,4],dirZ=result_Aalen[1,5],dirP=result_Aalen[1,6],
                      indirEff=result_Aalen[2,2],indirRSE=result_Aalen[2,4],indirZ=result_Aalen[2,5],indirP=result_Aalen[2,6])

  result
}

write.csv(Result_CVD_cmed2, file='Result_causalMediation_CVD_M2_v2.csv', row.names = F)


