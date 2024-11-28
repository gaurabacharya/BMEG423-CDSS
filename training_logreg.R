## Load some libraries I need to manipulate the data more easily
library('dplyr')
library('reshape2') 
library('zoo')
library('pROC')

## Clear workspace
rm(list=ls())

## Load all the data so we can quickly combine it and explore it. 
load(file.path("..","Data","training_2024-11-04.RData"))
SEPSISdat<-sepsis_data
load(file.path("..","Data","testing_2024-11-04.RData"))
nrow(SEPSISdat) # should be n=200112
SEPSISdat_test<-sepsis_data
rm("sepsis_data")

## Forward-fill missing values
SEPSISdat <- SEPSISdat %>% group_by(patient) %>% mutate_all(funs(na.locf(., na.rm = FALSE)))
SEPSISdat_test <- SEPSISdat_test %>% group_by(patient) %>% mutate_all(funs(na.locf(., na.rm = FALSE)))

## Get reference ranges for variables using 
## only non-sepsis patients as 'normal'
SEPSISdat_NOsepsis <- SEPSISdat[!SEPSISdat$patient %in% unique(SEPSISdat$patient[SEPSISdat$SepsisLabel==1]),]
SEPSISdat_NOsepsis <- SEPSISdat_NOsepsis[SEPSISdat_NOsepsis$ICULOS>1,seq(2,ncol(SEPSISdat_NOsepsis)-2)]
meanSEPSISdat <- round(sapply(SEPSISdat_NOsepsis,mean,na.rm=T),2)
sdSEPSISdat <- round(sapply(SEPSISdat_NOsepsis,sd,na.rm=T),2)

## Obtain the z-scores for all the variables
cols <- colnames(SEPSISdat)
cols <- cols[!(cols %in% c("patient","SepsisLabel","Sex"))]
SEPSISdat_zScores <- SEPSISdat[,2:ncol(SEPSISdat)]
SEPSISdat_test_zScores <- SEPSISdat_test[,2:ncol(SEPSISdat)]
for (c in cols){
  SEPSISdat_zScores[[c]] <- (SEPSISdat[[c]]-meanSEPSISdat[[c]])/sdSEPSISdat[[c]]
  SEPSISdat_test_zScores[[c]] <- (SEPSISdat_test[[c]]-meanSEPSISdat[[c]])/sdSEPSISdat[[c]]
}

## Replace values still missing with the mean
for (c in cols){
  SEPSISdat_zScores[[c]][is.na(SEPSISdat_zScores[[c]])]<-0
  SEPSISdat_zScores[[c]][is.infinite(SEPSISdat_zScores[[c]])]<-0
  SEPSISdat_test_zScores[[c]][is.na(SEPSISdat_test_zScores[[c]])]<-0
  SEPSISdat_test_zScores[[c]][is.infinite(SEPSISdat_test_zScores[[c]])]<-0
}

## Build a linear regression model using all the training data
cnames<-colnames(SEPSISdat_zScores)
form <- as.formula(paste0("SepsisLabel ~ ",paste0(cols,sep="",collapse="+")))
logReg <- glm(form,data=SEPSISdat_zScores,family=binomial(link='logit'))
summary(logReg$coefficients)
logReg_const <- logReg$coefficients[1]
logReg_coeffs <- logReg$coefficients[2:length(logReg$coefficients)]

## Quick but not necessarily great way to find a threshold
SEPSISdat_zScores$probSepsis <- predict(logReg,newdata=SEPSISdat_zScores,type=c("response"))
SEPSISdat_test_zScores$probSepsis <- predict(logReg,newdata=SEPSISdat_test_zScores,type=c("response"))
roc_logReg <- roc(SepsisLabel ~ probSepsis,data=SEPSISdat_zScores)
thresh<-coords(roc_logReg, "b", best.method="youden", input = "threshold", transpose = T,
               ret = c("threshold", "sensitivity","specificity","ppv","npv","fp","tp","fn","tn"))
thresh
# Plot the AUC
plot(roc_logReg,main=paste0('AUC=',round(roc_logReg$auc,3)))
roc_test <- roc(SepsisLabel ~ probSepsis,data=SEPSISdat_test_zScores)
plot(roc_test,add=T,col='red')
text(0.3,0.3,paste0('AUC_test=',round(roc_test$auc,3)),col="red")

## Report the values to put into my get_sepsis_score's load_sepsis_model function
myModel<- NULL
myModel$x_mean <- as.vector(meanSEPSISdat)
myModel$x_std <- as.vector(sdSEPSISdat)
myModel$const <- round(logReg_const,5)
myModel$coeffs <- round(as.vector(logReg_coeffs),5)
myModel$thresh <- round(thresh[1],3)
dput(myModel)

## Quick validation of score
x_norm <- SEPSISdat_zScores[200,1:23]
correct <- predict(logReg,x_norm,type="response")
score <- plogis(myModel$const + sum(x_norm * myModel$coeffs))
print(c(predict=correct,calculated=score,difference=abs(correct-score)<1e-5))

## Not particularly fast way to get utility score but it works
# This has no causality guarantee so please check with the official tool
source('evaluate_sepsis_score.R')
SEPSISdat_zScores$patient <- SEPSISdat$patient
SEPSISdat_zScores$PredictedProbability <- SEPSISdat_zScores$probSepsis
SEPSISdat_zScores$PredictedLabel <- SEPSISdat_zScores$probSepsis > myModel$thresh
SEPSISdat_test_zScores$patient <- SEPSISdat_test$patient
SEPSISdat_test_zScores$PredictedProbability <- SEPSISdat_test_zScores$probSepsis
SEPSISdat_test_zScores$PredictedLabel <- SEPSISdat_test_zScores$probSepsis > myModel$thresh
evaluate_sepsis_score(SEPSISdat_zScores)[5]
evaluate_sepsis_score(SEPSISdat_test_zScores)[5]
