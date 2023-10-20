#####################################################################
####			                SINGLE IMPUTATION		    	              ###
#####################################################################

#####################################################################
####		 TESTING VALIDITY OF IMPUTATION ON TRAINING DATA		    	###
#####################################################################

setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")

library(readxl)
data <- read_excel("HazelnutsAllData.xlsx")
mod.dat<-subset(data, Date %in% c("2021","2022"))
mod.dat$LAI_bin<-factor(mod.dat$LAI_bin,ordered=TRUE,levels=c("Open","Semi-open","Closed"))

#make 1000 reps training on 7/8th's of data and testing on 1/8
#keep track of proportion of successes imputed.ML=actual.ML

work.dat=mod.dat

reps=1000
nm=dim(work.dat)[1]
ai=1:nm 
test.frac=1/8;
test.num=round(length(ai)*test.frac)
pc=dm=rep(NA,reps)
conf=matrix(0,3,3) #3 because there are 3 levels open, semi-open, closed

for (r in 1:reps) {

  #form the test and training data subsets
  tei=sample(ai,test.num,replace=F)
  train.dat=work.dat[-tei,]
  test.dat=work.dat[tei,]

  #fit the model using the training data

    model<-lme(D13C~LAI_bin, random=~1|Location, data=train.dat, method='REML')
    beta=fixed.effects(model)

  #we will overwrite the true openness levels so save them
  Can.actual=test.dat$LAI_bin
  
  #impute levels by choosing the canopy density that puts the fitted D13C closest to the observed D13C
  #note that we don't use the random effects for Location in prediction as we will be using these
  #models to predict on archaeological data with no Location values.

  test.dat$LAI_bin[]="Open"
  Xm<-model.matrix(D13C~1+LAI_bin,data=test.dat)
  yhat.open<-Xm%*%beta                     #predict(model,newdata=test.dat)
  test.dat$LAI_bin[]="Semi-open"
  Xm<-model.matrix(D13C~LAI_bin,data=test.dat)
  yhat.semi<-Xm%*%beta                      #predict(model,newdata=test.dat)
  test.dat$LAI_bin[]="Closed"
  Xm<-model.matrix(D13C~LAI_bin,data=test.dat)
  yhat.closed<-Xm%*%beta                      #predict(model,newdata=test.dat)
  yhat=cbind(yhat.open,yhat.semi,yhat.closed)
  colnames(yhat)<-c("Open","Semi-open","Closed")
  imp<-apply((yhat-test.dat$D13C)^2,1,which.min)
  test.dat$Canopy.imputed<-factor(colnames(yhat)[imp],levels=c("Open","Semi-open","Closed"))

  conf=conf+ftable(test.dat$Canopy.imputed~Can.actual)

}

#confusion matrix and overall score
round(conf/matrix(rep(apply(conf,1,sum),3),3,3),2) ##This gives Table 2
round(sum(diag(conf))/sum(conf),2)

#####################################################################
####		 IMPUTATION OF OPENNESS USING ARCHAEOLOGICAL DATA		    	###
#####################################################################

setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")

library(readxl)
data <- read_excel("HazelnutsAllData.xlsx")

library(nlme)

## Subsetting modern data
mod.dat<-subset(batch, Date %in% c("2021","2022"))
nr=dim(mod.dat)[1]
mod.dat$LAI_bin<-factor(mod.dat$LAI_bin,levels=c("Open","Semi-open","Closed"))

## Model
model<-lme(D13C ~ 1 + LAI_bin, random=~1|Location, data=mod.dat, method="REML")
beta=fixed.effects(model)

## Subsetting archaeological data
batchA<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' )
batchB<-subset(data, Period!="Mesolithic")
batchC<-subset(batchB, Period=="Modern")
batchD<-subset(batchB, Period!="Modern")
## Correct charred nutshell d13C values for charring (BatchD)
batchD$normd13C<-batchD$normd13C-0.51
batchD$D13C<-(batchD$d13CCO2-batchD$normd13C)/(1+(batchD$normd13C/1000))
batch<-rbind(batchA, batchC, batchD)
old.dat<-subset(batch, Period!="Modern")

nold=dim(old.dat)[1]
guess.dat=old.dat; guess.dat$LAI_bin=factor(rep(NA,nold),
                                            levels=c("Open","Semi-open","Closed"))
guess.dat$LAI_bin[]="Open"
Xm<-model.matrix(D13C~1+LAI_bin,data=guess.dat)
yhat.open<-Xm%*%beta                      
guess.dat$LAI_bin[]="Semi-open"
Xm<-model.matrix(D13C~1+LAI_bin,data=guess.dat)
yhat.semi<-Xm%*%beta                      
guess.dat$LAI_bin[]="Closed"
Xm<-model.matrix(D13C~1+LAI_bin,data=guess.dat)
yhat.closed<-Xm%*%beta                   
yhat=cbind(yhat.open,yhat.semi,yhat.closed)
colnames(yhat)<-c("Open","Semi-open","Closed")
imp<-apply((yhat-guess.dat$D13C)^2,1,which.min)
old.dat$Openness.imputed<-factor(colnames(yhat)[imp],
                                 ordered=TRUE,levels=c("Open","Semi-open","Closed"))

old.dat <- data.frame(old.dat, old.dat$Openness.imputed)
