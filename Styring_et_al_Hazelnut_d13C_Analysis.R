########################################################################
### R script for statistical analysis of hazelnut shell d13C values  ###
###   Styring et al. 2024 Frontiers in Environmental Archaeology     ###
### Title: "Carbon isotope values of hazelnut shells: a new proxy    ### 
###                       for canopy density"                        ###
###   This R code runs all analyses presented in the manuscript      ###
###      Authors: A Styring, amy.styring@arch.ox.ac.uk               ###
########################################################################

# load the libraries into R
library(readxl)
library(plyr)
library(nlme)
library(MuMIn)
library(car)
library(lsmeans)
library(pgirmess)
library(clinfun)

# load data
data <- read_excel("Styring_et_al_HazelnutData.xlsx")

#####################################################################
###	3.1 Variability in d13C values within single hazelnut shells	###
#####################################################################

# Select data for nuts from which multiple samples were analysed 
nutdat<-subset(data, `Intra-tree/nut` %in% c(1,2,3))

# Get summary data for intra-nut variability
summary <- ddply(nutdat, c("SampleID"), 
                 function(x) c(d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), n=nrow(x)))

mean(summary$diffd13C)

#########################################################################
### 3.2 Variability in hazelnut shell d13C values within single trees ###
#########################################################################

# Select data for a sub-set of trees from each site with multiple nuts (n = 5)
tree.dat<-subset(data, `Intra-tree/nut` %in% "Y")
 
# Get summary data for intra-tree variability
se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(tree.dat$normd13C, INDEX = tree.dat$Location, FUN = se)
error <- data.frame(error)
CI<-error*1.96

summary <- ddply(tree.dat, c("Tree"), 
                 function(x) c(d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), n=nrow(x)))

summary <- cbind(summary, error, CI)
colnames(summary)[7]<-"CI"

max(summary$diffd13C)
mean(summary$diffd13C)
min(summary$CI)
max(summary$CI)
mean(summary$CI)

min(summary$sdC)
max(summary$sdC)
mean(summary$sdC)

############################################################################
### 3.3 The effect of canopy density on hazelnut shell d13C values (LAI) ###
############################################################################

#Select modern nutshell data
batch<-subset(data, Date %in% c("2021","2022"))

# Test for normality
shapiro.test(batch$normd13C)

# For Fig. 5 regression line
mod <-lm(normd13C ~ LAI, data=batch)
LAI<-seq(min(batch$LAI)-2, max(batch$LAI)+2, length.out=100)
preds<-predict(mod, newdata=data.frame(LAI), interval='confidence')
polygon(c(rev(LAI),LAI), c(rev(preds[,3]), preds[,2]), col=c(makeTransparent("grey80",120)), border=NA)
abline(mod, lty=2)

# Choose a model by AIC

interceptOnly <- gls(normd13C ~ 1, data=batch, method="ML")

randomInterceptOnly <- lme(normd13C ~ 1, data=batch, random=~1|Site, method="ML")

randomInterceptOnlySiteLoc <- lme(normd13C ~ 1, data=batch, random=~1|Site/Location, method="ML")

interceptLAI <- gls(normd13C ~ LAI, data=batch, method="ML")

interceptLAIDate <- gls(normd13C ~ LAI + Date, data=batch, method="ML")

randomInterceptLAI_Site <- lme(normd13C ~ LAI, data=batch, random=~1|Site, method="ML")

randomInterceptLAI_Loc <- lme(normd13C ~ LAI, data=batch, random=~1|Location, method="ML")

randomInterceptLAI_SiteLoc <- lme(normd13C ~ LAI, data=batch, random=~1|Site/Location, method="ML")

randomInterceptLAIDate_Site <- lme(normd13C ~ LAI + Date, data=batch, random=~1|Site, method="ML")

randomInterceptLAIDate_SiteLoc <- lme(normd13C ~ LAI + Date, data=batch, random=~1|Site/Location, method="ML")

anova(interceptOnly, randomInterceptOnly, randomInterceptOnlySiteLoc, interceptLAI, interceptLAIDate, 
      randomInterceptLAI_Site, randomInterceptLAI_Loc, randomInterceptLAI_SiteLoc, randomInterceptLAIDate_Site, 
      randomInterceptLAIDate_SiteLoc)

# model found: randomInterceptLAI_Loc
# normd13C ~ LAI, data=batch, random=~1|Location
summary(randomInterceptLAI_Loc)
intervals(randomInterceptLAI_Loc)

# Check normality of residuals
plot(randomInterceptLAI_Loc)
qqnorm(residuals(randomInterceptLAI_Loc)) # Residuals are normal

# Get conditional R2 
r.squaredGLMM(randomInterceptLAI_Loc) # Conditional R2 = 0.52


###################################################################################
### 3.3 The effect of canopy density on hazelnut shell d13C values (Categories) ###
###################################################################################

# Select modern nutshell data
batch<-subset(data, Date %in% c("2021","2022"))

# Check for normality
by(batch$normd13C, batch$LAI_bin, shapiro.test) # d13C values are normally distributed by openness 
quartz(6,6)
qqPlot(batch$normd13C[batch$LAI_bin=="Open"])
qqPlot(batch$normd13C[batch$LAI_bin=="Semi-open"])
qqPlot(batch$normd13C[batch$LAI_bin=="Closed"])
# Data are normally distributed

# Check for equality of variance
leveneTest(batch$normd13C, batch$LAI_bin) # Variances are similar between categories of canopy

# Perform nested anova
model<-lme(normd13C ~ LAI_bin, random=~1|Location, data=batch, method="REML")
anova.lme(model, type="sequential", adjustSigma=F) # There is a significant effect of canopy

# Calculate least squares mean of nutshell d13C values
leastsquare = lsmeans(model,
                      pairwise ~ LAI_bin,
                      adjust="tukey")       #  Tukey-adjusted comparisons

leastsquare

#########################################################################################
###	Imputing (predicting) canopy density from d13C: testing validity on training data	###
#########################################################################################

# Select modern nutshell data
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
  
  model<-lme(normd13C~LAI_bin, random=~1|Location, data=train.dat, method='REML')
  beta=fixed.effects(model)
  
  #we will overwrite the true openness levels so save them
  Can.actual=test.dat$LAI_bin
  
  #impute levels by choosing the canopy density that puts the fitted D13C closest to the observed D13C
  #note that we don't use the random effects for Location in prediction as we will be using these
  #models to predict on archaeological data with no Location values.
  
  test.dat$LAI_bin[]="Open"
  Xm<-model.matrix(normd13C~1+LAI_bin,data=test.dat)
  yhat.open<-Xm%*%beta                     #predict(model,newdata=test.dat)
  test.dat$LAI_bin[]="Semi-open"
  Xm<-model.matrix(normd13C~LAI_bin,data=test.dat)
  yhat.semi<-Xm%*%beta                      #predict(model,newdata=test.dat)
  test.dat$LAI_bin[]="Closed"
  Xm<-model.matrix(normd13C~LAI_bin,data=test.dat)
  yhat.closed<-Xm%*%beta                      #predict(model,newdata=test.dat)
  yhat=cbind(yhat.open,yhat.semi,yhat.closed)
  colnames(yhat)<-c("Open","Semi-open","Closed")
  imp<-apply((yhat-test.dat$normd13C)^2,1,which.min)
  test.dat$Canopy.imputed<-factor(colnames(yhat)[imp],levels=c("Open","Semi-open","Closed"))
  
  conf=conf+ftable(test.dat$Canopy.imputed~Can.actual)
  
}

#confusion matrix and overall score
round(conf/matrix(rep(apply(conf,1,sum),3),3,3),2) ##This gives Table 3
round(sum(diag(conf))/sum(conf),2)

#####################################################################################
###	    Imputing (predicting) canopy density from D13C: archaeological data	      ###
#####################################################################################

# Select modern nutshell data
mod.dat<-subset(data, Date %in% c("2021","2022"))
nr=dim(mod.dat)[1]
mod.dat$LAI_bin<-factor(mod.dat$LAI_bin,levels=c("Open","Semi-open","Closed"))

#fit the model using the modern data
model<-lme(D13C ~ 1 + LAI_bin, random=~1|Location, data=mod.dat, method="REML")
beta=fixed.effects(model)

# Select archaeological data
batchA<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' ) # Mesolithic data
batchB<-subset(data, Period!="Mesolithic")
batchC<-subset(batchB, Period!="Modern") # Other archaeological data

# Correct charred nutshell d13C values for charring (BatchC). This value is from Holguin et al. (in prep.)
batchC$normd13C<-batchC$normd13C-0.51
batchC$D13C<-(batchC$d13CCO2-batchC$normd13C)/(1+(batchC$normd13C/1000))
old.dat<-rbind(batchA, batchC)

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

#####################################################################################
###	3.4 Change in archaeological hazelnut shell d13C values during the Mesolithic ###
#####################################################################################

#Select Mesolithic (non-acid treated) nutshell data
batch<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' )

# Get summary data for intra-site variability
se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(batch$normd13C, INDEX = batch$Site, FUN = se)
error <- data.frame(error)
CI<-error*1.96

summary <- ddply(batch, c("Site"), 
                 function(x) c(d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), n=nrow(x),D13C=mean(x$D13C)))

summary <- cbind(summary, error, CI)
colnames(summary)[8]<-"CI"

max(summary$diffd13C)
mean(summary$diffd13C)
min(summary$CI)
max(summary$CI)
mean(summary$CI)

# Check for normality
by(batch$normd13C, batch$Site, shapiro.test) # d13C values are normally distributed within sites 
quartz(6,6)
qqPlot(batch$normd13C[batch$Site=="Slabälta 1"])
qqPlot(batch$normd13C[batch$Site=="Ringsjöholm"])
qqPlot(batch$normd13C[batch$Site=="Rönneholm 6_1"])
qqPlot(batch$normd13C[batch$Site=="Rönneholm 10_3"])
# Data are normally distributed

# Check for equality of variance
leveneTest(batch$normd13C, batch$Site) # Variances are similar between sites

# Analysis of variance comparing d13C values among sites
anova<-aov(batch$normd13C ~ batch$Site)
summary(anova)
sqrt((24.54-(3*1.655))/(24.54+59.58+1.655)) # This is the effect size, w 

# Pairwise comparisons using Tukey HSD
TukeyHSD(anova)

###################################################################################################
###	3.4 Change in archaeological hazelnut shell d13C values between the Mesolithic and Iron Age ###
###################################################################################################

# Select Mesolithic (non-acid treated) nutshell data
batchA<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' )

# Select other archaeological charred nutshell data
batchB<-subset(data, Period %in% c("Neolithic","Bronze Age", "Iron Age"))

# Correct charred nutshell d13C values for charring (batchB). This value is from Holguin et al. (in prep.)
batchB$normd13C<-batchB$normd13C-0.51
batchB$D13C<-(batchB$d13CCO2-batchB$normd13C)/(1+(batchB$normd13C/1000))

# Combine data from all archaeological nutshell samples
batch<-rbind(batchA, batchB)

# Get summary data for intra-period (Mesolithic, Neolithic, Bronze Age, Iron Age) variability
se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(batch$normd13C, INDEX = batch$Period, FUN = se)
error <- data.frame(error)
CI<-error*1.96

summary <- ddply(batch, c("Period"), 
                 function(x) c(d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), 
                               n=nrow(x),D13C=mean(x$D13C)))

summary <- cbind(summary, error, CI)
colnames(summary)[8]<-"CI"

# Remove Neolithic data because the samples are too few
batch<-subset(batch, Period!='Neolithic' ) 

# Check for normality
by(batch$D13C, batch$Period, shapiro.test) # d13C values are not normally distributed within periods 
quartz(6,6)
qqPlot(batch$D13C[batch$Period=="Mesolithic"])
qqPlot(batch$D13C[batch$Period=="Bronze Age"])
qqPlot(batch$D13C[batch$Period=="Iron Age"])
# Mesolithic data are not normally distributed

leveneTest(batch$D13C, batch$Period) ##Variances are similar between sites

kruskal.test(D13C ~ Period, data=batch)
batch$Period<-factor(batch$Period, levels=c("Mesolithic", "Bronze Age", "Iron Age"))
jonckheere.test(batch$D13C, as.numeric(batch$Period))
batch$Ranks<-rank(batch$D13C)
by(batch$Ranks, batch$Period, mean)

kruskalmc(D13C ~ Period, data=batch)
