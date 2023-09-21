
#####################################################################
####			       SUMMARY STATS - INTRA-TREE VARIATION		    	    ###
#####################################################################
setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")

library(readxl)
data <- read_excel("HazelnutsAllData.xlsx")

tree.dat<-subset(data, `Intra-tree/nut` %in% "Y")
## This is selecting a sub-set of trees from each site with multiple nuts 

se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(tree.dat$normd13C, INDEX = tree.dat$Tree, FUN = se)
error <- data.frame(error)
CI<-error*1.96

library(plyr)
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

#####################################################################
####			                  d13C v LAI		    	                  ###
#####################################################################
setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")

library(readxl)
data <- read_excel("HazelnutsAllData.xlsx")

batch<-subset(data, Date %in% c("2021","2022"))

shapiro.test(batch$normd13C)

## For Fig. 5 regression line
mod <-lm(normd13C ~ LAI, data=batch)
LAI<-seq(min(batch$LAI)-2, max(batch$LAI)+2, length.out=100)
preds<-predict(mod, newdata=data.frame(LAI), interval='confidence')
polygon(c(rev(LAI),LAI), c(rev(preds[,3]), preds[,2]), col=c(makeTransparent("grey80",120)), border=NA)
abline(mod, lty=2)

## Finding the best model
library(nlme)
library(lme4)
library(MuMIn)
interceptOnly <- gls(normd13C ~ 1, data=batch, method="ML")

randomInterceptOnly <- lme(normd13C ~ 1, data=batch, random=~1|Location, method="ML")

interceptLAI <- gls(normd13C ~ LAI, data=batch, method="ML")

interceptLAISite <- gls(normd13C ~ LAI + Site, data=batch, method="ML")

interceptLAIDate <- gls(normd13C ~ LAI + Date, data=batch, method="ML")

randomInterceptLAIT <- lme(normd13C ~ LAI, data=batch, random=~1|Location, method="ML")
summary(randomInterceptLAIT)
intervals(randomInterceptLAIT)
plot(randomInterceptLAIT)
qqnorm(residuals(randomInterceptLAIT)) ## Residuals are normal
r.squaredGLMM(randomInterceptLAIT) 
## R2GLMM = 0.52

randomInterceptLAISiteTree <- lme(normd13C ~ LAI + Site, data=batch, random=~1|Location, method="ML")
summary(randomInterceptLAISiteTree)
plot(randomInterceptLAISiteTree) ## Residuals are normal

randomInterceptLAIDateTree <- lme(normd13C ~ LAI + Date, data=batch, random=~1|Tree, method="ML")

randomInterceptLAITree <- lme(normd13C ~ LAI, data=batch, random=~1|Site/Tree, method="ML")

randomInterceptLAISiteDateTree <- lme(normd13C ~ LAI + Site + Date, data=batch, random=~1|Tree, method="ML")

anova(interceptOnly, randomInterceptOnly, interceptLAI, interceptLAISite, interceptLAIDate, 
      randomInterceptLAIT, randomInterceptLAISiteTree, randomInterceptLAIDateTree,
      randomInterceptLAISiteDateTree)
     
#####################################################################
##						      STATISTICS	- d13C and Openness				         ##
#####################################################################
require(dplyr)
require(tidyr)
require(ordinal)
require(lme4)
require(ggplot2)
require(car)
require(MuMIn)
require(nlme)
library(car)
library(readxl)

setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")
data <- read_excel("HazelnutsAllData.xlsx")

batch<-subset(data, Date %in% c("2021","2022"))

by(batch$normd13C, batch$LAI_bin, shapiro.test) ##d13C values are normally distributed by openness 
quartz(6,6)
qqPlot(batch$normd13C[batch$LAI_bin=="Open"])
qqPlot(batch$normd13C[batch$LAI_bin=="Semi-open"])
qqPlot(batch$normd13C[batch$LAI_bin=="Closed"])
## Data are normally distributed

leveneTest(batch$normd13C, batch$LAI_bin) ##Variances are similar between types of canopy

## Nested anova
model<-lme(normd13C ~ LAI_bin, random=~1|Location, data=batch, method="REML")
anova.lme(model, type="sequential", adjustSigma=F) ##Significant effect of canopy

library(multcompView)
library(lsmeans)
library(multcomp)
library(lmerTest)

leastsquare = lsmeans(model,
                      pairwise ~ LAI_bin,
                      adjust="tukey")       ###  Tukey-adjusted comparisons

leastsquare

model<-lmer(normd13C ~ LAI_bin + (1|Location), data=batch, REML=T)
anova(model)
rand(model)
difflsmeans(model, test.effs="LAI_bin")

## Linear model
model<-lme(normd13C ~ LAI_bin, random=~1|Location, data=batch, method="REML")
summary(model)

#####################################################################
##						      STATISTICS	- Mesolithic sites				         ##
#####################################################################
library(readxl)
setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")
data <- read_excel("HazelnutsAllData.xlsx")

batch<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' )

se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(batch$normd13C, INDEX = batch$Site, FUN = se)
error <- data.frame(error)
CI<-error*1.96

library(plyr)
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

library(car)

by(batch$normd13C, batch$Site, shapiro.test) ##d13C values are normally distributed by openness 
quartz(6,6)
qqPlot(batch$normd13C[batch$Site=="Slabälta 1"])
qqPlot(batch$normd13C[batch$Site=="Ringsjöholm"])
qqPlot(batch$normd13C[batch$Site=="Rönneholm 6_1"])
qqPlot(batch$normd13C[batch$Site=="Rönneholm 10_3"])
## Data are normally distributed

leveneTest(batch$normd13C, batch$Site) ##Variances are similar between sites

anova<-aov(batch$normd13C ~ batch$Site)
summary(anova)
sqrt((24.54-(3*1.655))/(24.54+59.58+1.655)) ## This is the effect size, w 

TukeyHSD(anova)

#####################################################################
##						  STATISTICS	- Archaeological sites				         ##
#####################################################################
library(readxl)
setwd("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Results")
data <- read_excel("HazelnutsAllData.xlsx")

batchA<-subset(data, `Acid-treated`=="N" & Period=='Mesolithic' )
batchB<-subset(data, Period %in% c("Neolithic","Bronze Age", "Iron Age"))
## Correct charred nutshell d13C values for charring (BatchB)
batchB$normd13C<-batchB$normd13C-0.51
batchB$D13C<-(batchB$d13CCO2-batchB$normd13C)/(1+(batchB$normd13C/1000))
batch<-rbind(batchA, batchB)

se <- function(x){          # x is a dummy variable for the function
  s <- sd(x, na.rm = TRUE)  # calculate the standard deviation
  n <- length(x[!is.na(x)]) # calculate the sample size
  se <- s/sqrt(n)           # standard error
  se                        # what the function will return
}

error <- tapply(batch$normd13C, INDEX = batch$Period, FUN = se)
error <- data.frame(error)
CI<-error*1.96

library(plyr)
summary <- ddply(batch, c("Period"), 
                 function(x) c(d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), 
                               n=nrow(x),D13C=mean(x$D13C)))

summary <- cbind(summary, error, CI)
colnames(summary)[8]<-"CI"


library(car)
library(clinfun)
library(pgirmess)

batch<-subset(batch, Period!='Neolithic' )

by(batch$normd13C, batch$Period, shapiro.test) ##d13C values are normally distributed by openness 
quartz(6,6)
qqPlot(batch$normd13C[batch$Period=="Mesolithic"])
qqPlot(batch$normd13C[batch$Period=="Bronze Age"])
qqPlot(batch$normd13C[batch$Period=="Iron Age"])
## Mesolithic data are not normally distributed

leveneTest(batch$normd13C, batch$Period) ##Variances are similar between sites

kruskal.test(normd13C ~ Period, data=batch)
batch$Period<-factor(batch$Period, levels=c("Mesolithic", "Bronze Age", "Iron Age"))
jonckheere.test(batch$normd13C, as.numeric(batch$Period))
batch$Ranks<-rank(batch$normd13C)
by(batch$Ranks, batch$Period, mean)

kruskalmc(normd13C ~ Period, data=batch)
