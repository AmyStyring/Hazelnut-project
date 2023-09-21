
###################################################################
####					        EXTRACTING STANDARDS					            ###
###################################################################
library(readxl)
library(dplyr)
library(multiway)

data1 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 1")
data2 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 2")
data3 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 3")
data4 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 4")
data5 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 5")
data6 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 6")
data7 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 7")
data8 <- read_excel("~/Library/CloudStorage/OneDrive-Nexus365/RESEARCH PROJECTS/Hazelnut Pilot/Paper/Supplementary Table 3.xlsx", sheet = "Data-Session 8")

dat<-rbind(data1, data2, data3, data4, data5, data6, data7, data8)
RawStandards <- dat[grep("SPRUCE|SEAL|CAFF|SALANINE", dat$ID), ]
RawStandards <- RawStandards[ ,c(2,15,19)]
names(RawStandards) <- c("ID", "Runfile","normd13C")

###################################################################
####					        EXTRACTING REPLICATES					            ###
###################################################################
for (i in c(2, 4:14, 16:20)){data[,i] <- as.numeric(as.character(data[,i]))}
data$Runfile <- as.character(data$Runfile)

RepCA<-dat[grep("DA$", dat$ID), ]
RepCA <- RepCA[order(as.character(RepCA$ID)), ]
RepCA <- RepCA[ ,c(2,15, 19)]
names(RepCA) <- c("ID", "Runfile","normd13C_DuplA")
RepCA$ID<-gsub('DA$','',RepCA$ID)

RepCB<-dat[grep("DB$", dat$ID), ]
RepCB <- RepCB[order(as.character(RepCB$ID)), ]
RepCB <- RepCB[ ,c(2,19)]
names(RepCB) <- c("ID", "normd13C_DuplB")
RepCB$ID<-gsub('DB$','',RepCB$ID)

RepCar<-merge(RepCA, RepCB, "ID")

###################   Difference between modern duplicates     #################

RepCar$Diff<-RepCar$normd13C_DuplA-RepCar$normd13C_DuplB
RepCar$Diff<-abs(RepCar$Diff)

################### Mean and Stdev for all analytical sessions #################

Cmean<-aggregate(RawStandards$normd13C, list(RawStandards$ID), mean)
Csd<-aggregate(RawStandards$normd13C, list(RawStandards$ID), sd)
N<-count(RawStandards, "ID")

all.standards<-cbind(N,Cmean, Csd)
all.standards<-all.standards[,c(3,2,4,6)]
names(all.standards)<-c("RM","Number", "d13Cmean","d13Csd")
all.standards$srm<-(all.standards$Number-1)*(all.standards$d13Csd^2)

################################################################################        
##                    Check and calibration standards                         ##
##          Pooled standard deviation of each standard (Ssrm) 						    ##
##                and the degrees of freedom (dfsrm)				   						    ##
################################################################################
dfsrm<-sum(all.standards$Number)-nrow(all.standards)
Ssrm<-sqrt(sum(all.standards$srm)/dfsrm)

################################################################################        
##                            Check standards                                 ##
################################################################################
checkC<-subset(all.standards,all.standards$RM=="SALANINE"|all.standards$RM=="SPRUCE")

CheckS.1<--25.44#### SPRUCE
CheckS.2<--27.18 #### SALANINE

y="SPRUCE" #### change if using different check standards - this should reflect the name of CheckS.1
fun1<-function(x,y) if(x==y) {CheckS.1} else {CheckS.2}
checkC$known<-mapply(fun1, checkC$RM, y)
CheckS.1sd<-0.02 #### SPRUCE
CheckS.2sd<-0.16 #### SALANINE

fun1<-function(x,y) if(x==y) {CheckS.1sd} else {CheckS.2sd}
checkC$knownsd<-mapply(fun1, checkC$RM, y)

checkC$Diff_measured_known<-checkC$d13Cmean-checkC$known
################################################################################
# RMS bias = the root mean square of the difference between the observed mean	 #
# and the known values of standard reference materials which are treated as    #
# unknowns during analysis (Szpak's "check standards").                 			 #
# u_cref = the root mean square of the known standard deviations of the RMs    #
# used as check standards.                                                		 # 
################################################################################

RMSbias<-sqrt(sumsq(checkC$Diff_measured_known)/nrow(checkC))
u_cref<-sqrt(sumsq(checkC$knownsd)/nrow(checkC))

################################################################################		 
##                                ACCURACY	                                  ##
################################################################################

x<-list(RMSbias, u_cref)
u_bias<-sqrt(sumsq(x))

################################################################################		 
##                               SAMPLE REPLICATES	                          ##
################################################################################
RepCar$Sd<-apply(subset(RepCar,select = c("normd13C_DuplA","normd13C_DuplB")),1,sd)

RepCar$Mean<-apply(subset(RepCar,select = c("normd13C_DuplA","normd13C_DuplB")),1,mean)

RepCar$Number<-2## this needs to change if you run a replicate more than twice in a run

RepCar$RepSsrm<-1*(RepCar$Sd^2)
dfrep<-sum(RepCar$Number)-nrow(RepCar)
Srep<-sqrt((sum(RepCar$RepSsrm))/dfrep)

################################################################################
##                                 PRECISION                                  ##
## Random errors within the laboratory = precision ( Menditto et al 2007)     ##
## or repeatability (Carter and Fry 2013) - summed standard deviation of      ##
## all repeated measurements during relevant analytical sessions - including  ##
##  check and calibration RMs (Ssrum) and the replicates (Srep)	              ##
################################################################################

uRw<-sqrt((Ssrm^2)+(Srep^2)/2)

################################################################################		 
##                      Standard uncertainty (Szpak 2017)	                    ##
################################################################################

y<-list(u_bias, uRw)
Uc<-sqrt(sumsq(y))

################################################################################		 
##                      VALUES TO REPORT (cf. Szpak 2017)                     ##
################################################################################
print(uRw) ## precision (u(Rw))
print(u_bias) ## accuracy, or systematic error (u(bias))
print(Uc) ## total analytical uncertainty (Uc)
