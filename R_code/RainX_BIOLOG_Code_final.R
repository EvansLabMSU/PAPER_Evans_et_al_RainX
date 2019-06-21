##Code for Community level physiological profiles (as measured by BIOLOG/ECOLOG plates). 
##First examines diversity (univariate, as reported in Fig 2b) and shifts in substrate use (Fig 5, PCAs).
##Note that figures were generated in R, but tweaked and finished in Illustrator
##Associated file: "RainX ECO metafile_se_forR.csv"
#Code assembled by: Sarah Evans evanssa6@msu.edu
#Associated with publication Evans et al. 2019. "Dispersal alters soil microbial community response to drought" in Env Microbiology

#Load relevant libraries
library(reshape2)
library(vegan)
library(plyr)

#Choose directory - user-defined
#setwd("~/Dropbox/_HOME TRANSFER/Active projects/RainX/Data/BIOLOG")

#Load in file
Data<-read.csv("RainX ECO metafile_se_forR.csv",  header=TRUE) #Current name: 

## Data manipulation in preparation for a separate univariate and multivariate analysis

#"Swing" CarbonSource into columns along the top, so we have a  data table amenable to analysis with multivariate tools
#(analogous to community data, OTUs go along the top, samples along the side or in a separate mapping file)
DataMultivar<-dcast(Data, Sample+Label.ID + Rain.Trt+Disp.Trt+Soil ~ CarbonSource, value.var="Final")
DataOnly<-DataMultivar[,6:37] #Separate predictor and response variables
pMap<-DataMultivar[,1:5]
DataOnly<-cbind(pMap[,1],DataOnly) 

############ Diversity ######################
Diversity<-diversity(DataOnly,index="shannon",MARGIN=1)
#these data shown in Fig2b. See other code for Fig 2 for plotting of these data

############ Multivariate analysis ##################
clay<-DataMultivar[DataMultivar$Soil=="clay",] #Only analyzing clay results in this paper
map<-clay[,1:5] #Converting clay data to two matrices: Mapping and substrate (Comm)
comm<-clay[,6:36]
comm<-comm/rowSums(comm) #Relative abundance
All<-paste(map$Soil, map$Disp.Trt, map$Rain.Trt, sep="_") #Create a name that includes all trts
map<-cbind(map,All) #add the new name
dist<-vegdist(comm, method='bray', na.rm=TRUE)

#Principal component analysis
apcoa<-cmdscale(dist, k=2, eig=FALSE, add=FALSE, x.ret=FALSE)
plot<-ordiplot(apcoa, type="n")
vecpcoa<-envfit(apcoa,comm)
plot(vecpcoa)
vecpcoa

points(plot, "sites", pch = 15, col = "red", map$All=="clay_Nonsterile_reduced") #filled, red, square
points(plot, "sites", pch = 15, col = "blue",map$All=="clay_Nonsterile_ambient") #filled, blue, square
points(plot, "sites", pch = 0, col = "red",map$All=="clay_Sterile _reduced") #open, red, square
points(plot, "sites", pch = 0, col = "blue",map$All=="clay_Sterile _ambient") #open, red, square
points(plot, "sites", pch = 16, col = "black",map$All=="clay_t0_ambient") #open, red, square
#Note- finished this figure in illustrator
#Convert to black and white symbols in Illustrator, eliminate some arrows

######### Statistics ############
#First remove T0 to examine only treatment effects
clay.comm.trt<-comm[map$Disp.Trt!="t0",]
clay.map.trt<-map[map$Disp.Trt!="t0",]
clay.comm.trt<-clay.comm.trt/rowSums(clay.comm.trt) #Rel abun
#Statistical test
adonis(clay.comm.trt ~ Disp.Trt*Rain.Trt, clay.map.trt) 

######### Mean distance analysis ###########
#This converts distance matrices (triangular) into rectangular matrices so I can calculate the mean distance between different groups
#Code does not include populating a table- did this manually
dist.b<-vegdist(clay.comm.trt, method='bray')
dist.b.M<-as.matrix(dist.b)

### Within Dispersal, distance between Ambient and Drought
colnames(dist.b.M)<-clay.map.trt$Label.ID
rownames(dist.b.M)<-clay.map.trt$Label.ID
Disp<-dist.b.M[clay.map.trt$Disp.Trt=="Nonsterile",clay.map.trt$Disp.Trt=="Nonsterile"]
clay.map.trt.Disp<-clay.map.trt[clay.map.trt$Disp.Trt=="Nonsterile",]
DispRed<-Disp[clay.map.trt.Disp$Rain.Trt=="reduced",]
DispRed<-DispRed[,1:8]
NoDisp<-dist.b.M[17:32, 17:32]
clay.map.trt.NoDisp<-clay.map.trt[17:32,]
NoDispRed<-NoDisp[clay.map.trt.NoDisp$Rain.Trt=="reduced",]
NoDispRed<-NoDispRed[,1:8]
groups<-cbind(c(NoDispRed),c(DispRed)) #Put in groups:
colnames(groups)<-c("NoDisp","Disp")
ttest<-t.test(groups)
mean(NoDispRed) #output means- record in different doc
mean(DispRed)
NoDispSE<-sd(groups[,1])/sqrt(sum(!is.na(groups[,1])))
DispSE<-sd(groups[,2])/sqrt(sum(!is.na(groups[,2])))

### Within Ambient or Reduced, distance between Disp and No
colnames(dist.b.M)<-clay.map.trt$Label.ID
rownames(dist.b.M)<-clay.map.trt$Label.ID
Reduced<-dist.b.M[clay.map.trt$Rain.Trt=="reduced",clay.map.trt$Rain.Trt=="reduced"]
dim(Reduced)
clay.map.trt.Red<-clay.map.trt[clay.map.trt$Rain.Trt=="reduced",]
ReducedNonsterile<-Reduced[clay.map.trt.Red$Disp.Trt=="Nonsterile",]
ReducedNonsterile<-ReducedNonsterile[,1:8]
Ambient<-dist.b.M[clay.map.trt$Rain.Trt=="ambient",clay.map.trt$Rain.Trt=="ambient"]
clay.map.trt.Amb<-clay.map.trt[clay.map.trt$Rain.Trt=="ambient",]
AmbientNonsterile<-Ambient[clay.map.trt.Amb$Rain.Trt=="ambient",]
AmbientNonsterile<-AmbientNonsterile[,1:8]
groupsAmbRed<-cbind(c(AmbientNonsterile),c(ReducedNonsterile)) #Put in groups
colnames(groupsAmbRed)<-c("Amb","Red")
ttestAmbRed<-t.test(groupsAmbRed)
mean(AmbientNonsterile)
mean(ReducedNonsterile)
AmbientNonsterileSE<-sd(groupsAmbRed[,1])/sqrt(sum(!is.na(groupsAmbRed[,1])))
ReducedNonsterileSE<-sd(groupsAmbRed[,2])/sqrt(sum(!is.na(groupsAmbRed[,2])))

