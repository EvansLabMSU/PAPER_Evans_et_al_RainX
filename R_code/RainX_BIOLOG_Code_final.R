## SIR data and other univariate similar data for RainX Project
rm(list=ls(all=TRUE))

setwd("~/Dropbox/_HOME TRANSFER/Active projects/RainX/Data/BIOLOG")

#Load in file
Data<-read.csv(file.choose(),  header=TRUE) #Current name: "RainX ECO metafile_se_forR.csv"

head(Data)
dim(Data)
library(reshape2)
library(vegan)
library(plyr)

#"Swing" CarbonSource into columns along the top, so we have a  data table amenable to analysis with multivariate tools
#(analogous to community data, OTUs go along the top, samples along the side or in a separate mapping file)
DataMultivar<-dcast(Data, Sample+Label.ID + Rain.Trt+Disp.Trt+Soil ~ CarbonSource, value.var="Final")
head(DataMultivar)

DataOnly<-cbind(pMap[,1],DataOnly)
pMap<-DataMultivar[,1:5]
DataOnly<-DataMultivar[,6:37]
clay.comm.trt<-clay.comm[clay.map$Disp.Trt!="t0",]
clay.map.trt<-clay.map[clay.map$Disp.Trt!="t0",]
clay.comm.trt<-clay.comm.trt/rowSums(clay.comm.trt) #Rel abun

############ Diversity ###########
Diversity<-diversity(DataOnly,index="shannon",MARGIN=1)

############ Multivariate analysis ##################
attach(DataMultivar)
clay<-DataMultivar[DataMultivar$Soil=="clay",]
DataPA<-decostand(DataOnly,method='pa')
#relative abundance
DataOnlyRelAbun<-DataOnly/rowSums(DataOnly)

#Converting to two matrices: Mapping and substrate (Comm)
clay.map<-clay[,1:5]
clay.comm<-clay[,6:36]

#Going forward keeping T0s in. 
#Name simpler names:
comm<-clay.comm
map<-clay.map
comm<-comm/rowSums(comm)
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
#Convert to black and white symbols in Illustrator

plot(vec) #Put vectors on top

######### Statistics ############
#Using version with no T0
adonis(clay.comm.trt ~ Disp.Trt*Rain.Trt, clay.map.trt) 

######### Mean distance analysis ###########
#Names 
clay.comm.trt
clay.map.trt
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
dim(NoDisp)
clay.map.trt.NoDisp<-clay.map.trt[17:32,]
NoDispRed<-NoDisp[clay.map.trt.NoDisp$Rain.Trt=="reduced",]
NoDispRed<-NoDispRed[,1:8]

#Put in groups:
groups<-cbind(c(NoDispRed),c(DispRed))
colnames(groups)<-c("NoDisp","Disp")
ttest<-t.test(groups)
mean(NoDispRed)
mean(DispRed)

NoDispSE<-sd(groups[,1])/sqrt(sum(!is.na(groups[,1])))
DispSE<-sd(groups[,2])/sqrt(sum(!is.na(groups[,2])))

### Within Ambient or Reduced, distance between Disp and No
colnames(dist.b.M)<-clay.map.trt$Label.ID
rownames(dist.b.M)<-clay.map.trt$Label.ID
#Reduced
Reduced<-dist.b.M[clay.map.trt$Rain.Trt=="reduced",clay.map.trt$Rain.Trt=="reduced"]
dim(Reduced)
clay.map.trt.Red<-clay.map.trt[clay.map.trt$Rain.Trt=="reduced",]
ReducedNonsterile<-Reduced[clay.map.trt.Red$Disp.Trt=="Nonsterile",]
ReducedNonsterile<-ReducedNonsterile[,1:8]

#Ambient
dim(dist.b.M)
Ambient<-dist.b.M[clay.map.trt$Rain.Trt=="ambient",clay.map.trt$Rain.Trt=="ambient"]
dim(Ambient)
clay.map.trt.Amb<-clay.map.trt[clay.map.trt$Rain.Trt=="ambient",]
AmbientNonsterile<-Ambient[clay.map.trt.Amb$Rain.Trt=="ambient",]
AmbientNonsterile<-AmbientNonsterile[,1:8]

#Put in groups:
groupsAmbRed<-cbind(c(AmbientNonsterile),c(ReducedNonsterile))
colnames(groupsAmbRed)<-c("Amb","Red")
ttestAmbRed<-t.test(groupsAmbRed)
mean(AmbientNonsterile)
mean(ReducedNonsterile)

AmbientNonsterileSE<-sd(groupsAmbRed[,1])/sqrt(sum(!is.na(groupsAmbRed[,1])))
ReducedNonsterileSE<-sd(groupsAmbRed[,2])/sqrt(sum(!is.na(groupsAmbRed[,2])))

