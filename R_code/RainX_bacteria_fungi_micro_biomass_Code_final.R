##Code for all annalyses other than community level physiological profiles (see RainX_BIOLOG_Code_final for those analyses)
##Note that figures were generated in R, but tweaked and finished in Illustrator
#Associated files for bacterial analyses: in R_intput_files/Bacteria
#"16s_map_rainx.txt" - treatments and metadata (see READ.me)
#"16S_RainX_80c_simple.txt" - taxonomic classifications of OTUs based off of SILVA v123
#"16S_RainX_OTU_table.txt" - OTU matrix constructed in USEARCH pipline
#"bact.rainx_tree.nwk" - phylogenetic tree of OTUS constructed using PASTA
#Associated files for fungal analyses: in R_intput_files/Fungi
#"ITS_map_rainx.txt" - treatments and metadata (see READ.me)
#"ITS_RainX_80c_simple.txt" - taxonomic classifications of OTUs based off of SILVA v123
#"ITS_RainX_OTU_table.txt" - OTU matrix constructed in USEARCH pipline
#Code assembled by: LUkas Bell-Dereske belldere@msu.edu
#Associated with publication Evans et al. 2019. "Dispersal alters soil microbial community response to drought" in Env Microbiology

library("phyloseq")
library("ggplot2")
library(gridExtra)
library(plyr); library(dplyr)
library(car)
library(MASS)
library(emmeans)
library(reshape2)
library(vegan)
library(ape)
library(usedist)
library(gdata)
library(RVAideMemoire) # for tests after PERMANOVA
#library(cowplot)
options(contrasts=c("contr.sum", "contr.poly"))
"%w/o%" <- function(x,y)!('%in%'(x,y))


#####Begin Bacterial Pre-analyses filtering######

#Load associated files and covert to phyloseq obj
setwd("D:/RainX/github/PAPER_Evans_et_al_RainX-master/R_input_files/Bacteria")
Rainx_16S<-read.table("16S_RainX_OTU_table.txt", header=T,row.names = 1)
bact.OTU = otu_table(Rainx_16S, taxa_are_rows = TRUE)
sum(otu_table(bact.OTU))
#6595777
bact.map=read.table("16s_map_rainx.txt", header=T,row.names = 1)

bact.taxa=as.matrix(read.table("16S_RainX_80c_simple.txt", header=T))
bact.TAX = tax_table(bact.taxa)
bact_tree=read_tree("bact.rainx_tree.nwk")
bact.data=phyloseq(bact.OTU, bact.TAX, sample_data(bact.map2),phy_tree(bact_tree))

ntaxa(bact.data)
#19378

sum(otu_table(bact.data))
#6595777


#remove Chloroplast Reads

bact.data<-subset_taxa(bact.data,Class!="c:Chloroplast")
ntaxa(bact.data)
#19227

sum(otu_table(bact.data))
#6526343

#remove Mitochondria Reads
bact.data<-subset_taxa(bact.data,Family!="f:Mitochondria")
ntaxa(bact.data)
#19076

sum(otu_table(bact.data))
#6471278




#Filter taxa with lees than 3 reads in RainX
bact.rainx<-subset_samples(bact.data, ProjectName=="RainX")

bact.rainx3<-prune_taxa(taxa_sums(bact.rainx) > 2, bact.rainx)
ntaxa(bact.rainx3)
#16974


min(taxa_sums(bact.rainx3))
sum(otu_table(bact.rainx3))
#6453022






#soil: no blanks removed but more than 3 reads across the entire exp
bact.soil<-subset_samples(bact.rainx3, SampleType=="Soil")
bact.soil<-prune_taxa(taxa_sums(bact.soil)>0,bact.soil)
ntaxa(bact.soil)
#14893


sum(otu_table(bact.soil))
#5247626
bact.soil.tree=bact.soil

#root the tree with a random node
phy_tree(bact.soil.tree)<-ape::root(phy_tree(bact.soil.tree), "OTU617", resolve.root=TRUE)




#rain with powerwater removed and more than 3 reads across the entire exp
#remove Power Water from Rain
PWblanks<-subset_samples(bact.data, PW=="PW")

PWblanks.tab <- as(otu_table(PWblanks), "matrix") # Taxa are rows
head(PWblanks.tab)
PWblank <- (PWblanks.tab > 0)
head(PWblank)
numsamplesblankPW <- apply(PWblank, 1, sum)
head(numsamplesblankPW)
datafilPW<- prune_taxa(numsamplesblankPW==0, bact.data)
ntaxa(datafilPW)
# 18885
#19076- 18885=191


sum(otu_table(datafilPW))
#5928142
#6471278-5928142=543136


bact.rainxPW<-subset_samples(datafilPW, ProjectName=="RainX")
bact.rainxPW_2=phyloseq(otu_table(bact.rainxPW), tax_table(bact.rainxPW), sample_data(bact.rainxPW),phy_tree(bact_tree))
ntaxa(bact.rainxPW_2)
#29662


sum(otu_table(bact.rainxPW_2))
#5928142





bact.rainxPW3<-prune_taxa(taxa_sums(bact.rainxPW_2) > 2, bact.rainxPW_2)
ntaxa(bact.rainxPW3)
#16852
#16974-16852=122


min(taxa_sums(bact.rainxPW3))
sum(otu_table(bact.rainxPW3))
#5924717
#6453022-5924717=528305



bact.rain<-subset_samples(bact.rainxPW3, SampleType=="Rain")
bact.rain<-prune_taxa(taxa_sums(bact.rain)>0,bact.rain)
ntaxa(bact.rain)
#4844
#4954-4844=110 #Number of OTUs removed with PW filtering


sum(otu_table(bact.rain))
#852909
#1205396-852909=352487 #Number of Reads removed with PW filtering

bact.rain.tree=phyloseq(otu_table(bact.rain), tax_table(bact.rain), sample_data(bact.map2),phy_tree(bact_tree)) 


#####END Bacterial Pre-analyses filtering######

#####Begin Fungal Pre-analyses filtering######

#Load associated files and covert to phyloseq obj
setwd("D:/RainX/github/PAPER_Evans_et_al_RainX-master/R_input_files/Fungi")
fung.tax=as.matrix(read.table("ITS_RainX_80c_simple.txt",header=T))
fung.otu=read.table("ITS_RainX_OTU_table.txt",header=T,row.names = 1)

fung.OTU = otu_table(fung.otu, taxa_are_rows = TRUE)
ntaxa(fung.OTU)
#4241

fung.TAX = tax_table(fung.tax)
#Read into phyloseq
fung.data=phyloseq(fung.OTU,fung.TAX)
fung.map=read.table("ITS_map_rainx.txt",header=T,row.names = 1)
fung.map$com.rain=with(fung.map, interaction(RainCom, RainLevel))
fung.data=phyloseq(fung.OTU,fung.TAX, sample_data(fung.map))
ntaxa(fung.data)
#4241

sum(otu_table(fung.data))
#2924860





fung.rainx<-subset_samples(fung.data, ProjectName=="RainX")
fung.rainx<-prune_taxa(taxa_sums(fung.rainx) > 0, fung.rainx)
ntaxa(fung.rainx)
#4234

sum(otu_table(fung.rainx))
#2903502



#filtered out anything below 3
fung.rainx3<-prune_taxa(taxa_sums(fung.rainx) > 2, fung.rainx)
ntaxa(fung.rainx3)
#3938

min(taxa_sums(fung.rainx3))
sum(otu_table(fung.rainx3))
#2903000






## Look at overlap with soil and rain
#soil: no blanks removed but more than 3 reads across the entire exp
fung.soil<-subset_samples(fung.rainx3, SampleType=="Soil")
fung.soil<-prune_taxa(taxa_sums(fung.soil)>0,fung.soil)
ntaxa(fung.soil)
#1839

sum(otu_table(fung.soil))
#2150541


#save(fung.soil, file = "fung.soil.RData")
#load("fung.soil.RData")

#rain with powerwater removed and more than 3 reads across the entire exp
#remove Power Water from Rain
F.PWblanks<-subset_samples(fung.data, PW=="PW")

F.PWblanks.tab <- as(otu_table(F.PWblanks), "matrix") # Taxa are rows
head(F.PWblanks.tab)
F.PWblank <- (F.PWblanks.tab > 0)
head(F.PWblank)
F.numsamplesblankPW <- apply(F.PWblank, 1, sum)
head(F.numsamplesblankPW)
F.datafilPW<- prune_taxa(F.numsamplesblankPW==0, fung.data)
ntaxa(F.datafilPW)
#4202
#4241-4202=39


sum(otu_table(F.datafilPW))
#2643797
#2924860-2643797=281063


fung.rainxPW<-subset_samples(F.datafilPW, ProjectName=="RainX")

ntaxa(fung.rainxPW)
#4202


sum(otu_table(fung.rainxPW))
#2643797





fung.rainxPW3<-prune_taxa(taxa_sums(fung.rainxPW) > 2, fung.rainxPW)
ntaxa(fung.rainxPW3)
#3908
#3938-3908=30


min(taxa_sums(fung.rainxPW3))
sum(otu_table(fung.rainxPW3))
#2643298
#2903000-2643298=259702



fung.rain<-subset_samples(fung.rainxPW3, SampleType=="Rain")
fung.rain<-prune_taxa(taxa_sums(fung.rain)>0,fung.rain)
ntaxa(fung.rain)
#2451
#2480-2451=29 #Number of OTUs removed with PW filtering


sum(otu_table(fung.rain))
#524165
#752459-524165=228294 #Number of Reads removed with PW filtering



#####END Fungal Pre-analyses filtering######





#####Begin Fig 1 Diversity Graphs and Analyses#####
#Figure 1 . Community properties after exposure to Dispersal and Drought treatments (N=8/treatment) for 6 months. 
#Microbial biomass carbon (a), metabolic diversity (b), (Shannon diversity index based on substrate utilization rates 
#of 32 substrates), Bacterial (c) and Fungal (d) richness (Chao1), bacterial (e) and fungal (f) evenness (Pielou's). 
#See Table 1 for corresponding ANOVA results. Solid line represents mean values of Time 0 soil cores (N=3). 

#Table 1. Summary of ANOVA p-values for microbial properties exposed to 6-months of Dispersal and 
#No Dispersal and Drought and Ambient conditions. 


#Fig 1a Microbial biomass carbon 

setwd("D:/RainX/github/PAPER_Evans_et_al_RainX-master/R_input_files/Funct_Data")

Zdata<-read.csv("Biomass_AssignZerotoNeg.csv", header=T, na.strings=c("","NA"))

Zdata<-Zdata[1:82,]
ZdataOther<-Zdata[65:82,]
Zdatatrt<-Zdata[1:64,]


ZClaytrt<-Zdatatrt[Zdatatrt$Soil=="clay",]

#Takeing out T0s - getting mean for a T0 line on boxplot
t0<-Zdata[Zdata$Trt=="t0",]
t0clay<-t0[t0$Soil=="clay",]

meant0clay<-mean(t0clay$MB.POC.ug_g_dry)



#Clay
ZClaytrt<-drop.levels(ZClaytrt)


levels(ZClaytrt$Disp.Trt)
head(ZClaytrt)

ZClaytrt$com.rain.time=with(ZClaytrt, interaction(Rain.Trt, Disp.Trt))
levels(ZClaytrt$com.rain.time)
ZClaytrt_POC=subset(ZClaytrt, MB.POC.ug_g_dry!="NA")

ZClaytrt_POC=ZClaytrt_POC %>% group_by(Rain.Trt, Disp.Trt)
ZClaytrt_sum=summarise_at(ZClaytrt_POC, vars(MB.TN.ug_g_dry,MB.POC.ug_g_dry),funs(mean,se=sd(.)/sqrt(n()),sd))




bact.clayE.micro_bio_g=ggplot(ZClaytrt_sum, aes(Rain.Trt, MB.POC.ug_g_dry_mean, ymin = MB.POC.ug_g_dry_mean-MB.POC.ug_g_dry_se, ymax = MB.POC.ug_g_dry_mean+MB.POC.ug_g_dry_se))

(bact.clayE_p_MB=bact.clayE.micro_bio_g+geom_errorbar(width=0.25)+geom_line(aes(group = Disp.Trt))+
    geom_point(size=5,aes(shape=Disp.Trt))+scale_shape_manual(values=c(16,1), name=NULL,
                                                              labels=c("Dispersal",
                                                                       "No Dispersal"))+
    scale_y_continuous(name = "Microbial biomass \n(C ug/g soil)")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = meant0clay, size=1.2))

#Analyses in Table 1

#microbiomass
sort(ZClaytrt$MB.POC.ug_g_dry)
micro_bio.clayE=lm(sqrt(MB.POC.ug_g_dry)~Rain.Trt * Disp.Trt, data=ZClaytrt)
qqPlot(studres(micro_bio.clayE))
hist(studres(micro_bio.clayE))
shapiro.test(studres(micro_bio.clayE))
#p-value = 0.8759

#looks okay....
Anova(micro_bio.clayE,type = 3)
#Rain.Trt            52.75  1   5.3622  0.028414 *  
#Disp.Trt           115.15  1  11.7048  0.001999 ** 
#Rain.Trt:Disp.Trt    6.89  1   0.7006  0.409918   


lsmeans(micro_bio.clayE, pairwise~Rain.Trt * Disp.Trt)
#contrast                                estimate       SE df t.ratio p.value
#ambient,nonsterile - reduced,nonsterile 1.668683 1.568260 27   1.064  0.7138
#ambient,nonsterile - ambient,sterile    2.916406 1.568260 27   1.860  0.2688
#ambient,nonsterile - reduced,sterile    6.474384 1.623304 27   3.988  0.0024
#reduced,nonsterile - ambient,sterile    1.247723 1.568260 27   0.796  0.8558
#reduced,nonsterile - reduced,sterile    4.805701 1.623304 27   2.960  0.0303
#ambient,sterile  - reduced,sterile      3.557978 1.623304 27   2.192  0.1511



#Fig 1b metabolic diversity
#Load in file
setwd("D:/RainX/github/PAPER_Evans_et_al_RainX-master/R_input_files/Funct_Data")
Data<-read.csv("RainX ECO metafile_se_forR.csv",  header=TRUE)
#Current name: "RainX ECO metafile_se_forR.csv"

head(Data)
dim(Data)



#Need to first rearrange data
#"Swing" CarbonSource into columns along the top, so we have a  data table amenable to analysis with multivariate tools
#(analogous to community data, OTUs go along the top, samples along the side or in a separate mapping file)
DataMultivar<-dcast(Data, Sample+Label.ID + Rain.Trt+Disp.Trt+Soil ~ CarbonSource, value.var="Final")
head(DataMultivar)

pMap<-DataMultivar[,1:5]
DataOnly<-DataMultivar[,6:37]
DataOnly<-cbind(pMap[,1],DataOnly)


Diversity<-diversity(DataOnly,index="shannon",MARGIN=1)


pMap$Soil
Div<-cbind(pMap,Diversity)
DiversityClay<-Div[Div$Soil=="clay",]
DiversityClay<-DiversityClay[DiversityClay$Disp.Trt!="t0",]
DiversityClay<-drop.levels(DiversityClay)



T0<-Div[Div$Disp.Trt=="t0",]
T0Clay<-T0[T0$Soil=="clay",]


DiversityClay$com.rain.time=with(DiversityClay, interaction(Rain.Trt, Disp.Trt))
levels(DiversityClay$com.rain.time)

DiversityClay=DiversityClay %>% group_by(Rain.Trt, Disp.Trt)
DiversityClay_sum=summarise_at(DiversityClay, vars(Diversity),funs(mean,se=sd(.)/sqrt(n()),sd))

scaleFUN2 <- function(x) sprintf("%.2f", x)


bact.clayE.bio_log_g=ggplot(DiversityClay_sum, aes(Rain.Trt, mean, ymin = mean-se, ymax = mean+se))

(bact.clayE_p_bio=bact.clayE.bio_log_g+geom_errorbar(width=0.25)+geom_line(aes(group = Disp.Trt))+
    geom_point(size=5,aes(shape=Disp.Trt))+scale_shape_manual(values=c(16,1), name=NULL,
                                                              labels=c("Dispersal",
                                                                       "No Dispersal"))+
    scale_y_continuous(name = "Metabolic diversity")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(T0Clay$Diversity), size=1.2))



#Analyses in Table 1
#Biolog diversity
#sort(ZClaytrt$MB.POC.ug_g_dry)
biolog.clayE=lm((Diversity)~Rain.Trt * Disp.Trt, data=DiversityClay)
qqPlot(studres(biolog.clayE))
hist(studres(biolog.clayE))
shapiro.test(studres(biolog.clayE))
#p-value = 0.6834
boxCox(biolog.clayE)
#looks okay....
Anova(biolog.clayE,type = 3)
#Rain.Trt           5.317  1  34.637 2.493e-06 ***
#Disp.Trt           5.061  1  32.969 3.673e-06 ***
#Rain.Trt:Disp.Trt  1.869  1  12.175  0.001621 ** 


lsmeans(biolog.clayE, pairwise~Rain.Trt * Disp.Trt)
#contrast                                   estimate        SE df t.ratio p.value
#ambient,Nonsterile - reduced,Nonsterile  1.29860006 0.1959035 28   6.629  <.0001
#ambient,Nonsterile - ambient,Sterile     1.27873253 0.1959035 28   6.527  <.0001
#ambient,Nonsterile - reduced,Sterile     1.61065369 0.1959035 28   8.222  <.0001
#reduced,Nonsterile - ambient,Sterile    -0.01986753 0.1959035 28  -0.101  0.9996
#reduced,Nonsterile - reduced,Sterile     0.31205363 0.1959035 28   1.593  0.3988
#ambient,Sterile  - reduced,Sterile       0.33192116 0.1959035 28   1.694  0.3455






#####Bacterial Diversities#####

#subset to only clay and end samples
bact.soil.end=subset_samples(bact.soil.tree, SampleTime=="End")
bact.clay.end=subset_samples(bact.soil.end, SoilType=="Clay")

alpha_meas = c("Observed", "Chao1", "Shannon")




scaleFUN3 <- function(x) sprintf("%.3f", x)

bact.clay.end.map=sample_data(bact.clay.end)
bact.clayE.t.divfil=estimate_richness(bact.clay.end,measures=alpha_meas)

bact.clayE.t.divfil=merge(bact.clayE.t.divfil, bact.clay.end.map, by ="row.names")
bact.clayE.t.divfil=mutate(bact.clayE.t.divfil, pielou=Shannon*(1/log(Observed)))


#subset so intial community is included in diversity metric
bact.clay.tree=subset_samples(bact.soil.tree, SoilType=="Clay")
bact.clayES.tree=subset_samples(bact.clay.tree, SampleTime!="Stl")
bact.clayES.tree<-prune_taxa(taxa_sums(bact.clayES.tree)>0,bact.clayES.tree)

#starting diversity included
#Diversities
alpha_meas = c("Observed", "Chao1", "Shannon")



bact.clayES.t.map=sample_data(bact.clayES.tree)
bact.clayES.t.divfil=estimate_richness(bact.clayES.tree,measures=alpha_meas)

bact.clayES.t.divfil=merge(bact.clayES.t.divfil, bact.clayES.t.map, by ="row.names")
bact.clayES.t.divfil=mutate(bact.clayES.t.divfil, pielou=Shannon*(1/log(Observed)))
bClay_t0<-bact.clayES.t.divfil[bact.clayES.t.divfil$SampleTime=="Start",]

b_meant0clay_choa1<-mean(bClay_t0$Chao1)
b_meant0clay_piel<-mean(bClay_t0$pielou)

#####Fungal Diversities#####

#Fungi

#subset to only clay and end samples
fung.soil.end=subset_samples(fung.soil, SampleTime=="End")
fung.clay.end=subset_samples(fung.soil.end, SoilType=="Clay")

alpha_meas = c("Observed", "Chao1", "Shannon")
fung.clay.end.map=sample_data(fung.clay.end)
fung.clayE.divfil=estimate_richness(fung.clay.end,measures=alpha_meas)

fung.clayE.divfil=merge(fung.clayE.divfil, fung.clay.end.map, by ="row.names")
fung.clayE.divfil=mutate(fung.clayE.divfil, pielou=Shannon*(1/log(Observed)))


#starting diversity included
fung.clay=subset_samples(fung.soil, SoilType=="Clay")
fung.clayES=subset_samples(fung.clay, SampleTime!="STL")
fung.clayES<-prune_taxa(taxa_sums(fung.clayES)>0,fung.clayES)
fung.clayES.map=sample_data(fung.clayES)
#Diversities
alpha_meas = c("Observed", "Chao1", "Shannon")



fung.clayES.divfil=estimate_richness(fung.clayES,measures=alpha_meas)

fung.clayES.divfil=merge(fung.clayES.divfil, fung.clayES.map, by ="row.names")
fung.clayES.divfil=mutate(fung.clayES.divfil, pielou=Shannon*(1/log(Observed)))

fClay_t0<-fung.clayES.divfil[fung.clayES.divfil$SampleTime=="Start",]
f_meant0fung_choa1<-mean(fClay_t0$Chao1)
f_meant0fung_piel<-mean(fClay_t0$pielou)


scaleFUN1 <- function(x) sprintf("%.1f", x)

scaleFUN3 <- function(x) sprintf("%.3f", x)



#Fig 1c Bacterial Richness Chao1

bact.clayE.t.divfil=bact.clayE.t.divfil %>% group_by(RainCom, RainLevel)
bact.clayE.t.divfil_sum=summarise_at(bact.clayE.t.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


bact.clayE.t.chao1_g=ggplot(bact.clayE.t.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(bact.clayE_p1=bact.clayE.t.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Bacterial richness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = b_meant0clay_choa1, size=1.2))

#Table 1 Analyses

#Chao1
bact.clayE.t.chao1=lm((Chao1)~RainCom*RainLevel, data=bact.clayE.t.divfil)
qqPlot(studres(bact.clayE.t.chao1))
hist(studres(bact.clayE.t.chao1))
shapiro.test(studres(bact.clayE.t.chao1))
#p-value = 0.1187
boxCox(bact.clayE.t.chao1)
#looks okay....
Anova(bact.clayE.t.chao1,type = 3)
#RainCom             2298348  1   29.6014 8.325e-06 ***
#RainLevel           1119429  1   14.4176 0.0007221 ***

lsmeans(bact.clayE.t.chao1, pairwise~RainLevel*RainCom)

#$contrasts
#contrast                                 estimate       SE df t.ratio p.value
#Ambient,NonSterile - Reduced,NonSterile  371.6430 139.3226 28   2.667  0.0573
#Ambient,NonSterile - Ambient,Sterile    -538.4250 139.3226 28  -3.865  0.0032
#Ambient,NonSterile - Reduced,Sterile    -161.9274 139.3226 28  -1.162  0.6550
#Reduced,NonSterile - Ambient,Sterile    -910.0680 139.3226 28  -6.532  <.0001
#Reduced,NonSterile - Reduced,Sterile    -533.5704 139.3226 28  -3.830  0.0035
#Ambient,Sterile - Reduced,Sterile        376.4976 139.3226 28   2.702  0.0532



#Fig 1d Fungal Richness Chao1

fung.clayE.divfil=fung.clayE.divfil %>% group_by(RainCom, RainLevel)
fung.clayE.divfil_sum=summarise_at(fung.clayE.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


fung.clayE.chao1_g=ggplot(fung.clayE.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(fung.clayE_p1=fung.clayE.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Fungal richness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y =element_text(size=18), axis.text.x=element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = f_meant0fung_choa1, size=1.2))


#Table 1 Analyses

#Chao1
fung.clayE.chao1=lm((Chao1)~RainCom*RainLevel, data=fung.clayE.divfil)
qqPlot(studres(fung.clayE.chao1))
hist(studres(fung.clayE.chao1))
shapiro.test(studres(fung.clayE.chao1))
#p-value = 0.1745
boxCox(fung.clayE.chao1)
#looks okay....
Anova(fung.clayE.chao1,type = 3)
#RainCom             11221  1   11.1246  0.002412 ** 

lsmeans(fung.clayE.chao1, pairwise~RainLevel*RainCom)

#$contrasts
#contrast                                 estimate       SE df t.ratio p.value
#Ambient,NonSterile - Reduced,NonSterile -13.09266 15.88003 28  -0.824  0.8424
#Ambient,NonSterile - Ambient,Sterile    -42.32857 15.88003 28  -2.666  0.0576
#Ambient,NonSterile - Reduced,Sterile    -45.66864 15.88003 28  -2.876  0.0361
#Reduced,NonSterile - Ambient,Sterile    -29.23591 15.88003 28  -1.841  0.2761
#Reduced,NonSterile - Reduced,Sterile    -32.57598 15.88003 28  -2.051  0.1939
#Ambient,Sterile - Reduced,Sterile        -3.34007 15.88003 28  -0.210  0.9966



#Fig 1e Bacterial evenness (Pielou's)

bact.clayE.t.piel_g=ggplot(bact.clayE.t.divfil_sum, aes(RainLevel, pielou_mean, ymin = pielou_mean-pielou_se, ymax = pielou_mean+pielou_se))

(bact.clayE_p2=bact.clayE.t.piel_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    scale_y_continuous(name = "Bacterial evenness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = b_meant0clay_piel, size=1.2))






#Table 1 Analyses

#pielou
bact.clayE.t.piel=lm((pielou)~RainCom*RainLevel, data=bact.clayE.t.divfil)
qqPlot(studres(bact.clayE.t.piel))
hist(studres(bact.clayE.t.piel))
shapiro.test(studres(bact.clayE.t.piel))
#0.521
boxCox(bact.clayE.t.piel)
#looks okay....
Anova(bact.clayE.t.piel,type = 3)
#RainCom            0.0001  1 3.7598e+00 0.0626378 .  
#RainLevel          0.0005  1 1.7625e+01 0.0002464 ***
#RainCom:RainLevel  0.0001  1 2.9577e+00 0.0965064 .  

lsmeans(bact.clayE.t.piel, pairwise~RainLevel*RainCom)
#$contrasts
#contrast                                    estimate         SE df t.ratio p.value
#Ambient,NonSterile - Reduced,NonSterile  0.004803702 0.00274103 28   1.753  0.3168
#Ambient,NonSterile - Ambient,Sterile    -0.007091511 0.00274103 28  -2.587  0.0681
#Ambient,NonSterile - Reduced,Sterile     0.004378812 0.00274103 28   1.598  0.3963
#Reduced,NonSterile - Ambient,Sterile    -0.011895212 0.00274103 28  -4.340  0.0009
#Reduced,NonSterile - Reduced,Sterile    -0.000424890 0.00274103 28  -0.155  0.9986
#Ambient,Sterile - Reduced,Sterile        0.011470322 0.00274103 28   4.185  0.0014



#Fig 1f Fungal evenness (Pielou's)

fung.clayE.piel_g=ggplot(fung.clayE.divfil_sum, aes(RainLevel, pielou_mean, ymin = pielou_mean-pielou_se, ymax = pielou_mean+pielou_se))

(fung.clayE_p2=fung.clayE.piel_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    scale_y_continuous(name = "Fungal evenness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text=element_text(size=18),axis.title.x=element_blank(), 
                     axis.title.y=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = f_meant0fung_piel, size=1.2))



#Table 1 Analyses


#pielou
fung.clayE.piel=lm((pielou)^2~RainCom*RainLevel, data=fung.clayE.divfil)
qqPlot(studres(fung.clayE.piel))
hist(studres(fung.clayE.piel))
shapiro.test(studres(fung.clayE.piel))
#p-value = 0.1147
boxCox(fung.clayE.piel)
#looks okay....
Anova(fung.clayE.piel,type = 3)
#nada sig

lsmeans(fung.clayE.piel, pairwise~RainLevel*RainCom)
#$contrasts
#contrast                                    estimate         SE df t.ratio p.value
#Ambient,NonSterile - Reduced,NonSterile -0.069333386 0.03717026 28  -1.865  0.2656
#Ambient,NonSterile - Ambient,Sterile    -0.058940853 0.03717026 28  -1.586  0.4027
#Ambient,NonSterile - Reduced,Sterile    -0.064148781 0.03717026 28  -1.726  0.3298
#Reduced,NonSterile - Ambient,Sterile     0.010392533 0.03717026 28   0.280  0.9922
#Reduced,NonSterile - Reduced,Sterile     0.005184606 0.03717026 28   0.139  0.9990
#Ambient,Sterile - Reduced,Sterile       -0.005207928 0.03717026 28  -0.140  0.9990



#Combined the graphs into a multipanel


plot_grid(bact.clayE_p_MB, bact.clayE_p_bio,bact.clayE_p1,fung.clayE_p1,bact.clayE_p2,fung.clayE_p2, 
          ncol = 2,align = "v", rel_heights = c(1,1,1.1))


#####END Fig 1 Diversity Graphs and Analyses#####

#####Begin Fig 2 Stacked Bar graphs####
#Taxon distribution

#Bacteria
#subset to only clay and end samples
bact.soil.end=subset_samples(bact.soil.tree, SampleTime=="End")
bact.clay.end=subset_samples(bact.soil.end, SoilType=="Clay")
bact.clay.end<-prune_taxa(taxa_sums(bact.clay.end)>0,bact.clay.end)

#extract OTU tables
bact.clay.end_NO_TREE=otu_table(bact.clay.end)
bact.rain_NO_TREE=otu_table(bact.rain.tree)
#combine files to make one phyloseq file with both rain and soil samples
bact.clay_rain.end=merge_phyloseq(bact.clay.end_NO_TREE,bact.rain_NO_TREE)
bact.clay_rain.end<-prune_taxa(taxa_sums(bact.clay_rain.end)>0,bact.clay_rain.end)
ntaxa(bact.clay_rain.end)
#14369

#make a new grouping variable with soil status and rain vs soil
bact.data_map=sample_data(bact.data)
bact.data_map$com.rain=with(bact.data_map, interaction(RainCom, RainLevel))

#Make a new Phyloseq obj with all of the associated data
bact.clay_rain.tree.end=phyloseq(otu_table(bact.clay_rain.end),tax_table(bact.data),bact.data_map,phy_tree(bact.data))

#merge OTUs by the soil and precipitation treatment
bact.clay_rain.tree.end_fact=merge_samples(bact.clay_rain.tree.end, "com.rain")
sample_names(bact.clay_rain.tree.end_fact)     

#combine the reads at Phylum level
get_taxa_unique(bact.clay_rain.tree.end, taxonomic.rank="Phylum")
#45
(bact.clay_rain.tree.end_fact.phylum<-tax_glom(bact.clay_rain.tree.end_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
TopPHYL = names(sort(taxa_sums(bact.clay_rain.tree.end_fact.phylum), TRUE)[1:10])
bact.clayT10 = prune_taxa(TopPHYL, bact.clay_rain.tree.end_fact.phylum)

Phyl_name_T10=get_taxa_unique(bact.clayT10, taxonomic.rank="Phylum")
PHyl_Name_T10 <- colsplit(Phyl_name_T10, ":", c("letter", "Phyl_name"))


#transform the read counts to prop of total reads

bact.clay_rain.tree.end_fact.phylum.prop=transform_sample_counts(bact.clay_rain.tree.end_fact.phylum, function(x)x/sum(x))

taxon_positions=c("Rain.Rain","NonSterile.Ambient","NonSterile.Reduced","Sterile.Ambient","Sterile.Reduced")
bact.clayT10.prop = prune_taxa(TopPHYL, bact.clay_rain.tree.end_fact.phylum.prop)

bact.clayT10.prop_otu=as.data.frame(t(otu_table(bact.clayT10.prop)))

#create an other taxa category
taxon_sums=c(sum(bact.clayT10.prop_otu$NonSterile.Ambient),sum(bact.clayT10.prop_otu$Sterile.Ambient),
             sum(bact.clayT10.prop_otu$Rain.Rain),sum(bact.clayT10.prop_otu$NonSterile.Reduced),
             sum(bact.clayT10.prop_otu$Sterile.Reduced))
other_spp=c(as.numeric(1-taxon_sums))
bact.clayT10.prop_OTU=rbind(bact.clayT10.prop_otu,other_spp)
summary(bact.clayT10.prop_OTU)
bact.clayT10.prop_OTU[,"Phylum"]=c(PHyl_Name_T10$Phyl_name,"Other Phyla")
summary(bact.clayT10.prop_OTU)
bact.clayT10.prop_OTU_M=melt(bact.clayT10.prop_OTU,id="Phylum")
phyl_order=c(sort(PHyl_Name_T10$Phyl_name),"Other Phyla")
summary(bact.clayT10.prop_OTU_M)



(p_bact_T10_v2.1_color=ggplot(bact.clayT10.prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_blank(),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Rain","Dispersal\nAmbient", 
                                                       "Dispersal\nDrought",
                                                       "No Dispersal\nAmbient", 
                                                       "No Dispersal\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))



#######fungi####



#merge the clay and rain phyloseq obj 
fung.soil.end=subset_samples(fung.soil, SampleTime=="End")
fung.clay.end=subset_samples(fung.soil.end, SoilType=="Clay")
fung.clay_rain.end=merge_phyloseq(fung.clay.end,fung.rain)


#Group by soil status and rain versus soil

fung.clay_rain.end_fact=merge_samples(fung.clay_rain.end, "com.rain")
sample_names(fung.clay_rain.end_fact)     

get_taxa_unique(fung.clay_rain.end, taxonomic.rank="Phylum")
#11
(fung.clay_rain.end_fact.phylum<-tax_glom(fung.clay_rain.end_fact, taxrank="Phylum"))
sort(taxa_sums(fung.clay_rain.end_fact.phylum))


#Sum read numbers by Phylum
get_taxa_unique(fung.clay_rain.end_fact.phylum, taxonomic.rank="Phylum")

#Take the top 10 phyla
TopPHYL = names(sort(taxa_sums(fung.clay_rain.end_fact.phylum), TRUE)[1:10])
fung.clayT10 = prune_taxa(TopPHYL, fung.clay_rain.end_fact.phylum)





#convert read numbers to prop of total reads
fung.clay_rain.end_fact.phylum.prop=transform_sample_counts(fung.clayT10, function(x)x/sum(x))
sort(get_taxa_unique(fung.clay_rain.end_fact.phylum.prop, taxonomic.rank="Phylum"))
Fphyl_order=c("p:Ascomycota","p:Basidiomycota","p:Chytridiomycota",
              "p:Entomophthoromycota","p:Glomeromycota","p:Kickxellomycota","p:Mortierellomycota",
              "p:Mucoromycota","p:Rozellomycota","UNKNOWN")
Fphyl_order_names=c("Ascomycota","Basidiomycota","Chytridiomycota",
                    "Entomophthoromycota","Glomeromycota","Kickxellomycota","Mortierellomycota",
                    "Mucoromycota","Rozellomycota","UNKNOWN")

taxon_positions=c("Rain.Rain","NonSterile.Ambient","NonSterile.Reduced","Sterile.Ambient","Sterile.Reduced")
fung.clayT10.prop = prune_taxa(TopPHYL, fung.clay_rain.end_fact.phylum.prop)


(p_fung_T10_v2.1_color=plot_bar(fung.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum, levels=Fphyl_order)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Paired",breaks=Fphyl_order,labels=Fphyl_order_names)+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),axis.title=element_text(size=20),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20))+
    guides(fill=guide_legend(title="Phyla"))+
    scale_x_discrete(limits = taxon_positions,labels=c("Rain","Dispersal\nAmbient", 
                                                       "Dispersal\nDrought",
                                                       "No Dispersal\nAmbient", 
                                                       "No Dispersal\nDrought")))



#Combine graphs
grid.arrange(p_bact_T10_v2.1_color,p_fung_T10_v2.1_color, nrow=1,ncol=2,widths=c(1,1))




#####End Fig 2 Stacked Bar graphs####

#####Begin Fig 3 Non-metric multidimensional scaling plots of bacterial####

#Figure 4.Non-metric multidimensional scaling plots of bacterial (a,c) and fungal (b,d) communities. 
#Plots a,b, show communities in rain samples (+, includes air deposition) collected over the 6-month 
#experiment, and c,d show soil communities. Soil mesocosms were experimentally altered by removing 
#immigrants (Dispersal and No Dispersal), and/or reducing the water input (Ambient and Drought). 
#Time 0 represents soils before treatments occurred. PerMANOVA statistics are shown in Table 1

#Bacterial Dataset
#including the intial community
#clay end and T0 only
bact.clay.tree=subset_samples(bact.soil.tree, SoilType=="Clay")
bact.clayES.tree=subset_samples(bact.clay.tree, SampleTime!="Stl")
bact.clayES.tree<-prune_taxa(taxa_sums(bact.clayES.tree)>0,bact.clayES.tree)

bact.clayES_NO_TREE=otu_table(bact.clayES.tree)
bact.rain_NO_TREE=otu_table(bact.rain.tree)
bact.clay_rain.ES=merge_phyloseq(bact.clayES_NO_TREE,bact.rain_NO_TREE)

ntaxa(bact.clay_rain.ES)
#15026
phy_tree(bact.data)<-ape::root(phy_tree(bact.data), "OTU617", resolve.root=TRUE)
bact.clay_rain.ES=phyloseq(otu_table(bact.clay_rain.ES),tax_table(bact.data),sample_data(bact.data),phy_tree(bact.data))

#ordinate
bact.clay_rain.ES.ord <- ordinate(bact.clay_rain.ES, method="NMDS",distance = "wunifrac")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1125909




(fig1_P1_t0=plot_ordination(bact.clay_rain.ES, bact.clay_rain.ES.ord, shape = "SampleType")+
    geom_point(aes(shape = SampleType), size=4, fill="gray",color="black")+scale_shape_manual(values=c(3,21), name=NULL)+
    geom_text(aes(label=SampleTimeTrunc),hjust=.5, vjust=1.5)+scale_y_continuous(labels=scaleFUN2)+
    theme_bw()+theme(axis.text=element_text(size=10), legend.position ="none", 
                     axis.title.y=element_text(size=12),axis.title.x=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    ggtitle(label = "Bacteria")+theme(plot.title = element_text(hjust = 0.5, size=24)))


#ordinate
bact.clayES.tree.ord <- ordinate(bact.clayES.tree, method="NMDS",distance = "wunifrac")
#*** Solution reached
#0.07064679
positions_t0 <- c("NonSterile.Ambient.Start", "NonSterile.Ambient.End", "NonSterile.Reduced.End", 
                  "Sterile.Ambient.End", "Sterile.Reduced.End")
(fig1_P2_t0=plot_ordination(bact.clayES.tree, bact.clayES.tree.ord)+
    geom_point(shape = 21,aes(fill=com.rain.time), size=4, stroke = 1)+
    scale_fill_manual(values=c("grey","grey", "white", "white","black"), name="Treatment",labels=c("Time 0","Dispersal Ambient", 
                                                                                                   "Dispersal Drought",
                                                                                                   "No Dispersal Ambient", 
                                                                                                   "No Dispersal Drought"),
                      breaks=positions_t0)+
    geom_point(shape = 20,aes(color=com.rain.time), size=2, stroke = 1)+
    scale_color_manual(values=c("black","grey","black","white","black"), name="Treatment",labels=c("Time 0","Dispersal Ambient", 
                                                                                                   "Dispersal Drought",
                                                                                                   "No Dispersal Ambient", 
                                                                                                   "No Dispersal Drought"),
                       breaks=positions_t0)+
    theme_bw()+theme(axis.text=element_text(size=10),legend.position ="none", 
                     axis.title=element_text(size=12),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))


#Results for table 1

#Stats examining the effect of drought and dispersal on community composition

bact.clay.end.UF.dis=distance(bact.clay.end,method="wunifrac", type = "samples")
bact.clay.end.map=sample_data(bact.clay.end)

adonis(bact.clay.end.UF.dis~bact.clay.end.map$RainLevel*bact.clay.end.map$RainCom, perm=9999)
#bact.clay.end.map$RainLevel                            1  0.020800 0.0207998 22.9571 0.35999 0.0001 ***
#bact.clay.end.map$RainCom                              1  0.007354 0.0073538  8.1165 0.12728 0.0002 ***
#bact.clay.end.map$RainLevel:bact.clay.end.map$RainCom  1  0.004256 0.0042562  4.6976 0.07366 0.0024 ** 



#Fungal Dataset
#clay end and T0 only
fung.clay=subset_samples(fung.soil, SoilType=="Clay")
fung.clayES=subset_samples(fung.clay, SampleTime!="STL")
fung.clayES<-prune_taxa(taxa_sums(fung.clayES)>0,fung.clayES)
fung.clay_rain.ES=merge_phyloseq(fung.clayES, fung.rain)


#ordinate
fung.clay_rain.ES.ord <- ordinate(fung.clay_rain.ES, method="NMDS",distance = "bray")
#*** Solution reached
#0.08376851

scaleFUN1 <- function(x) sprintf("%.1f", x)
(fig1_P3_t0=plot_ordination(fung.clay_rain.ES, fung.clay_rain.ES.ord, shape = "SampleType")+
    geom_point(aes(shape = SampleType), size=4, fill="gray",color="black")+scale_shape_manual(values=c(21,3), name=NULL,
                                                                                              labels=c("Soil                             ",
                                                                                                       "Rain                              "))+
    geom_text(aes(label=SampleTimeTrunc),hjust=.5, vjust=1.5)+scale_y_continuous(labels=scaleFUN1, limits = c(-2.5,3))+
    theme_bw()+theme(axis.text=element_text(size=10), legend.text = element_text(size =14),
                     axis.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    ggtitle(label = "Fungi")+theme(plot.title = element_text(hjust = 0.5, size=24)))


#ordinate
fung.clayES.ord <- ordinate(fung.clayES, method="NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.2530532
#save(fung.clayES.ord, file = "fung.clayES.ord_PW.RData")
#load("fung.clayES.ord_PW.RData")
positions_t0 <- c("NonSterile.Ambient.Start", "NonSterile.Ambient.End", "NonSterile.Reduced.End", 
                  "Sterile.Ambient.End", "Sterile.Reduced.End")
(fig1_P4_t0=plot_ordination(fung.clayES, fung.clayES.ord)+
    geom_point(shape = 21,aes(fill=com.rain.time), size=4, stroke = 1)+
    scale_fill_manual(values=c("grey","grey", "white", "white","black"), name=NULL,labels=c("Time 0","Dispersal Ambient", 
                                                                                            "Dispersal Drought",
                                                                                            "No Dispersal Ambient", 
                                                                                            "No Dispersal Drought"),
                      breaks=positions_t0)+
    geom_point(shape = 20,aes(color=com.rain.time), size=2, stroke = 1)+
    scale_color_manual(values=c("black","grey","black","white","black"), name=NULL,labels=c("Time 0","Dispersal Ambient", 
                                                                                    "Dispersal Drought",
                                                                                    "No Dispersal Ambient", 
                                                                                    "No Dispersal Drought"),
                       breaks=positions_t0)+
    theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                     axis.title.x=element_text(size=12),axis.title.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

#Results for table 1

#Stats examining the effect of drought and dispersal on community composition
fung.clay.end.dis=distance(fung.clay.end,method="bray", type = "samples")
fung.clay.end.map=sample_data(fung.clay.end)

adonis(fung.clay.end.dis~fung.clay.end.map$RainLevel*fung.clay.end.map$RainCom, perm=9999)
#fung.clay.end.map$RainLevel                            1    0.2921 0.29208 1.61311 0.05082 0.0093 **
#fung.clay.end.map$RainCom                              1    0.2332 0.23320 1.28795 0.04058 0.0838 . 


#combine the graphs into a panel
grid.arrange(fig1_P1_t0,fig1_P3_t0,fig1_P2_t0, fig1_P4_t0, nrow=2,ncol=2,heights=c(1.01,1), widths=c(1,1.5))


#####End Fig 3 Non-metric multidimensional scaling plots of bacterial####

#####Begin Mean Distance between Communty Centriods####

######Table 2 Analyses #####

#Bacteria
#Mean Pairwise distance

#Clay only
bact.clay.tree=subset_samples(bact.soil.tree, SoilType=="Clay")
bact.clay.tree<-prune_taxa(taxa_sums(bact.clay.tree) > 0, bact.clay.tree)
ntaxa(bact.clay.tree)
#12863
sum(otu_table(bact.clay.tree))
#2575563

#Calculate Weighted unifrac distance
bact.clay.UF=distance(bact.clay.tree,method="wunifrac", type = "samples")
#root as -- OTU617
bact.clay.map=sample_data(bact.clay.tree)


#pairwise distance using the usedist package

summary(bact.clay.map)

#Calculated the pairwise the distance between and within each treatment
bact.clay.UF_pair_dist=dist_groups(bact.clay.UF,bact.clay.map$com.rain.time)
head(bact.clay.UF_pair_dist)
nrow(bact.clay.UF_pair_dist)
#820
summary(bact.clay.UF_pair_dist)

#The below steps are for creating an easier way to parse the pairwise distance beetween treatments
bact.clay.UF_pair_dist_2=bact.clay.UF_pair_dist
colnames(bact.clay.UF_pair_dist_2)[c(1,2,3,4)]<-c("Item2","Item1","Group2","Group1")
nrow(bact.clay.UF_pair_dist_2)

bact.clay.UF_pair_dist_comb=rbind(bact.clay.UF_pair_dist,bact.clay.UF_pair_dist_2)
nrow(bact.clay.UF_pair_dist_comb)
summary(bact.clay.UF_pair_dist_comb)

#Calculated the mean distances to remove the pseudo-replication
bact.clay.UF_pair_dist_comb=bact.clay.UF_pair_dist_comb %>% group_by(Item1,Group1,Label)

bact.clay.UF_pair_dist_comb_sum=summarise_at(bact.clay.UF_pair_dist_comb,"Distance", funs(mean,n()))
nrow(bact.clay.UF_pair_dist_comb_sum)
#287
head(bact.clay.UF_pair_dist_comb_sum)
summary(bact.clay.UF_pair_dist_comb_sum)





#Distance between Dispersal Ambient and Dispersal Drought
NonSterile.Ambient.End_NonSterile.Reduced.End=subset(bact.clay.UF_pair_dist_comb_sum, 
                                                     Label=="Between NonSterile.Ambient.End and NonSterile.Reduced.End"&
                                                       Group1=="NonSterile.Ambient.End")
nrow(NonSterile.Ambient.End_NonSterile.Reduced.End)
#8
mean(NonSterile.Ambient.End_NonSterile.Reduced.End$mean)
#0.1175801
se(NonSterile.Ambient.End_NonSterile.Reduced.End$mean)
#0.002089044


#Distance between No Dispersal Ambient and No Dispersal Drought
Sterile.Ambient.End_Sterile.Reduced.End=subset(bact.clay.UF_pair_dist_comb_sum, 
                                               Label=="Between Sterile.Ambient.End and Sterile.Reduced.End"&
                                                 Group1=="Sterile.Ambient.End")
nrow(Sterile.Ambient.End_Sterile.Reduced.End)
#8
mean(Sterile.Ambient.End_Sterile.Reduced.End$mean)
#0.08555561
se(Sterile.Ambient.End_Sterile.Reduced.End$mean)
#0.002319076

Reduced_pair_dist=droplevels(rbind(NonSterile.Ambient.End_NonSterile.Reduced.End,Sterile.Ambient.End_Sterile.Reduced.End))


qqPlot(Reduced_pair_dist$mean)
hist(Reduced_pair_dist$mean)
shapiro.test((Reduced_pair_dist$mean))
#0.04857
CommC_m1=aov(mean~Group1, data=Reduced_pair_dist)
anova(CommC_m1)
#Group1     1 0.0041023 0.0041023  105.27 6.79e-08 ***


#Fungi
#Mean Pairwise distance

#Clay only
fung.clay=subset_samples(fung.soil, SoilType=="Clay")
fung.clay<-prune_taxa(taxa_sums(fung.clay) > 0, fung.clay)
ntaxa(fung.clay)
#1374
sum(otu_table(fung.clay))
#1084332

fung.clay.bray=distance(fung.clay,method="bray", type = "samples")
#
fung.clay.map=sample_data(fung.clay)

#pairwise distance using the usedist package

summary(fung.clay.map)

#Calculate the pairwise distance between and within each treatment
fung.clay_pair_dist=dist_groups(fung.clay.bray,fung.clay.map$com.rain.time)
head(fung.clay_pair_dist)
nrow(fung.clay_pair_dist)
#820
summary(fung.clay_pair_dist)

#The below section of code is to make it easier to parse the pairwise distances
fung.clay_pair_dist_2=fung.clay_pair_dist
colnames(fung.clay_pair_dist_2)[c(1,2,3,4)]<-c("Item2","Item1","Group2","Group1")
nrow(fung.clay_pair_dist_2)

fung.clay_pair_dist_comb=rbind(fung.clay_pair_dist,fung.clay_pair_dist_2)
nrow(fung.clay_pair_dist_comb)
summary(fung.clay_pair_dist_comb)

#Calculate the mean distance to remove the psuedo-replication 
fung.clay_pair_dist_comb=fung.clay_pair_dist_comb %>% group_by(Item1,Group1,Label)

fung.clay_pair_dist_comb_sum=summarise_at(fung.clay_pair_dist_comb,"Distance", funs(mean,n()))
nrow(fung.clay_pair_dist_comb_sum)
#287
head(fung.clay_pair_dist_comb_sum)
summary(fung.clay_pair_dist_comb_sum)

#Distance between Dispersal Ambient and Dispersal Drought
NonSterile.Ambient.End_NonSterile.Reduced.End=subset(fung.clay_pair_dist_comb_sum, 
                                                     Label=="Between NonSterile.Ambient.End and NonSterile.Reduced.End"&
                                                       Group1=="NonSterile.Ambient.End")
nrow(NonSterile.Ambient.End_NonSterile.Reduced.End)
#8
mean(NonSterile.Ambient.End_NonSterile.Reduced.End$mean)
#0.6296333
se(NonSterile.Ambient.End_NonSterile.Reduced.End$mean)
#0.01927116


#Distance between No Dispersal Ambient and No Dispersal Drought
Sterile.Ambient.End_Sterile.Reduced.End=subset(fung.clay_pair_dist_comb_sum, 
                                               Label=="Between Sterile.Ambient.End and Sterile.Reduced.End"&
                                                 Group1=="Sterile.Ambient.End")
nrow(Sterile.Ambient.End_Sterile.Reduced.End)
#8
mean(Sterile.Ambient.End_Sterile.Reduced.End$mean)
#0.5839926
se(Sterile.Ambient.End_Sterile.Reduced.End$mean)
#0.01034007


Reduced_pair_dist=droplevels(rbind(NonSterile.Ambient.End_NonSterile.Reduced.End,Sterile.Ambient.End_Sterile.Reduced.End))




qqPlot(Reduced_pair_dist$mean)
hist(Reduced_pair_dist$mean)
shapiro.test((Reduced_pair_dist$mean))
#0.3259
CommC_m1=aov(mean~Group1, data=Reduced_pair_dist)
anova(CommC_m1)
#Group1     1 0.0083323 0.0083323  4.3552 0.05566 .



######Table 2 Analyses #####




#####End Mean Distance between Communty Centriods####




#####Supplementary Data Analayses#####

####Table S1. ####

#Table S1. Mean distance (based on Weighted Unifrac for bacteria and Bray-Curtis for fungi and functional profiles) 
#between soil treatment groups and the initial soil community (see Fig. 4c, 4d in main text for ordination visual). 
#P-value tests the null hypothesis that pairwise comparison within factor is equal. Pairwise-permANOVA test if the 
#centroid of each soil treatment is significantly different than the centroid of the initial community (FDR correction).
#Bacteria
#Mean Pairwise distance

#Clay only
bact.clay.tree=subset_samples(bact.soil.tree, SoilType=="Clay")
bact.clay.tree<-prune_taxa(taxa_sums(bact.clay.tree) > 0, bact.clay.tree)
ntaxa(bact.clay.tree)
#12863
sum(otu_table(bact.clay.tree))
#2575563

#Calculate Weighted unifrac distance
bact.clay.UF=distance(bact.clay.tree,method="wunifrac", type = "samples")
#root as -- OTU617
bact.clay.map=sample_data(bact.clay.tree)


#pairwise distance using the usedist package

summary(bact.clay.map)
#Calculated the pairwise the distance between and within each treatment
bact.clay.UF_pair_dist=dist_groups(bact.clay.UF,bact.clay.map$com.rain.time)
head(bact.clay.UF_pair_dist)
nrow(bact.clay.UF_pair_dist)
#820
summary(bact.clay.UF_pair_dist)

#The below steps are for creating an easier way to parse the pairwise distance beetween treatments
bact.clay.UF_pair_dist_2=bact.clay.UF_pair_dist
colnames(bact.clay.UF_pair_dist_2)[c(1,2,3,4)]<-c("Item2","Item1","Group2","Group1")
nrow(bact.clay.UF_pair_dist_2)

bact.clay.UF_pair_dist_comb=rbind(bact.clay.UF_pair_dist,bact.clay.UF_pair_dist_2)
nrow(bact.clay.UF_pair_dist_comb)
summary(bact.clay.UF_pair_dist_comb)

#Calculated the mean distances to remove the pseudo-replication
bact.clay.UF_pair_dist_comb=bact.clay.UF_pair_dist_comb %>% group_by(Item1,Group1,Label)

bact.clay.UF_pair_dist_comb_sum=summarise_at(bact.clay.UF_pair_dist_comb,"Distance", funs(mean,n()))
nrow(bact.clay.UF_pair_dist_comb_sum)
#287
head(bact.clay.UF_pair_dist_comb_sum)
summary(bact.clay.UF_pair_dist_comb_sum)



#distance from intitial community
unique(bact.clay.UF_pair_dist_comb_sum$Label)


#Distance between Starting community and No Dispersal Ambient
Sterile.Ambient.End_Start=subset(bact.clay.UF_pair_dist_comb_sum, 
                                 Label=="Between Sterile.Ambient.End and NonSterile.Ambient.Start"&
                                   Group1=="Sterile.Ambient.End")
nrow(Sterile.Ambient.End_Start)
#8
mean(Sterile.Ambient.End_Start$mean)
#0.1643661
se(Sterile.Ambient.End_Start$mean)
#0.002640841


#Distance between Starting community and Dispersal Ambient
NonSterile.Ambient.End_Start=subset(bact.clay.UF_pair_dist_comb_sum, 
                                    Label=="Between NonSterile.Ambient.End and NonSterile.Ambient.Start"&
                                      Group1=="NonSterile.Ambient.End")
nrow(NonSterile.Ambient.End_Start)
#8
mean(NonSterile.Ambient.End_Start$mean)
#0.1402255
se(NonSterile.Ambient.End_Start$mean)
#0.002511714

#Distance between Starting community and No Dispersal Drought
Sterile.Reduced.End_Start=subset(bact.clay.UF_pair_dist_comb_sum, 
                                 Label=="Between Sterile.Reduced.End and NonSterile.Ambient.Start"&
                                   Group1=="Sterile.Reduced.End")
nrow(Sterile.Reduced.End_Start)
#8
mean(Sterile.Reduced.End_Start$mean)
#0.1876968
se(Sterile.Reduced.End_Start$mean)
#0.002649007

#Distance between Starting community and Dispersal Drought
NonSterile.Reduced.End_Start=subset(bact.clay.UF_pair_dist_comb_sum, 
                                    Label=="Between NonSterile.Reduced.End and NonSterile.Ambient.Start"&
                                      Group1=="NonSterile.Reduced.End")
nrow(NonSterile.Reduced.End_Start)
#8
mean(NonSterile.Reduced.End_Start$mean)
#0.182508
se(NonSterile.Reduced.End_Start$mean)
#0.002330879


Start_pair_dist=droplevels(rbind(Sterile.Ambient.End_Start,NonSterile.Ambient.End_Start,Sterile.Reduced.End_Start,NonSterile.Reduced.End_Start))


boxplot(mean~Group1, data=Start_pair_dist, ylab="Dissimilarity or distance from Initial Community")

qqPlot(Start_pair_dist$mean)
hist(Start_pair_dist$mean)
shapiro.test((Start_pair_dist$mean))
#0.05268
Start_m1=aov(mean~Group1, data=Start_pair_dist)
anova(Start_m1)
#Group1     3 0.011049 0.0036830  71.561 3.032e-13 ***

TukeyHSD(Start_m1)

#Distance between Starting community and the combined treatments
mean(Start_pair_dist$mean)
#0.1686991
se(Start_pair_dist$mean)
#0.00354834

#Fungi
#Mean Pairwise distance

#Clay only
fung.clay=subset_samples(fung.soil, SoilType=="Clay")
fung.clay<-prune_taxa(taxa_sums(fung.clay) > 0, fung.clay)
ntaxa(fung.clay)
#1374
sum(otu_table(fung.clay))
#1084332

fung.clay.bray=distance(fung.clay,method="bray", type = "samples")
#
fung.clay.map=sample_data(fung.clay)

#pairwise distance using the usedist package

summary(fung.clay.map)

#Calculate the pairwise distance between and within each treatment
fung.clay_pair_dist=dist_groups(fung.clay.bray,fung.clay.map$com.rain.time)
head(fung.clay_pair_dist)
nrow(fung.clay_pair_dist)
#820
summary(fung.clay_pair_dist)

#The below section of code is to make it easier to parse the pairwise distances
fung.clay_pair_dist_2=fung.clay_pair_dist
colnames(fung.clay_pair_dist_2)[c(1,2,3,4)]<-c("Item2","Item1","Group2","Group1")
nrow(fung.clay_pair_dist_2)

fung.clay_pair_dist_comb=rbind(fung.clay_pair_dist,fung.clay_pair_dist_2)
nrow(fung.clay_pair_dist_comb)
summary(fung.clay_pair_dist_comb)

#Calculate the mean distance to remove the psuedo-replication 
fung.clay_pair_dist_comb=fung.clay_pair_dist_comb %>% group_by(Item1,Group1,Label)

fung.clay_pair_dist_comb_sum=summarise_at(fung.clay_pair_dist_comb,"Distance", funs(mean,n()))
nrow(fung.clay_pair_dist_comb_sum)
#287
head(fung.clay_pair_dist_comb_sum)
summary(fung.clay_pair_dist_comb_sum)

#distance from intitial community
unique(fung.clay_pair_dist_comb_sum$Label)


#Distance from the Starting community and No Dispersal Ambient
Sterile.Ambient.End_Start=subset(fung.clay_pair_dist_comb_sum, 
                                 Label=="Between NonSterile.Ambient.Start and Sterile.Ambient.End"&
                                   Group1=="Sterile.Ambient.End")
nrow(Sterile.Ambient.End_Start)
#8
mean(Sterile.Ambient.End_Start$mean)
#0.6339347
se(Sterile.Ambient.End_Start$mean)
#0.01127929

#Distance from the Starting community and Dispersal Ambient
NonSterile.Ambient.End_Start=subset(fung.clay_pair_dist_comb_sum, 
                                    Label=="Between NonSterile.Ambient.End and NonSterile.Ambient.Start"&
                                      Group1=="NonSterile.Ambient.End")
nrow(NonSterile.Ambient.End_Start)
#8
mean(NonSterile.Ambient.End_Start$mean)
#0.6658635
se(NonSterile.Ambient.End_Start$mean)
#0.02263813

#Distance from the Starting community and No Dispersal Drought
Sterile.Reduced.End_Start=subset(fung.clay_pair_dist_comb_sum, 
                                 Label=="Between NonSterile.Ambient.Start and Sterile.Reduced.End"&
                                   Group1=="Sterile.Reduced.End")
nrow(Sterile.Reduced.End_Start)
#8
mean(Sterile.Reduced.End_Start$mean)
#0.6212315
se(Sterile.Reduced.End_Start$mean)
#0.01197376

#Distance from the Starting community and Dispersal Drought
NonSterile.Reduced.End_Start=subset(fung.clay_pair_dist_comb_sum, 
                                    Label=="Between NonSterile.Ambient.Start and NonSterile.Reduced.End"&
                                      Group1=="NonSterile.Reduced.End")
nrow(NonSterile.Reduced.End_Start)
#8
mean(NonSterile.Reduced.End_Start$mean)
#0.6276795
se(NonSterile.Reduced.End_Start$mean)
#0.008304026

Start_pair_dist=droplevels(rbind(Sterile.Ambient.End_Start,NonSterile.Ambient.End_Start,Sterile.Reduced.End_Start,NonSterile.Reduced.End_Start))



qqPlot(log(Start_pair_dist$mean))
hist(log(Start_pair_dist$mean))
shapiro.test(log(Start_pair_dist$mean))
#0.125
Start_m1=aov(log(mean)~Group1, data=Start_pair_dist)
anova(Start_m1)
#Group1     3 0.02060 0.0068666  1.7001 0.1897

TukeyHSD(Start_m1)

#all treatments to start
mean(Start_pair_dist$mean)
#0.6371773
se(Start_pair_dist$mean)
#0.007589323

####Table S1. ####



#####DOC Taxa Analyses####


#Table S2. To assess whether dispersal response was explained by response to the addition of labile DOC (dead cells), 
#we assessed changes in taxonomic groups that Cleveland et al. (2007) report to change under a DOC addition. 
#Cleveland et al. (2007), (DOC addition: 25g soils received 5mL of 650 mg/L C leachage, 15 g dry weight), found significant 
#differences in: Acidobacteria (decreased), Firmicutes (increased), Gammaproteobacteria (increased), and Nitrospora (decreased). 
#We found no significant changes in these groups under dispersal manipulation, suggesting shifts were not explained by response 
#to DOC from dead cells. Transformed for normality: Gammaproteobacteria (1/x), Firmicutes (cube root(x))

#clay focal taxa


#end only
bact.clayE.tree=subset_samples(bact.clay.tree, SampleTime=="End")


(bact.clayE.phylum<-tax_glom(bact.clayE.tree, taxrank="Phylum"))
TopPHYL = names(sort(taxa_sums(bact.clayE.phylum), TRUE)[1:10])
bact.clayET10 = prune_taxa(TopPHYL, bact.clayE.phylum)

plot_bar(bact.clayET10, x= "SampleID2", fill="Phylum")

Top20PHYL = names(sort(taxa_sums(bact.clayE.phylum), TRUE)[1:20])
bact.clayET20 = prune_taxa(Top20PHYL, bact.clayE.phylum)

plot_bar(bact.clayET20, x= "SampleID2", fill="Phylum")

bact.clayE.phylum.prop=transform_sample_counts(bact.clayE.phylum, function(x)x/sum(x))

bact.clayET10.prop = prune_taxa(TopPHYL, bact.clayE.phylum.prop)
plot_bar(bact.clayET10.prop, x= "SampleID2", fill="Phylum")+ylab("Proportion")

bact.clayET20.prop = prune_taxa(Top20PHYL, bact.clayE.phylum.prop)
plot_bar(bact.clayET20.prop, x= "SampleID2", fill="Phylum")

map.bact.soil=sample_data(bact.soil.tree)
Acido.soil.t <- subset_taxa(bact.soil.tree, Phylum=="p:Acidobacteria")
get_taxa_unique(Acido.soil.t, taxonomic.rank="Phylum")
Acido.soil.t.reads=sample_sums(Acido.soil.t)
soil.reads.t=sample_sums(bact.soil.tree)
Acido.soil.t.reads=cbind(Acido.soil.t.reads,soil.reads.t)
colnames(Acido.soil.t.reads)<-c("Acido.reads","total.reads")
Acido.soil.t.reads=merge(Acido.soil.t.reads, map.bact.soil, by ="row.names")
head(Acido.soil.t.reads)
Acido.soil.t.reads=mutate(Acido.soil.t.reads, Acido.prop=Acido.reads/total.reads)


Acido.clay.t.reads=subset(Acido.soil.t.reads, SoilType=="Clay")
Acido.clayE.t.reads=subset(Acido.clay.t.reads, SampleTime=="End")


#NUmber of Acidobacteria
bact.clayE.Acreads.t=lm((Acido.reads)~RainCom*RainLevel, data=Acido.clayE.t.reads)
qqPlot(studres(bact.clayE.Acreads.t))
hist(studres(bact.clayE.Acreads.t))
shapiro.test(studres(bact.clayE.Acreads.t))
#p-value = 0.8421
boxCox(bact.clayE.Acreads.t)
#looks okay....
Anova(bact.clayE.Acreads.t,type = 3)
#RainLevel           51432618  1    8.0472  0.008377 ** 
lsmeans(bact.clayE.Acreads.t, pairwise~RainLevel*RainCom)




#NUmber of Firmicutes
Firm.soil.t <- subset_taxa(bact.soil.tree, Phylum=="p:Firmicutes")
get_taxa_unique(Firm.soil.t, taxonomic.rank="Phylum")
Firm.soil.t.reads=sample_sums(Firm.soil.t)
soil.t.reads=sample_sums(bact.soil.tree)
Firm.soil.t.reads=cbind(Firm.soil.t.reads,soil.t.reads)
colnames(Firm.soil.t.reads)<-c("Firm.reads","total.reads")
Firm.soil.t.reads=merge(Firm.soil.t.reads, map.bact.soil, by ="row.names")
head(Firm.soil.t.reads)
Firm.soil.t.reads=mutate(Firm.soil.t.reads, Firm.prop=Firm.reads/total.reads)

Firm.clay.t.reads=subset(Firm.soil.t.reads, SoilType=="Clay")
Firm.clayE.t.reads=subset(Firm.clay.t.reads, SampleTime=="End")

bact.clayE.Fireads.t=lm(sqrt(Firm.reads)~RainCom*RainLevel, data=Firm.clayE.t.reads)
qqPlot(studres(bact.clayE.Fireads.t))
hist(studres(bact.clayE.Fireads.t))
shapiro.test(studres(bact.clayE.Fireads.t))
#p-value = 0.07459
boxCox(bact.clayE.Fireads.t)
#looks okay....
Anova(bact.clayE.Fireads.t,type = 3)
#RainLevel            391  1   13.0064  0.001194 ** 
#RainCom:RainLevel    102  1    3.3934  0.076071 . 
lsmeans(bact.clayE.Fireads.t, pairwise~RainLevel*RainCom)

#proteobacteria composition
#NUmber of gammaproteobacteria
Gproto.soil.t <- subset_taxa(bact.soil.tree, Class=="c:Gammaproteobacteria")
get_taxa_unique(Gproto.soil.t, taxonomic.rank="Class")
Gproto.soil.t.reads=sample_sums(Gproto.soil.t)
soil.t.reads=sample_sums(bact.soil.tree)
Gproto.soil.t.reads=cbind(Gproto.soil.t.reads,soil.t.reads)
colnames(Gproto.soil.t.reads)<-c("Gamma.reads","total.reads")
Gproto.soil.t.reads=merge(Gproto.soil.t.reads, map.bact.soil, by ="row.names")
head(Gproto.soil.t.reads)
Gproto.soil.t.reads=mutate(Gproto.soil.t.reads, Gamma.prop=Gamma.reads/total.reads)
Gproto.clay.t.reads=subset(Gproto.soil.t.reads, SoilType=="Clay")
Gproto.clayE.t.reads=subset(Gproto.clay.t.reads,SampleTime=="End")


bact.clayE.Greads.t=lm((Gamma.reads)^-1~RainCom*RainLevel, data=Gproto.clayE.t.reads)
qqPlot(studres(bact.clayE.Greads.t))
hist(studres(bact.clayE.Greads.t))
shapiro.test(studres(bact.clayE.Greads.t))
#p-value = 0.03308
boxCox(bact.clayE.Greads.t)
#looks okay....
Anova(bact.clayE.Greads.t,type = 3)
#RainLevel         3.5200e-07  1  20.3099 0.0001067 ***
lsmeans(bact.clayE.Greads.t, pairwise~RainLevel*RainCom)


#NUmber of Nitrospira
Nproto.soil.t <- subset_taxa(bact.soil.tree, Class=="c:Nitrospira")
get_taxa_unique(Nproto.soil.t, taxonomic.rank="Class")
Nproto.soil.t.reads=sample_sums(Nproto.soil.t)
soil.t.reads=sample_sums(bact.soil.tree)
Nproto.soil.t.reads=cbind(Nproto.soil.t.reads,soil.t.reads)
colnames(Nproto.soil.t.reads)<-c("Nitrospira.reads","total.reads")
Nproto.soil.t.reads=merge(Nproto.soil.t.reads, map.bact.soil, by ="row.names")
head(Nproto.soil.t.reads)
Nproto.soil.t.reads=mutate(Nproto.soil.t.reads, Nitrospira.prop=Nitrospira.reads/total.reads)



Nproto.clay.t.reads=subset(Nproto.soil.t.reads, SoilType=="Clay")
Nproto.clayE.t.reads=subset(Nproto.clay.t.reads,SampleTime=="End")

bact.clayE.Nreads.t=lm((Nitrospira.reads)~RainCom*RainLevel, data=Nproto.clayE.t.reads)
qqPlot(studres(bact.clayE.Nreads.t))
hist(studres(bact.clayE.Nreads.t))
shapiro.test(studres(bact.clayE.Nreads.t))
#p-value = 0.2473
boxCox(bact.clayE.Nreads.t)
#looks okay....
Anova(bact.clayE.Nreads.t,type = 3)
#nada sig 
lsmeans(bact.clayE.Nreads.t, pairwise~RainLevel*RainCom)


#####DOC Taxa Analyses####



#####Soil Taxa Shared with rain####
#Table S3. ANOVA p-values for differences in rain-OTUs (OTUs=number of taxa, reads=abundance) across soil treatments. 
#"Rain OTU" is an OTU found in any rain sample. Bacterial non-rain taxa were transformed for normality (x3). Reads are 
#the sum of taxa abundance from shared and unshared OTUs. Test statistics correspond to data shown in Fig. S2 and S3. 


#Figure S2. Number of rain-taxa (left column, 'shared with rain') and non-rain taxa (right column, 'Not shared with rain') 
#in soil Bacterial (a) and Fungal (b) communities. 


#Figure S3. Number of 'rain-taxa' (left column, 'shared with rain') and non-rain taxa (right column, 'Not shared with rain') 
#found in soil treatments for Bacteria (a) and Fungi (b). 

#Bacterial 
#Shared OTUs reads

bact.soil.n<-taxa_names(bact.soil.tree)
length(bact.soil.n)
#14893


bact.map.soil=sample_data(bact.soil.tree)

bact.rain.n<-taxa_names(bact.rain.tree)

length(bact.rain.n)
#4844



bact.Names_s_in_r<-bact.soil.n[bact.soil.n %in% bact.rain.n]
length(bact.Names_s_in_r) ## How many OTUs
#2785


length(bact.Names_s_in_r)/length(bact.soil.n)     ## Proportion OTUs in ocean, out of total in soil
#0.1870006


bact.filtered<-prune_taxa(bact.Names_s_in_r,bact.soil)

bact.sumFilteredsam=sample_sums(bact.filtered)
mean(bact.sumFilteredsam)
#42586.05


sum(sample_sums(bact.filtered))
#3492056


sum(bact.sumFilteredsam)/sum(sample_sums(bact.soil))
#0.6654544



bact.fil.sampl=as(bact.sumFilteredsam, "matrix")
colnames(bact.fil.sampl)="filt.sum"
mean(bact.fil.sampl)
#42586.05


bact.otabfil <- as(otu_table(bact.filtered), "matrix") # Taxa are rows
bact.present_absent.fil <- (bact.otabfil > 0)
bact.filt.rich <- apply(t(bact.present_absent.fil), 1, sum)
sh_bact.filt.rich=as.matrix(bact.filt.rich)
colnames(sh_bact.filt.rich)="bact.filt.rich"
mean(sh_bact.filt.rich)
#1063.598




bact.sampl=as(sample_sums(bact.soil.tree), "matrix")
colnames(bact.sampl)="all.sum"
mean(bact.sampl)
#63995.44


bact.otab <- as(otu_table(bact.soil.tree), "matrix") # Taxa are rows
bact.present_absent <- (bact.otab > 0)
bact.rich <- apply(t(bact.present_absent), 1, sum)
mean(bact.rich)
#3305.085




bact.soil.sum=cbind(bact.fil.sampl,sh_bact.filt.rich,bact.sampl,bact.rich)
bact.soil.sum=merge(bact.soil.sum, bact.map.soil, by="row.names")
bact.soil.sum$Prop_Reads_Shared=with(bact.soil.sum, bact.soil.sum$filt.sum/bact.soil.sum$all.sum)
mean(bact.soil.sum$Prop_Reads_Shared)
#0.6654133

bact.soil.sum$Prop_OTUs_Shared=with(bact.soil.sum, bact.soil.sum$bact.filt.rich/bact.soil.sum$bact.rich)
mean(bact.soil.sum$Prop_OTUs_Shared)
#0.3516787


bact.soil.sum$com.rain.soil.time=with(bact.soil.sum, interaction(RainCom, RainLevel,SoilType,SampleTime))



#end samples
bact.soilE.sum=subset(bact.soil.sum, SampleTime=="End")
mean(bact.soilE.sum$Prop_Reads_Shared)
#0.6542716

mean(bact.soilE.sum$Prop_OTUs_Shared)
#0.3181375


#non-overlapping OTUs
bact.soil_NO_TREE=otu_table(bact.soil.tree)
bact.rain_NO_TREE=otu_table(bact.rain.tree)
bact.rainx.merge=merge_phyloseq(bact.soil_NO_TREE,bact.rain_NO_TREE)

ntaxa(bact.rainx.merge)
#16952



bact.rainx.merge=phyloseq(otu_table(bact.rainx.merge),tax_table(bact.data),sample_data(bact.data),phy_tree(bact.data))
ntaxa(bact.rainx.merge)
#16952
bact.rainx.merge.tree=bact.rainx.merge


ntaxa(bact.rainx.merge.tree)
#16952



bact.Names_s_not_in_r<-bact.soil.n[bact.soil.n %w/o% bact.rain.n]
length(bact.Names_s_not_in_r)
#12108


bact.Names_r_not_in_s<-bact.rain.n[bact.rain.n %w/o% bact.soil.n]
length(bact.Names_r_not_in_s)
#2059



bact.Names_no_overlap=c(bact.Names_s_not_in_r,bact.Names_r_not_in_s)
length(bact.Names_no_overlap)
#14167


nonover.bact.filtered.tree<-prune_taxa(bact.Names_no_overlap,bact.rainx.merge.tree)
nonover.bact.filtered.tree=prune_taxa(taxa_sums(nonover.bact.filtered.tree) > 0, nonover.bact.filtered.tree)
ntaxa(nonover.bact.filtered.tree)
#14167


sum(otu_table(nonover.bact.filtered.tree))
#1980149



#look at soils
U.bact.soil.tree<-subset_samples(nonover.bact.filtered.tree, SampleType=="Soil")
U.bact.soil.tree<-prune_taxa(taxa_sums(U.bact.soil.tree)>0,U.bact.soil.tree)
ntaxa(U.bact.soil.tree)
#12108


sum(otu_table(U.bact.soil.tree))
#1755570




#prop not shared
length(bact.Names_s_not_in_r)/length(bact.soil.n)     ## Proportion OTUs in ocean, out of total in soil
#0.8129994


U.bact.soil.tree<-subset_samples(nonover.bact.filtered.tree, SampleType=="Soil")
U.bact.soil.tree<-prune_taxa(taxa_sums(U.bact.soil.tree)>0,U.bact.soil.tree)
ntaxa(U.bact.soil.tree)
#12108



U.bact.sumFilteredsam=sample_sums(U.bact.soil.tree)
mean(U.bact.sumFilteredsam)
#21409.39


sum(otu_table(U.bact.soil.tree))
#1755570


sum(U.bact.sumFilteredsam)/sum(sample_sums(bact.soil.tree))
#0.3345456


U.bact.fil.sampl=as(U.bact.sumFilteredsam, "matrix")
colnames(U.bact.fil.sampl)="filt.sum"
mean(U.bact.fil.sampl)
#21409.39


U.bact.otabfil <- as(otu_table(U.bact.soil.tree), "matrix") # Taxa are rows
U.bact.present_absent.fil <- (U.bact.otabfil > 0)
U.bact.filt.rich <- apply(t(U.bact.present_absent.fil), 1, sum)
U.bact.filt.rich=as.matrix(U.bact.filt.rich)
colnames(U.bact.filt.rich)="bact.filt.rich"
mean(U.bact.filt.rich)
#2241.488




bact.sampl=as(sample_sums(bact.soil.tree), "matrix")
colnames(bact.sampl)="all.sum"
mean(bact.sampl)
#63995.44


bact.otab <- as(otu_table(bact.soil.tree), "matrix") # Taxa are rows
bact.present_absent <- (bact.otab > 0)
bact.rich <- apply(t(bact.present_absent), 1, sum)
mean(bact.rich)
#3305.085




U.bact.soil.sum=cbind(U.bact.fil.sampl,U.bact.filt.rich,bact.sampl,bact.rich)
U.bact.soil.sum=merge(U.bact.soil.sum, bact.map.soil, by="row.names")
U.bact.soil.sum$Prop_Reads_Shared=with(U.bact.soil.sum, U.bact.soil.sum$filt.sum/U.bact.soil.sum$all.sum)
mean(U.bact.soil.sum$Prop_Reads_Shared)
#0.3345867



U.bact.soil.sum$Prop_OTUs_Shared=with(U.bact.soil.sum, U.bact.soil.sum$bact.filt.rich/U.bact.soil.sum$bact.rich)
mean(U.bact.soil.sum$Prop_OTUs_Shared)
#0.6483213


U.bact.soil.sum$com.rain.soil.time=with(U.bact.soil.sum, interaction(RainCom, RainLevel,SoilType,SampleTime))

#end samples
U.bact.soilE.sum=subset(U.bact.soil.sum, SampleTime=="End")
mean(U.bact.soilE.sum$Prop_Reads_Shared)
#0.3457284

mean(U.bact.soilE.sum$Prop_OTUs_Shared)
#0.6818625



row.names(bact.soil.sum)<-bact.soil.sum$Row.names
overl=c(rep("Shared",length(bact.soil.sum$Row.names)))
bact.soil.sum1=cbind(overl,bact.soil.sum)



row.names(U.bact.soil.sum)<-U.bact.soil.sum$Row.names
overl=c(rep("NotShared",length(U.bact.soil.sum$Row.names)))
U.bact.soil.sum1=cbind(overl,U.bact.soil.sum)

colnames(bact.soil.sum1)
colnames(U.bact.soil.sum1)
com.bact.soil.sum=rbind(bact.soil.sum1,U.bact.soil.sum1)
com.bact.soil.sum$com.rain.ovrl=with(com.bact.soil.sum, interaction(overl, com.rain.time))
length(com.bact.soil.sum$com.rain.ovrl)
#164
com.bact.soil.sum.end=subset(com.bact.soil.sum, SampleTime=="End")
length(com.bact.soil.sum.end$com.rain.ovrl)
#128
com.bact.clay.sum.end=subset(com.bact.soil.sum.end, SoilType=="Clay")
length(com.bact.clay.sum.end$com.rain.ovrl)
# 64

com.bact.clay.sum.end=com.bact.clay.sum.end %>% group_by(overl, RainLevel,RainCom)
com.bact.clay.sum.end_sum=summarise_at(com.bact.clay.sum.end, vars(bact.filt.rich,filt.sum),funs(mean,se=sd(.)/sqrt(n()),sd))
com.bact.clay.sum.end_sum$overl_RainCom=with(com.bact.clay.sum.end_sum,interaction(overl,RainCom))
com.bact.clay.sum.end_sum$overl_RainLevel=with(com.bact.clay.sum.end_sum,interaction(RainLevel,overl))

#####Bact Rain Overlap Graphs####




bact.clayE.OTUs_g2=ggplot(com.bact.clay.sum.end_sum, aes(overl_RainLevel, bact.filt.rich_mean, 
                                                         ymin = bact.filt.rich_mean-bact.filt.rich_se, ymax = bact.filt.rich_mean+bact.filt.rich_se))

(p_bact.clayE.OTUs_g2=bact.clayE.OTUs_g2+geom_errorbar(width=0.25)+geom_line(aes(group = overl_RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "# of OTUs")+xlab(NULL)+scale_x_discrete(labels=c("Ambient","Drought","Ambient","Drought"))+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_vline(xintercept = 2.5, size=1.2)+ggtitle(label = "Bacteria")+
    theme(plot.title = element_text(hjust = 0.5, size=28)))





bact.clayE.reads_g2=ggplot(com.bact.clay.sum.end_sum, aes(overl_RainLevel, filt.sum_mean, 
                                                          ymin = filt.sum_mean-filt.sum_se, ymax = filt.sum_mean+filt.sum_se))

(p_bact.clayE.reads_g2=bact.clayE.reads_g2+geom_errorbar(width=0.25)+geom_line(aes(group = overl_RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Reads")+xlab(NULL)+scale_x_discrete(labels=c("Ambient","Drought","Ambient","Drought"))+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_vline(xintercept = 2.5, size=1.2)+ggtitle(label = "Bacteria")+
    theme(plot.title = element_text(hjust = 0.5, size=28)))

#SHARED
bact.clayE.sum=subset(com.bact.clay.sum.end, overl=="Shared")
mean(bact.clayE.sum$Prop_Reads_Shared)
#0.6654412
mean(bact.clayE.sum$Prop_OTUs_Shared)
#0.3178995
#number of the OTUs shared with rain

sh.bact.clayE.t.OTU=lm((bact.filt.rich)~RainCom*RainLevel, data=bact.clayE.sum)
qqPlot(studres(sh.bact.clayE.t.OTU))
hist(studres(sh.bact.clayE.t.OTU))
shapiro.test(studres(sh.bact.clayE.t.OTU))
#p-value = 0.6587
boxCox(sh.bact.clayE.t.OTU)
#looks okay....
Anova(sh.bact.clayE.t.OTU,type = 3)
#RainCom              52488  1    14.8369 0.0006242 ***

lsmeans(sh.bact.clayE.t.OTU, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile  -11.375 29.73916 28  -0.382  0.9806
Ambient,NonSterile - Ambient,Sterile     -85.375 29.73916 28  -2.871  0.0365
Ambient,NonSterile - Reduced,Sterile     -88.000 29.73916 28  -2.959  0.0298
Reduced,NonSterile - Ambient,Sterile     -74.000 29.73916 28  -2.488  0.0838
Reduced,NonSterile - Reduced,Sterile     -76.625 29.73916 28  -2.577  0.0697
Ambient,Sterile - Reduced,Sterile         -2.625 29.73916 28  -0.088  0.9997"


#number of the reads shared with rain

sh.bact.clayE.t.reads=lm((filt.sum)~RainCom*RainLevel, data=bact.clayE.sum)
qqPlot(studres(sh.bact.clayE.t.reads))
hist(studres(sh.bact.clayE.t.reads))
shapiro.test(studres(sh.bact.clayE.t.reads))
#p-value = 0.147
boxCox(sh.bact.clayE.t.reads)
#looks okay....
Anova(sh.bact.clayE.t.reads,type = 3)
#RainLevel         6.7484e+08  1   13.6159 0.0009585 *** 


lsmeans(sh.bact.clayE.t.reads, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                  estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile -10088.000 3520.037 28  -2.866  0.0369
Ambient,NonSterile - Ambient,Sterile     -2206.125 3520.037 28  -0.627  0.9226
Ambient,NonSterile - Reduced,Sterile    -10487.125 3520.037 28  -2.979  0.0285
Reduced,NonSterile - Ambient,Sterile      7881.875 3520.037 28   2.239  0.1374
Reduced,NonSterile - Reduced,Sterile      -399.125 3520.037 28  -0.113  0.9995
Ambient,Sterile - Reduced,Sterile        -8281.000 3520.037 28  -2.353  0.1103"
#NOT SHARED
U.bact.clayE.sum=subset(com.bact.clay.sum.end, overl=="NotShared")
#number of the OTUs not shared with rain

U.bact.clayE.t.OTU=lm((bact.filt.rich)~RainCom*RainLevel, data=U.bact.clayE.sum)
qqPlot(studres(U.bact.clayE.t.OTU))
hist(studres(U.bact.clayE.t.OTU))
shapiro.test(studres(U.bact.clayE.t.OTU))
#p-value = 0.1015
boxCox(U.bact.clayE.t.OTU)
#looks okay....
Anova(U.bact.clayE.t.OTU,type = 3)
#RainCom             1002174  1   23.4474 4.269e-05 ***
#RainLevel            500750  1   11.7158  0.001926 ** 

lsmeans(U.bact.clayE.t.OTU, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile  255.750 103.3698 28   2.474  0.0863
Ambient,NonSterile - Ambient,Sterile    -348.375 103.3698 28  -3.370  0.0112
Ambient,NonSterile - Reduced,Sterile    -103.750 103.3698 28  -1.004  0.7485
Reduced,NonSterile - Ambient,Sterile    -604.125 103.3698 28  -5.844  <.0001
Reduced,NonSterile - Reduced,Sterile    -359.500 103.3698 28  -3.478  0.0086
Ambient,Sterile - Reduced,Sterile        244.625 103.3698 28   2.367  0.1073"


#number of the reads not shared with rain

U.bact.clayE.t.reads=lm((filt.sum)~RainCom*RainLevel, data=U.bact.clayE.sum)
qqPlot(studres(U.bact.clayE.t.reads))
hist(studres(U.bact.clayE.t.reads))
shapiro.test(studres(U.bact.clayE.t.reads))
#p-value = 0.07022
boxCox(U.bact.clayE.t.reads)
#looks okay....
Anova(U.bact.clayE.t.reads,type = 3)
#nada sig

lsmeans(U.bact.clayE.t.reads, pairwise~RainLevel*RainCom)



bact.filtered_map=sample_data(bact.filtered)
bact.filtered.divfil=estimate_richness(bact.filtered,measures=alpha_meas)

bact.filtered.divfil=merge(bact.filtered.divfil, bact.filtered_map, by ="row.names")
bact.filtered.divfil=mutate(bact.filtered.divfil, pielou=Shannon*(1/log(Observed)))


bact.filteredE.divfil=subset(bact.filtered.divfil, SampleTime=="End")
sh.bact.clayE.t.divfil=subset(bact.filteredE.divfil, SoilType=="Clay")

#Chao1
sh.bact.clayE.t.chao1=lm((Chao1)~RainCom*RainLevel, data=sh.bact.clayE.t.divfil)
qqPlot(studres(sh.bact.clayE.t.chao1))
hist(studres(sh.bact.clayE.t.chao1))
shapiro.test(studres(sh.bact.clayE.t.chao1))
#p-value = 0.7132
boxCox(sh.bact.clayE.t.chao1)
#looks okay....
Anova(sh.bact.clayE.t.chao1,type = 3)
#RainCom              64905  1    17.7013 0.0002405 ***

lsmeans(sh.bact.clayE.t.chao1, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                   estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile    3.475999 30.27662 28   0.115  0.9994
Ambient,NonSterile - Ambient,Sterile    -102.418374 30.27662 28  -3.383  0.0108
Ambient,NonSterile - Reduced,Sterile     -74.251845 30.27662 28  -2.452  0.0902
Reduced,NonSterile - Ambient,Sterile    -105.894372 30.27662 28  -3.498  0.0081
Reduced,NonSterile - Reduced,Sterile     -77.727844 30.27662 28  -2.567  0.0711
Ambient,Sterile - Reduced,Sterile         28.166529 30.27662 28   0.930  0.7889"


#t0 Samples Clay

bact.filtered.t0.divfil=subset(bact.filtered.divfil, SampleTime=="Start")
sh.bact.clay.t0.divfil=subset(bact.filtered.t0.divfil, SoilType=="Clay")
mean(sh.bact.clay.t0.divfil$Chao1)

sh.bact.clayE.t.divfil=sh.bact.clayE.t.divfil %>% group_by(RainCom, RainLevel)
sh.bact.clayE.t.divfil_sum=summarise_at(sh.bact.clayE.t.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


sh.bact.clayE.t.chao1_g=ggplot(sh.bact.clayE.t.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(sh.bact.clayE_p1=sh.bact.clayE.t.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Bacterial richness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(sh.bact.clay.t0.divfil$Chao1), size=1.2))


#Not shared



nonover.bact.filtered_map=sample_data(nonover.bact.filtered.tree)
nonover.bact.filtered.divfil=estimate_richness(nonover.bact.filtered.tree,measures=alpha_meas)

nonover.bact.filtered.divfil=merge(nonover.bact.filtered.divfil, nonover.bact.filtered_map, by ="row.names")
nonover.bact.filtered.divfil=mutate(nonover.bact.filtered.divfil, pielou=Shannon*(1/log(Observed)))

U.bact.filteredE.divfil=subset(nonover.bact.filtered.divfil, SampleTime=="End")
U.bact.clayE.t.divfil=subset(U.bact.filteredE.divfil, SoilType=="Clay")

#Chao1
U.bact.clayE.t.chao1=lm((Chao1)~RainCom*RainLevel, data=U.bact.clayE.t.divfil)
qqPlot(studres(U.bact.clayE.t.chao1))
hist(studres(U.bact.clayE.t.chao1))
shapiro.test(studres(U.bact.clayE.t.chao1))
#p-value = 0.1221
boxCox(U.bact.clayE.t.chao1)
#looks okay....
Anova(U.bact.clayE.t.chao1,type = 3)
#RainCom             1599656  1   29.7461 8.029e-06 ***
#RainLevel           1023221  1   19.0271 0.0001581 ***

lsmeans(U.bact.clayE.t.chao1, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                  estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile  368.58855 115.9494 28   3.179  0.0178
Ambient,NonSterile - Ambient,Sterile    -436.21163 115.9494 28  -3.762  0.0042
Ambient,NonSterile - Reduced,Sterile     -89.53075 115.9494 28  -0.772  0.8663
Reduced,NonSterile - Ambient,Sterile    -804.80018 115.9494 28  -6.941  <.0001
Reduced,NonSterile - Reduced,Sterile    -458.11930 115.9494 28  -3.951  0.0026
Ambient,Sterile - Reduced,Sterile        346.68088 115.9494 28   2.990  0.0278
"
U.bact.filtered.t0.divfil=subset(nonover.bact.filtered.divfil, SampleTime=="Start")
U.bact.clay.t0.divfil=subset(U.bact.filtered.t0.divfil, SoilType=="Clay")
mean(U.bact.clay.t0.divfil$Chao1)

U.bact.clayE.t.divfil=U.bact.clayE.t.divfil %>% group_by(RainCom, RainLevel)
U.bact.clayE.t.divfil_sum=summarise_at(U.bact.clayE.t.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


U.bact.clayE.t.chao1_g=ggplot(U.bact.clayE.t.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(U.bact.clayE_p1=U.bact.clayE.t.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Bacterial richness")+xlab(NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    theme_bw()+theme(legend.position = "none",axis.text=element_text(size=18),axis.title.x=element_blank(), 
                     axis.title.y=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(U.bact.clay.t0.divfil$Chao1), size=1.2))



#Fungi
#shared with soil
fung.soil.n<-taxa_names(fung.soil)
fung.map.soil=sample_data(fung.soil)

fung.rain.n<-taxa_names(fung.rain)

fung.Names_s_in_r<-fung.soil.n[fung.soil.n %in% fung.rain.n]
length(fung.Names_s_in_r) ## How many OTUs
#360


length(fung.Names_s_in_r)/length(fung.soil.n)     
#0.1957586


fung.filtered<-prune_taxa(fung.Names_s_in_r,fung.soil)
ntaxa(fung.filtered)
#360



sum(otu_table(fung.filtered))
#1520492


fung.sumFilteredsam=sample_sums(fung.filtered)
mean(fung.sumFilteredsam)
#18542.59


sum(fung.sumFilteredsam)/sum(sample_sums(fung.soil))
#0.7070277



sh.fung.fil.sampl=as(fung.sumFilteredsam, "matrix")
colnames(sh.fung.fil.sampl)="filt.sum"
mean(sh.fung.fil.sampl)
#18542.59


sh.fung.otabfil <- as(otu_table(fung.filtered), "matrix") # Taxa are rows
sh.fung.present_absent.fil <- (sh.fung.otabfil > 0)
sh.fung.filt.rich <- as(apply(t(sh.fung.present_absent.fil), 1, sum), "matrix")
colnames(sh.fung.filt.rich)="fung.filt.rich"
mean(sh.fung.filt.rich)
#73.52439


fung.sampl=as(sample_sums(fung.soil), "matrix")
colnames(fung.sampl)="all.sum"
mean(fung.sampl)
#26226.11


fung.otab <- as(otu_table(fung.soil), "matrix") # Taxa are rows
fung.present_absent <- (fung.otab > 0)
fung.rich <- apply(t(fung.present_absent), 1, sum)
mean(fung.rich)
#194.3293




fung.soil.sum=cbind(sh.fung.fil.sampl,sh.fung.filt.rich,fung.sampl,fung.rich)
fung.soil.sum=merge(fung.soil.sum, fung.map.soil, by="row.names")
fung.soil.sum$Prop_Reads_Shared=with(fung.soil.sum, fung.soil.sum$filt.sum/fung.soil.sum$all.sum)
mean(fung.soil.sum$Prop_Reads_Shared)
#0.7055521


fung.soil.sum$Prop_OTUs_Shared=with(fung.soil.sum, fung.soil.sum$fung.filt.rich/fung.soil.sum$fung.rich)
mean(fung.soil.sum$Prop_OTUs_Shared)
#0.391438


fung.soil.sum$com.rain.soil.time=with(fung.soil.sum, interaction(RainCom, RainLevel,SoilType,SampleTime))



fung.soilE.sum=subset(fung.soil.sum, SampleTime=="End")
mean(fung.soilE.sum$Prop_Reads_Shared)
#0.6843251

mean(fung.soilE.sum$Prop_OTUs_Shared)
#0.371035

#non-overlapping OTUs
fung.rainx.merge=merge_phyloseq(fung.soil,fung.rain)
ntaxa(fung.rainx.merge)
#3930


fung.Names_s_not_in_r<-fung.soil.n[fung.soil.n %w/o% fung.rain.n]
length(fung.Names_s_not_in_r)
#1479

fung.Names_r_not_in_s<-fung.rain.n[fung.rain.n %w/o% fung.soil.n]
length(fung.Names_r_not_in_s)
#2091
fung.Names_no_overlap=c(fung.Names_s_not_in_r,fung.Names_r_not_in_s)
length(fung.Names_no_overlap)
#3570


#soil
U.fung.rainx<-prune_taxa(fung.Names_no_overlap,fung.rainx.merge)
ntaxa(U.fung.rainx)
#3570

sum(otu_table(U.fung.rainx))
#899233


#look at soils
U.fung.soil<-subset_samples(U.fung.rainx, SampleType=="Soil")
U.fung.soil<-prune_taxa(taxa_sums(U.fung.soil)>0,U.fung.soil)
ntaxa(U.fung.soil)
#1479

sum(otu_table(U.fung.soil))
#630049



#prop not shared
length(fung.Names_s_not_in_r)/length(fung.soil.n)     ## Proportion OTUs in ocean, out of total in soil
#0.8042414

U.fung.sumFilteredsam=sample_sums(U.fung.soil)
mean(U.fung.sumFilteredsam)
#7683.524

sum(otu_table(U.fung.soil))
#630049

sum(U.fung.sumFilteredsam)/sum(sample_sums(fung.soil))
#0.2929723

U.fung.fil.sampl=as(U.fung.sumFilteredsam, "matrix")
colnames(U.fung.fil.sampl)="filt.sum"
mean(U.fung.fil.sampl)
#7683.524

U.fung.otabfil <- as(otu_table(U.fung.soil), "matrix") # Taxa are rows
U.fung.present_absent.fil <- (U.fung.otabfil > 0)
U.fung.filt.rich <- apply(t(U.fung.present_absent.fil), 1, sum)
U.fung.filt.rich=as.matrix(U.fung.filt.rich)
colnames(U.fung.filt.rich)="fung.filt.rich"
mean(U.fung.filt.rich)
#120.8049



fung.sampl=as(sample_sums(fung.soil), "matrix")
colnames(fung.sampl)="all.sum"
mean(fung.sampl)
#26226.11

fung.otab <- as(otu_table(fung.soil), "matrix") # Taxa are rows
fung.present_absent <- (fung.otab > 0)
fung.rich <- apply(t(fung.present_absent), 1, sum)
mean(fung.rich)
#194.3293

U.fung.soil.sum=cbind(U.fung.fil.sampl,U.fung.filt.rich,fung.sampl,fung.rich)
U.fung.soil.sum=merge(U.fung.soil.sum, fung.map.soil, by="row.names")
U.fung.soil.sum$Prop_Reads_Shared=with(U.fung.soil.sum, U.fung.soil.sum$filt.sum/U.fung.soil.sum$all.sum)
mean(U.fung.soil.sum$Prop_Reads_Shared)
#0.2944479

U.fung.soil.sum$Prop_OTUs_Shared=with(U.fung.soil.sum, U.fung.soil.sum$fung.filt.rich/U.fung.soil.sum$fung.rich)
mean(U.fung.soil.sum$Prop_OTUs_Shared)
#0.608562

U.fung.soil.sum$com.rain.soil.time=with(U.fung.soil.sum, interaction(RainCom, RainLevel,SoilType,SampleTime))

U.fung.soilE.sum=subset(U.fung.soil.sum, SampleTime=="End")
mean(U.fung.soilE.sum$Prop_Reads_Shared)
#0.3156749

mean(U.fung.soilE.sum$Prop_OTUs_Shared)
#0.628965




row.names(fung.soil.sum)<-fung.soil.sum$Row.names
overl=c(rep("Shared",length(fung.soil.sum$Row.names)))
fung.soil.sum1=cbind(overl,fung.soil.sum)



row.names(U.fung.soil.sum)<-U.fung.soil.sum$Row.names
overl=c(rep("NotShared",length(U.fung.soil.sum$Row.names)))
U.fung.soil.sum1=cbind(overl,U.fung.soil.sum)

colnames(fung.soil.sum1)
colnames(U.fung.soil.sum1)
com.fung.soil.sum=rbind(fung.soil.sum1,U.fung.soil.sum1)
com.fung.soil.sum$com.rain.ovrl=with(com.fung.soil.sum, interaction(overl, com.rain.time))
length(com.fung.soil.sum$com.rain.ovrl)
#164
com.fung.soil.sum.end=subset(com.fung.soil.sum, SampleTime=="End")
length(com.fung.soil.sum.end$com.rain.ovrl)
#128
com.fung.clay.sum.end=subset(com.fung.soil.sum.end, SoilType=="Clay")
length(com.fung.clay.sum.end$com.rain.ovrl)
# 64

com.fung.clay.sum.end=com.fung.clay.sum.end %>% group_by(overl, RainLevel,RainCom)
com.fung.clay.sum.end_sum=summarise_at(com.fung.clay.sum.end, vars(fung.filt.rich,filt.sum),funs(mean,se=sd(.)/sqrt(n()),sd))
com.fung.clay.sum.end_sum$overl_RainCom=with(com.fung.clay.sum.end_sum,interaction(overl,RainCom))
com.fung.clay.sum.end_sum$overl_RainLevel=with(com.fung.clay.sum.end_sum,interaction(RainLevel,overl))

#####FungiRain Overlap Graphs####




fung.clayE.OTUs_g2=ggplot(com.fung.clay.sum.end_sum, aes(overl_RainLevel, fung.filt.rich_mean, 
                                                         ymin = fung.filt.rich_mean-fung.filt.rich_se, ymax = fung.filt.rich_mean+fung.filt.rich_se))

(p_fung.clayE.OTUs_g2=fung.clayE.OTUs_g2+geom_errorbar(width=0.25)+geom_line(aes(group = overl_RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    ylab(NULL)+xlab(NULL)+scale_x_discrete(labels=c("Ambient","Drought","Ambient","Drought"))+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_vline(xintercept = 2.5, size=1.2)+ggtitle(label = "Fungi")+
    theme(plot.title = element_text(hjust = 0.5, size=28)))




fung.clayE.reads_g2=ggplot(com.fung.clay.sum.end_sum, aes(overl_RainLevel, filt.sum_mean, 
                                                          ymin = filt.sum_mean-filt.sum_se, ymax = filt.sum_mean+filt.sum_se))

(p_fung.clayE.reads_g2=fung.clayE.reads_g2+geom_errorbar(width=0.25)+geom_line(aes(group = overl_RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    ylab(NULL)+xlab(NULL)+scale_x_discrete(labels=c("Ambient","Drought","Ambient","Drought"))+
    scale_y_continuous(breaks = seq(0, 30000, by = 5000))+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_vline(xintercept = 2.5, size=1.2)+ggtitle(label = "Fungi")+
    theme(plot.title = element_text(hjust = 0.5, size=28)))


#OTU Number combined graph
grid.arrange(p_bact.clayE.OTUs_g2,p_fung.clayE.OTUs_g2, ncol=2,widths=c(1.1,1))

#Read number combined graph
grid.arrange(p_bact.clayE.reads_g2,p_fung.clayE.reads_g2, ncol=2,widths=c(1.1,1))


#SHARED
fung.clayE.sum=subset(com.fung.clay.sum.end, overl=="Shared")
mean(fung.clayE.sum$Prop_Reads_Shared)
#0.7068355

mean(fung.clayE.sum$Prop_OTUs_Shared)
#0.3781697

#number of the OTUs shared with rain

sh.fung.clayE.OTU=lm((fung.filt.rich)~RainCom*RainLevel, data=fung.clayE.sum)
qqPlot(studres(sh.fung.clayE.OTU))
hist(studres(sh.fung.clayE.OTU))
shapiro.test(studres(sh.fung.clayE.OTU))
#p-value = 0.2818
boxCox(sh.fung.clayE.OTU)
#looks okay....
Anova(sh.fung.clayE.OTU,type = 3)
#RainCom              780  1   10.1397  0.003543 ** 

lsmeans(sh.fung.clayE.OTU, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile   -4.375 4.385701 28  -0.998  0.7520
Ambient,NonSterile - Ambient,Sterile     -11.250 4.385701 28  -2.565  0.0714
Ambient,NonSterile - Reduced,Sterile     -12.875 4.385701 28  -2.936  0.0315
Reduced,NonSterile - Ambient,Sterile      -6.875 4.385701 28  -1.568  0.4127
Reduced,NonSterile - Reduced,Sterile      -8.500 4.385701 28  -1.938  0.2356
Ambient,Sterile - Reduced,Sterile         -1.625 4.385701 28  -0.371  0.9823"



#number of the reads shared with rain

sh.fung.clayE.reads=lm((filt.sum)~RainCom*RainLevel, data=fung.clayE.sum)
qqPlot(studres(sh.fung.clayE.reads))
hist(studres(sh.fung.clayE.reads))
shapiro.test(studres(sh.fung.clayE.reads))
#p-value = 0.1855
boxCox(sh.fung.clayE.reads)
#looks okay....
Anova(sh.fung.clayE.reads,type = 3)
#nada 


lsmeans(sh.fung.clayE.reads, pairwise~RainLevel*RainCom)

#SHARED
U.fung.clayE.sum=subset(com.fung.clay.sum.end, overl=="NotShared")
#number of the OTUs NOT shared with rain

U.fung.clayE.OTU=lm((fung.filt.rich)~RainCom*RainLevel, data=U.fung.clayE.sum)
qqPlot(studres(U.fung.clayE.OTU))
hist(studres(U.fung.clayE.OTU))
shapiro.test(studres(U.fung.clayE.OTU))
#p-value = 0.329
boxCox(U.fung.clayE.OTU)
#looks okay....
Anova(U.fung.clayE.OTU,type = 3)
#RainCom             5513  1   11.3725  0.002194 ** 

lsmeans(U.fung.clayE.OTU, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile   -8.125 11.00822 28  -0.738  0.8809
Ambient,NonSterile - Ambient,Sterile     -30.750 11.00822 28  -2.793  0.0435
Ambient,NonSterile - Reduced,Sterile     -29.875 11.00822 28  -2.714  0.0518
Reduced,NonSterile - Ambient,Sterile     -22.625 11.00822 28  -2.055  0.1926
Reduced,NonSterile - Reduced,Sterile     -21.750 11.00822 28  -1.976  0.2211
Ambient,Sterile - Reduced,Sterile          0.875 11.00822 28   0.079  0.9998"


#number of the reads NOT shared with rain

U.fung.clayE.reads=lm((filt.sum)~RainCom*RainLevel, data=U.fung.clayE.sum)
qqPlot(studres(U.fung.clayE.reads))
hist(studres(U.fung.clayE.reads))
shapiro.test(studres(U.fung.clayE.reads))
#p-value = 0.6274
boxCox(U.fung.clayE.reads)
#looks okay....
Anova(U.fung.clayE.reads,type = 3)
#RainLevel           31133941  1   4.6301 0.04019 *  
#RainCom:RainLevel   20457606  1   3.0423 0.09209 . 

lsmeans(U.fung.clayE.reads, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                 estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile -3571.875 1296.564 28  -2.755  0.0474
Ambient,NonSterile - Ambient,Sterile    -2076.750 1296.564 28  -1.602  0.3940
Ambient,NonSterile - Reduced,Sterile    -2450.375 1296.564 28  -1.890  0.2552
Reduced,NonSterile - Ambient,Sterile     1495.125 1296.564 28   1.153  0.6605
Reduced,NonSterile - Reduced,Sterile     1121.500 1296.564 28   0.865  0.8227
Ambient,Sterile - Reduced,Sterile        -373.625 1296.564 28  -0.288  0.9915"



fung.filtered_map=sample_data(fung.filtered)
fung.filtered.divfil=estimate_richness(fung.filtered,measures=alpha_meas)

fung.filtered.divfil=merge(fung.filtered.divfil, fung.filtered_map, by ="row.names")
fung.filtered.divfil=mutate(fung.filtered.divfil, pielou=Shannon*(1/log(Observed)))


fung.filteredE.divfil=subset(fung.filtered.divfil, SampleTime=="End")
sh.fung.clayE.divfil=subset(fung.filteredE.divfil, SoilType=="Clay")


#Chao1
sh.fung.clayE.chao1=lm((Chao1)~RainCom*RainLevel, data=sh.fung.clayE.divfil)
qqPlot(studres(sh.fung.clayE.chao1))
hist(studres(sh.fung.clayE.chao1))
shapiro.test(studres(sh.fung.clayE.chao1))
#p-value = 0.9229
boxCox(sh.fung.clayE.chao1)
#looks okay....
Anova(sh.fung.clayE.chao1,type = 3)
#RainCom             1016  1    8.7943  0.006116 **  

lsmeans(sh.fung.clayE.chao1, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                  estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile  -0.968750 5.373894 28  -0.180  0.9979
Ambient,NonSterile - Ambient,Sterile     -9.639583 5.373894 28  -1.794  0.2974
Ambient,NonSterile - Reduced,Sterile    -13.866667 5.373894 28  -2.580  0.0691
Reduced,NonSterile - Ambient,Sterile     -8.670833 5.373894 28  -1.614  0.3876
Reduced,NonSterile - Reduced,Sterile    -12.897917 5.373894 28  -2.400  0.1003
Ambient,Sterile - Reduced,Sterile        -4.227083 5.373894 28  -0.787  0.8599"

fung.filtered.t0.divfil=subset(fung.filtered.divfil, SampleTime=="Start")
sh.fung.clay.t0.divfil=subset(fung.filtered.t0.divfil, SoilType=="Clay")
mean(sh.fung.clay.t0.divfil$Chao1)

sh.fung.clayE.divfil=sh.fung.clayE.divfil %>% group_by(RainCom, RainLevel)
sh.fung.clayE.divfil_sum=summarise_at(sh.fung.clayE.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


sh.fung.clayE.chao1_g=ggplot(sh.fung.clayE.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(sh.fung.clayE_p1=sh.fung.clayE.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Fungal richness")+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(sh.fung.clay.t0.divfil$Chao1), size=1.2))

#Not shared

U.fung.filtered_map=sample_data(U.fung.rainx)
U.fung.filtered.divfil=estimate_richness(U.fung.rainx,measures=alpha_meas)

U.fung.filtered.divfil=merge(U.fung.filtered.divfil, U.fung.filtered_map, by ="row.names")
U.fung.filtered.divfil=mutate(U.fung.filtered.divfil, pielou=Shannon*(1/log(Observed)))


U.fung.filteredE.divfil=subset(U.fung.filtered.divfil, SampleTime=="End")
U.fung.clayE.divfil=subset(U.fung.filteredE.divfil, SoilType=="Clay")


#Chao1
U.fung.clayE.chao1=lm((Chao1)~RainCom*RainLevel, data=U.fung.clayE.divfil)
qqPlot(studres(U.fung.clayE.chao1))
hist(studres(U.fung.clayE.chao1))
shapiro.test(studres(U.fung.clayE.chao1))
#p-value = 0.5472
boxCox(U.fung.clayE.chao1)
#looks okay....
Anova(U.fung.clayE.chao1,type = 3)
#RainCom             6211  1   10.7420  0.002797 **  

lsmeans(U.fung.clayE.chao1, pairwise~RainLevel*RainCom)
"$contrasts
contrast                                  estimate       SE df t.ratio p.value
Ambient,NonSterile - Reduced,NonSterile -11.555714 12.02313 28  -0.961  0.7723
Ambient,NonSterile - Ambient,Sterile    -33.037227 12.02313 28  -2.748  0.0481
Ambient,NonSterile - Reduced,Sterile    -34.246683 12.02313 28  -2.848  0.0384
Reduced,NonSterile - Ambient,Sterile    -21.481513 12.02313 28  -1.787  0.3007
Reduced,NonSterile - Reduced,Sterile    -22.690969 12.02313 28  -1.887  0.2563
Ambient,Sterile - Reduced,Sterile        -1.209456 12.02313 28  -0.101  0.9996"

U.fung.filtered.t0.divfil=subset(U.fung.filtered.divfil, SampleTime=="Start")
U.fung.clay.t0.divfil=subset(U.fung.filtered.t0.divfil, SoilType=="Clay")
mean(U.fung.clay.t0.divfil$Chao1)

U.fung.clayE.divfil=U.fung.clayE.divfil %>% group_by(RainCom, RainLevel)
U.fung.clayE.divfil_sum=summarise_at(U.fung.clayE.divfil, vars(Chao1,pielou),funs(mean,se=sd(.)/sqrt(n()),sd))


U.fung.clayE.chao1_g=ggplot(U.fung.clayE.divfil_sum, aes(RainLevel, Chao1_mean, ymin = Chao1_mean-Chao1_se, ymax = Chao1_mean+Chao1_se))

(U.fung.clayE_p1=U.fung.clayE.chao1_g+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_y_continuous(name = "Fungal richness")+xlab(NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    theme_bw()+theme(legend.position = "none",axis.text=element_text(size=18),axis.title.x=element_blank(), 
                     axis.title.y=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(U.fung.clay.t0.divfil$Chao1), size=1.2))


#combined richness of shared and not shared taxa
grid.arrange(sh.bact.clayE_p1, sh.fung.clayE_p1, U.bact.clayE_p1,U.fung.clayE_p1, ncol=2,widths=c(1,1),heights=c(1,1.1))



#####Soil Taxa Shared with rain####



#####Drought Tolerant Taxa Analyses####

#Figure S1. Change in relative abundance across soil treatments of Phyla (a: Actinobacteria and b: Planctomycetes) 
#reported to increase in abundance under drought in previous studies (Evans et al. 2014, Evans and Wallenstein 2014, Bouskill et al. 2016). 
#Tables show 2-factor ANOVA test statistics. Bold line indicates mean abundance at time 0 (N=3 cores). 




#Planctomycetes
Planct.soil.t <- subset_taxa(bact.soil.tree, Phylum=="p:Planctomycetes")
get_taxa_unique(Planct.soil.t, taxonomic.rank="Phylum")
Planct.soil.t.reads=sample_sums(Planct.soil.t)
soil.reads.t=sample_sums(bact.soil.tree)
Planct.soil.t.reads=cbind(Planct.soil.t.reads,soil.reads.t)
colnames(Planct.soil.t.reads)<-c("Planct.reads","total.reads")
Planct.soil.t.reads=merge(Planct.soil.t.reads, bact.map2, by ="row.names")
head(Planct.soil.t.reads)
Planct.soil.t.reads=mutate(Planct.soil.t.reads, Planct.prop=Planct.reads/total.reads)

Planct.soilE.t.reads=subset(Planct.soil.t.reads, SampleTime=="End")
Planct.clayE.t.reads=subset(Planct.soilE.t.reads, SoilType=="Clay")


#Prop of Planctomycetes
bact.clayE.Preads.prop.t=lm((Planct.prop)~RainCom*RainLevel, data=Planct.clayE.t.reads)
qqPlot(studres(bact.clayE.Preads.prop.t))
hist(studres(bact.clayE.Preads.prop.t))
shapiro.test(studres(bact.clayE.Preads.prop.t))
#p-value = 0.8628
boxCox(bact.clayE.Preads.prop.t)
#looks okay....
Anova(bact.clayE.Preads.prop.t,type = 3)
#RainCom           0.003995  1   61.9890 1.415e-08 ***
#RainLevel         0.000808  1   12.5427  0.001415 ** 
#RainCom:RainLevel 0.000528  1    8.1895  0.007885 ** 
lsmeans(bact.clayE.Preads.prop.t, pairwise~RainLevel*RainCom)



positions <- c("NonSterile.Ambient.End", "NonSterile.Reduced.End", 
               "Sterile.Ambient.End", "Sterile.Reduced.End")


Planct.soil.t.readsT0=subset(Planct.soil.t.reads, SampleTime=="Start")
Planct.clayE.t.readsT0=subset(Planct.soil.t.readsT0, SoilType=="Clay")
mean(Planct.clayE.t.readsT0$Planct.prop)


#
scaleFUN3 <- function(x) sprintf("%.3f", x)
Pprop_clay=ggplot(Planct.clayE.t.reads_sum, aes(RainLevel, Planct.prop_mean, ymin =Planct.prop_mean-Planct.prop_se,
                                                ymax = Planct.prop_mean+Planct.prop_se))

(bact.clayE_Planct=Pprop_clay+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    scale_y_continuous(name = "Proportion of reads Planctomycetes",labels=scaleFUN3)+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(Planct.clayE.t.readsT0$Planct.prop), size=1.2))


#Actinobacteria
Actin.soil.t <- subset_taxa(bact.soil.tree, Phylum=="p:Actinobacteria")
get_taxa_unique(Actin.soil.t, taxonomic.rank="Phylum")
Actin.soil.t.reads=sample_sums(Actin.soil.t)
soil.reads.t=sample_sums(bact.soil.tree)
Actin.soil.t.reads=cbind(Actin.soil.t.reads,soil.reads.t)
colnames(Actin.soil.t.reads)<-c("Actin.reads","total.reads")
Actin.soil.t.reads=merge(Actin.soil.t.reads, bact.map2, by ="row.names")
head(Actin.soil.t.reads)
Actin.soil.t.reads=mutate(Actin.soil.t.reads, Actin.prop=Actin.reads/total.reads)

Actin.soilE.t.reads=subset(Actin.soil.t.reads, SampleTime=="End")
Actin.clayE.t.reads=subset(Actin.soilE.t.reads, SoilType=="Clay")


#Prop of pctinobacteria
bact.clayE.Actinreads.prop.t=lm((Actin.prop)~RainCom*RainLevel, data=Actin.clayE.t.reads)
qqPlot(studres(bact.clayE.Actinreads.prop.t))
hist(studres(bact.clayE.Actinreads.prop.t))
shapiro.test(studres(bact.clayE.Actinreads.prop.t))
#p-value = 0.1639
boxCox(bact.clayE.Actinreads.prop.t)
#looks okay....
Anova(bact.clayE.Actinreads.prop.t,type = 3)
#RainLevel         0.006387  1   56.5649 3.423e-08 ***
#RainCom:RainLevel 0.001648  1   14.5913 0.0006797 ***

lsmeans(bact.clayE.Actinreads.prop.t, pairwise~RainLevel*RainCom)
"$contrasts
 contrast                                estimate      SE df t.ratio p.value
 Ambient,NonSterile - Reduced,NonSterile  -0.0426 0.00531 28 -8.019  <.0001 
 Ambient,NonSterile - Ambient,Sterile     -0.0182 0.00531 28 -3.434  0.0095 
 Ambient,NonSterile - Reduced,Sterile     -0.0322 0.00531 28 -6.051  <.0001 
 Reduced,NonSterile - Ambient,Sterile      0.0244 0.00531 28  4.585  0.0005 
 Reduced,NonSterile - Reduced,Sterile      0.0105 0.00531 28  1.968  0.2240 
 Ambient,Sterile - Reduced,Sterile        -0.0139 0.00531 28 -2.617  0.0639 "



Actin.soil.t.readsT0=subset(Actin.soil.t.reads, SampleTime=="Start")
Actin.clayE.t.readsT0=subset(Actin.soil.t.readsT0, SoilType=="Clay")
mean(Actin.clayE.t.readsT0$Actin.prop)
positions <- c("NonSterile.Ambient.End", "NonSterile.Reduced.End", 
               "Sterile.Ambient.End", "Sterile.Reduced.End")

Actin.clayE.t.reads=Actin.clayE.t.reads %>% group_by(RainCom, RainLevel)
Actin.clayE.t.reads_sum=summarise_at(Actin.clayE.t.reads, vars(Actin.reads,Actin.prop),funs(mean,se=sd(.)/sqrt(n()),sd))

Actinprop_clay=ggplot(Actin.clayE.t.reads_sum, aes(RainLevel, Actin.prop_mean, ymin =Actin.prop_mean-Actin.prop_se,
                                                   ymax = Actin.prop_mean+Actin.prop_se))

(bact.clayE_Actin=Actinprop_clay+geom_errorbar(width=0.25)+geom_line(aes(group = RainCom))+
    geom_point(size=5,aes(shape=RainCom))+scale_shape_manual(values=c(16,1), name=NULL,
                                                             labels=c("Dispersal",
                                                                      "No Dispersal"))+
    scale_x_discrete(labels=c("Ambient","Drought"))+
    scale_y_continuous(name = "Proportion of reads Actinobacteria",labels=scaleFUN3)+xlab(NULL)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_hline(yintercept = mean(Actin.clayE.t.readsT0$Actin.prop), size=1.2))



#####Drought Tolerant Taxa Analyses####

