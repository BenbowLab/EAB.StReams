
#Load Packages
library(vegan)
library(indicspecies)

#Create dataset
EAB_16S_uni<-read.csv("~/Desktop/Core/EAB_exported-wunifrac/distance-matrix.tsv", sep="\t", header = T, row.names=1)
EAB_16S_map <- read.table("~/Desktop/Core/EAB_Map_Filtered_R.txt", header=T)
EAB_16S_map_r<-EAB_16S_map
row.names(EAB_16S_map_r)<-EAB_16S_map_r[,1]
EAB_16S_uni_map <-merge(EAB_16S_map_r, EAB_16S_uni, by=0)
row.names(EAB_16S_uni_map)<-EAB_16S_uni_map[,1]
EAB_16S_uni_map<-EAB_16S_uni_map[,-c(1)]
EAB_16S_uni<-EAB_16S_uni_map[15:ncol(EAB_16S_uni_map)]
EAB_16S_uni_com_samples<-as.vector(rownames(EAB_16S_uni))
#UNI-use output of names to subset columns into rows to make square matrix
EAB_16S_uni<-as.matrix(EAB_16S_uni)
EAB_16S_uni<-subset(EAB_16S_uni, select=c(EAB_16S_uni_com_samples))
str(EAB_16S_uni)
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env<-EAB_16S_uni_map
row.names(EAB_16S_uni_env)<-EAB_16S_uni_env[,1]
EAB_16S_uni_env<-EAB_16S_uni_env[,1:14]
EAB_16S_uni_env$Source<-as.factor(EAB_16S_uni_env$Source)
EAB_16S_uni_env$Stream<-as.factor(EAB_16S_uni_env$Stream)
EAB_16S_uni_env$Watershed<-as.factor(EAB_16S_uni_env$Watershed)
EAB_16S_uni_env$Time<-as.factor(EAB_16S_uni_env$Time)
EAB_16S_uni_env$Gap<-as.factor(EAB_16S_uni_env$Gap)
EAB_16S_uni_env$Location<-as.factor(EAB_16S_uni_env$Location)

#UNI-Overall permanova with unifrac distances
adonis(as.dist(EAB_16S_uni) ~ Source*Watershed*Gap*Stream+Time, data=EAB_16S_uni_env, permutations=999)

#upload phyla level info and run indicator analysis for gap/generate box plots
#Upload phyla level files for each run
EAB_16S_P<-read.table("~/Desktop/Core/EAB_otu_table_p.txt", sep="\t", header = T)
#Clasify as data.frame
EAB_16S_P<-data.frame(EAB_16S_P)
#Format data frame so the taxonomy is row name
row.names(EAB_16S_P)<-EAB_16S_P[,1]
#Delete taxonomy column
EAB_16S_P$OTUID<-NULL
#transpose
EAB_16S_P_t<-t(EAB_16S_P)
EAB_16S_P_t<-data.frame(EAB_16S_P_t)
str(EAB_16S_P_t)
names(EAB_16S_P_t)
#Merge metadata onto data table
EAB_16S_P_map <-merge(EAB_16S_map_r, EAB_16S_P_t, by=0)
EAB_16S_P_map[,1:15]<-sapply(EAB_16S_P_map[,1:15], as.factor)
EAB_16S_P_map_env<-EAB_16S_P_map[,1:15]
EAB_16S_P_map_com<-EAB_16S_P_map[,16:ncol(EAB_16S_P_map)]
EAB_p_indic<-signassoc(EAB_16S_P_map_com, cluster=EAB_16S_P_map_env$Gap,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_p_indic_sig<-subset(EAB_p_indic, psidak<=0.05)
#upload family level info and run indicator analysis for gap/generate box plots
#Upload family level files for each run
EAB_16S_F<-read.table("~/Desktop/Core/EAB_otu_table_f.txt", sep="\t", header = T)
#Clasify as data.frame
EAB_16S_F<-data.frame(EAB_16S_F)
#Format data frame so the taxonomy is row name
row.names(EAB_16S_F)<-EAB_16S_F[,1]
#Delete taxonomy column
EAB_16S_F$OTUID<-NULL
#transpose
EAB_16S_F_t<-t(EAB_16S_F)
EAB_16S_F_t<-data.frame(EAB_16S_F_t)
str(EAB_16S_F_t)
names(EAB_16S_F_t)
#Merge metadata onto data table
EAB_16S_F_map <-merge(EAB_16S_map_r, EAB_16S_F_t, by=0)
EAB_16S_F_map[,1:15]<-sapply(EAB_16S_F_map[,1:15], as.factor)
EAB_16S_F_map_env<-EAB_16S_F_map[,1:15]
EAB_16S_F_map_com<-EAB_16S_F_map[,16:ncol(EAB_16S_F_map)]
EAB_f_indic<-signassoc(EAB_16S_F_map_com, cluster=EAB_16S_F_map_env$Gap,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_f_indic_sig<-subset(EAB_f_indic, psidak<=0.05)
EAB_16S_F_map_agg<-aggregate(EAB_16S_F_map[16:ncol(EAB_16S_F_map)], by=list(Gap=EAB_16S_F_map$Gap, Source=EAB_16S_F_map$Source), FUN=sum)
EAB_16S_F_map_agg$k__Bacteria.p__Chloroflexi.c__Anaerolineae.o__Caldilineales.f__Caldilineaceae
#Found in gap terrestrial and aquatic litter, but not in non gap or live leaves in gap
EAB_16S_F_map_agg$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Aurantimonadaceae
EAB_16S_F_map_agg$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Methylobacteriaceae

#Upload phylogenetic diversity dataset
EAB_16S_fpd<-read.table("~/Desktop/Core/EAB_exported-faithpd/alpha-diversity.tsv", sep="\t", header = T, row.names=1)
EAB_16S_fpd_map <-merge(EAB_16S_map_r, EAB_16S_fpd, by=0)
#paste into paired sites datasheet
