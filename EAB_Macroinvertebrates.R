#Analyze macroinvertebrate community data from EAB MI survey

###################################################
#Load packages and functions and get color vectors
##################################################

#load packages
library(vegan)
library(reshape)

#Functions
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Make color vectors

##############################################
#Upload data and make "r friendly"
###############################################

#Upload .csv file
EAB_Macroinvertebrates<-read.csv("~/Documents/MSU/Research/Surveying/Macroinvertebrates/EAB_Survey_Macroinvertebrates.csv", sep = ",", header = T )
#Confirm header names
names(EAB_Macroinvertebrates)
#Create sample ID name by combining environmental variables into one name in new column
EAB_Macroinvertebrates$SampleID<-factor(paste(EAB_Macroinvertebrates$Stream_Name, EAB_Macroinvertebrates$Gap, EAB_Macroinvertebrates$Basket, EAB_Macroinvertebrates$Date, EAB_Macroinvertebrates$Sample_Type))
#Create new column that combines taxonomic variables to family level
EAB_Macroinvertebrates$Taxonomy<-(paste(EAB_Macroinvertebrates$Class, EAB_Macroinvertebrates$Order, EAB_Macroinvertebrates$Family))
#Delete pupae (named Hexapod NA NA in taxonomy column)
EAB_Inverts_Subset<-subset(EAB_Macroinvertebrates, Count >= 0 & Taxonomy != "Insecta NA NA" & Taxonomy != "Insecta Diptera NA" & Taxonomy != "Insecta Trichoptera NA" & Taxonomy != "Insecta Odonata NA")
#Convert Taxonomy charater to factor
EAB_Inverts_Subset$Taxonomy<-as.factor(EAB_Inverts_Subset$Taxonomy)
#Check levels of SampleID and Taxonomy variables
str(EAB_Inverts_Subset)
levels(EAB_Inverts_Subset$SampleID)
levels(EAB_Inverts_Subset$Taxonomy)

#Convert to community table for community analysis
EAB_M_Matrix <- cast(EAB_Inverts_Subset, Watershed + Stream_Name + Date + Gap + Basket ~ Taxonomy, value = "Count", fun.aggregate =sum)
#create file with only community data
str(EAB_M_Matrix)
EAB_Inverts_Community <- EAB_M_Matrix[,6:ncol(EAB_M_Matrix)]
#Replace "na" with 0
EAB_Inverts_Community[is.na(EAB_Inverts_Community)]<-0
str(EAB_Inverts_Community)
#Separate environmental data
EAB_Inverts_Env<-EAB_M_Matrix[,1:5]
str(EAB_Inverts_Env)

#Permanova for Stream, Date, Basket
adonis(EAB_Inverts_Community ~ Date*Basket*Watershed, data=EAB_Inverts_Env, method="bray", permutations=999)
#Date and watershed significant