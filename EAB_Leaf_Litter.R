
#Leaf litter analysis for EAB survey
#load libraries
library(plyr)
library(ggplot2)
library(reshape)
library(vegan)
library("RColorBrewer")
library(wesanderson)
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


EAB_LL<-read.csv("~/Documents/MSU/Research/Surveying/Leaf_Litter/EAB_Aq_Leaf_Litter.csv", sep=",", header = T)
#visualize with stacked bar plots
#Create new sample id label
EAB_LL$SampleID<-as.factor(paste(EAB_LL$Stream_name, EAB_LL$Date, EAB_LL$Gap_number, EAB_LL$Gap_location))
llg<- ddply(EAB_LL, "SampleID", transform, percent_count=Dry_mass_g/sum(Dry_mass_g)*100)
ggplot(llg, aes(x=Gap_location, y=percent_count, fill=Taxa)) + geom_bar(stat="identity")

#calculate alpha diversity values for each sample
EAB_LL_CM<-cast(EAB_LL, SampleID+Watershed+Stream_name+Gap_number+Gap_location+Date~Taxa, value="Dry_mass_g")
EAB_LL_CM[is.na(EAB_LL_CM)] <- 0
EAB_LL_com<-EAB_LL_CM[,7:ncol(EAB_LL_CM)]
EAB_LL_env<-EAB_LL_CM[,1:6]
adonis(EAB_LL_com ~ Gap_location, data=EAB_LL_env, method="bray", permutations=999)
#not significant
EAB_LL_CM$Alpha<-diversity (EAB_LL_com, index = "shannon")
EAB_LL_CM$Watershed<-factor(EAB_LL_CM$Watershed, levels=c("Kalamazoo", "Grand", "Clinton"))
EAB_LL_CM$Gap_location<-factor(EAB_LL_CM$Gap_location, levels=c("Upstream", "Gap", "Downstream"))

ggplot(EAB_LL_CM, aes(x=Gap_location, y=Alpha, colour=Watershed)) + 
  geom_boxplot() +
  ylab("Shannon Diversity") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))
anova(lmer(Alpha ~ Gap_location + (1|Gap_number/Stream_name), EAB_LL_CM))

