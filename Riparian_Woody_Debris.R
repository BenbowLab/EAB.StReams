# library
library(ggplot2)

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

#Upload Woody debris riparian count data
EAB_Rip_WD<-read.csv("~/Documents/MSU/Research/Surveying/Riparian/Riparian_Woody_Debris.csv", sep = ",", header = T )
#make percent ash add to 1
EAB_Rip_WD$Percent_ash<-EAB_Rip_WD$Percent_ash/100
#visualize differences between forest and gap overall
ggplot(EAB_Rip_WD, aes(x=Stream, y=Percent_ash, fill=Location)) + 
  geom_boxplot() +
  ylab("Percent Ash") +
  xlab("Watershed")+
  scale_colour_hue(name="Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))

#check for normality
shapiro.test(EAB_Rip_WD$Percent_ash)
#not normal
kruskal.test(Percent_ash ~ Stream, EAB_Rip_WD)
#Identical populations


#Paired t test
T_EAB_Rip<- cast(EAB_Rip_WD, Watershed ~ Location, value = "Percent_ash")
t.test(T_EAB_Rip$Forest,T_EAB_Rip$Gap,paired=T)


