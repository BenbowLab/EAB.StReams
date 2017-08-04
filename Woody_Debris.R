#Analysis of woody debris survey

#Upload libraries
library(ggplot2)
#Upload data
EAB_LWD<-read.csv("~/Documents/MSU/Research/Surveying/CWD_Survey.csv", sep=",", header = T)
#Create new variable for location
EAB_LWD$LocationID<-factor(paste(EAB_LWD$Gap_number, EAB_LWD$Gap_location))
#specify order
EAB_LWD$LocationID <- factor(EAB_LWD$LocationID, c("1 Upstream", "1 At", "1 Downstream", "2 Upstream", "2 At", "2 Downstream"))
#For total large woody debris
#Base plot
p<-ggplot(EAB_LWD, aes(x=Stream_name, y=LocationID, fill=Total_LWD))
#using geom_tile
p+geom_tile()
#For ratio of large woody debris
EAB_LWD$Ash_Ratio<-as.numeric(EAB_LWD$Number_Ash_LWD/EAB_LWD$Total_LWD)
#Base plot
r<-ggplot(EAB_LWD, aes(x=Stream_name, y=LocationID, fill=Ash_Ratio))
#using geom_tile
r+geom_tile()