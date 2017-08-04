#Analysis of woody debris survey

#Upload libraries
library(ggplot2)
#Upload data
EAB_LWD<-read.csv("~/Documents/MSU/Research/Surveying/CWD_Survey.csv", sep=",", header = T)
#Create new variable for location
EAB_LWD$LocationID<-factor(paste(EAB_LWD$Gap_number, EAB_LWD$Gap_location))
#specify order
EAB_LWD$LocationID <- factor(EAB_LWD$LocationID, c("1 Upstream", "1 Gap", "1 Downstream", "2 Upstream", "2 Gap", "2 Downstream"))
#Create new variable for total woody debris density
EAB_LWD$TWD_Density<-as.numeric(EAB_LWD$Total_LWD/EAB_LWD$Transect_Length_m)
#Create new variable for ash woody debris density
EAB_LWD$AWD_Density<-as.numeric(EAB_LWD$Number_Ash_LWD/EAB_LWD$Transect_Length_m)
#Create new variable for ratio of ash to total woody debris
EAB_LWD$Ash_Ratio<-as.numeric(EAB_LWD$Number_Ash_LWD/EAB_LWD$Total_LWD)

#Create heat maps
#Base plot
TWD<-ggplot(EAB_LWD, aes(x=Stream_name, y=LocationID, fill=TWD_Density))
#using geom_tile
TWD+
  geom_tile()+
  xlab("Stream name")+
  scale_x_discrete(limits=c("Augusta Creek", "Seven Mile Creek", "Sessions Creek", "Frayer Creek", "Spring Creek"))+
  ylab("Gap Location")+
  labs(fill="Woody Debris Density")+
  scale_fill_gradient(low="#f7fcf5", high="#00441b")

#Base plot
AWD<-ggplot(EAB_LWD, aes(x=Stream_name, y=LocationID, fill=AWD_Density))
#using geom_tile
AWD+
  geom_tile()+
  xlab("Stream name")+
  scale_x_discrete(limits=c("Augusta Creek", "Seven Mile Creek", "Sessions Creek", "Frayer Creek", "Spring Creek"))+
  ylab("Gap Location")+
  labs(fill="Ash WD Density")+
  scale_fill_gradient(low="#f7fcf5", high="#00441b")

#Base plot
DL<-ggplot(EAB_LWD, aes(x=Stream_name, y=LocationID, fill=Average_Decay_Level))
#using geom_tile
DL+
  geom_tile()+
  xlab("Stream name")+
  scale_x_discrete(limits=c("Augusta Creek", "Seven Mile Creek", "Sessions Creek", "Frayer Creek", "Spring Creek"))+
  ylab("Gap Location")+
  labs(fill="Average Decay Level")+
  scale_fill_gradient(low="#f7fcf5", high="#00441b")

