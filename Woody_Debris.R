#Analysis of woody debris survey

#Upload libraries
library(ggplot2)
library(nlme)
library(nlme)
library(lme4)
library(lmerTest)
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
#order for upstream, gap, downstream
str(EAB_LWD)
EAB_LWD$Gap_location<- factor(EAB_LWD$Gap_location, levels=c("Upstream", "Gap", "Downstream"))
#make gap number a factor
EAB_LWD$Gap_number<-as.factor(EAB_LWD$Gap_number)

#Plots for percent ash
#Create line plots
ggplot(EAB_LWD, aes(x=Watershed, y=Ash_Ratio, colour=Gap_location)) + 
  geom_boxplot() +
  ylab("Percent Ash") +
  scale_colour_hue(name="Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))

#Create heat maps
#using geom_tile
#Base plot
AWD<-ggplot(EAB_LWD, aes(x=Watershed, y=Gap_location, fill=Ash_Ratio))
AWD+
  geom_tile()+
  xlab("Watershed")+
  scale_x_discrete(limits=c("Kalamazoo", "Grand", "Clinton"))+
  ylab("Gap Location")+
  labs(fill="Percent Ash
Woody Debris")+
  scale_fill_gradient(low="#f7fcf5", high="#00441b")+
  theme_classic()+
  theme(axis.title = element_text(size = rel(2)))+
  theme(axis.text = element_text(size = rel(1.5)))+
  theme(legend.title = element_text(size=rel(1.5)))+
  theme(legend.text = element_text(size=rel(1)))

#Stats for percent ash in stream

#test for normality of distribution
shapiro.test(EAB_LWD$Ash_Ratio)
#We reject the null hypothesis and determine not normally distributed
anova(lmer(SR_Ash_Ratio ~ Gap_location + (1|Gap_number/Stream_name), EAB_LWD))

#Now look at decay levels
ggplot(EAB_LWD, aes(x=Watershed, y=Average_Decay_Level, colour=Gap_location)) + 
  geom_boxplot() +
  ylab("Decay Level") +
  scale_colour_hue(name="Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))

#Stats for decay level
#test for normality of distribution
shapiro.test(EAB_LWD$Average_Decay_Level)
#Not normal
anova(lmer(Average_Decay_Level ~ Gap_location + (1|Gap_number/Stream_name), EAB_LWD))
