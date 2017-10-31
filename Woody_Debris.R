#Analysis of woody debris survey

#Upload libraries
library(ggplot2)
library(lme4)
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
#test different values for watersheds
kruskal.test(Ash_Ratio ~ Watershed, EAB_LWD)
#We reject the null hypothesis and determine each watershed has different ash ratios. 
#Therefore use each watershed for further tests
#subset dataset by each watershed
EAB_LWD_KZ<-subset(EAB_LWD, Watershed=="Kalamazoo")
EAB_LWD_GR<-subset(EAB_LWD, Watershed=="Grand")
EAB_LWD_CL<-subset(EAB_LWD, Watershed=="Clinton")
#Determine if differences for each stream in each watershed
shapiro.test(EAB_LWD_KZ$Ash_Ratio)
#Kalamazoo not normal
shapiro.test(EAB_LWD_GR$Ash_Ratio)
#Grand not normal
shapiro.test(EAB_LWD_CL$Ash_Ratio)
#Clinton normal
wilcox.test(Ash_Ratio ~ Stream_name, EAB_LWD_KZ)
#We fail to reject the null hypothesis, do not have different ratios
wilcox.test(Ash_Ratio ~ Stream_name, EAB_LWD_GR)
#We fail to reject the null hypothesis, do not have different ratios
#Only ony stream in Clinton watershed
#We do not need to run different tests for each stream
EAB_LWD_KZ$Block<-factor(paste(EAB_LWD_KZ$Stream_name,EAB_LWD_KZ$Gap_number))
friedman.test(EAB_LWD_KZ$Ash_Ratio, EAB_LWD_KZ$Gap_location, EAB_LWD_KZ$Block)
#Not significant
EAB_LWD_GR$Block<-factor(paste(EAB_LWD_GR$Stream_name,EAB_LWD_GR$Gap_number))
friedman.test(EAB_LWD_GR$Ash_Ratio, EAB_LWD_GR$Gap_location, EAB_LWD_GR$Block)
#Not significant
EAB_LWD_CL$Block<-factor(paste(EAB_LWD_CL$Stream_name,EAB_LWD_CL$Gap_number))
anova(lm(Ash_Ratio ~ Gap_location + Block, EAB_LWD_CL))
#Not significant

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
#test different values for watersheds
kruskal.test(Average_Decay_Level ~ Watershed, EAB_LWD)
#not identical, so split up by watershed
#Determine if differences for each stream in each watershed
shapiro.test(EAB_LWD_KZ$Average_Decay_Level)
#Kalamazoo normal
shapiro.test(EAB_LWD_GR$Average_Decay_Level)
#Grand not normal
shapiro.test(EAB_LWD_CL$Average_Decay_Level)
#Clinton normal
t.test(EAB_LWD_KZ$Average_Decay_Level ~ EAB_LWD_KZ$Stream_name)
#no difference in streams
wilcox.test(Average_Decay_Level ~ Stream_name, EAB_LWD_GR)
#no difference in streams
#clinton only one stream
#Continue analyses within each watershed
friedman.test(EAB_LWD_KZ$Average_Decay_Level, EAB_LWD_KZ$Gap_location, EAB_LWD_KZ$Block)
#Not significant
friedman.test(EAB_LWD_GR$Average_Decay_Level, EAB_LWD_GR$Gap_location, EAB_LWD_GR$Block)
#Not significant
friedman.test(EAB_LWD_CL$Average_Decay_Level, EAB_LWD_CL$Gap_location, EAB_LWD_CL$Block)
#Not significant

