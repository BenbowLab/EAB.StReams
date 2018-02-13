#Analysis of woody debris survey

#Upload libraries
library(ggplot2)
library(nlme)
library(nlme)
library(lme4)
library(lmerTest)
#Upload data
EAB_ILWD<-read.csv("~/Documents/MSU/Research/Surveying/Stream_Woody_Debris.csv", sep=",", header = T)
#Create new variable for location
EAB_ILWD$LocationID<-factor(paste(EAB_ILWD$Stream, EAB_ILWD$Gap, EAB_ILWD$Basket))
#Create new variable for Area
EAB_ILWD$Area<-as.numeric(2*pi*((EAB_ILWD$Diameter)/2)*EAB_ILWD$In_Stream_Length+2*pi*(EAB_ILWD$Diameter/2)^2)
#Cast data into new data sheet with counts for each location
EAB_LWD<-aggregate(EAB_ILWD$Area, by=list(EAB_ILWD$LocationID,EAB_ILWD$Species,EAB_ILWD$Basket,EAB_ILWD$Stream), FUN=sum)
#Create Watershed variable
EAB_LWD$Watershed<-as.character(EAB_LWD$Group.4)
EAB_LWD$Watershed[EAB_LWD$Watershed == "Augusta Creek"]<-"Kalamazoo"
EAB_LWD$Watershed[EAB_LWD$Watershed == "Seven Mile Creek"]<-"Kalamazoo"
EAB_LWD$Watershed[EAB_LWD$Watershed == "Frayer Creek"]<-"Grand"
EAB_LWD$Watershed[EAB_LWD$Watershed == "Sessions Creek"]<-"Grand"
EAB_LWD$Watershed[EAB_LWD$Watershed == "Spring Creek"]<-"Clinton"
EAB_LWD$Watershed[EAB_LWD$Watershed == "Stoney Creek"]<-"Clinton"
#Create transect length variable
EAB_LWD$Length<-as.character(EAB_LWD$Group.1)
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 1 Upstream"]<-"32"
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 1 Gap"]<-"71"
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 1 Downstream"]<-"30"
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 2 Upstream"]<-"38"
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 2 Gap"]<-"114"
EAB_LWD$Length[EAB_LWD$Length == "Augusta Creek 2 Downstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 1 Upstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 1 Gap"]<-"33"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 1 Downstream"]<-"21"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 2 Upstream"]<-"29"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 2 Gap"]<-"31"
EAB_LWD$Length[EAB_LWD$Length == "Seven Mile Creek 2 Downstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 1 Upstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 1 Gap"]<-"104"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 1 Downstream"]<-"87"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 2 Upstream"]<-"36"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 2 Gap"]<-"29"
EAB_LWD$Length[EAB_LWD$Length == "Frayer Creek 2 Downstream"]<-"45"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 1 Upstream"]<-"28"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 1 Gap"]<-"26"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 1 Downstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 2 Upstream"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 2 Gap"]<-"20"
EAB_LWD$Length[EAB_LWD$Length == "Sessions Creek 2 Downstream"]<-"40"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 1 Upstream"]<-"54"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 1 Gap"]<-"202"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 1 Downstream"]<-"73"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 2 Upstream"]<-"36"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 2 Gap"]<-"34"
EAB_LWD$Length[EAB_LWD$Length == "Spring Creek 2 Downstream"]<-"58"
#Delete stoney creek
EAB_LWDS<-EAB_LWD[-c(52,53,54,55,56),]
#Create new dataset for total woody debris
EAB_TWD<-aggregate(EAB_LWDS$x, by=list(EAB_LWDS$Group.4,EAB_LWDS$Group.1,EAB_LWDS$Group.3,EAB_LWDS$Watershed, EAB_LWDS$Length), FUN=sum)
EAB_TWD$TWDAREA<-as.numeric(EAB_TWD$Group.5)
#Create new variable for total woody debris density
EAB_TWD$TWD_Density<-as.numeric(EAB_TWD$x/EAB_TWD$TWDAREA)
#order for watershed
EAB_TWD$Watershed<-factor(EAB_TWD$Group.4, levels=c("Kalamazoo", "Grand", "Clinton"))
#order for stream location
EAB_TWD$GapLocation<-factor(EAB_TWD$Group.3, levels=c("Upstream", "Gap", "Downstream"))
#Plot for total woody debris density
ggplot(EAB_TWD, aes(x=GapLocation, y=TWD_Density, fill=Watershed)) + 
  geom_boxplot() +
  ylab("TWD area per m in stream") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#Some basic stats
anova(lmer(TWD_Density ~ GapLocation + (1|Watershed), EAB_TWD))


################################
#Create new variable for ash woody debris density
EAB_LWD$AWD_Density<-as.numeric(EAB_LWD$Number_Ash_LWD/EAB_LWD$Transect_Length_m)
#Create new variable for ratio of ash to total woody debris
EAB_LWD$Ash_Ratio<-as.numeric(EAB_LWD$Number_Ash_LWD/EAB_LWD$Total_LWD)
#Create new variable for ash density percent
EAB_LWD$AshPer<-as.numeric(EAB_LWD$AWD_Density/EAB_LWD$TWD_Density)
#order for upstream, gap, downstream
str(EAB_LWD)
EAB_LWD$Gap_location<- factor(EAB_LWD$Gap_location, levels=c("Upstream", "Gap", "Downstream"))
#make gap number a factor
EAB_LWD$Gap_number<-as.factor(EAB_LWD$Gap_number)
#order for watershed
EAB_LWD$Watershed<-factor(EAB_LWD$Watershed, levels=c("Kalamazoo", "Grand", "Clinton"))

#Plots for percent ash
#Create box plots
ggplot(EAB_LWD, aes(x=Gap_location, y=Ash_Ratio, fill=Watershed)) + 
  geom_boxplot() +
  ylab("Percent Ash") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


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
shapiro.test(EAB_LWD$AshPer)
#We reject the null hypothesis and determine not normally distributed
#Transform arcsin function
EAB_LWD$arcsinAshPer<-asin(EAB_LWD$AshPer)
shapiro.test(EAB_LWD$arcsinAshPer)
#test for outliers
dotchart(EAB_LWD$arcsinAshPer)
#appears to be one outlier
plot(EAB_LWD$arcsinAsh_Ratio)
identify(EAB_LWD$arcsinAsh_Ratio)
#Create new data frame with outlier deleted (augusta creek gap 2 downstream)
EAB_LWD_NO<-EAB_LWD[-c(12),]

AR<-glmer(Ash_Ratio ~ Gap_location + Watershed + Gap_number + (1|Stream_name), family=binomial, data = EAB_LWD_NO)
summary(AR)
#no significance, possibly too few data points

drop1(AR)
#run model again without gap_number
AR1<-glmer(Ash_Ratio ~ Gap_location + Watershed + (1|Stream_name), family=binomial, data = EAB_LWD_NO)
summary(AR1)

drop1(AR1)

#try model with transformed data
ART<-lm(arcsinAsh_Ratio ~ Gap_location + Watershed, data = EAB_LWD_NO)
summary(ART)
alias(ART)
anova(lmer(arcsinAsh_Ratio ~ Gap_location + (1|Gap_number/Stream_name), EAB_LWD))

#Now look at decay levels
ggplot(EAB_LWD, aes(x=Gap_location, y=Average_Decay_Level, fill=Watershed)) + 
  geom_boxplot() +
  ylab("Decay Level") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

#Stats for decay level
#test for normality of distribution
shapiro.test(EAB_LWD$Average_Decay_Level)
#Not normal
anova(lmer(Average_Decay_Level ~ Gap_location + (1|Gap_number/Stream_name), EAB_LWD))

#analyze ash density
#Create line plots
anova(lmer(AWD_Density ~ Gap_location + (1|Gap_number/Stream_name), EAB_LWD))
ggplot(EAB_LWD, aes(x=Gap_location, y=AWD_Density, colour=Watershed)) + 
  geom_boxplot() +
  ylab("Ash Density") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))
