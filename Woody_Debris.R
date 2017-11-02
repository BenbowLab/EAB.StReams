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
#Create new variable for ash density percent
EAB_LWD$AshPer<-as.numeric(EAB_LWD$AWD_Density/EAB_LWD$TWD_Density)
#order for upstream, gap, downstream
str(EAB_LWD)
EAB_LWD$Gap_location<- factor(EAB_LWD$Gap_location, levels=c("Upstream", "Gap", "Downstream"))
#make gap number a factor
EAB_LWD$Gap_number<-as.factor(EAB_LWD$Gap_number)

#Plots for percent ash
#Create line plots
ggplot(EAB_LWD, aes(x=Gap_location, y=Ash_Ratio, colour=Watershed)) + 
  geom_boxplot() +
  ylab("Percent Ash") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#Create line plots
ggplot(EAB_LWD, aes(x=Gap_location, y=Ash_Ratio)) + 
  geom_boxplot() +
  ylab("Percent Ash") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14))

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
ggplot(EAB_LWD, aes(x=Gap_location, y=Average_Decay_Level, colour=Watershed)) + 
  geom_boxplot() +
  ylab("Decay Level") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

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

