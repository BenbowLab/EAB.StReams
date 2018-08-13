########################
#Analysis of paired sites
########################

#load packages
library(ggplot2)
library(rjags)
library(coda)
library(runjags)
library(MASS)
library(MCMCpack)
library(reshape)
library(ggpubr)
library(nlme)
library(nlme)
library(lme4)
library(lmerTest)
library(vegan)
#Functions

#Color vectors
sampling_time_col_vec<-c("#a6cee3", "#1f78b4", "#b2df8a")

#########################
#Exploratory data analysis
#######################

#Upload data set
EAB_Paired<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Sites_EAB_Survey.csv", sep = ",", header = T )

#Create variables that are standardized by transect length
#First calculate transect area (units m^2)
EAB_Paired$Terrestrial_Transect_Area_m2<-EAB_Paired$Terrestrial_Transect_length*EAB_Paired$Terrestrial_Transect_width
EAB_Paired$Aquatic_Transect_Area_m2<-EAB_Paired$Aquatic_Transect_length*EAB_Paired$Aquatic_Transect_width

#Use transect area to standardize values observed for each transect (per m2)

EAB_Paired_Stand<-EAB_Paired

EAB_Paired_Stand$Terrestrial_Total_CWD_Density_m<-EAB_Paired$Terrestrial_Total_CWD_Volume_m3/EAB_Paired$Terrestrial_Transect_Area_m2
EAB_Paired_Stand$Terrestrial_Ash_CWD_Density_m<-EAB_Paired$Terrestrial_Ash_CWD_V_m3/EAB_Paired$Terrestrial_Transect_Area_m2
EAB_Paired_Stand$Terrestrial_1_CWD_Density_m<-EAB_Paired$Terrestrial_1_CWD_V_m3/EAB_Paired$Terrestrial_Transect_Area_m2
EAB_Paired_Stand$Aquatic_Total_CWD_Density_m<-EAB_Paired$Aquatic_Total_CWD_Volume_m3/EAB_Paired$Aquatic_Transect_Area_m2
EAB_Paired_Stand$Aquatic_Ash_CWD_Density_m<-EAB_Paired$Aquatic_Ash_CWD_V_m3/EAB_Paired$Aquatic_Transect_Area_m2
EAB_Paired_Stand$Aquatic_1_CWD_Density_m<-EAB_Paired$Aquatic_1_CWD_V_m3/EAB_Paired$Aquatic_Transect_Area_m2

#Delete variables that aren't of interest
EAB_Paired_Stand$Terrestrial_Transect_length<-NULL
EAB_Paired_Stand$Terrestrial_Transect_width<-NULL
EAB_Paired_Stand$Terrestrial_Transect_Area_m2<-NULL
EAB_Paired_Stand$Terrestrial_Total_CWD_Volume_m3<-NULL
EAB_Paired_Stand$Terrestrial_Ash_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_Maple_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_Beech_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_Elm_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_1_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_2_CWD_V_m3<-NULL
EAB_Paired_Stand$Terrestrial_3_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_Transect_length<-NULL
EAB_Paired_Stand$Aquatic_Transect_width<-NULL
EAB_Paired_Stand$Aquatic_Transect_Area_m2<-NULL
EAB_Paired_Stand$Aquatic_Total_CWD_Volume_m3<-NULL
EAB_Paired_Stand$Aquatic_Ash_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_Maple_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_Beech_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_1_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_2_CWD_V_m3<-NULL
EAB_Paired_Stand$Aquatic_3_CWD_V_m3<-NULL

names(EAB_Paired_Stand)
#Look for trends in the data using plots and correlations
pairs(EAB_Paired_Stand)

#WATERSHED
#Gap Location
#ALL Richness Early

#ALL prop unknown early
ggplot(EAB_Paired_Stand, aes(x=ALL_Richness_Early, y=ALL_prop_Unknown_Early, size=Watershed, color=Gap_location)) +
  geom_point() +
  theme_classic()

##########################
#Modelling and figures for ESA 2018 presentation
#########################
EAB_Paired_Model<-EAB_Paired_Stand
#Upload data set that uses each gap location as a different variable, rather than a single category
EAB_Paired_Exp<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Sites_EAB_Survey_Expanded.csv", sep = ",", header = T )
#Make area variable again
EAB_Paired_Exp$Ter_Woods_Transect_Area_m2<-EAB_Paired_Exp$Terrestrial_Woods_Transect_length*EAB_Paired_Exp$Terrestrial_Transect_width
EAB_Paired_Exp$Ter_Gap_Transect_Area_m2<-EAB_Paired_Exp$Terrestrial_Gap_Transect_length*EAB_Paired_Exp$Terrestrial_Transect_width
EAB_Paired_Exp$Aq_Transect_Area_m2_US<-EAB_Paired_Exp$Aq_Transect_length_US*EAB_Paired_Exp$Aq_Transect_width
EAB_Paired_Exp$Aq_Transect_Area_m2_G<-EAB_Paired_Exp$Aq_Transect_length_G*EAB_Paired_Exp$Aq_Transect_width
EAB_Paired_Exp$Aq_Transect_Area_m2_DS<-EAB_Paired_Exp$Aq_Transect_length_DS*EAB_Paired_Exp$Aq_Transect_width

#standardize woody debris values by area surveyed
EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US<-EAB_Paired_Exp$Ter_Total_CWD_V_m3_US/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2
EAB_Paired_Exp$Ter_Tot_CWD_Density_m_G<-EAB_Paired_Exp$Ter_Total_CWD_V_m3_G/EAB_Paired_Exp$Ter_Gap_Transect_Area_m2
EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS<-EAB_Paired_Exp$Ter_Total_CWD_V_m3_DS/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2

EAB_Paired_Exp$Ter_Ash_CWD_Density_m_US<-EAB_Paired_Exp$Ter_Ash_CWD_V_m3_US/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2
EAB_Paired_Exp$Ter_Ash_CWD_Density_m_G<-EAB_Paired_Exp$Ter_Ash_CWD_V_m3_G/EAB_Paired_Exp$Ter_Gap_Transect_Area_m2
EAB_Paired_Exp$Ter_Ash_CWD_Density_m_DS<-EAB_Paired_Exp$Ter_Ash_CWD_V_m3_DS/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2

EAB_Paired_Exp$Ter_1_CWD_Density_m_US<-EAB_Paired_Exp$Ter_1_CWD_V_m3_US/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2
EAB_Paired_Exp$Ter_1_CWD_Density_m_G<-EAB_Paired_Exp$Ter_1_CWD_V_m3_G/EAB_Paired_Exp$Ter_Gap_Transect_Area_m2
EAB_Paired_Exp$Ter_1_CWD_Density_m_DS<-EAB_Paired_Exp$Ter_1_CWD_V_m3_DS/EAB_Paired_Exp$Ter_Woods_Transect_Area_m2

EAB_Paired_Exp$Aq_Tot_CWD_Density_m_US<-EAB_Paired_Exp$Aq_Tot_CWD_V_m3_US/EAB_Paired_Exp$Aq_Transect_Area_m2_US
EAB_Paired_Exp$Aq_Tot_CWD_Density_m_G<-EAB_Paired_Exp$Aq_Tot_CWD_V_m3_G/EAB_Paired_Exp$Aq_Transect_Area_m2_G
EAB_Paired_Exp$Aq_Tot_CWD_Density_m_DS<-EAB_Paired_Exp$Aq_Tot_CWD_V_m3_DS/EAB_Paired_Exp$Aq_Transect_Area_m2_DS

EAB_Paired_Exp$Aq_Ash_CWD_Density_m_US<-EAB_Paired_Exp$Aq_Ash_CWD_V_m3_US/EAB_Paired_Exp$Aq_Transect_Area_m2_US
EAB_Paired_Exp$Aq_Ash_CWD_Density_m_G<-EAB_Paired_Exp$Aq_Ash_CWD_V_m3_G/EAB_Paired_Exp$Aq_Transect_Area_m2_G
EAB_Paired_Exp$Aq_Ash_CWD_Density_m_DS<-EAB_Paired_Exp$Aq_Ash_CWD_V_m3_DS/EAB_Paired_Exp$Aq_Transect_Area_m2_DS

EAB_Paired_Exp$Aq_1_CWD_Density_m_US<-EAB_Paired_Exp$Aq_1_CWD_V_m3_US/EAB_Paired_Exp$Aq_Transect_Area_m2_US
EAB_Paired_Exp$Aq_1_CWD_Density_m_G<-EAB_Paired_Exp$Aq_1_CWD_V_m3_G/EAB_Paired_Exp$Aq_Transect_Area_m2_G
EAB_Paired_Exp$Aq_1_CWD_Density_m_DS<-EAB_Paired_Exp$Aq_1_CWD_V_m3_DS/EAB_Paired_Exp$Aq_Transect_Area_m2_DS

#First let's see if there are differences in the riparian factors based on gap location

#Woody debris
#Because there are two zero values (out of 18), and the response variable is continuous, transform Terrestrial_Total_CWD_Density_m variable so that these zeros are very small positive values
#zeros will be 0.0003 (half the minimum value)
EAB_Paired_Model$Ter_To_CWD_D_m_nz<-EAB_Paired_Model$Terrestrial_Total_CWD_Density_m
EAB_Paired_Model <- within(EAB_Paired_Model, Ter_To_CWD_D_m_nz[Terrestrial_Total_CWD_Density_m==0] <- 0.0003)
#check for normality
shapiro.test(EAB_Paired_Model$Ter_To_CWD_D_m_nz)
#W=0.8,p=0.001; not normally distributed
#Check for outliers
ggdensity(EAB_Paired_Model$Ter_To_CWD_D_m_nz)
#Skewed right
ggqqplot(EAB_Paired_Model$Ter_To_CWD_D_m_nz)
#clear non normal distribution, possible outlier
plot(EAB_Paired_Model$Ter_To_CWD_D_m_nz)
#due to small sample size, use tests that do not assume normality
#use a gamma distribution with log link function
TTCWD<-glmer(Ter_To_CWD_D_m_nz ~ Gap_location*Watershed + (1|Stream), family=Gamma(link=log), data = EAB_Paired_Model)
summary(TTCWD)
#Significant factors: intercept, gap location, watershed, gap*watershed location. AIC=-154.2, BIC=-144.4
#Visualize dataset
#create new variable *1000 for easier visualization (volume CWD per hectare)
EAB_Paired_Model$Ter_To_CWD_D_m_nz_hect<-EAB_Paired_Model$Ter_To_CWD_D_m_nz*1000
EAB_Paired_Model$Gap_location <- factor(EAB_Paired_Model$Gap_location, levels = c("Upstream","Gap","Downstream"))
ggplot(EAB_Paired_Model, aes(x=Gap_location, y=Ter_To_CWD_D_m_nz_hect, fill=Watershed)) + 
  geom_boxplot() +
  ylab("Terrestrial coarse wood (V/ha)") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=22,margin=margin(50,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=16),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16), legend.position="bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Time of EAB Invasion", labels=c("Early", "Mid", "Late"))

#Riparian vegetation richness (genera level)
#check for normality
shapiro.test(EAB_Paired_Model$Terrestrial_Transect_Live_Richness)
#W=0.9,p=0.09; normally distributed
#Check for outliers
ggdensity(EAB_Paired_Model$Terrestrial_Transect_Live_Richness)
#bimodal distribution
ggqqplot(EAB_Paired_Model$Terrestrial_Transect_Live_Richness)
#clear non normal distribution, no outliers
plot(EAB_Paired_Model$Terrestrial_Transect_Live_Richness)
#Use a poisson distribution
TR<-glmer(Terrestrial_Transect_Live_Richness ~ Gap_location+Watershed + (1|Stream), family=poisson, data = EAB_Paired_Model)
summary(TR)
#Significant factors: intercept. AIC=69.8, BIC=75.1
#Use a normal distribution
TRG<-lmer(Terrestrial_Transect_Live_Richness ~ Gap_location+Watershed + (1|Stream), data = EAB_Paired_Model)
summary(TRG)
AIC(TRG)
BIC(TRG)
#Significant factors: intercept, gap location and watershed. AIC=56.18, BIC=62.41. better fit
#Visualize dataset
ggplot(EAB_Paired_Model, aes(x=Gap_location, y=Terrestrial_Transect_Live_Richness, fill=Watershed)) + 
  geom_boxplot() +
  ylab("Riparian vegetation richness") +
  xlab("Gap Location")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20,margin=margin(50,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=16),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Time of\nEAB Invasion", labels=c("Early", "Mid", "Late"))

#Now move onto stream habitat conditions

#Because there is 1 zero value (out of 6), and the response variable is continuous, transform Ter_Total_CWD_V_m3_US variable so that these zeros are very small positive values
#zeros will be 0.0003 (half the minimum value for all terr cwd)
EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US_nz<-EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US
EAB_Paired_Exp <- within(EAB_Paired_Exp, Ter_Tot_CWD_Density_m_US_nz[Ter_Tot_CWD_Density_m_US==0] <- 0.0003)
#Woody debris US
#check for normality for terrestrial woody debris factor US
shapiro.test(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US_nz)
#W=0.7,p=0.02; not normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US_nz)
#bimodal distribution
ggqqplot(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US_nz)
#clear non normal distribution, Possible outlier
plot(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_US_nz)
#Use a gamma distribution
#check for normality for aq woody debris factor US
shapiro.test(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_US)
#W=0.95,p=0.77; normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_US)
#bimodal distribution
ggqqplot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_US)
#clear non normal distribution, Possible outlier
plot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_US)
#Use a normal distribution
ATCWDUS<-glm(Aq_Tot_CWD_Density_m_US ~ Ter_Tot_CWD_Density_m_US_nz+Watershed, data = EAB_Paired_Exp)
summary(ATCWDUS)
#Significant factor: intercept AIC=-47.796
#try again without terrestrial
ATCWDUSNT<-glm(Aq_Tot_CWD_Density_m_US ~ Watershed, data = EAB_Paired_Exp)
summary(ATCWDUSNT)
#Significant factor: intercept AIC=-45.629 - worse model
#create new variable *1000 for easier visualization (volume CWD per hectare)
EAB_Paired_Model$Aq_To_CWD_D_m_hect<-EAB_Paired_Model$Aquatic_Total_CWD_Density_m*1000
ggplot(EAB_Paired_Model, aes(x=Ter_To_CWD_D_m_nz_hect, y=Aq_To_CWD_D_m_hect)) +
  geom_point(size=3, se=FALSE) +
  geom_smooth(method=lm)+
  xlab("Terrestrial Coarse Woody Material (V/ha)")+
  ylab("Aquatic Coarse Woody Material (V/ha)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=22), axis.text.x=element_text(size=16),
        axis.title.y = element_text(size=22), axis.text.y=element_text(size=16))
#Woody debris G
#check for normality for terrestrial woody debris factor G
shapiro.test(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_G)
#W=0.7, p=0.01, not normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_G)
#bimodal distribution
plot(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_G)
#gamma distributed
#check for normality for aq woody debris factor G
shapiro.test(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_G)
#W=0.9, p=0.2; normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_G)
#right skewed distribution
ggqqplot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_G)
#clear non normal distribution, Possible outlier
plot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_G)
#Use a gamma distribution
ATCWDG<-glm(Aq_Tot_CWD_Density_m_G ~ Watershed+Aq_Tot_CWD_Density_m_US+Ter_Tot_CWD_Density_m_G, family = Gamma(link="log"), maxit = 100, data = EAB_Paired_Exp)
summary(ATCWDG)
#No significant factors, AIC=-30.47
#Try again without terrestrial
ATCWDGNT<-glm(Aq_Tot_CWD_Density_m_G ~ Watershed+Aq_Tot_CWD_Density_m_US, family = Gamma(link="log"), maxit = 100, data = EAB_Paired_Exp)
summary(ATCWDGNT)
#No significant factors, AIC=-32.34 - better model
#Try again without aquatic
ATCWDGNTNA<-glm(Aq_Tot_CWD_Density_m_G ~ Watershed, family = Gamma(link="log"), maxit = 100, data = EAB_Paired_Exp)
summary(ATCWDGNTNA)
#Intercept significant factor, AIC=-33.756 - better model

#Woody debris DS
EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS_nz<-EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS
EAB_Paired_Exp <- within(EAB_Paired_Exp, Ter_Tot_CWD_Density_m_DS_nz[Ter_Tot_CWD_Density_m_DS==0] <- 0.0003)
#check for normality for terrestrial woody debris factor G
shapiro.test(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS_nz)
#W=0.84, p=0.12, normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS_nz)
#bimodal right skewed distribution
plot(EAB_Paired_Exp$Ter_Tot_CWD_Density_m_DS_nz)
#gamma distributed
#check for normality for aq woody debris factor G
shapiro.test(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_DS)
#W=0.9, p=0.3; normally distributed
#Check for outliers
ggdensity(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_DS)
#right skewed distribution
ggqqplot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_DS)
#normal distribution
plot(EAB_Paired_Exp$Aq_Tot_CWD_Density_m_DS)
#Use a normal distribution
ATCWDDS<-glm(Aq_Tot_CWD_Density_m_DS ~ Watershed+Aq_Tot_CWD_Density_m_US+Ter_Tot_CWD_Density_m_DS_nz, data = EAB_Paired_Exp)
summary(ATCWDDS)
#No significant factors, AIC=-48.647
#Try again without terrestrial
ATCWDDSNT<-glm(Aq_Tot_CWD_Density_m_DS ~ Watershed+Aq_Tot_CWD_Density_m_US, data = EAB_Paired_Exp)
summary(ATCWDDSNT)
#No significant factors, AIC=-47.61 = worse fit model
#Try again without aquatic
ATCWDDSNA<-glm(Aq_Tot_CWD_Density_m_DS ~ Watershed+Ter_Tot_CWD_Density_m_DS_nz, data = EAB_Paired_Exp)
summary(ATCWDDSNA)
#No significant factors, AIC=-47.455 worse fit
#Try again without terrestrial or aquatic
ATCWDDSNTNA<-glm(Aq_Tot_CWD_Density_m_DS ~ Watershed, data = EAB_Paired_Exp)
summary(ATCWDDSNTNA)
#Watershed significant factor, and AIC=-48.647, better model (slightly)

#do specific types of woody debris later, now moving on to leaf litter

#US leaf litter richness
#Import expanded dataset that does not account for time (each leaf litter sample is a replicate)
EAB_Paired_Exp_NT<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Sites_EAB_Survey_Expanded_NTime.csv", sep = ",", header = T )
#Import dataset that does not account for time or gap location
EAB_Paired_NT<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Sites_EAB_Survey_NTime.csv", sep = ",", header = T )
#Model aquatic leaf litter richness
#check for normality for aquatic leaf litter richness
shapiro.test(EAB_Paired_NT$ALL_Rich)
#W=0.94, p=0.01, not normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT$ALL_Rich)
#Looks evenly distributed
plot(EAB_Paired_NT$ALL_Rich)
#Poisson distributed
ALLR<-glmer(ALL_Rich ~ Ter_Transect_Live_Rich+Watershed+Gap_Location+Sampling_Time+(1|Stream), family=poisson, data = EAB_Paired_NT)
summary(ALLR)
#Intercept significant, AIC=214.5, BIC=240.3

#visualize data
ggplot(EAB_Paired_NT, aes(x=Ter_Transect_Live_Rich, y=ALL_Rich)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se=FALSE)+
  xlab("Terrestrial vegetation richness")+
  ylab("Aquatic leaf litter richness")+
  theme_classic()+
  theme(axis.title.x = element_text(size=22), axis.text.x=element_text(size=16),
        axis.title.y = element_text(size=22), axis.text.y=element_text(size=16))

#Proportion unknown
#check for normality for aquatic leaf litter proportion unknown
shapiro.test(EAB_Paired_NT$ALL_prop_Uk)
#W=0.95, p=0.04, not normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT$ALL_prop_Uk)
#Slightly right skewed
plot(EAB_Paired_NT$ALL_prop_Uk)
#Use gamma distribution
#Start with all potential factors
ALLU<-glmer(ALL_prop_Uk ~ Watershed+Sampling_Time+Gap_Location+(1|Stream), family=Gamma(link="log"), data = EAB_Paired_NT)
summary(ALLU)
#Gap location and sampling time significant, intercept and watershed approaching significance, AIC=-7.7 BIC=10.2
#Try without watershed
ALLUNW<-glmer(ALL_prop_Uk ~ Sampling_Time+Gap_Location+(1|Stream), family=Gamma(link="log"), data = EAB_Paired_NT)
summary(ALLUNW)
#Gap location, sampling time and intercept significant, AIC=-9.3, BIC=4.7, better model

#Visualize dataset
EAB_Paired_NT$Sampling_Time <- factor(EAB_Paired_NT$Sampling_Time, levels = c("Early","Mid","Late"))
EAB_Paired_NT$Gap_Location <- factor(EAB_Paired_NT$Gap_Location, levels = c("Upstream","Gap","Downstream"))
ggplot(EAB_Paired_NT, aes(x=Gap_Location, y=ALL_prop_Uk, fill=Sampling_Time)) +
  geom_boxplot() +
  xlab("Gap Location")+
  ylab("Aquatic litter mass decomposed")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20,margin=margin(50,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=16),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),legend.position="bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_fill_manual(values=sampling_time_col_vec, name="Sampling Time", labels=c("Pre Fall", "Leaf Fall","Post Fall"))

#Proportion Oak
#Transform data so that there are no zeros, just very small postive numbers
EAB_Paired_NT$ALL_prop_Oak_nz<-EAB_Paired_NT$ALL_prop_Oak
EAB_Paired_NT <- within(EAB_Paired_NT, ALL_prop_Oak_nz[ALL_prop_Oak_nz==0] <- 0.025) #half of min pos value
#check for normality for aquatic leaf litter proportion Oak
shapiro.test(EAB_Paired_NT$ALL_prop_Oak_nz)
#W=0.9, p=0.0001, not normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT$ALL_prop_Oak_nz)
#Very right skewed
plot(EAB_Paired_NT$ALL_prop_Oak_nz)
#Use gamma distribution
#Start with all potential factors
ALLO<-glmer(ALL_prop_Oak_nz ~ Watershed*Sampling_Time+Gap_Location+(1|Stream), family=Gamma(link="log"), data = EAB_Paired_NT)
summary(ALLO)
#Intercept, watershed and sampling time significant, AIC=-23.5 BIC=2.4
#Try without Gap_Location
ALLONG<-glmer(ALL_prop_Oak_nz ~ Watershed*Sampling_Time+(1|Stream), family=Gamma(link="log"), data = EAB_Paired_NT)
summary(ALLONG)
#Intercept, watershed and sampling time interaction significant, AIC=-26.7, BIC=-4.8

#Visualize dataset
ggplot(EAB_Paired_NT, aes(x=Sampling_Time, y=ALL_prop_Elm, fill=Gap_Location)) +
  geom_boxplot() +
  xlab("Sampling Time")+
  ylab("Aquatic Proportion Elm")+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

#Do other leaf litter variables later

#Now move onto microbes
#check for normality for all microbe shannon diversity
shapiro.test(EAB_Paired_NT$ALL_16S_Sh_Div)
#W=0.96, p=0.5, normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT$ALL_16S_Sh_Div)
#slight skewness
plot(EAB_Paired_NT$ALL_16S_Sh_Div)
#Use normal distribution
#Start with simple model
ALL16<-lmer(ALL_16S_Sh_Div ~ Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALL16)
AIC(ALL16)
BIC(ALL16)
#Intercept significant, with gap location approaching significance, AIC=64.4 BIC=70
#Add Watershed
ALL16W<-lmer(ALL_16S_Sh_Div ~ Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16W)
AIC(ALL16W)
BIC(ALL16W)
#Intercept and gap location significant, with watershed approaching significance, AIC=62.2 BIC=73.5 better model
#Add ALL Richness
ALL16WR<-lmer(ALL_16S_Sh_Div ~ ALL_Rich+Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WR)
AIC(ALL16WR)
BIC(ALL16WR)
#Intercept and gap location significant, with watershed approaching significance, AIC=66.2 BIC=78.7 worse model
#Add Sampling Time
ALL16WT<-lmer(ALL_16S_Sh_Div ~ Sampling_Time+Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WT)
AIC(ALL16WT)
BIC(ALL16WT)
#Intercept and gap location significant, with watershed approaching significance, AIC=64.8 BIC=78.5 worse model
#Add prop unknown
ALL16WU<-lmer(ALL_16S_Sh_Div ~ ALL_prop_Uk*Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WU)
AIC(ALL16WU)
BIC(ALL16WU)
#Intercept significant AIC=42.2 BIC=60.4 better model
#Add ALL Richness
ALL16WUR<-lmer(ALL_16S_Sh_Div ~ ALL_Rich+ALL_prop_Uk*Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WUR)
AIC(ALL16WUR)
BIC(ALL16WUR)
#nothing significant AIC=44.4 BIC=63.7 worse model
#Add sampling day
ALL16WUT<-lmer(ALL_16S_Sh_Div ~ Sampling_Time+ALL_prop_Uk*Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WUT)
AIC(ALL16WUT)
BIC(ALL16WUT)
#intercept significant AIC=43.22 BIC=63.6 worse model
#Add proportion oak
ALL16WUO<-lmer(ALL_16S_Sh_Div ~ ALL_prop_Oak*ALL_prop_Uk*Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT)
summary(ALL16WUO)
AIC(ALL16WUO)
BIC(ALL16WUO)
#intercept significant, All prop oak, all prop Uk, gap location, watershed, and interactions significant AIC=-7.6 BIC=17.4 better model

#Visualize dataset
ggplot(EAB_Paired_NT, aes(x=ALL_prop_Oak, y=ALL_16S_Sh_Div, color=Gap_Location)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se=FALSE)+
  xlab("Aquatic leaf litter proportion oak")+
  ylab("Microbial Community Shannon Diversity")+
  scale_color_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a"), name="Gap Location")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title=element_text(size=18),
        legend.text = element_text(size=12))

#Use pd whole tree
#check for normality for all pd whole tree alpha diversity
shapiro.test(EAB_Paired_NT$ALL_16S_pd_Div)
#W=0.92, p=0.055, normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT$ALL_16S_pd_Div)
#slight right skewness
plot(EAB_Paired_NT$ALL_16S_pd_Div)
#Use normal distribution
#Start with simple model
ALLpd<-lmer(ALL_16S_pd_Div ~ Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALLpd)
AIC(ALLpd)
BIC(ALLpd)
#Intercept significant and gap location significant, AIC=157 BIC=163
#add watershed
ALLpdw<-lmer(ALL_16S_pd_Div ~ Watershed*Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALLpdw)
AIC(ALLpdw)
BIC(ALLpdw)
#Intercept significant and gap location approaching significant, AIC=132 BIC=143 better model
#add time
ALLpdwt<-lmer(ALL_16S_pd_Div ~ Sampling_Time*Watershed*Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALLpdwt)
AIC(ALLpdwt)
BIC(ALLpdwt)
#Intercept significant and gap location approachiing significant, AIC=102 BIC=120 better model
#add leaf litter richness
ALLpdwtr<-lmer(ALL_16S_pd_Div ~ ALL_Rich*Sampling_Time*Watershed*Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALLpdwtr)
AIC(ALLpdwtr)
BIC(ALLpdwtr)
#Intercept, all richness, watershed, gap location plus interactions significant, AIC=63 BIC=90 better model
#get rid of sampling time (can't include any more factors with sample size)
ALLpdwr<-lmer(ALL_16S_pd_Div ~ ALL_Rich*Watershed*Gap_Location+(1|Stream), data = EAB_Paired_NT)
summary(ALLpdwr)
AIC(ALLpdwr)
BIC(ALLpdwr)
#Intercept, all richness, watershed significant, AIC=110 BIC=128 worse model

#visualize data
ggplot(EAB_Paired_NT, aes(x=ALL_Rich, y=ALL_16S_pd_Div, color=Gap_Location)) +
  geom_point(size=5) +
  geom_smooth(method=lm, se=FALSE)+
  xlab("Aquatic leaf litter richness")+
  ylab("Microbial Phylogenetic Diversity")+
  scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"), name="Gap Location")+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.title=element_text(size=22),
        legend.text = element_text(size=16))

#determine whether downstream site is an outlier
plot(EAB_Paired_NT$ALL_Rich, EAB_Paired_NT$ALL_16S_pd_Div)
identify(EAB_Paired_NT$ALL_Rich, EAB_Paired_NT$ALL_16S_pd_Div)
#41, seven mile creek early upstream
EAB_Paired_NT_NO<-EAB_Paired_NT[-c(41),]
#check for normality for all pd whole tree alpha diversity
shapiro.test(EAB_Paired_NT_NO$ALL_16S_pd_Div)
#W=0.94, p=0.18, normally distributed
#Check for outliers
ggdensity(EAB_Paired_NT_NO$ALL_16S_pd_Div)
#slight right skewness, better than last run with outlier
plot(EAB_Paired_NT_NO$ALL_16S_pd_Div)
#Use normal distribution
#Start with simple model
ALLpdno<-lmer(ALL_16S_pd_Div ~ Gap_Location+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdno)
AIC(ALLpdno)
BIC(ALLpdno)
#Intercept significant, AIC=146 BIC=151
#Add watershed
ALLpdnow<-lmer(ALL_16S_pd_Div ~ Gap_Location*Watershed+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnow)
AIC(ALLpdnow)
BIC(ALLpdnow)
#Intercept significant and gap location approaching significance, AIC=130 BIC=140, better model
#Add time
ALLpdnowt<-lmer(ALL_16S_pd_Div ~ Gap_Location*Watershed*Sampling_Time+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnowt)
AIC(ALLpdnowt)
BIC(ALLpdnowt)
#Intercept significant and gap location approaching significance, AIC=100 BIC=117, better model
#Add leaf litter richness
ALLpdnowtr<-lmer(ALL_16S_pd_Div ~ Gap_Location*Watershed*Sampling_Time*ALL_Rich+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnowtr)
AIC(ALLpdnowtr)
BIC(ALLpdnowtr)
#Intercept gap location, sampling time, ALL_Rich and interactions significant , AIC=59 BIC=84, better model
#get rid of watershed (can't add more factors in)
ALLpdnotr<-lmer(ALL_16S_pd_Div ~ Gap_Location*Sampling_Time*ALL_Rich+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnotr)
AIC(ALLpdnotr)
BIC(ALLpdnotr)
#Intercept gap location, sampling time, ALL_Rich and interactions significant , AIC=77 BIC=97, worse model
#Add ALL prop decomposed
ALLpdnotru<-lmer(ALL_16S_pd_Div ~ Gap_Location*Sampling_Time*ALL_Rich*ALL_prop_Uk+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnotru)
AIC(ALLpdnotru)
BIC(ALLpdnotru)
#Intercept gap location, sampling time, ALL_Rich and interactions significant , AIC=72 BIC=94, worse model
#relplace ALL prop decomposed with prop Maple
ALLpdnotrm<-lmer(ALL_16S_pd_Div ~ Gap_Location*Sampling_Time*ALL_Rich*ALL_prop_Maple+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnotrm)
AIC(ALLpdnotrm)
BIC(ALLpdnotrm)
#Intercept gap location, sampling time, ALL_Rich and interactions significant , AIC=31 BIC=54, better model
#relplace ALL prop Maple with prop Elm
ALLpdnotre<-lmer(ALL_16S_pd_Div ~ Gap_Location*Sampling_Time*ALL_Rich*ALL_prop_Elm+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnotre)
AIC(ALLpdnotre)
BIC(ALLpdnotre)
#Intercept gap location, sampling time, ALL_Rich, ALL_Elm and interactions significant , AIC=41 BIC=62, worse model
#relplace ALL prop Elm with prop Beech
ALLpdnotrb<-lmer(ALL_16S_pd_Div ~ Gap_Location*Sampling_Time*ALL_Rich*ALL_prop_Beech+(1|Stream), data = EAB_Paired_NT_NO)
summary(ALLpdnotrb)
AIC(ALLpdnotrb)
BIC(ALLpdnotrb)
#Intercept gap location, sampling time, ALL_Rich and interactions significant , AIC=60 BIC=82, worse model


#visualize data
ggplot(EAB_Paired_NT_NO, aes(x=ALL_Rich, y=ALL_16S_pd_Div, color=Gap_Location)) +
  geom_point(size=5) +
  geom_smooth(method=lm, se=FALSE)+
  xlab("Aquatic leaf litter richness")+
  ylab("Microbial Faith's Phylogenetic Diversity")+
  scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"), name="Gap Location")+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.title=element_text(size=22),
        legend.text = element_text(size=16))
ggplot(EAB_Paired_NT_NO, aes(x=Gap_Location, y=ALL_16S_pd_Div, fill=Sampling_Time)) +
  geom_boxplot() +
  xlab("Gap Location")+
  ylab("Microbial Faith's Phylogenetic Diversity")+
  scale_fill_manual(values=sampling_time_col_vec, name="Sampling Time", labels=c("Pre Fall", "Leaf Fall", "Post Fall"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=20,margin=margin(50,0,0,0)),axis.title.y=element_text(size=19),
        axis.text.x=element_text(size=16),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))

#Now move to microbial community dataset
#upload unifrac distance matrix
EAB_16S_uni<-read.delim("~/Desktop/weighted_unifrac_otu_table_EAB_rare_meta.txt", sep = "\t", header = T )
row.names(EAB_16S_uni)<-EAB_16S_uni[,1]
#upload mapping file
EAB_16S_map <- read.table("~/Desktop/Demultiplexed/EAB_16S_Map_Filtered.txt", sep="\t", header=T)
row.names(EAB_16S_map)<-EAB_16S_map[,1]
#merge mapping file and unifrac distances
EAB_16S_uni_map <-merge(EAB_16S_map, EAB_16S_uni, by=0)
row.names(EAB_16S_uni_map)<-EAB_16S_uni_map[,1]
EAB_16S_uni_map<-EAB_16S_uni_map[,-c(1)]

#Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env<-EAB_16S_uni_map
row.names(EAB_16S_uni_env)<-EAB_16S_uni_env[,1]
EAB_16S_uni_env<-EAB_16S_uni_env[,1:21]

#subset to paired sites
EAB_16S_uni_map_pair<-subset(EAB_16S_uni_map, Stream=="StoneyCreek" & Gap=="1" | Stream=="SpringCreek" & Gap=="2" | Stream=="SevenMileCreek" & Gap=="1" | Stream=="FrayerCreek" & Gap=="1" | Stream=="SessionsCreek")
EAB_16S_uni_map_pair_com<-EAB_16S_uni_map_pair[,22:ncol(EAB_16S_uni_map_pair)]
EAB_16S_uni_map_pair_com_samples<-as.vector(rownames(EAB_16S_uni_map_pair_com))
#use output of names to subset columns into rows to make square matrix
EAB_16S_uni_map_pair_com<-as.matrix(EAB_16S_uni_map_pair_com)
EAB_16S_uni_map_pair_com<-subset(EAB_16S_uni_map_pair_com, select=c(EAB_16S_uni_map_pair_com_samples))
#Create environmental data matrix for community analysis with uni distances
EAB_16S_uni_map_pair_env<-EAB_16S_uni_map_pair[,1:21]
EAB_16S_uni_map_pair_ST<-factor(EAB_16S_uni_map_pair$SampleTime, levels=c("Pre Fall", "Leaf Fall", "Post Fall"))

#permanova test
adonis(as.dist(EAB_16S_uni_map_pair_com) ~ SampleTime*Location*Watershed, data=EAB_16S_uni_map_pair_env, permutations=999)
#visualize with nmds
EAB_NMDS_uni_pair<-metaMDS(as.dist(EAB_16S_uni_map_pair_com))
ordiplot(EAB_NMDS_uni_pair, type="n", main="Aquatic leaf litter 16S communities")
with(EAB_NMDS_uni_pair, points(EAB_NMDS_uni_pair, display="sites", col=sampling_time_col_vec[EAB_16S_uni_map_pair_ST], pch=19, pt.bg=sampling_time_col_vec))
with(EAB_NMDS_uni_pair, legend("topleft", legend=levels(EAB_16S_uni_map_pair_ST), bty="n", col=sampling_time_col_vec, pch=19, pt.bg=sampling_time_col_vec))
with(EAB_NMDS_uni_pair, ordiellipse(EAB_NMDS_uni_pair, EAB_16S_uni_map_pair_ST, kind="se", conf=0.95, lwd=2, col="#a6cee3", show.groups = "Pre Fall"))
with(EAB_NMDS_uni_pair, ordiellipse(EAB_NMDS_uni_pair, EAB_16S_uni_map_pair_ST, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "Leaf Fall"))
with(EAB_NMDS_uni_pair, ordiellipse(EAB_NMDS_uni_pair, EAB_16S_uni_map_pair_ST, kind="se", conf=0.95, lwd=2, col="#b2df8a", show.groups = "Post Fall"))


###########################################
#Modelling for ELME
#########################################

#Cast data into correct shape for model used in class
EAB_ELME<-EAB_Paired_Stand
EAB_ELME$ALL_prop_Unknown_Early<-NULL
EAB_ELME$ALL_prop_Unknown_Mid<-NULL
EAB_ELME$ALL_prop_Unknown_Late<-NULL
EAB_ELME$Terrestrial_Ash_CWD_Density_m<-NULL
EAB_ELME$Terrestrial_Total_CWD_Density_m<-NULL
EAB_ELME$Terrestrial_1_CWD_Density_m<-NULL
EAB_ELME$Aquatic_Ash_CWD_Density_m<-NULL
EAB_ELME$Aquatic_Total_CWD_Density_m<-NULL
EAB_ELME$Aquatic_1_CWD_Density_m<-NULL

EAB_ELME_ALL_R_E<-cast(EAB_ELME, Stream+Watershed+Gap_Size_m~Gap_location, value="ALL_Richness_Early") 
EAB_ELME_ALL_R_M<-cast(EAB_ELME, Stream+Watershed+Gap_Size_m~Gap_location, value="ALL_Richness_Mid") 
EAB_ELME_ALL_R_L<-cast(EAB_ELME, Stream+Watershed+Gap_Size_m~Gap_location, value="ALL_Richness_Late") 
EAB_ELME_ALL_R<-rbind(EAB_ELME_ALL_R_E,EAB_ELME_ALL_R_M,EAB_ELME_ALL_R_L)
EAB_ELME_RR<-cast(EAB_ELME, Stream+Watershed+Gap_Size_m~Gap_location, value="Terrestrial_Live_Rich_Stand") 
EAB_ELME_RR<-rbind(EAB_ELME_RR,EAB_ELME_RR,EAB_ELME_RR)

#model used in class
Watershed<-EAB_ELME_ALL_R$Watershed
gapsize<-EAB_ELME_ALL_R$Gap_Size_m
Stream<-EAB_ELME_ALL_R$Stream
ALLRUS<-EAB_ELME_ALL_R$Upstream
ALLRG<-EAB_ELME_ALL_R$Gap
ALLRDS<-EAB_ELME_ALL_R$Downstream
RRUS<-EAB_ELME_RR$Upstream
RRG<-EAB_ELME_RR$Gap
RRDS<-EAB_ELME_RR$Downstream

EAB_ELME_Model_Data<-data.frame(cbind(Watershed,gapsize,Stream,ALLRUS,ALLRG,ALLRDS,RRUS,RRG,RRDS))
EAB_ELME_Model_Data$Watershed<-as.factor(EAB_ELME_Model_Data$Watershed)

#data for predictions
ALLRUSP<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

inits<-list(variance = c(1,1), 
            M = c(0,0), 
            alpha1 = 0, alpha2 = 0,
            beta=0,
            theta=0,
            omega=0,
            mu1=0, mu2=0)

d <- read.jagsdata("~/Documents/MSU/Current_Coursework/ELME/ELME_EAB-Jagsdata.R")

m <- jags.model("~/Documents/MSU/Current_Coursework/ELME/enjoy/Courtney_JAGS.bug", d, inits=inits, n.chains=1)

update(m, 1000000)
x <- coda.samples(m, c("alpha1","alpha2","beta","theta", "omega", "mu1", "mu2", "M", "variance"), n.iter=1000000, thin = 1)
x0<-summary(x)
write.csv(x0[1],"parameters1.csv")
write.csv(x0[2],"parameters2.csv")
summary(x)  #quantiles
summary(x)$statistics #means and SDs

#chains and density plots
plot(x[,"alpha1"])
plot(x[,"alpha2"])
plot(x[,"beta"])
plot(x[,"theta"])
plot(x[,"omega"])
plot(x[,"M[1]"])
plot(x[,"M[2]"])
plot(x[,"variance[1]"])
plot(x[,"variance[2]"])
plot(x[,"mu1"])
plot(x[,"mu2"])

#or for all
png('plot1.png')
plot(x)
dev.off()

#Plot overall aquatic leaf litter upstream, downstream for each watershed
ggplot(EAB_ELME_Model_Data, aes(x=ALLRUS, y=ALLRG, color=Watershed)) +
  geom_point(size=5) +
  xlab("Leaf litter richness upstream")+
  ylab("Leaf litter richness gap")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title=element_text(size=18),
        legend.text = element_text(size=12))

ggplot(EAB_ELME_Model_Data, aes(x=ALLRG, y=ALLRDS, color=Watershed)) +
  geom_point(size=3) +
  xlab("Leaf litter richness gap")+
  ylab("Leaf litter richness downstream")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title=element_text(size=18),
        legend.text = element_text(size=12))

ggplot(EAB_ELME_Model_Data, aes(x=gapsize, y=ALLRG, color=Watershed)) +
  geom_point(size=3) +
  xlab("Gap size")+
  ylab("Leaf litter richness gap")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title=element_text(size=18),
        legend.text = element_text(size=12))

#Predictions
ALLRpred0 <- coda.samples(m, c("ALLRG.h","ALLRDS.h"), n.iter=1000000, thin = 1)
ALLRpred<-summary(ALLRpred0)$statistics
write.csv(ALLRpred,"ALLRGpred.csv")
ALLRpred

plot(ALLRpred[19:36,1],ALLRG, xlab = "predicted", ylab="observed")
xp<-c(0,10,20,30,40)
yp<-c(0,10,20,30,40)
lines(xp,yp)
cor(ALLRpred[19:36,1],ALLRG)

plot(ALLRpred[1:18,1],ALLRDS, xlab = "predicted", ylab="observed")
xp<-c(0,10,20,30,40)
yp<-c(0,10,20,30,40)
lines(xp,yp)
cor(ALLRpred[1:18,1],ALLRDS)
