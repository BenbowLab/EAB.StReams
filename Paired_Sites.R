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
library(lme4)
library(lmerTest)
library(vegan)
library(car)
library(DT)
library(moments)
library(EnvStats)
library(cAIC4)
library(e1071) 
library(forcats)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(polycor)
library(glmnet)
library(pastecs)
library(funrar)
library(plyr)
library(ecole)
library(indicspecies)
library(rcompanion)
library(caret)
library(ResourceSelection)
library(functional)
library(tidyverse)
library(devtools)
library(rpart)
library(party)
library(randomForest)
library(mvpart)
library(randomForestSRC)
library(rstatix)
library(broom)
library(leaps)
library(inTrees)
library(MuMIn)
library(reshape2)
library(data.table)
library(gvlma)
library(boxcoxmix)
library(VennDiagram)
library(psych)
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

#####################################################################
###
### Functions to calculate Predictive R-squared
###
#####################################################################


### This calculate the PRESS (predictive residual sum of squares), the lower, the better

#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#' 
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}

### This calculate the Predictive r-squared

#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
#'              the PRESS statistic.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}

### This calculate the R-squared, the adj-R-squared and the Predictive R-squared (from the functions above)

#' @title Model Fit Statistics
#' @description Returns lm model fit statistics R-squared, adjucted R-squared,
#'      predicted R-squared and PRESS.
#'      Thanks to John Mount for his 6-June-2014 blog post, R style tip: prefer functions that return data frames" for
#'      the idea \link{http://www.win-vector.com/blog/2014/06/r-style-tip-prefer-functions-that-return-data-frames}
#' @return Returns a data frame with one row and a column for each statistic
#' @param linear.model A \code{lm()} model.
model_fit_stats <- function(linear.model) {
  r.sqr <- summary(linear.model)$r.squared
  adj.r.sqr <- summary(linear.model)$adj.r.squared
  ratio.adjr2.to.r2 <- (adj.r.sqr/r.sqr)
  pre.r.sqr <- pred_r_squared(linear.model)
  press <- PRESS(linear.model)
  return.df <- data.frame("R-squared" = r.sqr, "Adj R-squared" = adj.r.sqr, 
                          "Ratio Adj.R2 to R2" = ratio.adjr2.to.r2, "Pred R-squared" = pre.r.sqr, PRESS = press)
  return(round(return.df,3))
}

# id: model id
# object: regsubsets object
# data: data used to fit regsubsets
# outcome: outcome variable
get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}

reg.output.search.with.test<- function (search_object) {  ## input an object from a regsubsets search
  ## First build a df listing model components and metrics of interest
  search_comp<-data.frame(R2=summary(search_object)$rsq,  
                          adjR2=summary(search_object)$adjr2,
                          BIC=summary(search_object)$bic,
                          CP=summary(search_object)$cp,
                          n_predictors=row.names(summary(search_object)$which),
                          summary(search_object)$which)
  ## Categorize different types of predictors based on whether '.' is present
  predictors<-colnames(search_comp)[(match("X.Intercept.",names(search_comp))+1):dim(search_comp)[2]]
  main_pred<-predictors[grep(pattern = ".", x = predictors, invert=T, fixed=T)]
  higher_pred<-predictors[grep(pattern = ".", x = predictors, fixed=T)]
  ##  Define a variable that indicates whether model should be reject, set to FALSE for all models initially.
  search_comp$reject_model<-FALSE  
  
  for(main_eff_n in 1:length(main_pred)){  ## iterate through main effects
    ## Find column numbers of higher level ters containing the main effect
    search_cols<-grep(pattern=main_pred[main_eff_n],x=higher_pred) 
    ## Subset models that are not yet flagged for rejection, only test these
    valid_model_subs<-search_comp[search_comp$reject_model==FALSE,]  
    ## Subset dfs with only main or higher level predictor columns
    main_pred_df<-valid_model_subs[,colnames(valid_model_subs)%in%main_pred]
    higher_pred_df<-valid_model_subs[,colnames(valid_model_subs)%in%higher_pred]
    
    if(length(search_cols)>0){  ## If there are higher level pred, test each one
      for(high_eff_n in search_cols){  ## iterate through higher level pred. 
        ##  Test if the intxn effect is present without main effect (working with whole column of models)
        test_responses<-((main_pred_df[,main_eff_n]==FALSE)&(higher_pred_df[,high_eff_n]==TRUE)) 
        valid_model_subs[test_responses,"reject_model"]<-TRUE  ## Set reject to TRUE where appropriate
      } ## End high_eff for
      ## Transfer changes in reject to primary df:
      search_comp[row.names(valid_model_subs),"reject_model"]<-valid_model_subs[,"reject_model"]
    } ## End if
  }  ## End main_eff for
  
  ## Output resulting table of all models named for original search object and current time/date
  current_time_date<-format(Sys.time(), "%m_%d_%y at %H_%M_%S")
  write.table(search_comp,paste(current_time_date,deparse(substitute(search_object)),
                                                                 "regSS_model_search.csv"),row.names=FALSE, col.names=TRUE, sep=",")
}  ## End reg.output.search.with.test fn

get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}

#Color vectors
Gap_location_col_vec<-c("#1b9e77", "#d95f02", "#7570b3")
sampling_time_col_vec<-c("#a6cee3", "#1f78b4", "#b2df8a")
Watershed_col_vec<-c("#af8dc3","#fdc086","#7fbf7b")
YGF_col_vec<-c("#67000d","#cb181d","#fb6a4a","#fcbba1")
GA_col_vec<-c("#c6dbef","#6baed6","#2171b5","#08306b","black")
ash_col_vec<-c("#b2df8a","#1f78b4")
ALL_Beech_col_vec<-c("#d8b365","#5ab4ac")
ALLR_col_vec<-c("#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d")
gapven_col_vec<-c("#1f78b4","#33a02c","#b15928","#999999")
woodsven_col_vec<-c("#a6cee3","#b2df8a","#ffff33","#d95f02")
WT_col_vec<-c("#99d8c9","#00441b")
GF_col_vec<-c("#66c2a5","#fdae61")
GFS_col_vec<-c("#00441b","#1b7837","#a6dba0","#543005","#bf812d","#dfc27d")
aspen_col_vec<-c("green","blue")
source_col_vec<-c("#1f78b4","#b15928","#33a02c")

#######################
#Calculate aquatic CWD metrics
###################
EAB_AqCWD<-read.csv("~/Documents/MSU/Research/Surveying/Stream_Woody_Debris_Modified_as_Costiganetal.csv", sep = ",", header = T )

#subset so only greater than 7.6 in diameter
EAB_AqCWD<-subset(EAB_AqCWD, Diameter>=7.6)

#Calculate volume of each piece of wood
EAB_AqCWD$V<-pi*((((EAB_AqCWD$Diameter/100)/2)^2)*(EAB_AqCWD$In_Stream_Length/100))
range(EAB_AqCWD$V)

#aggregate for count of wood in each reach
EAB_AqCWD_count<-aggregate(V~Stream+Gap+Basket, data=EAB_AqCWD, FUN=length)
#Add paired site dataset

#aggregate for volume sum of wood in each reach
EAB_AqCWD_Vol<-aggregate(V~Stream+Gap+Basket, data=EAB_AqCWD, FUN=sum)
#Add paired site dataset 

#Calculate volume of log jam wood per reach
EAB_AqCWD_logjamvol<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Log_jam=="Yes"), sum)
#Add to paired site dataset

#aggregate for count of log jam in each reach
EAB_AqCWD_logjamcount<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Log_jam=="Yes"), FUN=length)
#Add paired site dataset

#Create in situ variable
EAB_AqCWD$Insitu<-EAB_AqCWD$LW_Class%>% fct_collapse(Y = c("Bridge","Ramp"), N=c("Buried","Drift","Pinned"))

#aggregate for count of insitu in each reach
EAB_insituCWD_count<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Insitu=="Y"), FUN=length)
#Add paired site dataset

#aggregate for volume sum of insitu in each reach
EAB_insituCWD_Vol<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Insitu=="Y"), FUN=sum)
#Add paired site dataset 

#aggregate for count of ash in each reach
EAB_AshCWD_count<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Species=="Ash"), FUN=length)
#Add paired site dataset

#aggregate for volume sum of insitu in each reach
EAB_AshCWD_Vol<-aggregate(V~Stream+Gap+Basket, data=subset(EAB_AqCWD, Species=="Ash"), FUN=sum)
#Add paired site dataset 

#Upload paired site dataset
EAB_AqCWD_Paired<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_Aq_CWD.csv", sep = ",", header = T )

#Create densities and proportions
#Total amounts
#Number of logs per reach length
EAB_AqCWD_Paired$No_logs_per_length<-EAB_AqCWD_Paired$No_logs/EAB_AqCWD_Paired$Aquatic_Transect_length
#Number of logs per stream area
EAB_AqCWD_Paired$No_logs_per_area<-EAB_AqCWD_Paired$No_logs/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#VolCWD per reach length
EAB_AqCWD_Paired$VolCWD_per_length<-EAB_AqCWD_Paired$CWD_V_m3/EAB_AqCWD_Paired$Aquatic_Transect_length
#VolCWD per stream area
EAB_AqCWD_Paired$VolCWD_per_area<-EAB_AqCWD_Paired$CWD_V_m3/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#log jams
#Number of log jams per reach length
EAB_AqCWD_Paired$No_logjams_per_length<-EAB_AqCWD_Paired$No_logjams/EAB_AqCWD_Paired$Aquatic_Transect_length
#Number of log jams per stream area
EAB_AqCWD_Paired$No_logjams_per_area<-EAB_AqCWD_Paired$No_logjams/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#Number of log jam pieces per reach length
EAB_AqCWD_Paired$No_logjam_pieces_per_length<-EAB_AqCWD_Paired$No_logjam_CWD_pieces/EAB_AqCWD_Paired$Aquatic_Transect_length
#Number of log jam pieces per stream area
EAB_AqCWD_Paired$No_logjam_pieces_per_area<-EAB_AqCWD_Paired$No_logjam_CWD_pieces/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#VolCWD in log jams per reach length
EAB_AqCWD_Paired$VolCWDLJ_per_length<-EAB_AqCWD_Paired$Logjam_V_m3/EAB_AqCWD_Paired$Aquatic_Transect_length
#VolCWD in log jams per stream area
EAB_AqCWD_Paired$VolCWDLJ_per_area<-EAB_AqCWD_Paired$Logjam_V_m3/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#proportion pieces in logjams
EAB_AqCWD_Paired$prop_logjam_pieces<-EAB_AqCWD_Paired$No_logjam_CWD_pieces/EAB_AqCWD_Paired$No_logs
#proportion volume in logjams
EAB_AqCWD_Paired$prop_logjam_vol<-EAB_AqCWD_Paired$Logjam_V_m3/EAB_AqCWD_Paired$CWD_V_m3

#insitu
#Number of insitu pieces per reach length
EAB_AqCWD_Paired$No_insitu_per_length<-EAB_AqCWD_Paired$No_insitu/EAB_AqCWD_Paired$Aquatic_Transect_length
#Number of insitu pieces per stream area
EAB_AqCWD_Paired$No_insitu_per_area<-EAB_AqCWD_Paired$No_insitu/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#VolCWD in insitu per reach length
EAB_AqCWD_Paired$Vol_insitu_per_length<-EAB_AqCWD_Paired$insitu_V_m3/EAB_AqCWD_Paired$Aquatic_Transect_length
#VolCWD in insitu per stream area
EAB_AqCWD_Paired$Vol_insitu_per_area<-EAB_AqCWD_Paired$insitu_V_m3/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#proportion pieces insitu
EAB_AqCWD_Paired$prop_insitu_pieces<-EAB_AqCWD_Paired$No_insitu/EAB_AqCWD_Paired$No_logs
#proportion volume insitu
EAB_AqCWD_Paired$prop_insitu_vol<-EAB_AqCWD_Paired$insitu_V_m3/EAB_AqCWD_Paired$CWD_V_m3

#Ash
#Number of ash pieces per reach length
EAB_AqCWD_Paired$No_ash_per_length<-EAB_AqCWD_Paired$No_Ash/EAB_AqCWD_Paired$Aquatic_Transect_length
#Number of ash pieces per stream area
EAB_AqCWD_Paired$No_ash_per_area<-EAB_AqCWD_Paired$No_Ash/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#VolCWD in ash per reach length
EAB_AqCWD_Paired$Vol_ash_per_length<-EAB_AqCWD_Paired$Ash_V_m3/EAB_AqCWD_Paired$Aquatic_Transect_length
#VolCWD in ash per stream area
EAB_AqCWD_Paired$Vol_ash_per_area<-EAB_AqCWD_Paired$Ash_V_m3/(EAB_AqCWD_Paired$Aquatic_Transect_length*EAB_AqCWD_Paired$Aquatic_Transect_width)

#proportion pieces insitu
EAB_AqCWD_Paired$prop_ash_pieces<-EAB_AqCWD_Paired$No_Ash/EAB_AqCWD_Paired$No_logs
#proportion volume insitu
EAB_AqCWD_Paired$prop_ash_vol<-EAB_AqCWD_Paired$Ash_V_m3/EAB_AqCWD_Paired$CWD_V_m3

#aggregate into stream as experimental unit
EAB_AqCWD_Paired_ag<-aggregate(cbind(Aquatic_Transect_length,No_logs,Logjam_V_m3,No_logjam_CWD_pieces,
                                 CWD_V_m3,No_insitu,insitu_V_m3,No_logjams,No_Ash,Ash_V_m3)~Stream+
                                 Watershed+Aquatic_Transect_width+Year_Gapformation+Gap_Area+Gap_Diameter, data=EAB_AqCWD_Paired,
                               FUN=sum)
#Create densities and proportions

#Total amounts
#Number of logs per reach length
EAB_AqCWD_Paired_ag$No_logs_per_length<-EAB_AqCWD_Paired_ag$No_logs/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#Number of logs per stream area
EAB_AqCWD_Paired_ag$No_logs_per_area<-EAB_AqCWD_Paired_ag$No_logs/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#VolCWD per reach length
EAB_AqCWD_Paired_ag$VolCWD_per_length<-EAB_AqCWD_Paired_ag$CWD_V_m3/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#VolCWD per stream area
EAB_AqCWD_Paired_ag$VolCWD_per_area<-EAB_AqCWD_Paired_ag$CWD_V_m3/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#log jams
#Number of log jams per reach length
EAB_AqCWD_Paired_ag$No_logjams_per_length<-EAB_AqCWD_Paired_ag$No_logjams/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#Number of log jams per stream area
EAB_AqCWD_Paired_ag$No_logjams_per_area<-EAB_AqCWD_Paired_ag$No_logjams/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#Number of log jam pieces per reach length
EAB_AqCWD_Paired_ag$No_logjam_pieces_per_length<-EAB_AqCWD_Paired_ag$No_logjam_CWD_pieces/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#Number of log jam pieces per stream area
EAB_AqCWD_Paired_ag$No_logjam_pieces_per_area<-EAB_AqCWD_Paired_ag$No_logjam_CWD_pieces/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#VolCWD in log jams per reach length
EAB_AqCWD_Paired_ag$VolCWDLJ_per_length<-EAB_AqCWD_Paired_ag$Logjam_V_m3/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#VolCWD in log jams per stream area
EAB_AqCWD_Paired_ag$VolCWDLJ_per_area<-EAB_AqCWD_Paired_ag$Logjam_V_m3/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#proportion pieces in logjams
EAB_AqCWD_Paired_ag$prop_logjam_pieces<-EAB_AqCWD_Paired_ag$No_logjam_CWD_pieces/EAB_AqCWD_Paired_ag$No_logs
#proportion volume in logjams
EAB_AqCWD_Paired_ag$prop_logjam_vol<-EAB_AqCWD_Paired_ag$Logjam_V_m3/EAB_AqCWD_Paired_ag$CWD_V_m3

#insitu
#Number of insitu pieces per reach length
EAB_AqCWD_Paired_ag$No_insitu_per_length<-EAB_AqCWD_Paired_ag$No_insitu/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#Number of insitu pieces per stream area
EAB_AqCWD_Paired_ag$No_insitu_per_area<-EAB_AqCWD_Paired_ag$No_insitu/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#VolCWD in insitu per reach length
EAB_AqCWD_Paired_ag$Vol_insitu_per_length<-EAB_AqCWD_Paired_ag$insitu_V_m3/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#VolCWD in insitu per stream area
EAB_AqCWD_Paired_ag$Vol_insitu_per_area<-EAB_AqCWD_Paired_ag$insitu_V_m3/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#proportion pieces insitu
EAB_AqCWD_Paired_ag$prop_insitu_pieces<-EAB_AqCWD_Paired_ag$No_insitu/EAB_AqCWD_Paired_ag$No_logs
#proportion volume insitu
EAB_AqCWD_Paired_ag$prop_insitu_vol<-EAB_AqCWD_Paired_ag$insitu_V_m3/EAB_AqCWD_Paired_ag$CWD_V_m3

#Ash
#Number of ash pieces per reach length
EAB_AqCWD_Paired_ag$No_ash_per_length<-EAB_AqCWD_Paired_ag$No_Ash/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#Number of ash pieces per stream area
EAB_AqCWD_Paired_ag$No_ash_per_area<-EAB_AqCWD_Paired_ag$No_Ash/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#VolCWD in ash per reach length
EAB_AqCWD_Paired_ag$Vol_ash_per_length<-EAB_AqCWD_Paired_ag$Ash_V_m3/EAB_AqCWD_Paired_ag$Aquatic_Transect_length
#VolCWD in ash per stream area
EAB_AqCWD_Paired_ag$Vol_ash_per_area<-EAB_AqCWD_Paired_ag$Ash_V_m3/(EAB_AqCWD_Paired_ag$Aquatic_Transect_length*EAB_AqCWD_Paired_ag$Aquatic_Transect_width)

#proportion pieces insitu
EAB_AqCWD_Paired_ag$prop_ash_pieces<-EAB_AqCWD_Paired_ag$No_Ash/EAB_AqCWD_Paired_ag$No_logs
#proportion volume insitu
EAB_AqCWD_Paired_ag$prop_ash_vol<-EAB_AqCWD_Paired_ag$Ash_V_m3/EAB_AqCWD_Paired_ag$CWD_V_m3

####################
#Reduce output CWD variables by correlation analysis
###################
#Check normality assumptions

#Total CWD

#Number of logs per reach length
shapiro.test(EAB_AqCWD_Paired_ag$No_logs_per_length)
#W = 0.89675, p-value = 0.3069 normal

#VolCWD per reach length
shapiro.test(EAB_AqCWD_Paired_ag$VolCWD_per_length)
#W = 0.9322, p-value = 0.6047, normal

#log jams

#Number of log jams per reach length
shapiro.test(EAB_AqCWD_Paired_ag$No_logjams_per_length)
#W = 0.86348, p-value = 0.2014, normal
#Number of log jams per stream area
shapiro.test(EAB_AqCWD_Paired_ag$No_logjams_per_area)
#W = 0.81451, p-value = 0.07906, normal

#Number of log jam pieces per reach length
shapiro.test(EAB_AqCWD_Paired_ag$No_logjam_pieces_per_length)
#W = 0.94719, p-value = 0.7378, normal

#VolCWD in log jams per reach length
shapiro.test(EAB_AqCWD_Paired_ag$VolCWDLJ_per_length)
#W = 0.86728, p-value = 0.2107, normal

#proportion pieces in logjams
shapiro.test(EAB_AqCWD_Paired_ag$prop_logjam_pieces)
#W = 0.92063, p-value = 0.4581, normal
#proportion volume in logjams
shapiro.test(EAB_AqCWD_Paired_ag$prop_logjam_vol)
#W = 0.92333, p-value = 0.5539, normal

#insitu

#Number of insitu pieces per reach length
shapiro.test(EAB_AqCWD_Paired_ag$No_insitu_per_length)
#W = 0.92252, p-value = 0.3194, normal

#VolCWD in insitu per reach length
shapiro.test(EAB_AqCWD_Paired_ag$Vol_insitu_per_length)
#W = 0.897, p-value = 0.3521, normal

#proportion pieces insitu
shapiro.test(EAB_AqCWD_Paired_ag$prop_insitu_pieces)
#W = 0.90399, p-value = 0.1634, normal
#proportion volume insitu
shapiro.test(EAB_AqCWD_Paired_ag$prop_insitu_vol)
#W = 0.87537, p-value = 0.266, normal
EAB_AqCWD_Paired_ag %>%
  get_summary_stats(prop_insitu_vol, type = "mean_se")

#Ash

#Number of ash pieces per reach length
shapiro.test(EAB_AqCWD_Paired_ag$No_ash_per_length)
#W = 0.97721, p-value = 0.9369, normal
#Number of ash pieces per stream area
shapiro.test(EAB_AqCWD_Paired_ag$No_ash_per_area)
#W = 0.91198, p-value = 0.4495, normal

#VolCWD in ash per reach length
shapiro.test(EAB_AqCWD_Paired_ag$Vol_ash_per_length)
#W = 0.95043, p-value = 0.7438, normal
#VolCWD in ash per stream area
shapiro.test(EAB_AqCWD_Paired_ag$Vol_ash_per_area)
#W = 0.97827, p-value = 0.9426, normal

#proportion pieces ash
shapiro.test(EAB_AqCWD_Paired_ag$prop_ash_pieces)
#W = 0.96054, p-value = 0.8558, normal
#proportion volume insitu
shapiro.test(EAB_AqCWD_Paired_ag$prop_ash_vol)
#W = 0.95822, p-value = 0.8061, normal

#############################
#assumption of normality met on all output variables
#############################

#Get summary stats on variables for correlation analysis
AqCWDSum<-describeBy(EAB_AqCWD_Paired_ag[,17:40], EAB_AqCWD_Paired_ag$Year_Gapformation)
write.csv(AqCWDSum$Clinton, "AqCWDSumCl.csv")
write.csv(AqCWDSum$Grand, "AqCWDSumGr.csv")
write.csv(AqCWDSum$Kalamazoo, "AqCWDSumKa.csv")


#Run pearson correlation on all output variables to reduce variable number
names(EAB_AqCWD_Paired_ag)
AqCWDCor<-cor(EAB_AqCWD_Paired_ag[,17:40])
write.csv(AqCWDCor,'AqCWDCor.csv')

#exclude one variable at a time based on correlations |>0.6|
#Also, eliminate all "per area" variables

#prop ash pieces and no ash per length 0.87, delete prop ash pieces
#prop logjam pieces and no ash per length -0.94, delete prop logjam pieces
#no ash per length and vol insitu per length 0.86, delete no ash per length
#vol insitu per length and prop insitu vol 0.85, delete vol insitu per length
#no logjam pieces per length and prop logjam vol 0.98, delete no logjam pieces per length
#vol ash per length and vol cwd per length 0.9, delete vol ash per length
#prop logjam vol and volcwdlj per length 0.77, delete prop logjam vol
#prop insitu pieces and no insitu per length 0.97, delete prop insitu pieces
#no logs per length and vol cwd per length 0.81, delete vol cwd per length
#no logjams per length and volcwdlj per length 0.69, delete vol cwd lj per length
#no insitu per length and prop insitu vol 0.93, delete prop insitu vol

#4 variables left: no logs per length, no logjams per length, no insitu per length
#and prop ash vol

#Use these variables as output

##############
#Aquatic CWD watershed
###################

#model this way: one way anova aq CWD ~ watershed

#No_logs_per_length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logs_per_length, type = "mean_se")
EAB_AqCWD_Paired_ag %>%
  get_summary_stats(No_logs_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logs_per_length)
#no outliers
#build model
lm.No_logspm<-lm(No_logs_per_length ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.No_logspm))
#all within gray bar
shapiro_test(residuals(lm.No_logspm))
#p 0.964 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logs_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.No_logspm, 1)
#no relationship
#all assumptions  met
No_logspm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logs_per_length~Watershed)
No_logspm.aov
#not significant

#number of log jams per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjams_per_length, type = "mean_se")
EAB_AqCWD_Paired_ag %>%
  get_summary_stats(No_logjams_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logjams_per_length)
#no outliers
#build model
lm.NoLJpm<-lm(No_logjams_per_length ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.NoLJpm))
#all within gray bar
shapiro_test(residuals(lm.NoLJpm))
#p 0.999 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logjams_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.NoLJpm, 1)
#no relationship
#all assumptions  met
NoLJpm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logjams_per_length~Watershed)
NoLJpm.aov
#not significant: F=0.696 p=0.565 ges=0.317

#No insitu per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_insitu_per_length, type = "mean_se")
EAB_AqCWD_Paired_ag %>%
  get_summary_stats(No_insitu_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_insitu_per_length)
#no outliers
#build model
lm.No_insitupm<-lm(No_insitu_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.No_insitupm))
#points all inside gray bar
shapiro_test(residuals(lm.No_insitupm))
#p 0.923 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_insitu_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.No_insitupm, 1)
#no relationship
#all assumptions  met
No_insitupm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_insitu_per_length~Watershed)
No_insitupm.aov
#not significant: F=0.001, p=0.999, ges=0.000897

#prop ash volume
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_vol, type = "mean_se")
EAB_AqCWD_Paired_ag %>%
  get_summary_stats(prop_ash_vol, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_ash_vol)
#no outliers
#build model
lm.propashv<-lm(prop_ash_vol~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propashv))
#points all inside gray bar
shapiro_test(residuals(lm.propashv))
#p 0.849 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_ash_vol", facet.by = "Watershed")
#all within gray bar
plot(lm.propashv, 1)
#no relationship
#all assumptions  met
propashv.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_ash_vol~Watershed)
propashv.aov
#not significant: F=2.214 p=0.257, ges=0.596

##########################
#Summary Aquatic CWD by watershed
#nothing significant
##############################

#Get summary stats for dissertation
#logjam pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_pieces, type = "mean_se")
#logjam volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_vol, type = "mean_se")
#insitu pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_pieces, type = "mean_se")
#insitu volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_vol, type = "mean_se")
#ash pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_pieces, type = "mean_se")
#ash volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_vol, type = "mean_se")

#Try again using year since gap formation as predictor in linear model

#Change year of gap formation to years since gap formation
EAB_AqCWD_Paired_ag$YearSinceGF<-2016-EAB_AqCWD_Paired_ag$Year_Gapformation

#merge year since EAB invasion with 
#no logs per length
No_logspm.lm <- lm(No_logs_per_length~YearSinceGF, data=EAB_AqCWD_Paired_ag)
summary(No_logspm.lm)
#not significant p=0.269

#number of log jams per length
No_ljpm.lm <- lm(No_logjams_per_length~YearSinceGF, data=EAB_AqCWD_Paired_ag)
summary(No_ljpm.lm)
#not significant p=0.41

#No insitu per length
No_ipm.lm <- lm(No_insitu_per_length~YearSinceGF, data=EAB_AqCWD_Paired_ag)
summary(No_ipm.lm)
#not significant p=0.43

#prop ash volume
PA.lm <- lm(prop_ash_vol~YearSinceGF, data=EAB_AqCWD_Paired_ag)
summary(PA.lm)
#not significant p=0.39

############################
#Calculate terrestrial forest ecology metrics: standing
###########################
#Upload data set
EAB_Terrstand_Paired<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Riparian_Survey_Pat_r.csv", sep = ",", header = T )

#Create individual tree basal area (in meters squared) column
EAB_Terrstand_Paired$IndBA_m2<-pi*(EAB_Terrstand_Paired$DBH_m/2)^2

#aggregate dataset so that there are total sum of tree type for each transect/plot
EAB_Terrstand_Paired_vol_ag<-aggregate(IndBA_m2 ~ Stream+Gap_location+Species_Pat+Material,
                              data=EAB_Terrstand_Paired,FUN=sum)
#Calculate overall basal area for each plot
#Create variable with area surveyed
EAB_Terrstand_Paired_vol_ag$Area_Surveyed_ha<-ifelse(EAB_Terrstand_Paired_vol_ag$Gap_location=="Gap",
                        0.1116,0.0458) 
EAB_Terrstand_Paired_vol_ag$BA_m2perha<-EAB_Terrstand_Paired_vol_ag$IndBA_m2/EAB_Terrstand_Paired_vol_ag$Area_Surveyed_ha

#Cast datset so that stream and gap location are on y and species and material are on x
EAB_Terrstand_Paired_vol_cast<-dcast(EAB_Terrstand_Paired_vol_ag,Stream+Gap_location~Species_Pat+Material)
EAB_Terrstand_Paired_vol_cast[is.na(EAB_Terrstand_Paired_vol_cast)] <- 0
EAB_Terrstand_Paired_vol_cast$Gap_location<-revalue(EAB_Terrstand_Paired_vol_cast$Gap_location, c("US"="Upstream", "DS"="Downstream"))
colnames(EAB_Terrstand_Paired_vol_cast) <- c("Stream","Gap_location", paste0('BAperha_', colnames(EAB_Terrstand_Paired_vol_cast[3:ncol(EAB_Terrstand_Paired_vol_cast)])))

#aggregate dataset so that there are total count of tree type for each transect/plot
EAB_Terrstand_Paired_no_ag<-aggregate(IndBA_m2 ~ Stream+Gap_location+Species_Pat+Material,
                                       data=EAB_Terrstand_Paired,FUN=length)
#Calculate overall basal area for each plot
#Create variable with area surveyed
EAB_Terrstand_Paired_no_ag$Area_Surveyed_ha<-ifelse(EAB_Terrstand_Paired_no_ag$Gap_location=="Gap",
                                                     0.1116,0.0458) 
EAB_Terrstand_Paired_no_ag$No_m2perha<-EAB_Terrstand_Paired_no_ag$IndBA_m2/EAB_Terrstand_Paired_no_ag$Area_Surveyed_ha

#Cast datset so that stream and gap location are on y and species and material are on x
EAB_Terrstand_Paired_no_cast<-dcast(EAB_Terrstand_Paired_no_ag,Stream+Gap_location~Species_Pat+Material)
EAB_Terrstand_Paired_no_cast[is.na(EAB_Terrstand_Paired_no_cast)] <- 0
EAB_Terrstand_Paired_no_cast$Gap_location<-revalue(EAB_Terrstand_Paired_no_cast$Gap_location, c("US"="Upstream", "DS"="Downstream"))
colnames(EAB_Terrstand_Paired_no_cast) <- c("Stream","Gap_location", paste0('Noperha_', colnames(EAB_Terrstand_Paired_no_cast[3:ncol(EAB_Terrstand_Paired_no_cast)])))

#merge terrestrial standing datasets
EAB_Terrstand_Paired_stand<-merge(EAB_Terrstand_Paired_vol_cast,EAB_Terrstand_Paired_no_cast,by=c("Stream","Gap_location"))

#########################
#Calculate terrestrial forest ecology metrics: fallen cwd
#########################
EAB_TerrCWD_Paired<-read.csv("~/Documents/MSU/Research/Surveying/Riparian/CWM_EAB_Survey.csv", sep = ",", header = T )

#Aggregate so sum for each gap location
EAB_TerrCWD_Paired_ag_vol<-aggregate(Volume_m3 ~ Stream+Gap_location+Taxa,
                              data=subset(EAB_TerrCWD_Paired,Survey_Type=="Terrestrial"),
                              FUN=sum)
#Convert to volume per ha surveyed
EAB_TerrCWD_Paired_ag_vol$Area_Surveyed_ha<- ifelse(EAB_TerrCWD_Paired_ag_vol$Gap_location=="Gap",
                                             0.02,0.005) 
EAB_TerrCWD_Paired_ag_vol$VolCWD_m3perha<-EAB_TerrCWD_Paired_ag_vol$Volume_m3/EAB_TerrCWD_Paired_ag_vol$Area_Surveyed_ha
#Convert to volume per length surveyed
EAB_TerrCWD_Paired_ag_vol$Length_Surveyed_m<- ifelse(EAB_TerrCWD_Paired_ag_vol$Gap_location=="Gap",
                                                150,25) 
EAB_TerrCWD_Paired_ag_vol$VolCWD_m3perm<-EAB_TerrCWD_Paired_ag_vol$Volume_m3/EAB_TerrCWD_Paired_ag_vol$Length_Surveyed_m

#Cast datset so that stream and gap location are on y and species on x for vol per area
EAB_TerrCWD_Paired_volperha_cast<-dcast.data.table(setDT(EAB_TerrCWD_Paired_ag_vol),Stream+Gap_location~Taxa,
                                            value.var=c("VolCWD_m3perha", "VolCWD_m3perm"))
EAB_TerrCWD_Paired_volperha_cast[is.na(EAB_TerrCWD_Paired_volperha_cast)] <- 0

#Aggregate so count for each gap location
EAB_TerrCWD_Paired_ag_no<-aggregate(Volume_m3 ~ Stream+Gap_location+Taxa,
                                     data=subset(EAB_TerrCWD_Paired,Survey_Type=="Terrestrial"),
                                     FUN=length)
#Convert to number per ha surveyed
EAB_TerrCWD_Paired_ag_no$Area_Surveyed_ha<- ifelse(EAB_TerrCWD_Paired_ag_no$Gap_location=="Gap",
                                                    0.02,0.005) 
EAB_TerrCWD_Paired_ag_no$NoCWDperha<-EAB_TerrCWD_Paired_ag_no$Volume_m3/EAB_TerrCWD_Paired_ag_no$Area_Surveyed_ha
#Convert to number per length surveyed
EAB_TerrCWD_Paired_ag_no$Length_Surveyed_m<- ifelse(EAB_TerrCWD_Paired_ag_no$Gap_location=="Gap",
                                                     150,25) 
EAB_TerrCWD_Paired_ag_no$NoCWDperm<-EAB_TerrCWD_Paired_ag_no$Volume_m3/EAB_TerrCWD_Paired_ag_no$Length_Surveyed_m

#Cast datset so that stream and gap location are on y and species on x for vol per area
EAB_TerrCWD_Paired_noperha_cast<-dcast.data.table(setDT(EAB_TerrCWD_Paired_ag_no),Stream+Gap_location~Taxa,
                                                   value.var=c("NoCWDperha", "NoCWDperm"))
EAB_TerrCWD_Paired_noperha_cast[is.na(EAB_TerrCWD_Paired_noperha_cast)] <- 0

#merge no and vol datasets
EAB_TerrCWD_Paired_stand<-merge(EAB_TerrCWD_Paired_volperha_cast,EAB_TerrCWD_Paired_noperha_cast, by=c("Stream","Gap_location"))

##########################
#Merge terrestrial datasets and look for correlations
##########################
EAB_Terr_Paired<-merge(EAB_Terrstand_Paired_stand,EAB_TerrCWD_Paired_stand,by=c("Stream","Gap_location"),
                       all=TRUE)
EAB_Terr_Paired[is.na(EAB_Terr_Paired)] <- 0
#Calculate total amounts for snag and cwd
EAB_Terr_Paired_dead<-EAB_Terr_Paired[, -grep("_Live", colnames(EAB_Terr_Paired))]
EAB_Terr_Paired_dead$Total_Snag_BA_m2pha<-rowSums(EAB_Terr_Paired_dead[, grep("BAperha_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$Total_Snag_Nopha<-rowSums(EAB_Terr_Paired_dead[, grep("Noperha_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$Total_CWD_Vol_m3pha<-rowSums(EAB_Terr_Paired_dead[, grep("VolCWD_m3perha_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$Total_CWD_Vol_m3pm<-rowSums(EAB_Terr_Paired_dead[, grep("VolCWD_m3perm_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$Total_CWD_Nopha<-rowSums(EAB_Terr_Paired_dead[, grep("NoCWDperha_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$Total_CWD_Nopm<-rowSums(EAB_Terr_Paired_dead[, grep("NoCWDperm_", names(EAB_Terr_Paired_dead))])
EAB_Terr_Paired_dead$prop_Ash_snag_no<-EAB_Terr_Paired_dead$Noperha_Ash_Snag/EAB_Terr_Paired_dead$Total_Snag_Nopha
EAB_Terr_Paired_dead$prop_Ash_snag_BA<-EAB_Terr_Paired_dead$BAperha_Ash_Snag/EAB_Terr_Paired_dead$Total_Snag_BA_m2pha
EAB_Terr_Paired_dead$prop_Ash_CWD_no<-EAB_Terr_Paired_dead$NoCWDperm_Ash/EAB_Terr_Paired_dead$Total_CWD_Nopm
EAB_Terr_Paired_dead$prop_Ash_CWD_vol<-EAB_Terr_Paired_dead$VolCWD_m3perha_Ash/EAB_Terr_Paired_dead$Total_CWD_Vol_m3pha

#Create new dataset that subsets just gap
EAB_Terr_Paired_Gap<-subset(EAB_Terr_Paired_dead, Gap_location=="Gap")
#delete variables with 0 values
EAB_Terr_Paired_Gap<-EAB_Terr_Paired_Gap[, colSums(EAB_Terr_Paired_Gap != 0) > 0]
#Calculate summary stats for manuscript
EAB_Terr_Paired_Gap %>%  get_summary_stats(Total_CWD_Vol_m3pha, type = "mean_se")
EAB_Terr_Paired_Gap %>%  get_summary_stats(Total_Snag_BA_m2pha, type = "mean_se")

#Create new dataset that subsets just forest
EAB_Terr_Paired_F<-subset(EAB_Terr_Paired_dead, Gap_location!="Gap")
#delete variables with 0 values
EAB_Terr_Paired_F<-EAB_Terr_Paired_F[, colSums(EAB_Terr_Paired_F != 0) > 0]
#Calculate summary stats for manuscript
EAB_Terr_Paired_F %>%  get_summary_stats(Total_CWD_Vol_m3pha, type = "mean_se")
EAB_Terr_Paired_F %>%  get_summary_stats(Total_Snag_BA_m2pha, type = "mean_se")

#########################
#CWD linear regression with terrestrial datasets
#######################

#Merge with terrestrial dataset
EAB_CWD_Paired<-merge(EAB_AqCWD_Paired_ag, EAB_Terr_Paired_Gap, by="Stream")

#Test normality for terrestrial variables used:
#Number of terrestrial cwd logs per ha
shapiro.test(EAB_CWD_Paired$Total_CWD_Nopha)
#W = 0.78069, p-value = 0.03912, not normal
skewness(EAB_CWD_Paired$Total_CWD_Nopha)
#right skewed, log transform
EAB_CWD_Paired$log10Total_CWD_Nopha<-log10(EAB_CWD_Paired$Total_CWD_Nopha+1)
shapiro.test(EAB_CWD_Paired$log10Total_CWD_Nopha)
#W = 0.88412, p-value = 0.288, normal, use log10 transformation

#number of snags per ha
shapiro.test(EAB_CWD_Paired$Total_Snag_Nopha)
#W = 0.90225, p-value = 0.3874, normal
ggqqplot(EAB_CWD_Paired$Total_Snag_Nopha)
#all within gray bar - normally distributed

#Volume CWD per ha
shapiro.test(EAB_CWD_Paired$Total_CWD_Vol_m3pha)
#W = 0.72619, p-value = 0.01153, not normal
skewness(EAB_CWD_Paired$Total_CWD_Vol_m3pha)
#right skewed, log transform
EAB_CWD_Paired$log10Total_CWD_Vol_m3pha<-log10(EAB_CWD_Paired$Total_CWD_Vol_m3pha+1)
shapiro.test(EAB_CWD_Paired$log10Total_CWD_Vol_m3pha)
#W = 0.83712, p-value = 0.1234, normal
ggqqplot(EAB_CWD_Paired$log10Total_CWD_Vol_m3pha)
#all within gray bar - normally distributed, use log10 transformation

#volume snags per ha
shapiro.test(EAB_CWD_Paired$Total_Snag_BA_m2pha)
#W = 0.85384, p-value = 0.169, normal
ggqqplot(EAB_CWD_Paired$Total_Snag_BA_m2pha)
#all within gray bar - normally distributed

###################
#linear regression

#noCWD per stream length - CWD use log transformation
nCWDpl.nCWD.lm<-lm(No_logs_per_length~log10Total_CWD_Nopha, data=EAB_CWD_Paired)
summary(nCWDpl.nCWD.lm)
#not significant
#snag - normal
nCWDpl.ns.lm<-lm(No_logs_per_length~Total_Snag_Nopha, data=EAB_CWD_Paired)
summary(nCWDpl.ns.lm)
#not significant

#Number of log jams per stream length
#CWD - not normal,
nljpl.nCWD.lm<-lm(No_logjams_per_length~log10Total_CWD_Nopha, data=EAB_CWD_Paired)
summary(nljpl.nCWD.lm)
#not significant
#snag - normal Total_Snag_Nopha
nljpl.ns.lm<-lm(No_logjams_per_length~Total_Snag_Nopha, data=EAB_CWD_Paired)
summary(nljpl.ns.lm)
#not significant

#number pieces insitu
nipl.nCWD.lm<-lm(No_insitu_per_length~log10Total_CWD_Nopha, data=EAB_CWD_Paired)
summary(nipl.nCWD.lm)
#not significant
#snag - normal 
nipl.ns.lm<-lm(No_insitu_per_length~Total_Snag_Nopha, data=EAB_CWD_Paired)
summary(nipl.ns.lm)
#significant = p<0.01
#visualize
propinsituvollab<-expression(paste("Aquatic ", italic("In situ"), " (#/m)"))
ggplot(EAB_CWD_Paired, aes(y=No_insitu_per_length, x=Total_Snag_Nopha)) +
  geom_point(size=3) +
  xlab("Number Standing Dead Trees per ha")+
  ylab(propinsituvollab)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")

#prop CWD volume ash
#CWD - log transformation, log10Total_CWD_Vol_m3pha
pa.vCWD.lm<-lm(prop_Ash_CWD_vol~log10Total_CWD_Vol_m3pha, data=EAB_CWD_Paired)
summary(pa.vCWD.lm)
#not significant
#snag - normal Total_Snag_BA_m2pha
pa.vs.lm<-lm(prop_Ash_CWD_vol~Total_Snag_BA_m2pha, data=EAB_CWD_Paired)
summary(pa.vs.lm)
#not significant

#######################################
#Summary of linear regressions:
#only one sigificant - prop insitu vol ~ total CWD vol per ha
########################################

####################################
#Water chemistry
#################################
#Upload data
Paired_H2OEnv<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Env.csv", sep = ",", header = T )
#Merge with terrestrial dataset
EAB_Ter_H20<-merge(Paired_H2OEnv, EAB_Terr_Paired,by=c("Stream","Gap_location"))
#merge with aquatic cwd dataset
EAB_Ter_H2o_aqcwd<-merge(EAB_Ter_H20,EAB_AqCWD_Paired,by=c("Stream","Gap_location"))

#Want to see whether gap location, watershed, terr total basal area has impact
#first try to model like h2o~watershed+gap location+total basal area+1|stream mixed model
#assumptions of linearity, homogeneity of variance, residuals of model are normally distributed

#create total basal area variable
EAB_Ter_H2o_aqcwd$Total_Live_BA<-(rowSums(EAB_Ter_H2o_aqcwd[, grep("BAperha_",names(EAB_Ter_H2o_aqcwd))]))-(rowSums(EAB_Ter_H2o_aqcwd[, grep("_snag",names(EAB_Ter_H2o_aqcwd))]))

#summary stats for manuscript
#Create new dataset that subsets just gap
EAB_Ter_H2o_aqcwd_Gap<-subset(EAB_Ter_H2o_aqcwd, Gap_location=="Gap")
#Calculate summary stats for manuscript
EAB_Ter_H2o_aqcwd_Gap %>%  get_summary_stats(Total_Live_BA, type = "mean_se")

#Create new dataset that subsets just forest
EAB_Ter_H2o_aqcwd_F<-subset(EAB_Ter_H2o_aqcwd, Gap_location!="Gap")
#Calculate summary stats for manuscript
EAB_Ter_H2o_aqcwd %>%  get_summary_stats(Total_Live_BA, type = "mean_se")

#Change order of levels for gap location
EAB_Ter_H2o_aqcwd$Gap_location<- factor(EAB_Ter_H2o_aqcwd$Gap_location, levels = c("Upstream", "Gap", "Downstream"))
#first do mixed models for gap location

#Water_temp
#build model
lmer.wt= lme(Water_temp~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
            data=EAB_Ter_H2o_aqcwd,
            method="REML")
#check assumptions
hist(residuals(lmer.wt),col="darkgray")
#not normal
#Check for normality
shapiro.test(EAB_Ter_H2o_aqcwd$Water_temp)
#W = 0.83708, p-value = 3.517e-06, not normal
hist(EAB_Ter_H2o_aqcwd$Water_temp)
skewness(EAB_Ter_H2o_aqcwd$Water_temp)
#not very skewed, so don't transform
range(EAB_Ter_H2o_aqcwd$Water_temp) 
shapiro.test(EAB_Ter_H2o_aqcwd$Total_Live_BA)
#W = 0.94163, p-value = 0.0109, not normal
hist(EAB_Ter_H2o_aqcwd$Total_Live_BA)
skewness(EAB_Ter_H2o_aqcwd$Total_Live_BA)
#very small positive skew, log transform
EAB_Ter_H2o_aqcwd$log10Total_Live_BA<-log10(EAB_Ter_H2o_aqcwd$Total_Live_BA)
shapiro.test(EAB_Ter_H2o_aqcwd$log10Total_Live_BA)
#worse, so use untransformed
summary(lmer.wt)
#not significant
EAB_Ter_H2o_aqcwd %>%
  group_by(Watershed) %>%
  get_summary_stats(Water_temp, type = "mean_se")

#mS_cm
#build model
lmer.mScm= lme(mScm~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
             data=EAB_Ter_H2o_aqcwd,
             method="REML")
#check assumptions
hist(residuals(lmer.mScm),col="darkgray")
#normal
summary(lmer.mScm)
#not significant
EAB_Ter_H2o_aqcwd %>%
  group_by(Watershed) %>%
  get_summary_stats(mScm, type = "mean_se")

#X.DO
#build model
lmer.DO= lme(X.DO~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
               data=EAB_Ter_H2o_aqcwd,
               method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.DO),col="darkgray")
#normal
summary(lmer.DO)
#downstream significantly lower
#visualize
ggplot(EAB_Ter_H2o_aqcwd, aes(x=Total_Live_BA, y=X.DO, color=Gap_location)) +
  geom_point(size=3) +
  ylab("Percent Dissolved Oxygen")+
  xlab("Riparian Total Live BA")+
  labs(color = "Gap Location")+
  scale_color_manual(values=Gap_location_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#summary stats
EAB_Ter_H2o_aqcwd %>%
  group_by(Gap_location) %>%
  get_summary_stats(X.DO, type = "mean_se")
ggplot(EAB_Ter_H2o_aqcwd, aes(x=Gap_location, y=X.DO)) +
  geom_boxplot() +
  ylab("Percent Dissolved Oxygen")+
  xlab("Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
EAB_Ter_H2o_aqcwd$Watershed = forcats::fct_rev(factor(EAB_Ter_H2o_aqcwd$Watershed))
ggplot(EAB_Ter_H2o_aqcwd, aes(x=Gap_location, y=X.DO, fill=Watershed)) +
  geom_boxplot() +
  ylab("Percent Dissolved Oxygen")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))

#pH
#build model
lmer.pH= lme(pH~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
             data=EAB_Ter_H2o_aqcwd,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.pH),col="darkgray")
#normal
summary(lmer.pH)
#downstream and Grand significantly higher
#visualize
ggplot(EAB_Ter_H2o_aqcwd, aes(x=Watershed, y=pH, color=Gap_location)) +
  geom_point(size=3) +
  ylab("pH")+
  xlab("Watershed")+
  labs(color = "Gap Location")+
  scale_color_manual(values=Gap_location_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#summary
EAB_Ter_H2o_aqcwd %>%
  group_by(Gap_location) %>%
  get_summary_stats(pH, type = "mean_se")
EAB_Ter_H2o_aqcwd %>%
  group_by(Watershed) %>%
  get_summary_stats(pH, type = "mean_se")

#NTU
#build model
lmer.NTU= lme(NTU~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
             data=EAB_Ter_H2o_aqcwd,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.NTU),col="darkgray")
#not normal
shapiro.test(EAB_Ter_H2o_aqcwd$NTU)
#W = 0.65131, p-value = 2.017e-09, not normal
hist(EAB_Ter_H2o_aqcwd$NTU)
#log transform
EAB_Ter_H2o_aqcwd$log10NTU<-log10(EAB_Ter_H2o_aqcwd$NTU+1)
shapiro.test(EAB_Ter_H2o_aqcwd$log10NTU)
#W = 0.97511, p-value = 0.3948, normal use this transformation
lmer.log10NTU= lme(log10NTU~Watershed+Gap_location+Total_Live_BA, random=~1|Stream,
              data=EAB_Ter_H2o_aqcwd,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10NTU),col="darkgray")
#normal
summary(lmer.log10NTU)
#nothing significant
EAB_Ter_H2o_aqcwd %>%
  group_by(Watershed) %>%
  get_summary_stats(NTU, type = "mean_se")

#######################
#Water chemistry summary
#%DO significantly lower downstream
#pH significantly higher downstream and in grand river watershed
#total basal area does not influence
###########################

#Now try with year since gap formation as predictor rather than watershed
#Create time since gap formation variable
EAB_Ter_H2o_aqcwd$YearSinceGF<-2016-EAB_Ter_H2o_aqcwd$Year_Gapformation
#Make sure the input variables aren't correlated
hetcor(EAB_Ter_H2o_aqcwd$YearSinceGF,EAB_Ter_H2o_aqcwd$Stream,EAB_Ter_H2o_aqcwd$Gap_location,
       EAB_Ter_H2o_aqcwd$Total_Live_BA)
#Stream is highly correlated to year of gap formation. So only include year of gap formation in these models, don't include the mixed effect

#Water temp
lm.wty= lm(Water_temp~YearSinceGF+Gap_location+Total_Live_BA, data=EAB_Ter_H2o_aqcwd)
plot(lm.wty,1)
#no fitted pattern and red line approximately at 0
plot(lm.wty, 3)
#no fitted pattern and red line approximately horizontal
plot(lm.wty, 2)
#looks skewed at each end. test normality of water temp
shapiro.test(EAB_Ter_H2o_aqcwd$Water_temp)
range(EAB_Ter_H2o_aqcwd$Water_temp)
#W = 0.83708, p-value = 3.517e-06, not normal, log transform
#range is biologically significant
EAB_Ter_H2o_aqcwd$log10Water_temp<-log10(EAB_Ter_H2o_aqcwd$Water_temp+1)
lm.wtylog= lm(log10Water_temp~YearSinceGF+Gap_location+Total_Live_BA, data=EAB_Ter_H2o_aqcwd)
plot(lm.wtylog,1)
#no fitted pattern and red line approximately at 0
plot(lm.wtylog, 3)
#no fitted pattern and red line approximately horizontal
plot(lm.wtylog, 2)
#Didn't really improve. So use non transformed model for interpretability
summary(lm.wty)
summary(lm.wtylog)
#not significant

#mS_cm
range(EAB_Ter_H2o_aqcwd$mScm)
#0.129-1.1

#X.DO
sort(EAB_Ter_H2o_aqcwd$X.DO)
#ranges from 57.3 to 139.0, so biologically significant range
shapiro.test(EAB_Ter_H2o_aqcwd$X.DO)
#W = 0.96855, p-value = 0.1925, normal
#build model
lm.doy= lm(X.DO~YearSinceGF+Gap_location+Total_Live_BA, data=EAB_Ter_H2o_aqcwd)
plot(lm.doy,1)
#no fitted pattern and red line approximately at 0
plot(lm.doy, 3)
#no fitted pattern and red line approximately horizontal
plot(lm.doy, 2)
#normal
summary(lm.doy)
#antly lower downstream
#visualize
ggplot(EAB_Ter_H2o_aqcwd, aes(x=Year_Gapformation, y=X.DO, color=Gap_location)) +
  geom_point() +
  ylab("Percent Dissolved Oxygen")+
  xlab("Year of Gap Formation")+
  scale_color_manual(values=Gap_location_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))

#pH
sort(EAB_Ter_H2o_aqcwd$pH)
#range 7.95 to 8.68 - not expected to be biologically relevant

#NTU
sort(EAB_Ter_H2o_aqcwd$NTU)
#0.5-110 range
boxplot.stats(EAB_Ter_H2o_aqcwd$NTU)$out
#Three outliers - 46, 53, and 110
hist(EAB_Ter_H2o_aqcwd$NTU)
#skewed, so log tranform
EAB_Ter_H2o_aqcwd$log10NTU<-log10(EAB_Ter_H2o_aqcwd$NTU+1)
boxplot.stats(EAB_Ter_H2o_aqcwd$log10NTU)$out
#no outliers
#build model
lm.ntuy= lm(log10NTU~YearSinceGF+Gap_location+Total_Live_BA, data=EAB_Ter_H2o_aqcwd)
plot(lm.ntuy,1)
#no fitted pattern and red line approximately at 0
plot(lm.ntuy, 3)
#no fitted pattern and red line approximately horizontal
plot(lm.ntuy, 2)
#normal
summary(lm.ntuy)
#nothing significant
EAB_Ter_H2o_aqcwd %>%
  group_by(Watershed) %>%
  get_summary_stats(NTU, type = "mean_se")

#########################
#Aquatic Leaf Litter
#########################

#upload 
#get proportions
EAB_ALL<-read.csv("~/Documents/MSU/Research/Surveying/Leaf_Litter/EAB_Aq_Leaf_Litter.csv", sep = ",", header = T )
names(EAB_ALL)
EAB_ALL_cast<-cast(EAB_ALL, Stream_name + Gap_number + Gap_location + Date ~ Taxa)
EAB_ALL_cast[is.na(EAB_ALL_cast)] <- 0
EAB_ALL_abund<-decostand(EAB_ALL_cast, method = "total")
write.csv(EAB_ALL_abund,'ALL_abund.csv')
EAB_ALL_cast$Totalmass<-rowSums(EAB_ALL_cast[5:26])
write.csv(EAB_ALL_cast,'EAB_ALL_cast.csv')

#add results to paired sites dataset
#combine poplar and aspen
#rename "birch" "hophornbeam
#Fix here in case any sample labels flipped in field notes
#Sessions Creek 9/4/2016 GIBII and GIBIII switched
#Sessions Creek 10/23/2016 GIIBIII and GIBII switched
#Spring Creek 10/22/2016 GIBIII and GIBII switched

EAB_Paired_ALL<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/ALL_dataset.csv", sep = ",", header = T )
#Merge with other dataset
EAB_Ter_H2o_aqcwd_ALL<-merge(EAB_Ter_H2o_aqcwd,EAB_Paired_ALL, by=c("Stream","Gap_location","Date"))

#Create summary table
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_Richness, type = "mean_se")
#"ALL_prop_Ash"
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Ash, type = "mean_se")
#"ALL_prop_Aspen" 
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Aspen, type = "mean_se")
#"ALL_prop_Basswood"
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Basswood, type = "mean_se")
#"ALL_prop_Beech"  
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Beech, type = "mean_se")
#"ALL_prop_Hophornbeam"
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Hophornbeam, type = "mean_se")
#"ALL_prop_Black_Cherry"  
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Black_Cherry, type = "mean_se")
#"ALL_prop_Black_Walnut"
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Black_Walnut, type = "mean_se")
#"ALL_prop_Buckthorn" 
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Buckthorn, type = "mean_se")
#"ALL_prop_Dogwood"  
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Dogwood, type = "mean_se")
#"ALL_prop_Elm" 
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Elm, type = "mean_se")
#"ALL_prop_Grass"      
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Grass, type = "mean_se")
#"ALL_prop_Maple"     
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Maple, type = "mean_se")
#"ALL_prop_Oak" 
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Oak, type = "mean_se")
EAB_Ter_H2o_aqcwd_ALL %>%
  get_summary_stats(ALL_prop_Oak, type = "mean_se")
#"ALL_prop_Unknown" 
EAB_Ter_H2o_aqcwd_ALL %>%
  group_by(Watershed) %>%
  get_summary_stats(ALL_prop_Unknown, type = "mean_se")

#model two ways:
#first all~watershed+gaplocation+1|stream linear mixed effect
#if linear assumption isn't met, use Mixed effects logistic regression for binary outcome
#next correlation analysis with terrestrial live equivalent

#mixed models
#ALL richness
#build model
lmer.AR= lme(ALL_Richness~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Ter_H2o_aqcwd_ALL,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.AR),col="darkgray")
#normal
summary(lmer.AR)
#nothing significant
EAB_Ter_H2o_aqcwd_ALL$Watershed = forcats::fct_rev(factor(EAB_Ter_H2o_aqcwd_ALL$Watershed))
ggplot(EAB_Ter_H2o_aqcwd_ALL, aes(x=Gap_location, y=ALL_Richness, fill=Watershed)) +
  geom_boxplot() +
  ylab("Aquatic Leaf Litter Richness")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))

#all prop ash
#build model
lmer.ALLA= lme(ALL_prop_Ash~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Ter_H2o_aqcwd_ALL,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.ALLA),col="darkgray")
#not normal, try binomial
glmer.ALLA<-glmer(ALL_prop_Ash~Watershed+Gap_location+1|Stream,family=binomial,weights=ALL_Total_weight,
             data=EAB_Ter_H2o_aqcwd_ALL)
#fail to converge, so see if Stream is adding to model
glm.ALLA<-glm(ALL_prop_Ash~Watershed+Gap_location,family=binomial,weights=ALL_Total_weight,
              data=EAB_Ter_H2o_aqcwd_ALL)
anova(glmer.ALLA,glm.ALLA)
#0.9997, so adding stream doesn't impact model
#check assumptions
hist(residuals(glm.ALLA),col="darkgray")
#normal
summary(glm.ALLA)
#intercept significant
ggplot(EAB_Ter_H2o_aqcwd_ALL, aes(x=Gap_location, y=ALL_prop_Ash, fill=Watershed)) +
  geom_boxplot() +
  ylab("Aquatic Leaf Litter Ash")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#all prop oak
#build model
lmer.ALLO= lme(ALL_prop_Oak~Watershed+Gap_location, random=~1|Stream,
               data=EAB_Ter_H2o_aqcwd_ALL,
               method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.ALLO),col="darkgray")
#not normal, try binomial
glmer.ALLO<-glmer(ALL_prop_Oak~Watershed+Gap_location+1|Stream,family=binomial,weights=ALL_Total_weight,
                  data=EAB_Ter_H2o_aqcwd_ALL)
#fail to converge, so see if Stream is adding to model
glm.ALLO<-glm(ALL_prop_Oak~Watershed+Gap_location,family=binomial,weights=ALL_Total_weight,
              data=EAB_Ter_H2o_aqcwd_ALL)
anova(glmer.ALLO,glm.ALLO)
#0.158, so adding stream doesn't impact model
#check assumptions
hist(residuals(glm.ALLO),col="darkgray")
#normal
summary(glm.ALLO)
#intercept significant

##########################
#Summary: aquatic leaf litter composition not altered by watershed or gap location
########################

#Now try with year since gap formation rather than watershed
#CAn't use mixed model because stream and year are correlated, so use linar model

#ALL richness
#build model
lm.ARy= lm(ALL_Richness~YearSinceGF+Gap_location, data=EAB_Ter_H2o_aqcwd_ALL)
plot(lm.ARy,1)
#no fitted pattern and red line approximately at 0
plot(lm.ARy, 3)
#no fitted pattern and red line approximately horizontal
plot(lm.ARy, 2)
#looks normal
summary(lm.ARy)
#not significant
#visualize
EAB_Ter_H2o_aqcwd_ALL$Year_Gapformation_factor<-as.factor(EAB_Ter_H2o_aqcwd_ALL$Year_Gapformation)
ggplot(EAB_Ter_H2o_aqcwd_ALL, aes(x=Gap_location, y=ALL_Richness, fill=Year_Gapformation_factor)) +
  geom_boxplot() +
  ylab("Aquatic Leaf Litter Richness")+
  xlab("Gap Location")+
  scale_fill_manual(values=YGF_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))

#all prop ash
#use logistic model for proportional data
glm.ALLA<-glm(ALL_prop_Ash~YearSinceGF+Gap_location,binomial,data=EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLA)
#not significant, residenual deviance < df
plot(glm.ALLA,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLA,3)
#positive trend, try with log transformation
glm.ALLAlog<-glm(log10(ALL_prop_Ash+1)~YearSinceGF+Gap_location,binomial,data=EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLAlog)
#not significant, residenual deviance < df
plot(glm.ALLAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLAlog,3)
#positive trend, try with binary outcomes
#convert to presence absence of ash
EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Ash<-EAB_Ter_H2o_aqcwd_ALL$ALL_prop_Ash
EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Ash[EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Ash > 0] <- 1 
glm.ALLAbi<- glm(ALL_pres_Ash~YearSinceGF+Gap_location, family = binomial,data = EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLAbi)
#not significant
plot(glm.ALLAbi,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLAbi,3)
#positive skew
#since nothing is fixing it and nothing is significant, use the original porportional outcomes result

#all prop oak
#build model
#all prop ash
#use logistic model for proportional data
glm.ALLO<-glm(ALL_prop_Oak~YearSinceGF+Gap_location,binomial,data=EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLO)
#not significant, residenual deviance < df
plot(glm.ALLO,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLO,3)
#close to flat, a bit positive
plot(glm.ALLO,2)
#a bit skewed, try log transformation
glm.ALLOlog<-glm(log10(ALL_prop_Oak+1)~YearSinceGF+Gap_location,binomial,data=EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLOlog)
#not significant, residenual deviance < df
plot(glm.ALLOlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLOlog,3)
#slight positive, but not bad
plot(glm.ALLOlog,2)
#some skew, so convert to presence absence of oak
EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Oak<-EAB_Ter_H2o_aqcwd_ALL$ALL_prop_Oak
EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Oak[EAB_Ter_H2o_aqcwd_ALL$ALL_pres_Oak > 0] <- 1 
glm.ALLObi<- glm(ALL_pres_Oak~YearSinceGF+Gap_location, family = binomial,data = EAB_Ter_H2o_aqcwd_ALL)
summary(glm.ALLObi)
#not significant
plot(glm.ALLObi,1)
#no fitted pattern and red line approximately at 0
plot(glm.ALLObi,3)
#flat
plot(glm.ALLObi,2)
#very skewed
#since nothing is fixing it and nothing is significant, use the original porportional outcomes result

#Move to correlation analysis

##ALL richness
#test normality
shapiro.test(EAB_Ter_H2o_aqcwd_ALL$ALL_Richness)
#W = 0.94056, p-value = 0.009839, not normal use non parametric
#Terrestrial richness
names(EAB_Ter_H2o_aqcwd_ALL)
EAB_Ter_H2o_aqcwd_ALL$Terrestrial_Live_Richness<-rowSums(EAB_Ter_H2o_aqcwd_ALL[c(44,45,47,49,51,53:55,57,58,60,61,63,65,67,69,70)]> 0)

cor.AR.TerrRich<- cor.test(EAB_Ter_H2o_aqcwd_ALL$ALL_Richness,
                           EAB_Ter_H2o_aqcwd_ALL$Terrestrial_Live_Richness,
                                method = "spearman")
cor.AR.TerrRich
#not significant

#all prop ash
#not normal so use spearman
EAB_Ter_H2o_aqcwd_ALL$prop_Ash_live_ba<-EAB_Ter_H2o_aqcwd_ALL$BAperha_Ash_Live/EAB_Ter_H2o_aqcwd_ALL$Total_Live_BA
cor.Alla.Terrpropash<- cor.test(EAB_Ter_H2o_aqcwd_ALL$ALL_prop_Ash,
                           EAB_Ter_H2o_aqcwd_ALL$prop_Ash_live_ba,
                           method = "spearman")
cor.Alla.Terrpropash

#all prop oak
#not normal so use spearman
EAB_Ter_H2o_aqcwd_ALL$prop_Oak_live_ba<-EAB_Ter_H2o_aqcwd_ALL$BAperha_Oak_Live/EAB_Ter_H2o_aqcwd_ALL$Total_Live_BA
cor.Allo.Terrpropo<- cor.test(EAB_Ter_H2o_aqcwd_ALL$ALL_prop_Oak,
                                EAB_Ter_H2o_aqcwd_ALL$prop_Oak_live_ba,
                                method = "spearman")
cor.Allo.Terrpropo

#######################
#Summary: aquatic leaf litter not correlated with adjacent terrestrial litter
##########################

###################
#Microbes
##################

###Microbial communities
#Use QIIME2 to process sequences.
#In mapping file, use corrected gap/basket labels (see above with leaf litter)

#Upload map
EAB_16S_map <- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_Map_Filtered_R_rev.csv", header=T)
#Combine map with other env variables
EAB_16S_map$Source<-revalue(EAB_16S_map$Source, c("AquaticLeafLitter"="ALL", "LiveLeaves"="LL", "TerrestrialLitter"="TL"))
EAB_16S_map$YearSinceGF<-2016-EAB_16S_map$Year_Gapformation
#for just aquatic
EAB_Ter_H2o_cwd_ALL_mics<-merge(EAB_16S_map, EAB_Ter_H2o_aqcwd_ALL,by=c("Stream","Watershed","Gap_location","Date","Days_since_start","Year_Gapformation","YearSinceGF"),all.x=TRUE)
EAB_Ter_H2o_cwd_ALL_mics<-subset(EAB_Ter_H2o_cwd_ALL_mics, SampleID!="NA")
row.names(EAB_Ter_H2o_cwd_ALL_mics)<-EAB_Ter_H2o_cwd_ALL_mics$SampleID
#create new Gap variable in EAB_Ter_H2o_cwd_ALL_mics with Y,N to mirror terrestrial leaf litte rmicrobes
EAB_Ter_H2o_cwd_ALL_mics$Gap<-EAB_Ter_H2o_cwd_ALL_mics$Gap_location
EAB_Ter_H2o_cwd_ALL_mics$Gap<-revalue(EAB_Ter_H2o_cwd_ALL_mics$Gap, c("Downstream"="N", "Gap"="Y", "Upstream"="N", "Woods"="N"))

#####Do this later###
######

#Create rarefaction graph with alpha diversity values
EAB_Shannon<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/rev microbes/EAB_shannon.csv", header=T)
EAB_Shannon$sample.id<- gsub("\\-", "_", EAB_Shannon$sample.id)
rownames(EAB_Shannon)<-EAB_Shannon[,1]
EAB_Shannon[,1]<-NULL
EAB_Shannon<-data.frame(t(EAB_Shannon))
EAB_Shannon$iteration<-NULL
str(EAB_Shannon)
EAB_Sh_m<-melt(EAB_Shannon, id.vars="sequences.per.sample")
any(is.na(EAB_Sh_m))
EAB_Sh_c_m<-cast(EAB_Sh_m, variable ~ sequences.per.sample, fun.aggregate=mean, na.omit=TRUE)
EAB_Sh_c_m$Calculation<-rep("mean",147)
EAB_Sh_c_var<-cast(EAB_Sh_m, variable ~ sequences.per.sample, fun.aggregate=var)
EAB_Sh_c_var$Calculation<-rep("variance",147)
EAB_Sh_c_m_var<-rbind(EAB_Sh_c_m,EAB_Sh_c_var)
names(EAB_Sh_c_m_var)[names(EAB_Sh_c_m_var)=="variable"] <- "SampleID"
EAB_Sh_c_m_var$Calculation<-as.factor(EAB_Sh_c_m_var$Calculation)
EAB_Sh_c_m_var<-as.data.frame(EAB_Sh_c_m_var)
EAB_Sh_c_m_var_m<-melt(EAB_Sh_c_m_var, id.vars=c("Calculation","SampleID"))
EAB_Sh_c_m_var_c<-cast(EAB_Sh_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
EAB_Sh_map <-merge(EAB_16S_map, EAB_Sh_c_m_var_c, by="SampleID")
names(HC_Sh_map)[names(HC_Sh_map)=="variable"] <- "Sequences_per_sample"
any(is.na(HC_Sh_map))
HC_Sh_map$Sequences_per_sample<-as.numeric(as.character(HC_Sh_map$Sequences_per_sample))
HC_Sh_map$Source<-factor(revalue(HC_Sh_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_Sh_map_sum_m <- summarySE(HC_Sh_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_Sh_map_sum_v <- summarySE(HC_Sh_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_Sh_map_sum_v$StandDev<-sqrt(HC_Sh_map_sum_v$variance)
HC_Sh_map_sum_v$StandEr<-HC_Sh_map_sum_v$StandDev/sqrt(HC_Sh_map_sum_v$N)
HC_Sh_sum_m_sd<-merge(HC_Sh_map_sum_m,HC_Sh_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_Sh_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Shannon H' diversity") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Faith's diversity
HC_Faith_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Faith_q2.txt", sep="\t", header=T)
rownames(HC_Faith_q2)<-HC_Faith_q2[,1]
HC_Faith_q2[,1]<-NULL
HC_Faith<-data.frame(t(HC_Faith_q2))
HC_Faith$iteration<-NULL
str(HC_Faith)
HC_F_m<-melt(HC_Faith, id.vars="sequences.per.sample")
HC_F_c_m<-cast(HC_F_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_F_c_m$Calculation<-rep("mean",175)
HC_F_c_var<-cast(HC_F_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_F_c_var$Calculation<-rep("variance",175)
HC_F_c_m_var<-rbind(HC_F_c_m,HC_F_c_var)
names(HC_F_c_m_var)[names(HC_F_c_m_var)=="variable"] <- "SampleID"
HC_F_c_m_var$Calculation<-as.factor(HC_F_c_m_var$Calculation)
HC_F_c_m_var<-as.data.frame(HC_F_c_m_var)
HC_F_c_m_var_m<-melt(HC_F_c_m_var, id.vars=c("Calculation","SampleID"))
HC_F_c_m_var_c<-cast(HC_F_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_F_map <-merge(Hunt_Creek_16S_map, HC_F_c_m_var_c, by="SampleID")
names(HC_F_map)[names(HC_F_map)=="variable"] <- "Sequences_per_sample"
HC_F_map$Sequences_per_sample<-as.numeric(as.character(HC_F_map$Sequences_per_sample))
HC_F_map$Source<-factor(revalue(HC_F_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_F_map_sum_m <- summarySE(HC_F_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_F_map_sum_v <- summarySE(HC_F_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_F_map_sum_v$StandDev<-sqrt(HC_F_map_sum_v$variance)
HC_F_map_sum_v$StandEr<-HC_F_map_sum_v$StandDev/sqrt(HC_F_map_sum_v$N)
HC_F_sum_m_sd<-merge(HC_F_map_sum_m,HC_F_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_F_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Faith's phylogenetic diversity") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Observed OTUs
HC_Ob_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Obs_q2.txt", sep="\t", header=T)
rownames(HC_Ob_q2)<-HC_Ob_q2[,1]
HC_Ob_q2[,1]<-NULL
HC_Ob<-data.frame(t(HC_Ob_q2))
HC_Ob$iteration<-NULL
str(HC_Ob)
HC_o_m<-melt(HC_Ob, id.vars="sequences.per.sample")
HC_o_c_m<-cast(HC_o_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_o_c_m$Calculation<-rep("mean",175)
HC_o_c_var<-cast(HC_o_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_o_c_var$Calculation<-rep("variance",175)
HC_o_c_m_var<-rbind(HC_o_c_m,HC_o_c_var)
names(HC_o_c_m_var)[names(HC_o_c_m_var)=="variable"] <- "SampleID"
HC_o_c_m_var$Calculation<-as.factor(HC_o_c_m_var$Calculation)
HC_o_c_m_var<-as.data.frame(HC_o_c_m_var)
HC_o_c_m_var_m<-melt(HC_o_c_m_var, id.vars=c("Calculation","SampleID"))
HC_o_c_m_var_c<-cast(HC_o_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_o_map <-merge(Hunt_Creek_16S_map, HC_o_c_m_var_c, by="SampleID")
names(HC_o_map)[names(HC_o_map)=="variable"] <- "Sequences_per_sample"
HC_o_map$Sequences_per_sample<-as.numeric(as.character(HC_o_map$Sequences_per_sample))
HC_o_map$Source<-factor(revalue(HC_o_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_o_map_sum_m <- summarySE(HC_o_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_o_map_sum_v <- summarySE(HC_o_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_o_map_sum_v$StandDev<-sqrt(HC_o_map_sum_v$variance)
HC_o_map_sum_v$StandEr<-HC_o_map_sum_v$StandDev/sqrt(HC_o_map_sum_v$N)
HC_o_sum_m_sd<-merge(HC_o_map_sum_m,HC_o_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_o_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean observed OTUs") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Chao1
HC_Ch_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Chao1_q2.txt", sep="\t", header=T)
rownames(HC_Ch_q2)<-HC_Ch_q2[,1]
HC_Ch_q2[,1]<-NULL
HC_Ch<-data.frame(t(HC_Ch_q2))
HC_Ch$iteration<-NULL
str(HC_Ch)
HC_c_m<-melt(HC_Ch, id.vars="sequences.per.sample")
HC_c_c_m<-cast(HC_c_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_c_c_m$Calculation<-rep("mean",175)
HC_c_c_var<-cast(HC_c_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_c_c_var$Calculation<-rep("variance",175)
HC_c_c_m_var<-rbind(HC_c_c_m,HC_c_c_var)
names(HC_c_c_m_var)[names(HC_c_c_m_var)=="variable"] <- "SampleID"
HC_c_c_m_var$Calculation<-as.factor(HC_c_c_m_var$Calculation)
HC_c_c_m_var<-as.data.frame(HC_c_c_m_var)
HC_c_c_m_var_m<-melt(HC_c_c_m_var, id.vars=c("Calculation","SampleID"))
HC_c_c_m_var_c<-cast(HC_c_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_c_map <-merge(Hunt_Creek_16S_map, HC_c_c_m_var_c, by="SampleID")
names(HC_c_map)[names(HC_c_map)=="variable"] <- "Sequences_per_sample"
HC_c_map$Sequences_per_sample<-as.numeric(as.character(HC_c_map$Sequences_per_sample))
HC_c_map$Source<-factor(revalue(HC_c_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_c_map_sum_m <- summarySE(HC_c_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_c_map_sum_v <- summarySE(HC_c_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_c_map_sum_v$StandDev<-sqrt(HC_c_map_sum_v$variance)
HC_c_map_sum_v$StandEr<-HC_c_map_sum_v$StandDev/sqrt(HC_c_map_sum_v$N)
HC_c_sum_m_sd<-merge(HC_c_map_sum_m,HC_c_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_c_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Chao1 estimator") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#####Created facetted graph for manuscript
HC_Sh_sum_m_sd$alpha<-rep("A. Shannon Diveristy",50)
HC_F_sum_m_sd$alpha<-rep("B. Faith's Phylogenetic Diveristy",50)
HC_o_sum_m_sd$alpha<-rep("C. Observed OTUs",50)
HC_c_sum_m_sd$alpha<-rep("D. Chao1 Estimator",50)
HC_sum_m_sd<-rbind(HC_Sh_sum_m_sd,HC_F_sum_m_sd,HC_o_sum_m_sd,HC_c_sum_m_sd)
ggplot(HC_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of Reads Sampled") +
  ylab("Mean Alpha Diversity Value (+/- SEM)") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Times New Roman", size=36, margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(family="Times New Roman", size=36, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(family="Times New Roman", size=18),
        axis.text.y = element_text(family="Times New Roman", size=18),
        legend.position = 0, 
        strip.text.x = element_text(family="Times New Roman", face="bold", size = 24)) +
  scale_color_manual(values=five_col_vec_CBBaSH,                      
                     limits=c("Carcass","Biofilm","B. brunneicolor","S. mutata","H. flavescens"),
                     labels=c("Carcass","Biofilm",expression(paste(italic('B. brunneicolor'))),expression(paste(italic('S. mutata'))),expression(paste(italic('H. flavescens'))))) +
  facet_wrap(~alpha,scales="free_y")

###Back here again
#########

#Upload OTU table
EAB_16S_OTU<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_OTUs_rev.csv", header=T)
#Format data frame so the OTU is row name
row.names(EAB_16S_OTU)<-EAB_16S_OTU[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_OTU$ID<-NULL
EAB_16S_OTU_t<-t(EAB_16S_OTU)
EAB_16S_OTU_t<-data.frame(EAB_16S_OTU_t, check.names = FALSE)

#merge OTUs
EAB_Paired_TA_OTU<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_OTU_t, by=0, no.dups=T)

#Combine with environmental variables
EAB_16S_OTU_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_OTU_t, by=0)
#delete OTUs not found in any samples
names(EAB_16S_OTU_map)
EAB_16S_OTU<-EAB_16S_OTU_map[,182:ncol(EAB_16S_OTU_map)]
EAB_16S_OTU<-EAB_16S_OTU[, colSums(EAB_16S_OTU != 0) > 0]
#3512 OTUs found in all samples

#Upload family level taxonomy
EAB_16S_f<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_f_rev.csv", header=T)
#Format data frame so the family is row name
row.names(EAB_16S_f)<-EAB_16S_f[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_f$ID<-NULL
EAB_16S_f_t<-t(EAB_16S_f)
EAB_16S_f_t<-data.frame(EAB_16S_f_t, check.names = FALSE)

#Combine with environmental variables
EAB_16S_f_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_f_t, by=0)
#delete families not found in 
names(EAB_16S_f_map)
EAB_16S_f<-EAB_16S_f_map[,182:ncol(EAB_16S_f_map)]
EAB_16S_f<-EAB_16S_f_aq[, colSums(EAB_16S_f_aq != 0) > 0]
sort(colSums(EAB_16S_f))
#most common family is Bacillaceae
stat.desc(EAB_16S_f$`k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae`)
#260 +/- 26
EAB_16S_f_map_GF<-EAB_16S_f_map$Gap
EAB_16S_f_map_W<-EAB_16S_f_map$Watershed
EAB_16S_f_map_SW<-paste(EAB_16S_f_map$Source,EAB_16S_f_map$Watershed)
EAB_16S_f_map_YGF<-as.factor(EAB_16S_f_map$Year_Gapformation)

#Upload phylogenetic diversity
EAB_16S_fpd<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_fpd_rev.csv", header=T)
#Format data frame so the ID is row name
row.names(EAB_16S_fpd)<-EAB_16S_fpd[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_fpd$X<-NULL
#Combine with environmental variables
EAB_16S_fpd_map<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_fpd, by=0)

####Terrestrial connection microbes - all habitats

#Model watershed + gap + source linear model

#Faith's PD for all habitats
range(EAB_16S_fpd_map$faith_pd)
#1.96-34.04
#build model
lm.FPD= lm(faith_pd~YearSinceGF+Gap+Source, data=EAB_16S_fpd_map)
#check assumptions
plot(lm.FPD,1)
#flat
plot(lm.FPD,3)
#mostly flat, slight positive
plot(lm.FPD,2)
#skew towards top, try log transformed
lm.FPDlog= lm(log10(faith_pd)~YearSinceGF+Gap+Source, data=EAB_16S_fpd_map)
#check assumptions
plot(lm.FPDlog,1)
#flat
plot(lm.FPDlog,3)
#flat
plot(lm.FPDlog,2)
#much better
summary(lm.FPDlog)
EAB_16S_fpd_map %>% emmeans_test(log10faith_pd ~ Source, p.adjust.method = "bonferroni",
                                  model = lmer.FPD)
#TLL significantly greater than LL and ALL
#visualize
ggplot(EAB_16S_fpd_map, aes(x=Source, y=faith_pd)) +
  geom_point(size=3) +
  ylab("Microbial Community Phylogenetic Diversity")+
  scale_x_discrete(labels=c("ALL" = "Aquatic Leaf Litter", "LL" = "Live Leaves",
                              "TL" = "Terrestrial Leaf Litter"),
                   limits=c("LL","TL","ALL"))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#Summary
EAB_16S_fpd_map %>%
  group_by(Source) %>%
  get_summary_stats(faith_pd, type = "mean_se")

#Chao1
#Upload chao1
EAB_16S_cha<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_chao_rev.csv", header=T)
#Format data frame so the ID is row name
row.names(EAB_16S_cha)<-EAB_16S_cha[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_cha$X<-NULL
#Combine with environmental variables
EAB_16S_cha_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_cha, by=0)

#start modelling
range(EAB_16S_cha_map$chao1)
#27-646
#build model
lm.cha= lm(chao1~YearSinceGF+Gap+Source, data=EAB_16S_cha_map)
#check assumptions
plot(lm.cha,1)
#flat
plot(lm.cha,3)
#positive, try log transformed
lm.chalog= lm(log10(chao1)~YearSinceGF+Gap+Source, data=EAB_16S_cha_map)
#check assumptions
plot(lm.chalog,1)
#flat
plot(lm.chalog,3)
#flat, sight positive
plot(lm.chalog,2)
#looks okay, slight skew
summary(lm.chalog)
#source significant
EAB_16S_cha_map %>% emmeans_test(log10chao1 ~ Source, p.adjust.method = "bonferroni",
                                 model = lmer.logcha)
#TLL greater than ALL and LL
#visualize
ggplot(EAB_16S_cha_map, aes(x=Source, y=chao1)) +
  geom_point(size=3) +
  ylab("Microbial Community Chao 1 Richness")+
  scale_x_discrete(labels=c("ALL" = "Aquatic Leaf Litter", "LL" = "Live Leaves",
                            "TL" = "Terrestrial Leaf Litter"),
                   limits=c("LL","TL","ALL"))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#Summary
EAB_16S_cha_map %>%
  group_by(Source) %>%
  get_summary_stats(chao1, type = "mean_se")

#Now look for community level differences using PERMANOVA and adonis

#Upload unifrac distance matrix
EAB_16S_uni<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_WUnifrac_rev.csv", header = T, row.names=1)

#Combine distance matrix with environmental variables
row.names(EAB_16S_uni)
row.names(EAB_Ter_H2o_cwd_ALL_mics)
EAB_16S_uni_map<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_uni, by="row.names")
row.names(EAB_16S_uni_map)<-EAB_16S_uni_map[,1]
EAB_16S_uni_map<-EAB_16S_uni_map[,-c(1)]
names(EAB_16S_uni_map)
EAB_16S_uni<-as.matrix(EAB_16S_uni_map[,c(181:ncol(EAB_16S_uni_map))])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env<-EAB_16S_uni_map[,1:180]
EAB_16S_uni_env$Gap<-revalue(EAB_16S_uni_env$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_GF<-EAB_16S_uni_env$Gap
EAB_16S_uni_env_GFS<-as.factor(paste(EAB_16S_uni_env$Gap.y,EAB_16S_uni_env$Source))
EAB_16S_uni_env_S<-as.factor(EAB_16S_uni_env$Source)
EAB_16S_uni_env_S<-revalue(EAB_16S_uni_env_S, c("ALL"="Aquatic Leaf Litter", "LL"="Live Leaves",
                                                "TL"="Terrestrial Leaf Litter"))
EAB_16S_uni_env_S<-factor(EAB_16S_uni_env_S, levels = c("Aquatic Leaf Litter", "Terrestrial Leaf Litter", "Live Leaves"))
EAB_16S_uni_env_W<-as.factor(EAB_16S_uni_env$Watershed)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(EAB_16S_uni) ~ (Source+YearSinceGF+Gap)^2, data=EAB_16S_uni_env,
       permutations=9999)
#source significant

#Visualize via nmds
EAB_Paired_Mic_NMDS<-metaMDS(as.dist(EAB_16S_uni))
#stress 0.12

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mic_NMDS)

#NMDS plot for source
ordiplot(EAB_Paired_Mic_NMDS, type="n")
with(EAB_Paired_Mic_NMDS, points(EAB_Paired_Mic_NMDS, display="sites", col=source_col_vec[EAB_16S_uni_env_S], pch=19))
with(EAB_Paired_Mic_NMDS, legend("topleft", legend=levels(EAB_16S_uni_env_S), bty="n", col=source_col_vec, pch=19, pt.bg=source_col_vec))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "Aquatic Leaf Litter"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#33a02c", show.groups = "Live Leaves"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#b15928", show.groups = "Terrestrial Leaf Litter"))

#NMDSplot for watershed
ordiplot(EAB_Paired_Mic_NMDS, type="n")
with(EAB_Paired_Mic_NMDS, points(EAB_Paired_Mic_NMDS, display="sites", col=Watershed_col_vec[EAB_16S_uni_env_W], pch=19))
with(EAB_Paired_Mic_NMDS, legend("topleft", legend=levels(EAB_16S_uni_env_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_W, kind="se", conf=0.95, lwd=2, col="#af8dc3", show.groups = "Clinton"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_W, kind="se", conf=0.95, lwd=2, col="#f7f7f7", show.groups = "Grand"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#indicator species analysis for watershed
EAB_Mic_Com_W_indic<-signassoc(EAB_16S_f, cluster=EAB_16S_f_map_W,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_W_indic_sig<-subset(EAB_Mic_Com_W_indic, psidak<=0.05)
#7 indicator families for watershed

#indicator species analysis for Year since gap formation
EAB_Mic_Com_YGF_indic<-data.frame(signassoc(EAB_16S_f, cluster=EAB_16S_f_map_YGF,  mode=0, alternative = "two.sided",control = how(nperm=999)))
EAB_Mic_Com_YGF_indic_sig<-subset(EAB_Mic_Com_YGF_indic, psidak<=0.05)
#6 indicator families for year of gap formation

#indicator species analysis for watershedxleaf type
EAB_Mic_Com_SW_indic<-signassoc(EAB_16S_f, cluster=EAB_16S_f_map_SW,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_SW_indic_sig<-subset(EAB_Mic_Com_SW_indic, psidak<=0.05)
#30 indicator families for watershedxleaf type
#k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;__ indicates gaps in clinton
EAB_16S_f_map %>%
  group_by(Source,Watershed) %>%
  get_summary_stats("k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Frankiaceae", type = "mean_se")
EAB_16S_f_map %>%
  group_by(Source,Watershed) %>%
  get_summary_stats("k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae", type = "mean_se")

#now move on to venn diagrams to find unique OTUs
#create venn diagram for shared OTUs among terrestrial/aquatic gap/woods.
#first list OTU names found in each group
names(EAB_Paired_TA_OTU)
EAB_16S_OTU_map_ag_venn<-aggregate(EAB_Paired_TA_OTU[181:ncol(EAB_Paired_TA_OTU)], 
                                   by=list(Gap=EAB_Paired_TA_OTU$Gap,Source=EAB_Paired_TA_OTU$Source),
                                   FUN=sum)
str(EAB_16S_OTU_map_ag_venn)
EAB_16S_OTU_map_ag_venn$GapSource<-as.factor(paste(EAB_16S_OTU_map_ag_venn$Gap, EAB_16S_OTU_map_ag_venn$Source))
levels(EAB_16S_OTU_map_ag_venn$GapSource)
row.names(EAB_16S_OTU_map_ag_venn)<-EAB_16S_OTU_map_ag_venn$GapSource
EAB_16S_OTU_map_ag_venn$GapSource<-NULL
EAB_16S_OTU_map_ag_venn$Gap<-NULL
EAB_16S_OTU_map_ag_venn$Source<-NULL
row.names(EAB_16S_OTU_map_ag_venn)
NAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["N ALL",] > 0]
YAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["Y ALL",] > 0]
NLiveLeaves<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["N LL",] > 0]
YLiveLeaves<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["Y LL",] > 0]
NTerrestrialLitter<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["N TL",] > 0]
YTerrestrialLitter<- colnames(EAB_16S_OTU_map_ag_venn)[EAB_16S_OTU_map_ag_venn["Y TL",] > 0]
NoGap<-unique(c(NAquaticLeafLitter,NLiveLeaves,NTerrestrialLitter))
YesGap<-unique(c(YAquaticLeafLitter,YLiveLeaves,YTerrestrialLitter))
AquaticLeafLitter<-unique(c(NAquaticLeafLitter,YAquaticLeafLitter))
LiveLeaves<-unique(c(NLiveLeaves,YLiveLeaves))
TerrestrialLitter<-c(NTerrestrialLitter,YTerrestrialLitter)
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
quartz()
GapVen<-venn_diagram4(YAquaticLeafLitter, YLiveLeaves, YTerrestrialLitter, NoGap,
                      "AL Gap  ", "TLL Gap", "TL Gap", "For",
                      colors=gapven_col_vec)
#No OTU's common between all gap sample types and not found in gap
#find relative abundances
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
quartz()
GapVen3<-venn_diagram3(YAquaticLeafLitter,YLiveLeaves,YTerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
ForVen3<-venn_diagram3(NAquaticLeafLitter,NLiveLeaves,NTerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
quartz()
WoodsVen<-venn_diagram4(NAquaticLeafLitter, NLiveLeaves, NTerrestrialLitter, YesGap,
                        "AL For", "TLL For", "TL For", "Gap",
                        colors=woodsven_col_vec)
#No OTUs common between all gap sample types and not found in gap
quartz()
Ven3<-venn_diagram3(AquaticLeafLitter,LiveLeaves,TerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
venn.diagram(x = list(YLiveLeaves,NLiveLeaves),category.names = c("Gap" , "Forest"),
  filename = 'LiveLivesGapForest.png',output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,
  scaled=FALSE)
venn.diagram(x = list(YTerrestrialLitter,NTerrestrialLitter),
             category.names = c("Gap" , "Forest"),filename = 'TerrestrialLitterGapForest.png',
             output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,scaled=FALSE)
venn.diagram(x = list(YAquaticLeafLitter,NAquaticLeafLitter),
             category.names = c("Gap" , "Forest"),filename = 'AquaticLitterGapForest.png',
             output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,scaled=FALSE)
EAB_16S_OTU_map_ag_venn_ALL<-aggregate(EAB_Paired_TA_OTU[178:ncol(EAB_Paired_TA_OTU)], 
                                   by=list(Gap=EAB_Paired_TA_OTU$Gap_location,
                                           Source=EAB_Paired_TA_OTU$Source),FUN=sum)
str(EAB_16S_OTU_map_ag_venn_ALL)
EAB_16S_OTU_map_ag_venn_ALL<-subset(EAB_16S_OTU_map_ag_venn_ALL,Source=="ALL")
row.names(EAB_16S_OTU_map_ag_venn_ALL)<-EAB_16S_OTU_map_ag_venn_ALL$Gap
EAB_16S_OTU_map_ag_venn_ALL$Gap<-NULL
EAB_16S_OTU_map_ag_venn_ALL$Source<-NULL
row.names(EAB_16S_OTU_map_ag_venn_ALL)
USAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_ALL)[EAB_16S_OTU_map_ag_venn_ALL["Upstream",] > 0]
GAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_ALL)[EAB_16S_OTU_map_ag_venn_ALL["Gap",] > 0]
DSAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_ALL)[EAB_16S_OTU_map_ag_venn_ALL["Downstream",] > 0]

venn.diagram(x = list(USAquaticLeafLitter,GAquaticLeafLitter,DSAquaticLeafLitter),
             category.names=c("Upstream","Gap","Downstream"),
             filename = 'AquaticLitterGL.png',
             output=TRUE,fill=Gap_location_col_vec,cat.cex =1.5,cex=1.5,scaled=FALSE)
EAB_16S_OTU_map_ag_venn_w<-aggregate(EAB_Paired_TA_OTU[181:ncol(EAB_Paired_TA_OTU)], 
                                   by=list(Watershed=EAB_Paired_TA_OTU$Watershed,Source=EAB_Paired_TA_OTU$Source),
                                   FUN=sum)
str(EAB_16S_OTU_map_ag_venn_w)
EAB_16S_OTU_map_ag_venn_w$WSource<-as.factor(paste(EAB_16S_OTU_map_ag_venn_w$Watershed, EAB_16S_OTU_map_ag_venn_w$Source))
levels(EAB_16S_OTU_map_ag_venn_w$WSource)
row.names(EAB_16S_OTU_map_ag_venn_w)<-EAB_16S_OTU_map_ag_venn_w$WSource
EAB_16S_OTU_map_ag_venn_w$WSource<-NULL
EAB_16S_OTU_map_ag_venn_w$Watershed<-NULL
EAB_16S_OTU_map_ag_venn_w$Source<-NULL
row.names(EAB_16S_OTU_map_ag_venn_w)
CAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Clinton ALL",] > 0]
GAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Grand ALL",] > 0]
KAquaticLeafLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Kalamazoo ALL",] > 0]
CLiveLeaves<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Clinton LL",] > 0]
GLiveLeaves<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Grand LL",] > 0]
KLiveLeaves<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Kalamazoo LL",] > 0]
CTerrestrialLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Clinton TL",] > 0]
GTerrestrialLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Grand TL",] > 0]
KTerrestrialLitter<- colnames(EAB_16S_OTU_map_ag_venn_w)[EAB_16S_OTU_map_ag_venn_w["Kalamazoo TL",] > 0]
Clinton<-unique(c(CAquaticLeafLitter,CLiveLeaves,CTerrestrialLitter))
Grand<-unique(c(GAquaticLeafLitter,GLiveLeaves,GTerrestrialLitter))
Kalamazoo<-unique(c(KAquaticLeafLitter,KLiveLeaves,KTerrestrialLitter))
CG<-unique(c(Clinton,Grand))
CK<-unique(c(Clinton,Kalamazoo))
GK<-unique(c(Grand,Kalamazoo))
quartz()
ClintonVen<-venn_diagram4(CAquaticLeafLitter, CLiveLeaves, CTerrestrialLitter, GK,
                      "AL Cl  ", "TLL Cl", "TL Cl", "Gr+Ka")
#No OTU's common between all Clinton sample types and not found in Grand or Kalamazoo
quartz()
ClintonVen3<-venn_diagram3(CAquaticLeafLitter,CLiveLeaves,CTerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
quartz()
GrandVen<-venn_diagram4(GAquaticLeafLitter, GLiveLeaves, GTerrestrialLitter, CK,
                          "AL Gr  ", "TLL Gr", "TL Gr", "Cl+Ka")
#1 OTU common between all Grand sample types and not found in Clinton or Kalamazoo
#find relative abundance
GrandVen$`AL Gr  _TLL Gr_TL Gr`
#"0e1737f92c9ab36c66f868ee6f762e5c"
EAB_Paired_TA_OTU %>%
  group_by(Watershed,Source) %>%
  get_summary_stats("0e1737f92c9ab36c66f868ee6f762e5c", type = "mean_se")
#<1% relative abundance
#find taxonomy
#upload OTU taxonomy
EAB_paired_OTU_tax<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_taxonomy.csv", header=T)
EAB_paired_OTU_tax[EAB_paired_OTU_tax$FeatureID %in% GrandVen$`AL Gr  _TLL Gr_TL Gr`,]
#k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Kineosporiaceae
quartz()
GrandVen3<-venn_diagram3(GAquaticLeafLitter,GLiveLeaves,GTerrestrialLitter, "Aquatic Litter",
                           "Terrestrial Litter","Live Leaves",colors=source_col_vec)
quartz()
KalamazooVen<-venn_diagram4(KAquaticLeafLitter, KLiveLeaves, KTerrestrialLitter, CG,
                        "AL Ka  ", "TLL Ka", "TL Ka", "Cl+Gr")
#1 OTU common between all Kalamazoo sample types and not found in Clinton or Grand
#find relative abundance
KalamazooVen$`AL Ka  _TLL Ka_TL Ka`
#"720ffeffab6263dc2410b57601b70bcd"
EAB_Paired_TA_OTU %>%
  group_by(Watershed,Source) %>%
  get_summary_stats("720ffeffab6263dc2410b57601b70bcd", type = "mean_se")
#<1% relative abundance
EAB_paired_OTU_tax[EAB_paired_OTU_tax$FeatureID %in% KalamazooVen$`AL Ka  _TLL Ka_TL Ka`,]
#k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae; g__; s__
quartz()
KalamazooVen3<-venn_diagram3(KAquaticLeafLitter,KLiveLeaves,KTerrestrialLitter, "Aquatic Litter",
                         "Terrestrial Litter","Live Leaves",colors=source_col_vec)
quartz()
WatVen3<-venn_diagram3(Clinton,Grand,Kalamazoo, "Clinton",
                       "Grand","Kalamazoo",colors=Watershed_col_vec)
quartz()
WatLLVen3<-venn_diagram3(CLiveLeaves,GLiveLeaves,KLiveLeaves, "Clinton",
                       "Grand","Kalamazoo",colors=Watershed_col_vec)
quartz()
WatTLVen3<-venn_diagram3(CTerrestrialLitter,GTerrestrialLitter,KTerrestrialLitter, "Clinton",
                         "Grand","Kalamazoo",colors=Watershed_col_vec)
quartz()
WatALLVen3<-venn_diagram3(CAquaticLeafLitter,GAquaticLeafLitter,KAquaticLeafLitter, "Clinton",
                         "Grand","Kalamazoo",colors=Watershed_col_vec)
EAB_16S_OTU_map_ag_venn_wg<-aggregate(EAB_Paired_TA_OTU[178:ncol(EAB_Paired_TA_OTU)], 
                                     by=list(Watershed=EAB_Paired_TA_OTU$Watershed,Gap=EAB_Paired_TA_OTU$Gap),
                                     FUN=sum)
str(EAB_16S_OTU_map_ag_venn_wg)
EAB_16S_OTU_map_ag_venn_wg$WG<-as.factor(paste(EAB_16S_OTU_map_ag_venn_wg$Watershed, EAB_16S_OTU_map_ag_venn_wg$Gap))
levels(EAB_16S_OTU_map_ag_venn_wg$WG)
row.names(EAB_16S_OTU_map_ag_venn_wg)<-EAB_16S_OTU_map_ag_venn_wg$WG
EAB_16S_OTU_map_ag_venn_wg$WG<-NULL
EAB_16S_OTU_map_ag_venn_wg$Watershed<-NULL
EAB_16S_OTU_map_ag_venn_wg$Gap<-NULL
row.names(EAB_16S_OTU_map_ag_venn_wg)
CN<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Clinton N",] > 0]
GN<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Grand N",] > 0]
KN<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Kalamazoo N",] > 0]
CY<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Clinton Y",] > 0]
GY<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Grand Y",] > 0]
KY<- colnames(EAB_16S_OTU_map_ag_venn_wg)[EAB_16S_OTU_map_ag_venn_wg["Kalamazoo Y",] > 0]
quartz()
WatForVen<-venn_diagram4(CN, GN, KN, YesGap,
                          "Cl For  ", "Gr For", "Ka For", "Gap")
#9 OTU's common between all watersheds forests and not found in the gap
#find relative abundances
WatForUniqueASVs<-WatForVen$`Cl For  _Gr For_Ka For`
RAWatForUniqueASVs<-EAB_Paired_TA_OTU %>%
  group_by(Watershed,Gap) %>%
  get_summary_stats(c(all_of(WatForUniqueASVs)), type = "mean_se")
#none over 1%
#taxonomy
WatForUTax<-EAB_paired_OTU_tax[EAB_paired_OTU_tax$FeatureID %in% WatForUniqueASVs,]
WatForUTax$Taxon
quartz()
ForWatVen3<-venn_diagram3(CN,GN,KN, "Clinton",
                           "Grand","Kalamazoo",colors=Watershed_col_vec)
quartz()
WatGapVen<-venn_diagram4(CY, GY, KY, NoGap,
                        "Cl Gap  ", "Gr Gap", "Ka Gap", "For")
#8 OTUs common between all watershed gaps and not found in forest
#find relative abundance
WatGapUniqueASVs<-WatGapVen$`Cl Gap  _Gr Gap_Ka Gap`
RAWatGapUniqueASVs<-EAB_Paired_TA_OTU %>%
  group_by(Watershed,Gap) %>%
  get_summary_stats(c(all_of(WatGapUniqueASVs)), type = "mean_se")
#none over 1%
WatGapUTax<-EAB_paired_OTU_tax[EAB_paired_OTU_tax$FeatureID %in% WatGapUniqueASVs,]
WatGapUTax$Taxon
quartz()
GapWatVen3<-venn_diagram3(CY,GY,KY, "Clinton",
                         "Grand","Kalamazoo",colors=Watershed_col_vec)
venn.diagram(x = list(CY,CN),category.names = c("Gap" , "Forest"),
             filename = 'ClintonGapForest.png',output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,
             scaled=FALSE)
venn.diagram(x = list(GY,GN),category.names = c("Gap" , "Forest"),
             filename = 'GrandGapForest.png',output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,
             scaled=FALSE)
venn.diagram(x = list(KY,KN),category.names = c("Gap" , "Forest"),
             filename = 'KalamazooGapForest.png',output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,
             scaled=FALSE)

#################################
#Summary all sources microbes
#terrestrial leaf litter have significantly higher fpd than aquatic leaf litter and live leaves
#terrestrial leaf litter higher chao1 richness compared to aquatic leaf litter and live leaves
#source influence community structure
#No OTUs only found in gap or forest, not the other
##############################

#Now model each habitat separately

#Live leaves
#Faith's PD for live leaves
EAB_16S_fpd_map_LL<-subset(EAB_16S_fpd_map, Source=="LL")
range(EAB_16S_fpd_map_LL$faith_pd)
#5.5-15.3
#build model
lm.FPD.LL= lm(faith_pd~YearSinceGF+Gap, data=EAB_16S_fpd_map_LL)
#check assumptions
plot(lm.FPD.LL,1)
#no pattern
plot(lm.FPD.LL,3)
#fairly flat
plot(lm.FPD.LL,2)
#looks good, no outliers
summary(lm.FPD.LL)
#intercept significant

#Chao1 for live leaves
EAB_16S_cha_map_LL<-subset(EAB_16S_cha_map, Source=="LL")
#start modelling
range(EAB_16S_cha_map_LL$chao1)
#27-151
#build model
lm.cha.ll= lm(chao1~YearSinceGF+Gap,data=EAB_16S_cha_map_LL)
#check assumptions
plot(lm.cha.ll,1)
#flat
plot(lm.cha.ll,3)
#flat
plot(lm.cha.ll,2)
#skewed, try log transforming
lm.cha.lllog= lm(log10(chao1)~YearSinceGF+Gap,data=EAB_16S_cha_map_LL)
#check assumptions
plot(lm.cha.lllog,1)
#flat
plot(lm.cha.lllog,3)
#slight down
plot(lm.cha.lllog,2)
#better, use log tranformed
summary(lm.cha.lllog)
#nothing significant

#Now look for community level differences using PERMANOVA and adonis in LL

#subset unifrac map for live leaves
EAB_16S_uni_map_LL<-subset(EAB_16S_uni_map,Source=="LL")
names(EAB_16S_uni_map_LL)
EAB_16S_uni_LL<-as.matrix(EAB_16S_uni_map_LL[,c(233:241)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_LL<-EAB_16S_uni_map_LL[,1:180]
EAB_16S_uni_env_LL$Gap<-revalue(EAB_16S_uni_env_LL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_LL_GF<-EAB_16S_uni_env_LL$Gap
EAB_16S_uni_env_LL_W<-as.factor(EAB_16S_uni_env_LL$Watershed)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(EAB_16S_uni_LL) ~ (YearSinceGF+Gap)^2, data=EAB_16S_uni_env_LL,
       permutations=999)
#nothing significant

#Visualize via nmds
EAB_Paired_Mic_NMDS_LL<-metaMDS(as.dist(EAB_16S_uni_LL))
#stress=0.09
#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mic_NMDS_LL)

#NMDSplot for watershed
ordiplot(EAB_Paired_Mic_NMDS_LL, type="n")
with(EAB_Paired_Mic_NMDS_LL, points(EAB_Paired_Mic_NMDS_LL, display="sites", col=Watershed_col_vec[EAB_16S_uni_env_LL_W], pch=19))
with(EAB_Paired_Mic_NMDS_LL, legend("topleft", legend=levels(EAB_16S_uni_env_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_W, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Grand"))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#NMDSplot for gap
ordiplot(EAB_Paired_Mic_NMDS_LL, type="n")
with(EAB_Paired_Mic_NMDS_LL, points(EAB_Paired_Mic_NMDS_LL, display="sites", col=GF_col_vec[EAB_16S_uni_env_LL_GF], pch=19))
with(EAB_Paired_Mic_NMDS_LL, legend("topleft", legend=levels(EAB_16S_uni_env_GF), bty="n", col=GF_col_vec, pch=19, pt.bg=GF_col_vec))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_GF, kind="se", conf=0.95, lwd=2, col="#66c2a5", show.groups = "Forest"))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_GF, kind="se", conf=0.95, lwd=2, col="#fdae61", show.groups = "Gap"))

#################################
#Summary live leave microbes
#FPd and chao not altered
#nothing influences community structure
##############################

#Terrestrial litter
EAB_16S_f_map_TL<-subset(EAB_16S_f_map,Source=="TL")
names(EAB_16S_f_map_TL)
EAB_16S_f_TL<-EAB_16S_f_map_TL[182:ncol(EAB_16S_f_map_TL)]
sort(colSums(EAB_16S_f_TL))
EAB_16S_f_TL %>%
  get_summary_stats('k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae', type = "mean_se")

#Faith's PD for terrestrial litter
EAB_16S_fpd_map_TL<-subset(EAB_16S_fpd_map, Source=="TL")
range(EAB_16S_fpd_map_TL$faith_pd)
#9.8-34
#build model
lm.FPD.TL= lm(faith_pd~YearSinceGF+Gap,data=EAB_16S_fpd_map_TL)
#check assumptions
plot(lm.FPD.TL,1)
#flat
plot(lm.FPD.TL,3)
#flat
plot(lm.FPD.TL,2)
#looks good
summary(lm.FPD.TL)
#nothing significant

#Chao1 for terrestrial leaves
EAB_16S_cha_map_TL<-subset(EAB_16S_cha_map, Source=="TL")
#start modelling
range(EAB_16S_cha_map_TL$chao1)
#82-646
#build model
lm.cha.tl= lm(chao1~YearSinceGF+Gap,data=EAB_16S_cha_map_TL)
#check assumptions
plot(lm.cha.tl,1)
#flat
plot(lm.cha.tl,3)
#flat
plot(lm.cha.tl,2)
#not bad
summary(lm.cha.tl)
#nothing significant

#Now look for community level differences using PERMANOVA and adonis in TL

#subset unifrac map for terrestrial litter
EAB_16S_uni_map_TL<-subset(EAB_16S_uni_map,Source=="TL")
names(EAB_16S_uni_map_TL)
EAB_16S_uni_TL<-as.matrix(EAB_16S_uni_map_TL[,c(242:253)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_TL<-EAB_16S_uni_map_TL[,1:180]
EAB_16S_uni_env_TL$Gap<-revalue(EAB_16S_uni_env_TL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_TL_GF<-EAB_16S_uni_env_TL$Gap
EAB_16S_uni_env_TL_W<-as.factor(EAB_16S_uni_env_TL$Watershed)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(EAB_16S_uni_TL) ~ (YearSinceGF+Gap)^2, data=EAB_16S_uni_env_TL,
       permutations=9999)
#nothing significant

#################################
#Summary terrestrial leaf litter microbes
#FPd and chao not altered
#nothing significant influence community structure
##############################

#Aquatic Leaf Litter

#Faith's PD for aquatic leaf litter
EAB_16S_fpd_map_ALL<-subset(EAB_16S_fpd_map, Source=="ALL")
range(EAB_16S_fpd_map_ALL$faith_pd)
#2-22
#build model
str(EAB_16S_fpd_map_ALL)
lm.FPD.ALL= lm(faith_pd~YearSinceGF+Gap_location+ALL_Richness,data=EAB_16S_fpd_map_ALL)
#check assumptions
plot(lm.FPD.ALL,1)
#flat
plot(lm.FPD.ALL,3)
#slight positive
plot(lm.FPD.ALL,2)
#some skew, try log transformed
lm.FPD.ALLlog= lm(log10(faith_pd)~YearSinceGF+Gap_location+ALL_Richness,data=EAB_16S_fpd_map_ALL)
#check assumptions
plot(lm.FPD.ALLlog,1)
#flat
plot(lm.FPD.ALLlog,3)
#flat
plot(lm.FPD.ALLlog,2)
#much better
summary(lm.FPD.ALLlog)
#nothing significant
ggplot(EAB_16S_fpd_map_ALL, aes(x=Gap_location, y=faith_pd)) +
  geom_boxplot() +
  ylab("Microbial Community Phylogenetic Diversity")+
  xlab("Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
EAB_16S_fpd_map_ALL$Watershed = forcats::fct_rev(factor(EAB_16S_fpd_map_ALL$Watershed))
ggplot(EAB_16S_fpd_map_ALL, aes(x=Gap_location, y=faith_pd, fill=Watershed)) +
  geom_boxplot() +
  ylab("Microbial Community Phylogenetic Diversity")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))

#Chao1 for aquatic leaf litter
EAB_16S_cha_map_ALL<-subset(EAB_16S_cha_map, Source=="ALL")
#start modelling
range(EAB_16S_cha_map_ALL$chao1)
#32-345
#build model
lm.cha.all= lm(chao1~YearSinceGF+Gap_location+ALL_Richness,data=EAB_16S_cha_map_ALL)
#check assumptions
plot(lm.cha.all,1)
#flat
plot(lm.cha.all,3)
#flat
plot(lm.cha.all,2)
#some skew, log transform
lm.cha.alllog= lm(log10(chao1)~YearSinceGF+Gap_location+ALL_Richness,data=EAB_16S_cha_map_ALL)
#check assumptions
plot(lm.cha.alllog,1)
#flat
plot(lm.cha.alllog,3)
#flat
plot(lm.cha.alllog,2)
#better
summary(lm.cha.alllog)
#nothing significant

#Now look for community level differences using PERMANOVA and adonis in all

#subset unifrac map for terrestrial litter
EAB_16S_uni_map_ALL<-subset(EAB_16S_uni_map,Source=="ALL")
names(EAB_16S_uni_map_ALL)
EAB_16S_uni_ALL<-as.matrix(EAB_16S_uni_map_ALL[,c(181:232)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_ALL<-EAB_16S_uni_map_ALL[,1:180]
EAB_16S_uni_env_ALL$Gap<-revalue(EAB_16S_uni_env_ALL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_ALL_GF<-EAB_16S_uni_env_ALL$Gap
EAB_16S_uni_env_ALL_W<-as.factor(EAB_16S_uni_env_ALL$Watershed)
EAB_16S_uni_env_ALL_GL<-factor(EAB_16S_uni_env_ALL$Gap_location)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(EAB_16S_uni_ALL) ~ (YearSinceGF+Gap_location+ALL_Richness)^2, data=EAB_16S_uni_env_ALL,
       permutations=9999)

#################################
#Summary aquatic leaf litter microbes
#FPd and chao not altered
#all richness influences community structure
##############################

################
#Macroinvertebrates
###################

#Upload dataset
#Macroinvertebrate communities
#This uploaded dataset has already had dates switched (see above leaf litter)
Paired_Macroinvertebrates<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/Paired_Macroinvertebrates.csv", sep = ",", header = T )
#Confirm header names
names(Paired_Macroinvertebrates)
#Create sample ID name by combining environmental variables into one name in new column
Paired_Macroinvertebrates$SampleID<-factor(paste(Paired_Macroinvertebrates$Stream, Paired_Macroinvertebrates$Date, Paired_Macroinvertebrates$Gap, Paired_Macroinvertebrates$Gap_location))
#Create new column that combines taxonomic variables to family level
Paired_Macroinvertebrates$Taxonomy<-as.factor(paste(Paired_Macroinvertebrates$Class, Paired_Macroinvertebrates$Order, Paired_Macroinvertebrates$Family, sep=""))
#taxonomy levels
levels(Paired_Macroinvertebrates$Taxonomy)

#Convert to community matrix based on sample ID
EAB_M_Matrix <- dcast(as.data.table(Paired_Macroinvertebrates), Stream + Date + Gap + Gap_location ~ Taxonomy, value.var = "Count", fun.aggregate =sum)
#Find total invertebrate #
sum(EAB_M_Matrix[,5:ncol(EAB_M_Matrix)])
#Find most abundant taxa
sort(colSums(EAB_M_Matrix[,5:ncol(EAB_M_Matrix)]))
#most abundant taxa elmidae with 519 individuals
stat.desc(EAB_M_Matrix$`InsectaColeopteraElmidae`)
#make relative abundances
str(EAB_M_Matrix)
EAB_M_MatrixRA<-EAB_M_Matrix
EAB_M_MatrixRA[,5:ncol(EAB_M_MatrixRA)]<-data.frame(make_relative(as.matrix(EAB_M_MatrixRA[,5:ncol(EAB_M_MatrixRA)])))
EAB_M_MatrixRA[is.na(EAB_M_MatrixRA)] <- 0
str(EAB_M_MatrixRA)
sort(colMeans(EAB_M_MatrixRA))
stat.desc(EAB_M_MatrixRA$InsectaColeopteraElmidae)

#combine eab_Paired_ALL_Stand with macs
EAB_Paired_Macs_Stand<-merge(EAB_Ter_H2o_aqcwd_ALL, EAB_M_Matrix, by=c("Stream","Gap_location","Date"))
#split up into environmental and community
names(EAB_Paired_Macs_Stand)
EAB_Paired_Macs_Stand_com<-EAB_Paired_Macs_Stand[,172:ncol(EAB_Paired_Macs_Stand)]
EABMacroComtotbray0Stand<-as.matrix(bray0(EAB_Paired_Macs_Stand_com))
EAB_Paired_Macs_Stand_env<-EAB_Paired_Macs_Stand[,1:171]
EAB_Paired_Macs_Stand_GA<-EAB_Paired_Macs_Stand_env$Gap_Area
EAB_Paired_Macs_Stand_GA_cat<-as.factor(EAB_Paired_Macs_Stand_GA)
EAB_Paired_Macs_Stand_YGF<-EAB_Paired_Macs_Stand_env$Year_Gapformation
EAB_Paired_Macs_Stand_YGF_cat<-as.factor(EAB_Paired_Macs_Stand_YGF)
EAB_Paired_Macs_Stand_GL<-EAB_Paired_Macs_Stand_env$Gap_location
EAB_Paired_Macs_Stand_W<-EAB_Paired_Macs_Stand_env$Watershed

#analyze community
adonis2(as.dist(EABMacroComtotbray0Stand) ~ YearSinceGF*Gap_location,
       data=EAB_Paired_Macs_Stand_env,permutations=9999)
#Year of gap formation and stream significant

#Watershed
EAB_Inverts_Com_W_indic<-signassoc(EAB_Paired_Macs_Stand_com, cluster=EAB_Paired_Macs_Stand_W,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Inverts_Com_W_indic_sig<-subset(EAB_Inverts_Com_W_indic, psidak<=0.05)
#10 families indicate watershed:
#Clinton: elmidae, ephemerellidae, heptageniidae, nemouridae
#Grand: Gammaridae
#Kalamazoo: Psephenidae, Baetidae, Taeniopterygidae, hydropsychidae, rhacophilidae

#Year of gap formation
EAB_Inverts_Com_YGF_indic<-data.frame(signassoc(EAB_Paired_Macs_Stand_com, cluster=EAB_Paired_Macs_Stand_YGF_cat,  mode=0, alternative = "two.sided",control = how(nperm=9999)))
EAB_Inverts_Com_YGF_indic_sig<-subset(EAB_Inverts_Com_YGF_indic, psidak<=0.05)
#9 families indicate year of gap formation:
#Grand/2011: Gammaridae
#Kalamazoo/2014: Elmidae, Psephenidae, Chironomidae, Baetidae, Ephemerellidae, Heptageniidae, hydropsychidae, rhacophilidae

#Gap Location
EAB_Inverts_Com_GL_indic<-signassoc(EAB_Paired_Macs_Stand_com, cluster=EAB_Paired_Macs_Stand_GL,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Inverts_Com_GL_indic_sig<-subset(EAB_Inverts_Com_GL_indic, psidak<=0.05)
#0 families indicate gap location

#so work with these 9 families for quantitative responses:
#Determine which are in all watersheds
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaColeopteraElmidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaColeopteraElmidae", type = "mean_se")
#in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaColeopteraPsephenidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaColeopteraPsephenidae", type = "mean_se")
#missing in frayer, sessions and spring
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaDipteraChironomidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaDipteraChironomidae", type = "mean_se")
#in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraBaetidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaEphemeropteraBaetidae", type = "mean_se")
#missing in spring creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraEphemerellidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaEphemeropteraEphemerellidae", type = "mean_se")
#missing in stoney creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraHeptageniidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaEphemeropteraHeptageniidae", type = "mean_se")
#missing in spring creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaTrichopteraHydropsychidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaTrichopteraHydropsychidae", type = "mean_se")
#found in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaTrichopteraRhyacophilidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("InsectaTrichopteraRhyacophilidae", type = "mean_se")
#missing in all but sessions and seven mile
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("MalacostracaAmphipodaGammaridae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("MalacostracaAmphipodaGammaridae", type = "mean_se")
#in all streams
#Elmidae, Chironomidae, Hydropsychidae, and Gammaridae found in all streams

#so work with hydropsychidae, Gammaridae, Chironomidae and Elmidae

#NMDS analysis
EAB_Paired_Mac_NMDS<-metaMDS(as.dist(EABMacroComtotbray0Stand))

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mac_NMDS)

#NMDS plot for watershed
ordiplot(EAB_Paired_Mac_NMDS, type="n")
with(EAB_Paired_Mac_NMDS, points(EAB_Paired_Mac_NMDS, display="sites", col=Watershed_col_vec[EAB_Paired_Macs_Stand_W], pch=19))
with(EAB_Paired_Mac_NMDS, legend("topleft", legend=levels(EAB_Paired_Macs_Stand_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_W, kind="se", conf=0.95, lwd=2, col="#af8dc3", show.groups = "Clinton"))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_W, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Grand"))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#NMDS plot for Year of gap formation
ordiplot(EAB_Paired_Mac_NMDS, type="n")
with(EAB_Paired_Mac_NMDS, points(EAB_Paired_Mac_NMDS, display="sites", col=YGF_col_vec[EAB_Paired_Macs_Stand_YGF_cat], pch=19))
with(EAB_Paired_Mac_NMDS, legend("topleft", legend=levels(EAB_Paired_Macs_Stand_YGF_cat), bty="n", col=YGF_col_vec, pch=19, pt.bg=YGF_col_vec))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#67000d", show.groups = "2006"))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#cb181d", show.groups = "2008"))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#fb6a4a", show.groups = "2011"))
with(EAB_Paired_Mac_NMDS, ordiellipse(EAB_Paired_Mac_NMDS, EAB_Paired_Macs_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#fcbba1", show.groups = "2014"))

#Now test models on responses of richness, diversity, and populations of indicator taxa

#Richness
#Calculate richness
EAB_Paired_Macs_Stand$Mac_Richness<-rowSums(EAB_Paired_Macs_Stand_com > 0)
stat.desc(EAB_Paired_Macs_Stand$Mac_Richness)
#Check for normality
shapiro.test(EAB_Paired_Macs_Stand$Mac_Richness)
#W = 0.97191, p-value = 0.2336, normal
range(EAB_Paired_Macs_Stand$Mac_Richness)

#Model richness as mixed model like microbes
#build model
lm.mr= lm(Mac_Richness~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
#check assumptions
plot(lm.mr,1)
#flat
plot(lm.mr,3)
#flat to slightly positive
plot(lm.mr,2)
#looks okay
summary(lm.mr)
#year since gap formation and downstream lowest

#visualize
ggplot(EAB_Paired_Macs_Stand, aes(x=Year_Gapformation_factor, y=Mac_Richness, color=Gap_location)) +
  geom_point(position=position_jitter(h=0.2, w=0.2),size=3)+
  ylab("Macroinvertebrate Family Richness")+
  xlab("Year of Gap Formation")+
  scale_color_manual(values=Gap_location_col_vec)+
  labs(color = "Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#summarize
EAB_Paired_Macs_Stand %>%
  group_by(Gap_location) %>%
  get_summary_stats("Mac_Richness", type = "mean_se")
#1 Upstream     Mac_Richness    18  7.56 0.908
#2 Gap          Mac_Richness    18  7    0.932
#3 Downstream   Mac_Richness    18  5.28 0.956
ggplot(EAB_Paired_Macs_Stand, aes(x=Gap_location, y=Mac_Richness)) +
  geom_boxplot()+
  ylab("Macroinvertebrate Family Richness")+
  xlab("Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
EAB_Paired_Macs_Stand$Watershed= forcats::fct_rev(factor(EAB_Paired_Macs_Stand$Watershed))
ggplot(EAB_Paired_Macs_Stand, aes(x=Gap_location, y=Mac_Richness, fill=Watershed)) +
  geom_boxplot() +
  ylab("Macroinvertebrate Family Richness")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#change to early, mid late for presentations
ggplot(EAB_Paired_Macs_Stand, aes(x=Gap_location, y=Mac_Richness, fill=Watershed)) +
  geom_boxplot() +
  ylab("Macroinvertebrate Family Richness")+
  xlab("Gap Location")+
  scale_fill_manual(values=Watershed_col_vec, name="EAB Invasion", labels = c("Late", "Mid", "Early"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
ggplot(EAB_Paired_Macs_Stand, aes(x=Gap_location, y=Mac_Richness, fill=Year_Gapformation_factor)) +
  geom_boxplot() +
  ylab("Macroinvertebrate Family Richness")+
  xlab("Gap Location")+
  scale_fill_manual(values=YGF_col_vec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#Model diversity as mixed model like microbes
#Calculate diversity
EAB_Paired_Macs_Stand$Mac_Simp<-diversity(EAB_Paired_Macs_Stand_com, index="invsimpson")
#Change inf to 1
EAB_Paired_Macs_Stand$Mac_Simp[!is.finite(EAB_Paired_Macs_Stand$Mac_Simp)] <- 1
#build model
lm.md= lm(Mac_Simp~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
#check assumptions
plot(lm.md,1)
#flat
plot(lm.md,3)
#slight positive
plot(lm.md,2)
#looks good
summary(lm.md)
#intercept, years, and downstream significant
#visualize
ggplot(EAB_Paired_Macs_Stand, aes(x=Watershed, y=Mac_Simp, color=Gap_location)) +
  geom_point(position=position_jitter(h=0.1, w=0.1),size=3)+
  ylab("Macroinvertebrate Inverse Simpson's Diversity")+
  scale_color_manual(values=Gap_location_col_vec)+
  labs(color = "Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))
#summarize
EAB_Paired_Macs_Stand %>%
  group_by(Gap_location) %>%
  get_summary_stats("Mac_Simp", type = "mean_se")
#1 Upstream     Mac_Simp    18  4.26 0.42 
#2 Gap          Mac_Simp    18  3.71 0.392
#3 Downstream   Mac_Simp    18  3.09 0.416

#Model total abundance
#Create total macroinvertebrate abundance column
names(EAB_Paired_Macs_Stand)
EAB_Paired_Macs_Stand$Total_Macs<-rowSums(EAB_Paired_Macs_Stand[172:218])
lm.tma= lm(Total_Macs~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
#check assumptions
plot(lm.tma,1)
#flat
plot(lm.tma,3)
#slight positive
plot(lm.tma,2)
#pretty bad skew towards the end, try transforming
lm.tmalog= lm(log10(Total_Macs+1)~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
#check assumptions
plot(lm.tmalog,1)
#flat
plot(lm.tmalog,3)
#flat
plot(lm.tmalog,2)
#better
summary(lm.tmalog)
#intercept and time significant
EAB_Paired_Macs_Stand %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("Total_Macs", type = "mean_se")
#Because total abundances different for each stream, use relative abundances for taxa specific analyses

#so work with these 5 families for quantitative responses:
#Elmidae, Chironomidae, Hydropsychidae, and Gammaridae found in all streams
#Use the same relative response modeling as the leaf packs
#combine relative abundances with environmental data
#combine eab_Paired_ALL_Stand with macs
EAB_Paired_Macs_StandRA<-merge(EAB_Ter_H2o_aqcwd_ALL, EAB_M_MatrixRA, by=c("Stream","Gap_location","Date"))

#Elmidae
#use logistic model for proportional data
glm.ElmRA<-glm(InsectaColeopteraElmidae~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.ElmRA)
#not significant, residenual deviance < df
plot(glm.ElmRA,1)
#no fitted pattern and red line approximately at 0
plot(glm.ElmRA,3)
#flat
plot(glm.ElmRA,2)
#Skew at end, try with log tranformed
glm.ElmRAlog<-glm(log10(InsectaColeopteraElmidae+1)~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.ElmRAlog)
#not significant, residenual deviance < df
plot(glm.ElmRAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.ElmRAlog,3)
#flat
plot(glm.ElmRAlog,2)
#worse skew, go with non transformed

#Chironomidae
glm.ChRA<-glm(InsectaDipteraChironomidae~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.ChRA)
#not significant, residenual deviance < df
plot(glm.ChRA,1)
#no fitted pattern and red line approximately at 0
plot(glm.ChRA,3)
#slight positive
plot(glm.ChRA,2)
#Skew at end, try with log tranformed
glm.ChRAlog<-glm(log10(InsectaDipteraChironomidae+1)~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.ChRAlog)
#not significant, residenual deviance < df
plot(glm.ChRAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.ChRAlog,3)
#flat
plot(glm.ChRAlog,2)
#worse skew, go with non transformed

#Hydropsychidae
glm.HyRA<-glm(InsectaTrichopteraHydropsychidae~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.HyRA)
#not significant, residenual deviance < df
plot(glm.HyRA,1)
#no fitted pattern and red line approximately at 0
plot(glm.HyRA,3)
#flat
plot(glm.HyRA,2)
#Skew at both ends, try with log tranformed
glm.HyRAlog<-glm(log10(InsectaTrichopteraHydropsychidae+1)~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.HyRAlog)
#not significant, residenual deviance < df
plot(glm.HyRAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.HyRAlog,3)
#flat
plot(glm.HyRAlog,2)
#worse skew, go with non transformed

#Gammaridae
glm.GaRA<-glm(MalacostracaAmphipodaGammaridae~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.GaRA)
#not significant, residenual deviance < df
plot(glm.GaRA,1)
#no fitted pattern and red line approximately at 0
plot(glm.GaRA,3)
#flat
plot(glm.GaRA,2)
#Skew at both ends, try with log tranformed
glm.GaRAlog<-glm(log10(MalacostracaAmphipodaGammaridae+1)~YearSinceGF+Gap_location,binomial,data=EAB_Paired_Macs_StandRA)
summary(glm.GaRAlog)
#not significant, residenual deviance < df
plot(glm.GaRAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.GaRAlog,3)
#flat
plot(glm.GaRAlog,2)
#worse skew, go with non transformed

#Try again with absolute abundances
#use linear models
#Elmidae
lm.Elm<-lm(InsectaColeopteraElmidae~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Elm,1)
#no fitted pattern and red line approximately at 0
plot(lm.Elm,3)
#slight positive
plot(lm.Elm,2)
#Skew at end, try with log tranformed
lm.Elmlog<-lm(log10(InsectaColeopteraElmidae+1)~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
#not significant, residenual deviance < df
plot(lm.Elmlog,1)
#no fitted pattern and red line approximately at 0
plot(lm.Elmlog,3)
#slight positive
plot(lm.Elmlog,2)
#much better
summary(lm.Elmlog)
#Intercept and year since gap formation significant

#Chironomidae
lm.Ch<-lm(InsectaDipteraChironomidae~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Ch,1)
#no fitted pattern and red line approximately at 0
plot(lm.Ch,3)
#slight positive
plot(lm.Ch,2)
#Skew at end, try with log tranformed
lm.Chlog<-lm(log10(InsectaDipteraChironomidae+1)~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Chlog,1)
#no fitted pattern and red line approximately at 0
plot(lm.Chlog,3)
#slight positive
plot(lm.Chlog,2)
#good
summary(lm.Chlog)
#intercept and year since gf significant

#Hydropsychidae
lm.Hy<-lm(InsectaTrichopteraHydropsychidae~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Hy,1)
#no fitted pattern and red line approximately at 0
plot(lm.Hy,3)
#slight positive
plot(lm.Hy,2)
#Skew at end, try with log tranformed
lm.Hylog<-lm(log10(InsectaTrichopteraHydropsychidae+1)~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Hylog,1)
#no fitted pattern and red line approximately at 0
plot(lm.Hylog,3)
#slight positive
plot(lm.Hylog,2)
#good
summary(lm.Hylog)
#year since gap formation and intercept significant

#Gammaridae
lm.Ga<-lm(MalacostracaAmphipodaGammaridae~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Ga,1)
#no fitted pattern and red line approximately at 0
plot(lm.Ga,3)
#slight positive
plot(lm.Ga,2)
#Skew at end, try with log tranformed
lm.Galog<-lm(log10(MalacostracaAmphipodaGammaridae+1)~YearSinceGF+Gap_location,data=EAB_Paired_Macs_Stand)
plot(lm.Galog,1)
#no fitted pattern and red line approximately at 0
plot(lm.Galog,3)
#slight positive
plot(lm.Galog,2)
#better
summary(lm.Galog)
#not significant

###################################
#Functional Feeding Groups
##########################

#See names
names(EAB_M_Matrix)
EAB_M_Matrix$Gap<-NULL
Macrosmelt<- melt(EAB_M_Matrix, id=c("Stream","Date","Gap_location")) 
names(Macrosmelt)[names(Macrosmelt) == "variable"] <- "Taxa"
#Combine FFG with dataset
MacroFFGs<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_Macro_FFGs.csv", sep = ",", header = T)
MacrosMFFG<-merge(Macrosmelt, MacroFFGs, by="Taxa")
str(MacrosMFFG)
MacrosFFG<-dcast(MacrosMFFG, Stream+Date+Gap_location~FFG, sum)

#Find most abundant FFG
colSums(MacrosFFG[,4:ncol(MacrosFFG)])
#most abundant FFG collectorgatherer with 1272 individuals
stat.desc(MacrosFFG$CollectorGatherer)
#mean 23.556 se 5.095
#make relative abundances
str(MacrosFFG)
EAB_FFG_MatrixRA<-MacrosFFG
EAB_FFG_MatrixRA[,4:ncol(EAB_FFG_MatrixRA)]<-data.frame(make_relative(as.matrix(EAB_FFG_MatrixRA[,4:ncol(EAB_FFG_MatrixRA)])))
EAB_FFG_MatrixRA[is.na(EAB_FFG_MatrixRA)] <- 0
sort(colMeans(EAB_FFG_MatrixRA[,4:ncol(EAB_FFG_MatrixRA)]))
#collector gatherer
stat.desc(EAB_FFG_MatrixRA$CollectorGatherer)

#combine eab_Paired_ALL_Stand with macs
EAB_Paired_FFG_Stand<-merge(EAB_Ter_H2o_aqcwd_ALL,MacrosFFG, by=c("Stream","Gap_location","Date"))
#split up into environmental and community
names(EAB_Paired_FFG_Stand)
EAB_Paired_FFG_Stand_com<-EAB_Paired_FFG_Stand[,171:ncol(EAB_Paired_FFG_Stand)]
EABFFGComtotbray0Stand<-as.matrix(bray0(EAB_Paired_FFG_Stand_com))
EAB_Paired_FFG_Stand_env<-EAB_Paired_FFG_Stand[,1:170]
EAB_Paired_FFG_Stand_GA<-EAB_Paired_FFG_Stand_env$Gap_Area
EAB_Paired_FFG_Stand_GA_cat<-as.factor(EAB_Paired_FFG_Stand_GA)
EAB_Paired_FFG_Stand_YGF<-EAB_Paired_FFG_Stand_env$Year_Gapformation
EAB_Paired_FFG_Stand_YGF_cat<-as.factor(EAB_Paired_FFG_Stand_YGF)
EAB_Paired_FFG_Stand_GL<-EAB_Paired_FFG_Stand_env$Gap_location
EAB_Paired_FFG_Stand_W<-EAB_Paired_FFG_Stand_env$Watershed

#analyze community
adonis2(as.dist(EABFFGComtotbray0Stand) ~ Year_Gapformation*Gap_location,
       data=EAB_Paired_FFG_Stand_env,permutations=9999)
#Year of gap formation significant

EAB_FFG_Com_YGF_indic<-data.frame(signassoc(EAB_Paired_FFG_Stand_com, cluster=EAB_Paired_FFG_Stand_YGF,  mode=0, alternative = "two.sided",control = how(nperm=9999)))
EAB_FFG_Com_YGF_indic_sig<-subset(EAB_FFG_Com_YGF_indic, psidak<=0.05)
#All 5 families indicate YGF
#Collector Filterer Grazer collector gatherer and Predator Kalamazoo
#shredder Grand

#Gap Location
EAB_FFG_Com_GL_indic<-data.frame(signassoc(EAB_Paired_FFG_Stand_com, cluster=EAB_Paired_FFG_Stand_GL,  mode=0, alternative = "two.sided",control = how(nperm=9999)))
EAB_FFG_Com_GL_indic_sig<-subset(EAB_FFG_Com_GL_indic, psidak<=0.05)
#no indicator ffgs

#Year of gap formation with relative abundances
EAB_FFG_Com_YGF_indicRA<-data.frame(signassoc(EAB_FFG_MatrixRA[,4:ncol(EAB_FFG_MatrixRA)], cluster=EAB_Paired_FFG_Stand_YGF,  mode=0, alternative = "two.sided",control = how(nperm=9999)))
EAB_FFG_Com_YGF_indicRA_sig<-subset(EAB_FFG_Com_YGF_indicRA, psidak<=0.05)
#Collector gatherers significantly indicate 2014

#Gap Location
EAB_FFG_Com_GL_indicRA<-data.frame(signassoc(EAB_FFG_MatrixRA[,4:ncol(EAB_FFG_MatrixRA)], cluster=EAB_Paired_FFG_Stand_GL,  mode=0, alternative = "two.sided",control = how(nperm=9999)))
EAB_FFG_Com_GL_indicRA_sig<-subset(EAB_FFG_Com_GL_indicRA, psidak<=0.05)
#no indicator ffgs

#NMDS analysis
EAB_Paired_FFG_NMDS<-metaMDS(as.dist(EABFFGComtotbray0Stand))
#stress 0.0906 

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_FFG_NMDS)

#NMDS plot for watershed
ordiplot(EAB_Paired_FFG_NMDS, type="n")
with(EAB_Paired_FFG_NMDS, points(EAB_Paired_FFG_NMDS, display="sites", col=Watershed_col_vec[EAB_Paired_FFG_Stand_W], pch=19))
with(EAB_Paired_FFG_NMDS, legend("topleft", legend=levels(EAB_Paired_FFG_Stand_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_W, kind="se", conf=0.95, lwd=2, col="#af8dc3", show.groups = "Clinton"))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_W, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Grand"))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#NMDS plot for watershed
ordiplot(EAB_Paired_FFG_NMDS, type="n")
with(EAB_Paired_FFG_NMDS, points(EAB_Paired_FFG_NMDS, display="sites", col=YGF_col_vec[EAB_Paired_FFG_Stand_YGF_cat], pch=19))
with(EAB_Paired_FFG_NMDS, legend("topleft", legend=levels(EAB_Paired_FFG_Stand_YGF_cat), bty="n", col=YGF_col_vec, pch=19, pt.bg=YGF_col_vec))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#67000d", show.groups = "2006"))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#cb181d", show.groups = "2008"))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#fb6a4a", show.groups = "2011"))
with(EAB_Paired_FFG_NMDS, ordiellipse(EAB_Paired_FFG_NMDS, EAB_Paired_FFG_Stand_YGF_cat, kind="se", conf=0.95, lwd=2, col="#fcbba1", show.groups = "2014"))

#Model relative abundances of collector gatherers
#Combine FFGs with other metadata
EAB_Paired_FFG_Stand_RA<-merge(EAB_Ter_H2o_aqcwd_ALL,EAB_FFG_MatrixRA, by=c("Stream","Gap_location","Date"))
#collector gatherers
#use logistic model for proportional data
glm.CGRA<-glm(CollectorGatherer~YearSinceGF+Gap_location,binomial,data=EAB_Paired_FFG_Stand_RA)
summary(glm.CGRA)
#not significant, residenual deviance < df
plot(glm.CGRA,1)
#no fitted pattern and red line approximately at 0
plot(glm.CGRA,3)
#slight decrease
plot(glm.CGRA,2)
#looks okay, but try with log transformed
glm.CGRAlog<-glm(log10(CollectorGatherer+1)~YearSinceGF+Gap_location,binomial,data=EAB_Paired_FFG_Stand_RA)
summary(glm.CGRAlog)
#not significant, residenual deviance < df
plot(glm.CGRAlog,1)
#no fitted pattern and red line approximately at 0
plot(glm.CGRAlog,3)
#slight decrease
plot(glm.CGRAlog,2)
#Isn't really improved so go with original
EAB_Paired_FFG_Stand_RA %>%
  group_by(Year_Gapformation) %>%
  get_summary_stats("CollectorGatherer", type = "mean_se")

#Create a figure for the manuscript that panels several responses to gap and year of gap formation
#DO, ALL Richness, macroinvertebrate richness, Collector gatherer RA
CJFASFigdf<-data.frame(EAB_Paired_Macs_Stand$Gap_location,EAB_Paired_Macs_Stand$X.DO,
                       EAB_Paired_Macs_Stand$Year_Gapformation,EAB_Paired_Macs_Stand$ALL_Richness,
              EAB_Paired_Macs_Stand$Mac_Richness, EAB_Paired_Macs_Stand$InsectaDipteraChironomidae)
CJFASFigdf$EAB_Paired_Macs_Stand.InsectaDipteraChironomidae<-log10(CJFASFigdf$EAB_Paired_Macs_Stand.InsectaDipteraChironomidae+1)
names(CJFASFigdf)
CJFASFigdfmelt<- melt(CJFASFigdf, id = c("EAB_Paired_Macs_Stand.Gap_location","EAB_Paired_Macs_Stand.Year_Gapformation"))
names(CJFASFigdfmelt)
levels(CJFASFigdfmelt$variable)
levels(CJFASFigdfmelt$variable) <- list("A. % Dissolved Oxygen"="EAB_Paired_Macs_Stand.X.DO", 
                                        "B. Aquatic Litter Richness"="EAB_Paired_Macs_Stand.ALL_Richness",
                                        "C. Macroinvtebrate Richness"="EAB_Paired_Macs_Stand.Mac_Richness",
                                        "D. log(Chironomidae Abundance)"="EAB_Paired_Macs_Stand.InsectaDipteraChironomidae")
ggplot(CJFASFigdfmelt, aes(x=EAB_Paired_Macs_Stand.Year_Gapformation, y=value, 
                           color=EAB_Paired_Macs_Stand.Gap_location)) +
  geom_point(position=position_jitter(h=0.3, w=0.3),size=1)+
  xlab("Year of Gap Formation")+
  scale_x_continuous(breaks=c(2006,2008,2011,2014))+
  scale_color_manual(values=Gap_location_col_vec)+
  labs(color = "Gap Location")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_blank(),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=12),legend.text = element_text(size=10))+
  facet_grid(variable~.,scales = "free_y",labeller = label_wrap_gen(width=20))
