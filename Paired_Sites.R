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

#Calculate volume of each piece of wood
EAB_AqCWD$V<-pi*((((EAB_AqCWD$Diameter/100)/2)^2)*(EAB_AqCWD$In_Stream_Length/100))

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
                                 Watershed+Aquatic_Transect_width, data=EAB_AqCWD_Paired,
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

##############
#Aquatic CWD watershed
###################

#model this way: one way anova aq CWD ~ watershed

#number of logs per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logs_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logs_per_length)
#no outliers
#build model
lm.nologpm<-lm(No_logs_per_length ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.nologpm))
#all within gray bar
shapiro_test(residuals(lm.nologpm))
#p 0.799 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logs_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.nologpm, 1)
#no relationship
#all assumptions  met
nologpm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logs_per_length~Watershed)
nologpm.aov
#not significant

#number of logs per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logs_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logs_per_area)
#no outliers
#build model
lm.nologpha<-lm(No_logs_per_area ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.nologpha))
#all within gray bar
shapiro_test(residuals(lm.nologpha))
#p 0.231 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logs_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.nologpha, 1)
#no relationship
#all assumptions  met
nologpha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logs_per_area~Watershed)
nologpha.aov
#not significant

#volume CWD per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWD_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(VolCWD_per_length)
#no outliers
#build model
lm.volcwdpm<-lm(VolCWD_per_length ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volcwdpm))
#all within gray bar
shapiro_test(residuals(lm.volcwdpm))
#p 0.933 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "VolCWD_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.volcwdpm, 1)
#no relationship
#all assumptions  met
volcwdpm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(VolCWD_per_length~Watershed)
volcwdpm.aov
#not significant

#volume CWD per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWD_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(VolCWD_per_area)
#no outliers
#build model
lm.volcwdpha<-lm(VolCWD_per_area ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volcwdpha))
#all within gray bar
shapiro_test(residuals(lm.volcwdpha))
#p 0.948 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "VolCWD_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.volcwdpha, 1)
#no relationship
#all assumptions  met
volcwdpha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(VolCWD_per_area~Watershed)
volcwdpha.aov
#not significant

#number of log jams per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
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
#p 1.00 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logjams_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.NoLJpm, 1)
#no relationship
#all assumptions  met
NoLJpm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logjams_per_length~Watershed)
NoLJpm.aov
#not significant

#number of log jams per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjams_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logjams_per_area)
#no outliers
#build model
lm.NoLJpha<-lm(No_logjams_per_area ~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.NoLJpha))
#points at end just outside of gray bar
shapiro_test(residuals(lm.NoLJpha))
#p 0.747 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logjams_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.NoLJpha, 1)
#no relationship
#all assumptions  met
NoLJpha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logjams_per_area~Watershed)
NoLJpha.aov
#not significant

#number of log jam pieces per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjam_pieces_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logjam_pieces_per_length)
#no outliers
#build model
lm.NoLJppm<-lm(No_logjam_pieces_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.NoLJppm))
#points all inside gray bar
shapiro_test(residuals(lm.NoLJppm))
#p 0.968 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logjam_pieces_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.NoLJppm, 1)
#no relationship
#all assumptions  met
NoLJppm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logjam_pieces_per_length~Watershed)
NoLJppm.aov
#not significant

#number of log jam pieces per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjam_pieces_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_logjam_pieces_per_area)
#no outliers
#build model
lm.NoLJppha<-lm(No_logjam_pieces_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.NoLJppha))
#points all inside gray bar
shapiro_test(residuals(lm.NoLJppha))
#p 0.717 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_logjam_pieces_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.NoLJppha, 1)
#no relationship
#all assumptions  met
NoLJppha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_logjam_pieces_per_area~Watershed)
NoLJppha.aov
#not significant

#prop log jam pieces
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_pieces, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_logjam_pieces)
#no outliers
#build model
lm.propLJp<-lm(prop_logjam_pieces~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propLJp))
#points all inside gray bar
shapiro_test(residuals(lm.propLJp))
#p 0.634 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_logjam_pieces", facet.by = "Watershed")
#all within gray bar
plot(lm.propLJp, 1)
#no relationship
#all assumptions  met
propLJp.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_logjam_pieces~Watershed)
propLJp.aov
#not significant

#log jam volume per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWDLJ_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(VolCWDLJ_per_length)
#no outliers
#build model
lm.volLJpm<-lm(VolCWDLJ_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volLJpm))
#points all inside gray bar
shapiro_test(residuals(lm.volLJpm))
#p 0.103 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "VolCWDLJ_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.volLJpm, 1)
#no relationship
#all assumptions  met
volLJpm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(VolCWDLJ_per_length~Watershed)
volLJpm.aov
#not significant

#log jam volume per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWDLJ_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(VolCWDLJ_per_area)
#no outliers
#build model
lm.volLJpha<-lm(VolCWDLJ_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volLJpha))
#points all inside gray bar
shapiro_test(residuals(lm.volLJpha))
#p 0.445 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "VolCWDLJ_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.volLJpha, 1)
#no relationship
#all assumptions  met
volLJpha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(VolCWDLJ_per_area~Watershed)
volLJpha.aov
#not significant

#prop log jam volume
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_vol, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_logjam_vol)
#no outliers
#build model
lm.propLJvol<-lm(prop_logjam_vol~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propLJvol))
#points all inside gray bar
shapiro_test(residuals(lm.propLJvol))
#p 0.176 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_logjam_vol", facet.by = "Watershed")
#all within gray bar
plot(lm.propLJvol, 1)
#no relationship
#all assumptions  met
propLJvol.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_logjam_vol~Watershed)
propLJvol.aov
#not significant

#number of insitu per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_insitu_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_insitu_per_length)
#no outliers
#build model
lm.noinsitupm<-lm(No_insitu_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.noinsitupm))
#points all inside gray bar
shapiro_test(residuals(lm.noinsitupm))
#p 0.53 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_insitu_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.noinsitupm, 1)
#no relationship
#all assumptions  met
noinsitupm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_insitu_per_length~Watershed)
noinsitupm.aov
#not significant

#number of insitu per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_insitu_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_insitu_per_area)
#no outliers
#build model
lm.noinsitupha<-lm(No_insitu_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.noinsitupha))
#points all inside gray bar
shapiro_test(residuals(lm.noinsitupha))
#p 0.244 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_insitu_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.noinsitupha, 1)
#no relationship
#all assumptions  met
noinsitupha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_insitu_per_area~Watershed)
noinsitupha.aov
#not significant

#prop insitu pieces
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_pieces, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_insitu_pieces)
#no outliers
#build model
lm.propinsitup<-lm(prop_insitu_pieces~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propinsitup))
#points all inside gray bar
shapiro_test(residuals(lm.propinsitup))
#p 0.452 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_insitu_pieces", facet.by = "Watershed")
#all within gray bar
plot(lm.propinsitup, 1)
#no relationship
#all assumptions  met
propinsitup.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_insitu_pieces~Watershed)
propinsitup.aov
#not significant

#volume insitu per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_insitu_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(Vol_insitu_per_length)
#no outliers
#build model
lm.volinsitupm<-lm(Vol_insitu_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volinsitupm))
#points all inside gray bar
shapiro_test(residuals(lm.volinsitupm))
#p 0.917 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "Vol_insitu_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.volinsitupm, 1)
#no relationship
#all assumptions  met
volinsitupm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(Vol_insitu_per_length~Watershed)
volinsitupm.aov
#not significant

#volume insitu per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_insitu_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(Vol_insitu_per_area)
#no outliers
#build model
lm.volinsitupha<-lm(Vol_insitu_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volinsitupha))
#points all inside gray bar
shapiro_test(residuals(lm.volinsitupha))
#p 0.714 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "Vol_insitu_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.volinsitupha, 1)
#no relationship
#all assumptions  met
volinsitupha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(Vol_insitu_per_area~Watershed)
volinsitupha.aov
#not significant

#prop volume insitu
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_vol, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_insitu_vol)
#no outliers
#build model
lm.propinsituv<-lm(prop_insitu_vol~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propinsituv))
#points all inside gray bar
shapiro_test(residuals(lm.propinsituv))
#p 0.393 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_insitu_vol", facet.by = "Watershed")
#all within gray bar
plot(lm.propinsituv, 1)
#no relationship
#all assumptions  met
propinsituv.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_insitu_vol~Watershed)
propinsituv.aov
#not significant

#number of ash per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_ash_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_ash_per_length)
#no outliers
#build model
lm.noapm<-lm(No_ash_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.noapm))
#points all inside gray bar
shapiro_test(residuals(lm.noapm))
#p 0.390 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_ash_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.noapm, 1)
#no relationship
#all assumptions  met
noapm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_ash_per_length~Watershed)
noapm.aov
#not significant

#number of ash per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_ash_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(No_ash_per_area)
#no outliers
#build model
lm.noapha<-lm(No_ash_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.noapha))
#points all inside gray bar
shapiro_test(residuals(lm.noapha))
#p 0.647 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "No_ash_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.noapha, 1)
#no relationship
#all assumptions  met
noapha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(No_ash_per_area~Watershed)
noapha.aov
#not significant

#prop ash pieces
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_pieces, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_ash_pieces)
#no outliers
#build model
lm.propap<-lm(prop_ash_pieces~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propap))
#points all inside gray bar
shapiro_test(residuals(lm.propap))
#p 0.411 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_ash_pieces", facet.by = "Watershed")
#all within gray bar
plot(lm.propap, 1)
#no relationship
#all assumptions  met
propap.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_ash_pieces~Watershed)
propap.aov
#watershed significant, tukey's posthoc test
propap.ph<- EAB_AqCWD_Paired_ag%>% tukey_hsd(prop_ash_pieces~Watershed)
propap.ph
#clinton signficantly greater than grand
#visualize
propap.ph<-propap.ph%>% add_xy_position(x = "Watershed")
ggplot(EAB_AqCWD_Paired_ag, aes(x=Watershed, y=prop_ash_pieces)) +
  geom_point() +
  ylab("Ash CWD Pieces Relative Abundance")+
  xlab("Watershed") +
  stat_pvalue_manual(propap.ph, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(propap.aov, detailed = TRUE),
    caption = get_pwc_label(propap.ph)
  )+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        plot.subtitle = element_text(size = 14),plot.caption=element_text(size=14),
        legend.position = "none")

#volume ash per length
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_ash_per_length, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(Vol_ash_per_length)
#no outliers
#build model
lm.volapm<-lm(Vol_ash_per_length~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volapm))
#points all inside gray bar
shapiro_test(residuals(lm.volapm))
#p 0.526 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "Vol_ash_per_length", facet.by = "Watershed")
#all within gray bar
plot(lm.volapm, 1)
#no relationship
#all assumptions  met
volapm.aov <- EAB_AqCWD_Paired_ag%>% anova_test(Vol_ash_per_length~Watershed)
volapm.aov
#nothing significant

#volume ash per area
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_ash_per_area, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(Vol_ash_per_area)
#no outliers
#build model
lm.volapha<-lm(Vol_ash_per_area~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.volapha))
#points all inside gray bar
shapiro_test(residuals(lm.volapha))
#p 0.920 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "Vol_ash_per_area", facet.by = "Watershed")
#all within gray bar
plot(lm.volapha, 1)
#no relationship
#all assumptions  met
volapha.aov <- EAB_AqCWD_Paired_ag%>% anova_test(Vol_ash_per_area~Watershed)
volapha.aov
#nothing significant

#prop volume ash
#summarize data
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_vol, type = "mean_se")
#Check assumptions
#outliers
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  identify_outliers(prop_ash_vol)
#no outliers
#build model
lm.propav<-lm(prop_ash_vol~ Watershed,data=EAB_AqCWD_Paired_ag)
#qqplot
ggqqplot(residuals(lm.propav))
#points all inside gray bar
shapiro_test(residuals(lm.propav))
#p 0.131 normally distributed residuals
#sample size too low to see normality for each watershed
ggqqplot(EAB_AqCWD_Paired_ag, "prop_ash_vol", facet.by = "Watershed")
#all within gray bar
plot(lm.propav, 1)
#no relationship
#all assumptions  met
propav.aov <- EAB_AqCWD_Paired_ag%>% anova_test(prop_ash_vol~Watershed)
propav.aov
#nothing significant

##########################
#Summary Aquatic CWD by watershed
#nothing significant except
#the proportion of pieces of CWD that were ash was greater in clinton compared to grand
##############################

#Get summary stats for dissertation
#nubmer of logs per area
EAB_AqCWD_Paired_ag$No_logs_per_ha<-1000*EAB_AqCWD_Paired_ag$No_logs_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logs_per_ha, type = "mean_se")
#Volume CWD per area
EAB_AqCWD_Paired_ag$VolCWD_per_ha<-1000*EAB_AqCWD_Paired_ag$VolCWD_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWD_per_ha, type = "mean_se")
#no logjams per area
EAB_AqCWD_Paired_ag$No_logjams_per_ha<-1000*EAB_AqCWD_Paired_ag$No_logjams_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjams_per_ha, type = "mean_se")
#no logjam pieces per area
EAB_AqCWD_Paired_ag$No_logjam_pieces_per_ha<-1000*EAB_AqCWD_Paired_ag$No_logjam_pieces_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_logjam_pieces_per_ha, type = "mean_se")
#logjam pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_pieces, type = "mean_se")
#logjam volume per stream area
EAB_AqCWD_Paired_ag$VolCWDLJ_per_ha<-1000*EAB_AqCWD_Paired_ag$VolCWDLJ_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(VolCWDLJ_per_ha, type = "mean_se")
#logjam volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_logjam_vol, type = "mean_se")
#no insitu per area
EAB_AqCWD_Paired_ag$No_insitu_per_ha<-1000*EAB_AqCWD_Paired_ag$No_insitu_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_insitu_per_ha, type = "mean_se")
#insitu pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_pieces, type = "mean_se")
#insitu volume per stream area
EAB_AqCWD_Paired_ag$Vol_insitu_per_ha<-1000*EAB_AqCWD_Paired_ag$Vol_insitu_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_insitu_per_ha, type = "mean_se")
#insitu volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_insitu_vol, type = "mean_se")
#no ash per area
EAB_AqCWD_Paired_ag$No_ash_per_ha<-1000*EAB_AqCWD_Paired_ag$No_ash_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(No_ash_per_ha, type = "mean_se")
#ash pieces relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_pieces, type = "mean_se")
#ash volume per stream area
EAB_AqCWD_Paired_ag$Vol_ash_per_ha<-1000*EAB_AqCWD_Paired_ag$Vol_ash_per_area
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(Vol_ash_per_ha, type = "mean_se")
#ash volume relative abundance
EAB_AqCWD_Paired_ag %>%
  group_by(Watershed) %>%
  get_summary_stats(prop_ash_vol, type = "mean_se")
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

#########################
#CWD correlation analysis with terrestrial datasets
#######################

#Merge with terrestrial dataset
EAB_CWD_Paired<-merge(EAB_AqCWD_Paired_ag, EAB_Terr_Paired_Gap, by="Stream")

#look for correltaions between aquatic CWD, terrestrial cwd and terrestrial snag variables
#use pearson correlations if assumptions met

#Test normality for terrestrial variables used:
#Number of terrestrial cwd logs per ha
shapiro.test(EAB_CWD_Paired$Total_CWD_Nopha)
#W = 0.78069, p-value = 0.03912, not normal
skewness(EAB_CWD_Paired$Total_CWD_Nopha)
#right skewed, log transform
EAB_CWD_Paired$log10Total_CWD_Nopha<-log10(EAB_CWD_Paired$Total_CWD_Nopha+1)
shapiro.test(EAB_CWD_Paired$log10Total_CWD_Nopha)
#W = 0.79085, p-value = 0.04855, not normal but better, square root transform
EAB_CWD_Paired$sqrtTotal_CWD_Nopha<-sqrt(EAB_CWD_Paired$Total_CWD_Nopha+1)
shapiro.test(EAB_CWD_Paired$sqrtTotal_CWD_Nopha)
#W = 0.78567, p-value = 0.0435, not normal
#use nonparametric test like spearman correlation

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
#all within gray bar - normally distributed

#volume snags per ha
shapiro.test(EAB_CWD_Paired$Total_Snag_BA_m2pha)
#W = 0.85384, p-value = 0.169, normal
ggqqplot(EAB_CWD_Paired$Total_Snag_BA_m2pha)
#all within gray bar - normally distributed

#no ash cwd per ha
shapiro.test(EAB_CWD_Paired$NoCWDperha_Ash)
#W = 0.67658, p-value = 0.003469, not normal
skewness(EAB_CWD_Paired$NoCWDperha_Ash)
#0.5371711, log transform
EAB_CWD_Paired$log10NoCWDperha_Ash<-log10(EAB_CWD_Paired$NoCWDperha_Ash+1)
shapiro.test(EAB_CWD_Paired$log10NoCWDperha_Ash)
#W = 0.74478, p-value = 0.0177, still not normal, so use non-parametric test

#no ash snag per ha
shapiro.test(EAB_CWD_Paired$Noperha_Ash_Snag)
#W = 0.91384, p-value = 0.4621 normally distributed

#prop ash cwd pieces
shapiro.test(EAB_CWD_Paired$prop_Ash_CWD_no)
#W = 0.9758, p-value = 0.9289, normal

#prop ash snag no
shapiro.test(EAB_CWD_Paired$prop_Ash_snag_no)
#W = 0.81323, p-value = 0.07704, normal

#vol ash cwd per ha
shapiro.test(EAB_CWD_Paired$VolCWD_m3perha_Ash)
#W = 0.71615, p-value = 0.009105, not normal
skewness(EAB_CWD_Paired$VolCWD_m3perha_Ash)
#positive, so log transform
EAB_CWD_Paired$log10VolCWD_m3perha_Ash<-log10(EAB_CWD_Paired$VolCWD_m3perha_Ash+1)
shapiro.test(EAB_CWD_Paired$log10VolCWD_m3perha_Ash)
#W = 0.88808, p-value = 0.3083, normal so use this transformation

#BA ash snag per ha
shapiro.test(EAB_CWD_Paired$BAperha_Ash_Snag)
#W = 0.78571, p-value = 0.04355, not normal
skewness(EAB_CWD_Paired$BAperha_Ash_Snag)
#right, so log transform
EAB_CWD_Paired$log10BAperha_Ash_Snag<-log10(EAB_CWD_Paired$BAperha_Ash_Snag+1)
shapiro.test(EAB_CWD_Paired$log10BAperha_Ash_Snag)
#W = 0.84313, p-value = 0.1384, normally distributed

#prop ash cwd vol
shapiro.test(EAB_CWD_Paired$prop_Ash_CWD_vol)
#W = 0.83729, p-value = 0.1238, normal

#prop ash snag vol
shapiro.test(EAB_CWD_Paired$prop_Ash_snag_BA)
#W = 0.96623, p-value = 0.8662, normal

###################
#correlation tests
#Number of logs per stream length
#test normality
shapiro.test(EAB_CWD_Paired$No_logs_per_length)
#W = 0.89675, p-value = 0.3551, normally distributed
ggqqplot(EAB_CWD_Paired$No_logs_per_length)
#one point slightly higher than gray bar
#CWD - not normal
cor.nolpm.nocwdpha<- cor.test(EAB_CWD_Paired$No_logs_per_length, EAB_CWD_Paired$Total_CWD_Nopha, 
                method = "spearman")
cor.nolpm.nocwdpha
#not significant
#snag - normal so use pearson
cor.nolpm.nospha<-cor.test(EAB_CWD_Paired$No_logs_per_length, EAB_CWD_Paired$Total_Snag_Nopha,
                           method="pearson")
cor.nolpm.nospha
#not significant

#Number of logs per stream area
#test normality
shapiro.test(EAB_CWD_Paired$No_logs_per_area)
#W = 0.8761, p-value = 0.2516, normally distributed
ggqqplot(EAB_CWD_Paired$No_logs_per_area)
#all points within gray bar
#CWD - not normal
cor.nolpha.nocwdpha<- cor.test(EAB_CWD_Paired$No_logs_per_area, EAB_CWD_Paired$Total_CWD_Nopha, 
                             method = "spearman")
cor.nolpha.nocwdpha
#not significant
#snag - normal so use pearson
cor.nolpha.nospha<-cor.test(EAB_CWD_Paired$No_logs_per_area, EAB_CWD_Paired$Total_Snag_Nopha,
                           method="pearson")
cor.nolpha.nospha
#significant positive
#Visualize
#convert to pieces per ha
EAB_CWD_Paired$No_logs_per_ha<-EAB_CWD_Paired$No_logs_per_area/0.0001
ggplot(EAB_CWD_Paired, aes(y=No_logs_per_ha, x=Total_Snag_Nopha)) +
  geom_point(size=3) +
  xlab("Snag Frequency (Total Number per ha)")+
  ylab("Aquatic CWD Frequency (Total Number per ha)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")

#Volume CWD per stream length
#test normality
shapiro.test(EAB_CWD_Paired$VolCWD_per_length)
#W = 0.9322, p-value = 0.5972, normally distributed
ggqqplot(EAB_CWD_Paired$VolCWD_per_length)
#all points within gray bar
#CWD - use log transformation
cor.volCWDpm.volcwdpha<- cor.test(EAB_CWD_Paired$VolCWD_per_length,
                               EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,method = "pearson")
cor.volCWDpm.volcwdpha
#not significant
#snag - normal so use pearson
cor.volCWDpm.volspha<-cor.test(EAB_CWD_Paired$VolCWD_per_length,EAB_CWD_Paired$Total_Snag_BA_m2pha,
                            method="pearson")
cor.volCWDpm.volspha
#not significant

#Number of log jams per stream length
#test normality
shapiro.test(EAB_CWD_Paired$No_logjams_per_length)
#W = 0.86348, p-value = 0.2014, normally distributed
ggqqplot(EAB_CWD_Paired$No_logjams_per_length)
#all points within gray bar
#CWD - not normal
cor.noLJpm.nocwdpha<- cor.test(EAB_CWD_Paired$No_logjams_per_length,
                                  EAB_CWD_Paired$Total_CWD_Nopha,method = "spearman")
cor.noLJpm.nocwdpha
#significant negative
#Visualize
ggplot(EAB_CWD_Paired, aes(x=No_logjams_per_length, y=Total_CWD_Nopha)) +
  geom_point(size=3) +
  geom_smooth(method=lm, color="black", se=FALSE) +
  ylab("Riparian CWD Density (Total Number per ha)")+
  xlab("Logjam Frequency (Total Number per m)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")
#snag - normal so use pearson
cor.noLJpm.nospha<-cor.test(EAB_CWD_Paired$No_logjams_per_length,EAB_CWD_Paired$Total_Snag_Nopha,
                               method="pearson")
cor.noLJpm.nospha
#not significant

#Number of log jams per stream area
#test normality
shapiro.test(EAB_CWD_Paired$No_logjams_per_area)
#W = 0.81451, p-value = 0.07906, normally distributed
ggqqplot(EAB_CWD_Paired$No_logjams_per_area)
#one point slightly above gray area
#CWD - not normal
cor.noLJpha.nocwdpha<- cor.test(EAB_CWD_Paired$No_logjams_per_area,
                               EAB_CWD_Paired$Total_CWD_Nopha,method = "spearman")
cor.noLJpha.nocwdpha
#significant negative
#Visualize
#covert logjams per m to logjams per ha
EAB_CWD_Paired$No_logjams_per_ha<-1000*EAB_CWD_Paired$No_logjams_per_area
ggplot(EAB_CWD_Paired, aes(x=No_logjams_per_ha, y=Total_CWD_Nopha)) +
  geom_point(size=3) +
  geom_smooth(method=lm, color="black", se=FALSE) +
  ylab("Riparian CWD Density (Total Number per ha)")+
  xlab("Logjam Density (Total Number per ha)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")
#snag - normal so use pearson
cor.noLJpha.nospha<-cor.test(EAB_CWD_Paired$No_logjams_per_area,EAB_CWD_Paired$Total_Snag_Nopha,
                            method="pearson")
cor.noLJpha.nospha
#not significant

#Number of log jam pieces per stream length
#test normality
shapiro.test(EAB_CWD_Paired$No_logjam_pieces_per_length)
#W = 0.94719, p-value = 0.7175, normally distributed
ggqqplot(EAB_CWD_Paired$No_logjam_pieces_per_length)
#all points within gray area
#CWD - not normal
cor.noLJppm.nocwdpha<- cor.test(EAB_CWD_Paired$No_logjam_pieces_per_length,
                                EAB_CWD_Paired$Total_CWD_Nopha,method = "spearman")
cor.noLJppm.nocwdpha
#not significant
#snag - normal so use pearson
cor.noLJppm.nospha<-cor.test(EAB_CWD_Paired$No_logjam_pieces_per_length,EAB_CWD_Paired$Total_Snag_Nopha,
                             method="pearson")
cor.noLJppm.nospha
#significant negative
#Visualize
ggplot(EAB_CWD_Paired, aes(x=No_logjam_pieces_per_length, y=Total_Snag_Nopha)) +
  geom_point(size=3) +
  geom_smooth(method=lm, color="black", se=FALSE) +
  ylab("Snag Density (Total Number per ha)")+
  xlab("Logjam CWD Pieces Frequency (Total Number per m)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")

#Number of log jam pieces per stream area
#test normality
shapiro.test(EAB_CWD_Paired$No_logjam_pieces_per_area)
#W = 0.96328, p-value = 0.8446, normally distributed
ggqqplot(EAB_CWD_Paired$No_logjam_pieces_per_area)
#all points within gray area
#CWD - not normal
cor.noLJppha.nocwdpha<- cor.test(EAB_CWD_Paired$No_logjam_pieces_per_area,
                                EAB_CWD_Paired$Total_CWD_Nopha,method = "spearman")
cor.noLJppha.nocwdpha
#not significant
#snag - normal so use pearson
cor.noLJppha.nospha<-cor.test(EAB_CWD_Paired$No_logjam_pieces_per_area,EAB_CWD_Paired$Total_Snag_Nopha,
                             method="pearson")
cor.noLJppha.nospha
#not significant

#proportion log jam pieces
#test normality
shapiro.test(EAB_CWD_Paired$prop_logjam_pieces)
#W = 0.92063, p-value = 0.5099, normally distributed
ggqqplot(EAB_CWD_Paired$prop_logjam_pieces)
#all points within gray area
#CWD - not normal
cor.propLJp.nocwdpha<- cor.test(EAB_CWD_Paired$prop_logjam_pieces,
                                 EAB_CWD_Paired$Total_CWD_Nopha,method = "spearman")
cor.propLJp.nocwdpha
#not significant
#snag - normal so use pearson
cor.propLJp.nospha<-cor.test(EAB_CWD_Paired$prop_logjam_pieces,EAB_CWD_Paired$Total_Snag_Nopha,
                              method="pearson")
cor.propLJp.nospha
#significant negative
#Visualize
ggplot(EAB_CWD_Paired, aes(x=prop_logjam_pieces, y=Total_Snag_Nopha)) +
  geom_point(size=3) +
  geom_smooth(method=lm, color="black", se=FALSE) +
  ylab("Snag Density (Total Number per ha)")+
  xlab("Logjam CWD Pieces Relative Abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")

#log jam volume per stream length
#test normality
shapiro.test(EAB_CWD_Paired$VolCWDLJ_per_length)
#W = 0.86728, p-value = 0.2156, normally distributed
#CWD - log transformation
cor.volLJpm.volcwdpha<- cor.test(EAB_CWD_Paired$VolCWDLJ_per_length,
                                EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,method = "pearson")
cor.volLJpm.volcwdpha
#not significant
#snag - normal so use pearson
cor.volLJpm.baspha<-cor.test(EAB_CWD_Paired$VolCWDLJ_per_length,EAB_CWD_Paired$Total_Snag_BA_m2pha,
                             method="pearson")
cor.volLJpm.baspha
#not significant

#log jam volume per stream area
#test normality
shapiro.test(EAB_CWD_Paired$VolCWDLJ_per_area)
#W = 0.87469, p-value = 0.2455, normally distributed
#CWD - log transformation
cor.volLJpha.volcwdpha<- cor.test(EAB_CWD_Paired$VolCWDLJ_per_area,
                                 EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,method = "pearson")
cor.volLJpha.volcwdpha
#not significant
#snag - normal so use pearson
cor.volLJpha.baspha<-cor.test(EAB_CWD_Paired$VolCWDLJ_per_area,EAB_CWD_Paired$Total_Snag_BA_m2pha,
                             method="pearson")
cor.volLJpha.baspha
#not significant

#prop log jam volue
#test normality
shapiro.test(EAB_CWD_Paired$prop_logjam_vol)
#W = 0.92333, p-value = 0.5297, normally distributed
#CWD - log transformation
cor.propLJv.volcwdpha<- cor.test(EAB_CWD_Paired$prop_logjam_vo,
                                  EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,method = "pearson")
cor.propLJv.volcwdpha
#not significant
#snag - normal so use pearson
cor.propLJv.baspha<-cor.test(EAB_CWD_Paired$prop_logjam_vo,EAB_CWD_Paired$Total_Snag_BA_m2pha,
                              method="pearson")
cor.propLJv.baspha
#not significant

#no insitu per stream length
#test normality
shapiro.test(EAB_CWD_Paired$No_insitu_per_length)
#W = 0.92252, p-value = 0.5237, normally distributed
#CWD - not normal
cor.noinsitupm.nocwdpha<- cor.test(EAB_CWD_Paired$No_insitu_per_length,
                                 EAB_CWD_Paired$Total_CWD_Nopha,
                                 method = "spearman")
cor.noinsitupm.nocwdpha
#not significant
#snag - normal so use pearson
cor.noinsitupm.nospha<-cor.test(EAB_CWD_Paired$No_insitu_per_length,
                                EAB_CWD_Paired$Total_Snag_Nopha,
                                method="pearson")
cor.noinsitupm.nospha
#not significant

#no insitu per stream area
#test normality
shapiro.test(EAB_CWD_Paired$No_insitu_per_area)
#W = 0.8563, p-value = 0.1768, normally distributed
#CWD - not normal
cor.noinsitupha.nocwdpha<- cor.test(EAB_CWD_Paired$No_insitu_per_area,
                                   EAB_CWD_Paired$Total_CWD_Nopha,
                                   method = "spearman")
cor.noinsitupha.nocwdpha
#not significant
#snag - normal so use pearson
cor.noinsitupha.nospha<-cor.test(EAB_CWD_Paired$No_insitu_per_area,
                                EAB_CWD_Paired$Total_Snag_Nopha,
                                method="pearson")
cor.noinsitupha.nospha
#not significant

#prop insitu pieces
#test normality
shapiro.test(EAB_CWD_Paired$prop_insitu_pieces)
#W = 0.90399, p-value = 0.3981, normally distributed
#CWD - not normal
cor.propinsitup.nocwdpha<- cor.test(EAB_CWD_Paired$prop_insitu_pieces,
                                    EAB_CWD_Paired$Total_CWD_Nopha,
                                    method = "spearman")
cor.propinsitup.nocwdpha
#not significant
#snag - normal so use pearson
cor.propinsitup.nospha<-cor.test(EAB_CWD_Paired$prop_insitu_pieces,
                                 EAB_CWD_Paired$Total_Snag_Nopha,
                                 method="pearson")
cor.propinsitup.nospha
#not significant

#volume insitu per stream length
#test normality
shapiro.test(EAB_CWD_Paired$Vol_insitu_per_length)
#W = 0.897, p-value = 0.3565, normally distributed
#CWD - log transform
cor.volinsitupm.volcwdpha<- cor.test(EAB_CWD_Paired$Vol_insitu_per_length,
                                    EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,
                                    method = "pearson")
cor.volinsitupm.volcwdpha
#significant positive
#Visualize
ripCWDDenslab<-expression(paste("Riparian CWD Density (",m^3," per ha)"))
insituvolpllab<-expression(paste("Aquatic ", italic("In situ"), " Density (",m^3," per m)"))
ggplot(EAB_CWD_Paired, aes(y=Vol_insitu_per_length, x=Total_CWD_Vol_m3pha)) +
  geom_point(size=3) +
  xlab(ripCWDDenslab)+
  ylab(insituvolpllab)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")
#snag - normal so use pearson
cor.volinsitupm.baspha<-cor.test(EAB_CWD_Paired$Vol_insitu_per_length,
                                 EAB_CWD_Paired$Total_Snag_BA_m2pha,
                                 method="pearson")
cor.volinsitupm.baspha
#not significant

#volume insitu per stream area
#test normality
shapiro.test(EAB_CWD_Paired$Vol_insitu_per_area)
#W = 0.88006, p-value = 0.2694, normally distributed
#CWD - log transform
cor.volinsitupha.volcwdpha<- cor.test(EAB_CWD_Paired$Vol_insitu_per_area,
                                     EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,
                                     method = "pearson")
cor.volinsitupha.volcwdpha
#not significant
#snag - normal so use pearson
cor.volinsitupha.baspha<-cor.test(EAB_CWD_Paired$Vol_insitu_per_area,
                                 EAB_CWD_Paired$Total_Snag_BA_m2pha,
                                 method="pearson")
cor.volinsitupha.baspha
#not significant

#prop insitu volume
#test normality
shapiro.test(EAB_CWD_Paired$prop_insitu_vol)
#W = 0.87537, p-value = 0.2484, normally distributed
#CWD - log transform
cor.propinsituvol.volcwdpha<- cor.test(EAB_CWD_Paired$prop_insitu_vol,
                                      EAB_CWD_Paired$log10Total_CWD_Vol_m3pha,
                                      method = "pearson")
cor.propinsituvol.volcwdpha
#significant positive
#visualize
propinsituvollab<-expression(paste("Aquatic ", italic("In situ"), " Relative Volume"))
ggplot(EAB_CWD_Paired, aes(y=prop_insitu_vol, x=Total_CWD_Vol_m3pha)) +
  geom_point(size=3) +
  xlab(ripCWDDenslab)+
  ylab(propinsituvollab)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")
#snag - normal so use pearson
cor.propinsituvol.baspha<-cor.test(EAB_CWD_Paired$prop_insitu_vol,
                                  EAB_CWD_Paired$Total_Snag_BA_m2pha,
                                  method="pearson")
cor.propinsituvol.baspha
#not significant

#Number of ash per stream length
#test normality
shapiro.test(EAB_CWD_Paired$No_ash_per_length)
#W = 0.97721, p-value = 0.9369, normally distributed
#CWD - not normal
cor.noapm.noacwdpha<- cor.test(EAB_CWD_Paired$No_ash_per_length,
                                       EAB_CWD_Paired$NoCWDperha_Ash,
                                       method = "spearman")
cor.noapm.noacwdpha
#not significant
#snag - normal so use pearson
cor.noapm.noaspha<-cor.test(EAB_CWD_Paired$No_ash_per_length,
                                   EAB_CWD_Paired$Noperha_Ash_Snag,
                                   method="pearson")
cor.noapm.noaspha
#not significant

#Number of ash per stream area
#test normality
shapiro.test(EAB_CWD_Paired$No_ash_per_area)
#W = 0.91198, p-value = 0.4495, normally distributed
#CWD - not normal
cor.noapha.noacwdpha<- cor.test(EAB_CWD_Paired$No_ash_per_area,
                               EAB_CWD_Paired$NoCWDperha_Ash,
                               method = "spearman")
cor.noapha.noacwdpha
#not significant
#snag - normal so use pearson
cor.noapha.noaspha<-cor.test(EAB_CWD_Paired$No_ash_per_area,
                            EAB_CWD_Paired$Noperha_Ash_Snag,
                            method="pearson")
cor.noapha.noaspha
#significant postive
#Visualize, convert stream area to ha
EAB_CWD_Paired$No_ash_per_ha<-1000*EAB_CWD_Paired$No_ash_per_area
ggplot(EAB_CWD_Paired, aes(y=No_ash_per_ha, x=Noperha_Ash_Snag)) +
  geom_point(size=3) +
  xlab("Ash Snag Frequency (Number per ha)")+
  ylab("Aquatic Ash CWD Frequency (Number per ha)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none")

#prop ash pieces
#test normality
shapiro.test(EAB_CWD_Paired$prop_ash_pieces)
#W = 0.96054, p-value = 0.8239, normally distributed
#CWD - normal
cor.propap.propacwdp<- cor.test(EAB_CWD_Paired$prop_ash_pieces,
                                EAB_CWD_Paired$prop_Ash_CWD_no,
                                method = "pearson")
cor.propap.propacwdp
#not significant
#snag - normal so use pearson
cor.propap.propasp<-cor.test(EAB_CWD_Paired$prop_ash_pieces,
                             EAB_CWD_Paired$prop_Ash_snag_no,
                             method="pearson")
cor.propap.propasp
#not significant

#vol ash per stream length
#test normality
shapiro.test(EAB_CWD_Paired$Vol_ash_per_length)
#W = 0.95043, p-value = 0.7438, normally distributed
#CWD - log10 transform
cor.volapm.volacwdpha<- cor.test(EAB_CWD_Paired$Vol_ash_per_length,
                                EAB_CWD_Paired$log10VolCWD_m3perha_Ash,
                                method = "pearson")
cor.volapm.volacwdpha
#not significant
#snag - log10 transformation
cor.volapm.baaspha<-cor.test(EAB_CWD_Paired$Vol_ash_per_length,
                             EAB_CWD_Paired$log10BAperha_Ash_Snag,
                             method="pearson")
cor.volapm.baaspha
#not significant

#vol ash per stream area
#test normality
shapiro.test(EAB_CWD_Paired$Vol_ash_per_area)
#W = 0.97827, p-value = 0.9426, normally distributed
#CWD - log10 transform
cor.volapha.volacwdpha<- cor.test(EAB_CWD_Paired$Vol_ash_per_area,
                                 EAB_CWD_Paired$log10VolCWD_m3perha_Ash,
                                 method = "pearson")
cor.volapha.volacwdpha
#not significant
#snag - log10 transformation
cor.volapha.baaspha<-cor.test(EAB_CWD_Paired$Vol_ash_per_area,
                             EAB_CWD_Paired$log10BAperha_Ash_Snag,
                             method="pearson")
cor.volapha.baaspha
#not significant

#prop ash volume
#test normality
shapiro.test(EAB_CWD_Paired$prop_ash_vol)
#W = 0.95822, p-value = 0.8059, normally distributed
#CWD - normal
cor.propavol.propacwdvol<- cor.test(EAB_CWD_Paired$prop_ash_vol,
                                  EAB_CWD_Paired$prop_Ash_CWD_vol,
                                  method = "pearson")
cor.propavol.propacwdvol
#not significant
#snag - normal
cor.propavol.propasba<-cor.test(EAB_CWD_Paired$prop_ash_vol,
                              EAB_CWD_Paired$prop_Ash_snag_BA,
                              method="pearson")
cor.propavol.propasba
#not significant

#######################################
#Summary of correlation analyses:
#no snags pos cor no logs per area
#no cwd cor neg no log jams per length
#no cwd cor neg no log jams per area
#no snags neg cor no log jam pieces per length
#no snags neg cor prop pieces in log jams
#vol cwd pos cor vol insitu per length
#vol cwd pos cor prop insitu vol
#no ash snag pos cor no ash per area
#only include per area in dissertation
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

#Move to correlation analysis

##ALL richness
#test normality
shapiro.test(EAB_Ter_H2o_aqcwd_ALL$ALL_Richness)
#W = 0.94056, p-value = 0.009839, not normal use non parametric
#Terrestrial richness
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
EAB_16S_map <- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_Map_Filtered_R.csv", header=T)
#Combine map with other env variables
EAB_16S_map$Source<-revalue(EAB_16S_map$Source, c("Aquatic Leaf Litter"="ALL", "Live Leaves"="LL", "Terrestrial Litter"="TL"))
#for just aquatic
EAB_Ter_H2o_cwd_ALL_mics<-merge(EAB_Ter_H2o_aqcwd_ALL, EAB_16S_map, by=c("Stream","Watershed","Gap_location","Days_since_start"), all=T, no.dups=T)
EAB_Ter_H2o_cwd_ALL_mics<-subset(EAB_Ter_H2o_cwd_ALL_mics, SampleID!="NA")
row.names(EAB_Ter_H2o_cwd_ALL_mics)<-EAB_Ter_H2o_cwd_ALL_mics$SampleID
#create new Gap variable in EAB_Ter_H2o_cwd_ALL_mics with Y,N to mirror terrestrial leaf litte rmicrobes
EAB_Ter_H2o_cwd_ALL_mics$Gap<-EAB_Ter_H2o_cwd_ALL_mics$Gap_location
EAB_Ter_H2o_cwd_ALL_mics$Gap<-revalue(EAB_Ter_H2o_cwd_ALL_mics$Gap, c("Downstream"="N", "Gap"="Y", "Upstream"="N", "Woods"="N"))

#Upload OTU table
EAB_16S_OTU<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_OTUs.csv", header=T)
#Format data frame so the OTU is row name
row.names(EAB_16S_OTU)<-EAB_16S_OTU[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_OTU$ID<-NULL
#Delete sample 56
EAB_16S_OTU$EAB_ALL_2016_56<-NULL
EAB_16S_OTU_t<-t(EAB_16S_OTU)
EAB_16S_OTU_t<-data.frame(EAB_16S_OTU_t, check.names = FALSE)

#merge OTUs and delete 56
EAB_Paired_TA_OTU<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_OTU_t, by=0, no.dups=T)

#Combine with environmental variables
EAB_16S_OTU_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_OTU_t, by=0)
#delete OTUs not found in any samples
names(EAB_16S_OTU_map)
EAB_16S_OTU<-EAB_16S_OTU_map[,177:ncol(EAB_16S_OTU_map)]
EAB_16S_OTU<-EAB_16S_OTU[, colSums(EAB_16S_OTU != 0) > 0]
#3182 OTUs found in all samples

#Upload family level taxonomy
EAB_16S_f<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_f.csv", header=T)
#Format data frame so the family is row name
row.names(EAB_16S_f)<-EAB_16S_f[,1]
#Delete otu id column, now that otu id is row name
EAB_16S_f$ID<-NULL
#Delete sample 56
EAB_16S_f$EAB_ALL_2016_56<-NULL
EAB_16S_f_t<-t(EAB_16S_f)
EAB_16S_f_t<-data.frame(EAB_16S_f_t, check.names = FALSE)

#Combine with environmental variables
EAB_16S_f_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_f_t, by=0)
#delete families not found in aquatic samples
names(EAB_16S_f_map)
EAB_16S_f<-EAB_16S_f_map[,177:ncol(EAB_16S_f_map)]
EAB_16S_f<-EAB_16S_f_aq[, colSums(EAB_16S_f_aq != 0) > 0]
sort(colSums(EAB_16S_f))
#most common family is Bacillaceae
stat.desc(EAB_16S_f$`k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae`)
#220 +/- 23
EAB_16S_f_map_GF<-EAB_16S_f_map$Gap
EAB_16S_f_map_W<-EAB_16S_f_map$Watershed

#Upload phylogenetic diversity
EAB_16S_fpd<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_fpd.csv", header=T)
#Format data frame so the ID is row name
row.names(EAB_16S_fpd)<-EAB_16S_fpd[,1]
#Delete sample 56
EAB_16S_fpd<-subset(EAB_16S_fpd, X!="EAB_ALL_2016_56")
#Delete otu id column, now that otu id is row name
EAB_16S_fpd$X<-NULL
#Combine with environmental variables
EAB_16S_fpd_map<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_fpd, by=0)

####Terrestrial connection microbes - all habitats

#Model watershed + gap + source 1|stream mixed effect

#Faith's PD for all habitats
range(EAB_16S_fpd_map$faith_pd)
#2.25-32.78
#build model
lmer.FPD= lme(faith_pd~Watershed+Gap+Source, random=~1|Stream,
               data=EAB_16S_fpd_map,
               method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.FPD),col="darkgray")
shapiro.test(residuals(lmer.FPD))
#W = 0.90199, p-value = 5.231e-05,not normal, transform
EAB_16S_fpd_map$log10faith_pd<-log10(EAB_16S_fpd_map$faith_pd)
#build model
lmer.logFPD= lme(log10faith_pd~Watershed+Gap+Source, random=~1|Stream,
              data=EAB_16S_fpd_map,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.logFPD),col="darkgray")
shapiro.test(residuals(lmer.logFPD))
#W = 0.986, p-value = 0.6361 normal
summary(lmer.logFPD)
#intercept and source significant
#for interpretability, see summary for not transformed
summary(lmer.FPD)
anova(lmer.FPD)
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
EAB_16S_cha<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_chao.csv", header=T)
#Format data frame so the ID is row name
row.names(EAB_16S_cha)<-EAB_16S_cha[,1]
#Delete sample 56
EAB_16S_cha<-subset(EAB_16S_cha, X!="EAB_ALL_2016_56")
#Delete otu id column, now that otu id is row name
EAB_16S_cha$X<-NULL
#Combine with environmental variables
EAB_16S_cha_map <-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_cha, by=0)

#start modelling
range(EAB_16S_cha_map$chao1)
#27-728
#build model
lmer.cha= lme(chao1~Watershed+Gap+Source, random=~1|Stream,
              data=EAB_16S_cha_map,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cha),col="darkgray")
shapiro.test(residuals(lmer.cha))
#W = 0.86739, p-value = 2.859e-06,not normal, transform
EAB_16S_cha_map$log10chao1<-log10(EAB_16S_cha_map$chao1)
#build model
lmer.logcha= lme(log10chao1~Watershed+Gap+Source, random=~1|Stream,
                 data=EAB_16S_cha_map,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.logcha),col="darkgray")
shapiro.test(residuals(lmer.logcha))
#W = 0.98882, p-value = 0.7987 normal
summary(lmer.logcha)
#intercept and source significant
summary(lmer.cha)
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
EAB_16S_uni<-read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_WUnifrac.tsv", sep="\t", header = T, row.names=1)

#Combine distance matrix with environmental variables
row.names(EAB_16S_uni)
row.names(EAB_Ter_H2o_cwd_ALL_mics)
EAB_16S_uni_map<-merge(EAB_Ter_H2o_cwd_ALL_mics, EAB_16S_uni, by="row.names")
row.names(EAB_16S_uni_map)<-EAB_16S_uni_map[,1]
EAB_16S_uni_map<-EAB_16S_uni_map[,-c(1)]
#delete EAB_ALL_56_2016 because error in sequencing
names(EAB_16S_uni_map)
EAB_16S_uni<-as.matrix(EAB_16S_uni_map[,c(176:205,207:246)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env<-EAB_16S_uni_map[,1:175]
EAB_16S_uni_env$Gap<-revalue(EAB_16S_uni_env$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_GF<-EAB_16S_uni_env$Gap
EAB_16S_uni_env_GFS<-as.factor(paste(EAB_16S_uni_env$Gap.y,EAB_16S_uni_env$Source))
EAB_16S_uni_env_S<-as.factor(EAB_16S_uni_env$Source)
EAB_16S_uni_env_W<-as.factor(EAB_16S_uni_env$Watershed)

#UNI-Overall permanova with unifrac distances
adonis(as.dist(EAB_16S_uni) ~ (Source+Watershed+Gap)^2+Stream, data=EAB_16S_uni_env,
       permutations=999)
#source and source:watershed interaction significant

#Visualize via nmds
EAB_Paired_Mic_NMDS<-metaMDS(as.dist(EAB_16S_uni))

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mic_NMDS)

#NMDS plot for source
ordiplot(EAB_Paired_Mic_NMDS, type="n")
with(EAB_Paired_Mic_NMDS, points(EAB_Paired_Mic_NMDS, display="sites", col=source_col_vec[EAB_16S_uni_env_S], pch=19))
with(EAB_Paired_Mic_NMDS, legend("topleft", legend=levels(EAB_16S_uni_env_S), bty="n", col=source_col_vec, pch=19, pt.bg=source_col_vec))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "ALL"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#b15928", show.groups = "LL"))
with(EAB_Paired_Mic_NMDS, ordiellipse(EAB_Paired_Mic_NMDS, EAB_16S_uni_env_S, kind="se", conf=0.95, lwd=2, col="#33a02c", show.groups = "TL"))

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
#11 indicator families for watershed

#now move on to venn diagrams to find unique OTUs
#create venn diagram for shared OTUs among terrestrial/aquatic gap/woods.
#first list OTU names found in each group
names(EAB_Paired_TA_OTU)
EAB_16S_OTU_map_ag_venn<-aggregate(EAB_Paired_TA_OTU[177:ncol(EAB_Paired_TA_OTU)], 
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
#only one otu is common between all gap sample types and not found in gap
GapVen$`AL Gap  _TLL Gap_TL Gap`
#"d9bf12bfe7c9ef750e953c02b9244ea1"
#find relative abundances
stat.desc(EAB_Paired_TA_OTU$d9bf12bfe7c9ef750e953c02b9244ea1)
#mean 0.39
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
GapVen3<-venn_diagram3(YAquaticLeafLitter,YLiveLeaves,YTerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
ForVen3<-venn_diagram3(NAquaticLeafLitter,NLiveLeaves,NTerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
GapVenOTU<- summarySE(EAB_Paired_TA_OTU, measurevar="d9bf12bfe7c9ef750e953c02b9244ea1", groupvars=c("Source","Days_since_start"))
ggplot(GapVenOTU, aes(x=Days_since_start, y=d9bf12bfe7c9ef750e953c02b9244ea1, color=Source)) +
  geom_errorbar(aes(ymin=d9bf12bfe7c9ef750e953c02b9244ea1-se, ymax=d9bf12bfe7c9ef750e953c02b9244ea1+se), width=.4) +
  geom_point(size=1.5) +
  geom_smooth(method=lm, se=F)+
  xlab("Time")+
  ylab("Gap OTU Abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=18),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
#ONly found once in litter for each sample
EAB_Paired_Taxonomy<- read.csv("~/Documents/MSU/Research/Surveying/Paired_Sites/EAB_Paired_taxonomy.csv", header=T)
subset(EAB_Paired_Taxonomy, FeatureID=="d9bf12bfe7c9ef750e953c02b9244ea1")
#k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Methylobacteriaceae; g__; s__

quartz()
WoodsVen<-venn_diagram4(NAquaticLeafLitter, NLiveLeaves, NTerrestrialLitter, YesGap,
                        "AL For", "TLL For", "TL For", "Gap",
                        colors=woodsven_col_vec)
#only one otu is common between all gap sample types and not found in gap
WoodsVen$`AL For_TLL For_TL For`
#"a1fe4f3c8d68400364e4b68adfc2a786"
WoodsVenOTU<- summarySE(EAB_Paired_TA_OTU, measurevar="a1fe4f3c8d68400364e4b68adfc2a786", groupvars=c("Source","Days_since_start"))
ggplot(WoodsVenOTU, aes(x=Days_since_start, y=a1fe4f3c8d68400364e4b68adfc2a786, color=Source)) +
  geom_errorbar(aes(ymin=a1fe4f3c8d68400364e4b68adfc2a786-se, ymax=a1fe4f3c8d68400364e4b68adfc2a786+se), width=.4) +
  geom_point(size=1.5) +
  geom_smooth(method=lm, se=F)+
  xlab("Time")+
  ylab("Woods OTU Abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=18),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
#similarly, only in low abundance
subset(EAB_Paired_Taxonomy, FeatureID=="a1fe4f3c8d68400364e4b68adfc2a786")
#k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Novosphingobium
Ven3<-venn_diagram3(AquaticLeafLitter,LiveLeaves,TerrestrialLitter, "Aquatic Litter",
                       "Terrestrial Litter","Live Leaves",colors=source_col_vec)
venn.diagram(x = list(YLiveLeaves,NLiveLeaves),category.names = c("Forest" , "Gap"),
  filename = 'LiveLivesGapForest.png',output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,
  scaled=FALSE)
venn.diagram(x = list(YTerrestrialLitter,NTerrestrialLitter),
             category.names = c("Forest" , "Gap"),filename = 'TerrestrialLitterGapForest.png',
             output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,scaled=FALSE)
venn.diagram(x = list(YAquaticLeafLitter,NAquaticLeafLitter),
             category.names = c("Forest" , "Gap"),filename = 'AquaticLitterGapForest.png',
             output=TRUE,fill=GF_col_vec,cat.cex =1.5,cex=1.5,scaled=FALSE)
EAB_16S_OTU_map_ag_venn_ALL<-aggregate(EAB_Paired_TA_OTU[177:ncol(EAB_Paired_TA_OTU)], 
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
#################################
#Summary all sources microbes
#live leaves and terrestrial leaf litter have significantly higher fpd than aquatic leaf litter
#terrestrial leaf litter higher chao1 richness compared to aquatic leaf litter
#source and source:watershed interaction influence community structure
#11 indicator families for watershed
#only 2 otus only found in gap or forest, not the other
##############################

#Now model each habitat separately

#Live leaves
#Faith's PD for live leaves
EAB_16S_fpd_map_LL<-subset(EAB_16S_fpd_map, Source=="LL")
range(EAB_16S_fpd_map_LL$faith_pd)
#5.365991-14.963764
#build model
lmer.FPD.LL= lme(faith_pd~Watershed+Gap, random=~1|Stream,
              data=EAB_16S_fpd_map_LL,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.FPD.LL),col="darkgray")
shapiro.test(residuals(lmer.FPD.LL))
#W = 0.97321, p-value = 0.9207,normal
summary(lmer.FPD.LL)
#intercept significant

#Chao1 for live leaves
EAB_16S_cha_map_LL<-subset(EAB_16S_cha_map, Source=="LL")
#start modelling
range(EAB_16S_cha_map_LL$chao1)
#27.33333-134.53846
#build model
lmer.cha.ll= lme(chao1~Watershed+Gap, random=~1|Stream,
              data=EAB_16S_cha_map_LL,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cha.ll),col="darkgray")
shapiro.test(residuals(lmer.cha.ll))
#W = 0.93988, p-value = 0.5806,normal
summary(lmer.cha.ll)
#intercept significant

#Now look for community level differences using PERMANOVA and adonis in LL

#subset unifrac map for live leaves
EAB_16S_uni_map_LL<-subset(EAB_16S_uni_map,Source=="LL")
names(EAB_16S_uni_map_LL)
EAB_16S_uni_LL<-as.matrix(EAB_16S_uni_map_LL[,c(227:235)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_LL<-EAB_16S_uni_map_LL[,1:175]
EAB_16S_uni_env_LL$Gap<-revalue(EAB_16S_uni_env_LL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_LL_GF<-EAB_16S_uni_env_LL$Gap
EAB_16S_uni_env_LL_W<-as.factor(EAB_16S_uni_env_LL$Watershed)

#UNI-Overall permanova with unifrac distances
adonis(as.dist(EAB_16S_uni_LL) ~ (Watershed+Gap)^2+Stream, data=EAB_16S_uni_env_LL,
       permutations=999)
#watershed significant

#Visualize via nmds
EAB_Paired_Mic_NMDS_LL<-metaMDS(as.dist(EAB_16S_uni_LL))

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mic_NMDS_LL)

#NMDSplot for watershed
ordiplot(EAB_Paired_Mic_NMDS_LL, type="n", xlim=c(-0.4,0.5), ylim=c(-0.2,0.15))
with(EAB_Paired_Mic_NMDS_LL, points(EAB_Paired_Mic_NMDS_LL, display="sites", col=Watershed_col_vec[EAB_16S_uni_env_LL_W], pch=19))
with(EAB_Paired_Mic_NMDS_LL, legend("topleft", legend=levels(EAB_16S_uni_env_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_W, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Grand"))
with(EAB_Paired_Mic_NMDS_LL, ordiellipse(EAB_Paired_Mic_NMDS_LL, EAB_16S_uni_env_LL_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#indicator species analysis for watershed
EAB_16S_f_map_LL<-subset(EAB_16S_f_map,Source=="LL")
names(EAB_16S_f_map_LL)
EAB_16S_f_LL<-EAB_16S_f_map_LL[177:ncol(EAB_16S_f_map_LL)]
EAB_16S_f_map_LL_W<-EAB_16S_f_map_LL$Watershed
EAB_Mic_Com_LL_W_indic<-signassoc(EAB_16S_f_LL, cluster=EAB_16S_f_map_LL_W,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_LL_W_indic_sig<-subset(EAB_Mic_Com_LL_W_indic, psidak<=0.05)
#2 indicator families for watershed
rownames(EAB_Mic_Com_LL_W_indic_sig)
#k__Bacteria;p__Armatimonadetes;c__[Fimbriimonadia];o__[Fimbriimonadales];f__[Fimbriimonadaceae] - Kalamazoo
#k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae - Kalamazoo
#summary stats
EAB_16S_f_map_LL %>%
  group_by(Watershed) %>%
  get_summary_stats('k__Bacteria;p__Armatimonadetes;c__[Fimbriimonadia];o__[Fimbriimonadales];f__[Fimbriimonadaceae]', type = "mean_se")
#Clinton 0(0),Grand 0.3(0.3),Kalamazoo 4.25(0,25)
EAB_16S_f_map_LL %>%
  group_by(Watershed) %>%
  get_summary_stats('k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae', type = "mean_se")
#Clinton 0.5(0.5),Grand 1.3(0.3),Kalamazoo 56(21)

#################################
#Summary live leave microbes
#FPd and chao not altered
#watershed significant influence community structure
#2 indicator families greatest in kalamazoo watershed
##############################

#Terrestrial litter
#Faith's PD for terrestrial litter
EAB_16S_fpd_map_TL<-subset(EAB_16S_fpd_map, Source=="TL")
range(EAB_16S_fpd_map_TL$faith_pd)
#9.331483-32.777676
#build model
lmer.FPD.TL= lme(faith_pd~Watershed+Gap, random=~1|Stream,
                 data=EAB_16S_fpd_map_TL,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.FPD.TL),col="darkgray")
shapiro.test(residuals(lmer.FPD.TL))
#W = 0.96691, p-value = 0.8537,normal
summary(lmer.FPD.TL)
#intercept significant

#Chao1 for terrestrial leaves
EAB_16S_cha_map_TL<-subset(EAB_16S_cha_map, Source=="TL")
#start modelling
range(EAB_16S_cha_map_TL$chao1)
#81.5000-728.3229
#build model
lmer.cha.tl= lme(chao1~Watershed+Gap, random=~1|Stream,
                 data=EAB_16S_cha_map_TL,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cha.tl),col="darkgray")
shapiro.test(residuals(lmer.cha.tl))
#W = 0.96653, p-value = 0.8496,normal
summary(lmer.cha.tl)
#intercept significant

#Now look for community level differences using PERMANOVA and adonis in TL

#subset unifrac map for terrestrial litter
EAB_16S_uni_map_TL<-subset(EAB_16S_uni_map,Source=="TL")
names(EAB_16S_uni_map_TL)
EAB_16S_uni_TL<-as.matrix(EAB_16S_uni_map_TL[,c(236:246)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_TL<-EAB_16S_uni_map_TL[,1:175]
EAB_16S_uni_env_TL$Gap<-revalue(EAB_16S_uni_env_TL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_TL_GF<-EAB_16S_uni_env_TL$Gap
EAB_16S_uni_env_TL_W<-as.factor(EAB_16S_uni_env_TL$Watershed)

#UNI-Overall permanova with unifrac distances
adonis(as.dist(EAB_16S_uni_TL) ~ (Watershed+Gap)^2+Stream, data=EAB_16S_uni_env_TL,
       permutations=999)
#watershed*gap interaction significant

#Visualize via nmds
EAB_Paired_Mic_NMDS_TL<-metaMDS(as.dist(EAB_16S_uni_TL))

#Stressplot Nmds
stressplot(EAB_Paired_Mic_NMDS_TL)

#NMDSplot for watershed
ordiplot(EAB_Paired_Mic_NMDS_TL, type="n",ylim=c(-0.3,0.2))
with(EAB_Paired_Mic_NMDS_TL, points(EAB_Paired_Mic_NMDS_TL, display="sites", col=Watershed_col_vec[EAB_16S_uni_env_TL_W], pch=19))
with(EAB_Paired_Mic_NMDS_TL, legend("topleft", legend=levels(EAB_16S_uni_env_TL_W), bty="n", col=Watershed_col_vec, pch=19, pt.bg=Watershed_col_vec))
with(EAB_Paired_Mic_NMDS_TL, ordiellipse(EAB_Paired_Mic_NMDS_TL, EAB_16S_uni_env_TL_W, kind="se", conf=0.95, lwd=2, col="#af8dc3", show.groups = "Clinton"))
with(EAB_Paired_Mic_NMDS_TL, ordiellipse(EAB_Paired_Mic_NMDS_TL, EAB_16S_uni_env_TL_W, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Grand"))
with(EAB_Paired_Mic_NMDS_TL, ordiellipse(EAB_Paired_Mic_NMDS_TL, EAB_16S_uni_env_TL_W, kind="se", conf=0.95, lwd=2, col="#7fbf7b", show.groups = "Kalamazoo"))

#NMDSplot for gap/forest
ordiplot(EAB_Paired_Mic_NMDS_TL, type="n")
with(EAB_Paired_Mic_NMDS_TL, points(EAB_Paired_Mic_NMDS_TL, display="sites", col=GF_col_vec[EAB_16S_uni_env_TL_GF], pch=19))
with(EAB_Paired_Mic_NMDS_TL, legend("topleft", legend=levels(EAB_16S_uni_env_TL_GF), bty="n", col=GF_col_vec, pch=19, pt.bg=GF_col_vec))
with(EAB_Paired_Mic_NMDS_TL, ordiellipse(EAB_Paired_Mic_NMDS_TL, EAB_16S_uni_env_TL_GF, kind="se", conf=0.95, lwd=2, col="#66c2a5", show.groups = "Forest"))
with(EAB_Paired_Mic_NMDS_TL, ordiellipse(EAB_Paired_Mic_NMDS_TL, EAB_16S_uni_env_TL_GF, kind="se", conf=0.95, lwd=2, col="#fdae61", show.groups = "Gap"))

#indicator species analysis for watershed
EAB_16S_f_map_TL<-subset(EAB_16S_f_map,Source=="TL")
names(EAB_16S_f_map_TL)
EAB_16S_f_TL<-EAB_16S_f_map_TL[177:ncol(EAB_16S_f_map_TL)]
EAB_16S_f_map_TL_W<-EAB_16S_f_map_TL$Watershed
EAB_Mic_Com_TL_W_indic<-signassoc(EAB_16S_f_TL, cluster=EAB_16S_f_map_TL_W,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_TL_W_indic_sig<-subset(EAB_Mic_Com_TL_W_indic, psidak<=0.05)
#2 indicator families for watershed both grand
rownames(EAB_Mic_Com_TL_W_indic_sig)
#"k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Patulibacteraceae" - grand
#"k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae" - grand
#summary stats
EAB_16S_f_map_TL %>%
  group_by(Watershed) %>%
  get_summary_stats("k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Patulibacteraceae", type = "mean_se")
#Clinton 7.8(1.6),Grand 0(0),Kalamazoo 12(6.3)
EAB_16S_f_map_TL %>%
  group_by(Watershed) %>%
  get_summary_stats("k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae", type = "mean_se")
#Clinton 42(17),Grand 1(1),Kalamazoo 21(7)

#indicator species analysis for gapforest
EAB_16S_f_map_TL_GF<-EAB_16S_f_map_TL$Gap
EAB_Mic_Com_TL_GF_indic<-signassoc(EAB_16S_f_TL, cluster=EAB_16S_f_map_TL_GF,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_TL_GF_indic_sig<-subset(EAB_Mic_Com_TL_GF_indic, psidak<=0.05)
#No indicator families

#################################
#Summary terrestrial leaf litter microbes
#FPd and chao not altered
#watershed*gap interaction significant influence community structure
#2 indicator families abscent in grand watershed
##############################

#Aquatic Leaf Litter

#Faith's PD for aquatic leaf litter
EAB_16S_fpd_map_ALL<-subset(EAB_16S_fpd_map, Source=="ALL")
range(EAB_16S_fpd_map_ALL$faith_pd)
#2.244885-19.462688
#build model
lmer.FPD.ALL= lme(faith_pd~Watershed+Gap_location+ALL_Richness, random=~1|Stream,
                 data=EAB_16S_fpd_map_ALL,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.FPD.ALL),col="darkgray")
shapiro.test(residuals(lmer.FPD.ALL))
#W = 0.87069, p-value = 6.908e-05,not normal, try with log transformed
#build model
lmer.FPD.log10ALL= lme(log10faith_pd~Watershed+Gap_location+ALL_Richness, random=~1|Stream,
                  data=EAB_16S_fpd_map_ALL,
                  method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.FPD.log10ALL),col="darkgray")
shapiro.test(residuals(lmer.FPD.log10ALL))
#W = 0.9773, p-value = 0.458,normal
summary(lmer.FPD.log10ALL)
#intercept significant

#Chao1 for aquatic leaf litter
EAB_16S_cha_map_ALL<-subset(EAB_16S_cha_map, Source=="ALL")
#start modelling
range(EAB_16S_cha_map_ALL$chao1)
#37.2500-338.1667
#build model
lmer.cha.all= lme(chao1~Watershed+Gap_location+ALL_Richness, random=~1|Stream,
                 data=EAB_16S_cha_map_ALL,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cha.all),col="darkgray")
shapiro.test(residuals(lmer.cha.all))
#W = 0.87517, p-value = 9.335e-05, not normal, try log transformed
#build model
lmer.log10cha.all= lme(log10chao1~Watershed+Gap_location+ALL_Richness, random=~1|Stream,
                  data=EAB_16S_cha_map_ALL,
                  method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10cha.all),col="darkgray")
shapiro.test(residuals(lmer.log10cha.all))
#W = 0.9806, p-value = 0.5908, normal
summary(lmer.log10cha.all)
#intercept significant

#Now look for community level differences using PERMANOVA and adonis in all

#subset unifrac map for terrestrial litter
EAB_16S_uni_map_ALL<-subset(EAB_16S_uni_map,Source=="ALL")
names(EAB_16S_uni_map_ALL)
EAB_16S_uni_ALL<-as.matrix(EAB_16S_uni_map_ALL[,c(176:205,207:226)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
EAB_16S_uni_env_ALL<-EAB_16S_uni_map_ALL[,1:175]
EAB_16S_uni_env_ALL$Gap<-revalue(EAB_16S_uni_env_ALL$Gap, c("N"="Forest", "Y"="Gap"))
EAB_16S_uni_env_ALL_GF<-EAB_16S_uni_env_ALL$Gap
EAB_16S_uni_env_ALL_W<-as.factor(EAB_16S_uni_env_ALL$Watershed)
EAB_16S_uni_env_ALL_GL<-factor(EAB_16S_uni_env_ALL$Gap_location)

#UNI-Overall permanova with unifrac distances
adonis(as.dist(EAB_16S_uni_ALL) ~ (Watershed+Gap_location+ALL_Richness)^2+Stream, data=EAB_16S_uni_env_ALL,
       permutations=999)
#all richness and all richness*gaplocation interaction significant

#Visualize via nmds
EAB_Paired_Mic_NMDS_ALL<-metaMDS(as.dist(EAB_16S_uni_ALL))

#Stressplot macroinvertebrate Nmds
stressplot(EAB_Paired_Mic_NMDS_ALL)

#NMDSplot for gap location
ordiplot(EAB_Paired_Mic_NMDS_ALL, type="n")
with(EAB_Paired_Mic_NMDS_ALL, points(EAB_Paired_Mic_NMDS_ALL, display="sites", col=Gap_location_col_vec[EAB_16S_uni_env_ALL_GL], pch=19))
with(EAB_Paired_Mic_NMDS_ALL, legend("topleft", legend=levels(EAB_16S_uni_env_ALL_GL), bty="n", col=Gap_location_col_vec, pch=19, pt.bg=Gap_location_col_vec))
with(EAB_Paired_Mic_NMDS_ALL, ordiellipse(EAB_Paired_Mic_NMDS_ALL, EAB_16S_uni_env_ALL_GL, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(EAB_Paired_Mic_NMDS_ALL, ordiellipse(EAB_Paired_Mic_NMDS_ALL, EAB_16S_uni_env_ALL_GL, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(EAB_Paired_Mic_NMDS_ALL, ordiellipse(EAB_Paired_Mic_NMDS_ALL, EAB_16S_uni_env_ALL_GL, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))

#indicator species analysis for gap location
EAB_16S_f_map_ALL<-subset(EAB_16S_f_map,Source=="ALL")
names(EAB_16S_f_map_ALL)
EAB_16S_f_ALL<-EAB_16S_f_map_ALL[177:ncol(EAB_16S_f_map_ALL)]
EAB_16S_f_map_ALL_GL<-factor(EAB_16S_f_map_ALL$Gap_location)
EAB_Mic_Com_ALL_GL_indic<-signassoc(EAB_16S_f_ALL, cluster=EAB_16S_f_map_ALL_GL,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Mic_Com_ALL_GL_indic_sig<-subset(EAB_Mic_Com_ALL_GL_indic, psidak<=0.05)
#1 indicator family for gap locaiton downstream
rownames(EAB_Mic_Com_ALL_GL_indic_sig)
#"k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__Gemmatales;f__Gemmataceae" - downstream
#summary stats
EAB_16S_f_map_ALL %>%
  group_by(Gap_location) %>%
  get_summary_stats("k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__Gemmatales;f__Gemmataceae", type = "mean_se")
#US 0, g 0, ds 0.375 (0.155)

#################################
#Summary aquatic leaf litter microbes
#FPd and chao not altered
#all richness and all richness*gaplocation influence community structure
#1 indicator family for gap location Gemmataceae which only found in ds locations
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
EAB_M_Matrix <- dcast(Paired_Macroinvertebrates, Stream + Date + Gap + Gap_location ~ Taxonomy, value.var = "Count", fun.aggregate =sum)
#Find total invertebrate #
sum(EAB_M_Matrix[,5:ncol(EAB_M_Matrix)])
#Find most abundant taxa
sort(colSums(EAB_M_Matrix[,5:ncol(EAB_M_Matrix)]))
#most abundant taxa elmidae with 519 individuals
stat.desc(EAB_M_Matrix$`InsectaColeopteraElmidae`)
#make relative abundances
str(EAB_M_Matrix)
EAB_M_MatrixRA<-data.frame(make_relative(as.matrix(EAB_M_Matrix[,5:ncol(EAB_M_Matrix)])))
EAB_M_MatrixRA[is.na(EAB_M_MatrixRA)] <- 0
str(EAB_M_MatrixRA)
sort(colMeans(EAB_M_MatrixRA))
stat.desc(EAB_M_MatrixRA$InsectaColeopteraElmidae)

#combine eab_Paired_ALL_Stand with macs
EAB_Paired_Macs_Stand<-merge(EAB_Ter_H2o_aqcwd_ALL, EAB_M_Matrix, by=c("Stream","Gap_location","Date"))
#split up into environmental and community
names(EAB_Paired_Macs_Stand)
EAB_Paired_Macs_Stand_com<-EAB_Paired_Macs_Stand[,167:ncol(EAB_Paired_Macs_Stand)]
EABMacroComtotbray0Stand<-as.matrix(bray0(EAB_Paired_Macs_Stand_com))
EAB_Paired_Macs_Stand_env<-EAB_Paired_Macs_Stand[,1:166]
EAB_Paired_Macs_Stand_GA<-EAB_Paired_Macs_Stand_env$Gap_Area
EAB_Paired_Macs_Stand_GA_cat<-as.factor(EAB_Paired_Macs_Stand_GA)
EAB_Paired_Macs_Stand_YGF<-EAB_Paired_Macs_Stand_env$Year_Gapformation
EAB_Paired_Macs_Stand_YGF_cat<-as.factor(EAB_Paired_Macs_Stand_YGF)
EAB_Paired_Macs_Stand_GL<-EAB_Paired_Macs_Stand_env$Gap_location
EAB_Paired_Macs_Stand_W<-EAB_Paired_Macs_Stand_env$Watershed

#analyze community
adonis(as.dist(EABMacroComtotbray0Stand) ~ Watershed*Gap_location+Stream,
       data=EAB_Paired_Macs_Stand_env,permutations=9999)
#Watershed and stream significant

#Watershed
EAB_Inverts_Com_W_indic<-signassoc(EAB_Paired_Macs_Stand_com, cluster=EAB_Paired_Macs_Stand_W,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Inverts_Com_W_indic_sig<-subset(EAB_Inverts_Com_W_indic, psidak<=0.05)
#10 families indicate watershed:
#Clinton: elmidae, ephemerellidae, heptageniidae, nemouridae
#Grand: Gammaridae
#Kalamazoo: Psephenidae, Baetidae, Taeniopterygidae, hydropsychidae, rhacophilidae

#Gap Location
EAB_Inverts_Com_GL_indic<-signassoc(EAB_Paired_Macs_Stand_com, cluster=EAB_Paired_Macs_Stand_GL,  mode=0, alternative = "two.sided",control = how(nperm=999))
EAB_Inverts_Com_GL_indic_sig<-subset(EAB_Inverts_Com_GL_indic, psidak<=0.05)
#0 families indicate gap location

#so work with these 10 families for quantitative responses:
#Determine which are in all watersheds
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaColeopteraElmidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaColeopteraElmidae", type = "mean_se")
#in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraEphemerellidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaEphemeropteraEphemerellidae", type = "mean_se")
#missing in stoney creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraHeptageniidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaEphemeropteraHeptageniidae", type = "mean_se")
#missing in spring creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaPlecopteraNemouridae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaPlecopteraNemouridae", type = "mean_se")
#missing in spring creek and stoney creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("MalacostracaAmphipodaGammaridae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("MalacostracaAmphipodaGammaridae", type = "mean_se")
#in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaColeopteraPsephenidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaColeopteraPsephenidae", type = "mean_se")
#missing in frayer, sessions and spring
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaEphemeropteraBaetidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaEphemeropteraBaetidae", type = "mean_se")
#missing in spring creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaPlecopteraTaeniopterygidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaPlecopteraTaeniopterygidae", type = "mean_se")
#only found in augusta and seven mile creek
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaTrichopteraHydropsychidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaTrichopteraHydropsychidae", type = "mean_se")
#found in all streams
EAB_Paired_Macs_Stand %>%
  group_by(Stream) %>%
  get_summary_stats("InsectaTrichopteraRhyacophilidae", type = "mean_se")
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("InsectaTrichopteraRhyacophilidae", type = "mean_se")
#missing in all but sessions and seven mile

#so work with hydropsychidae, Gammaridae and Elmidae

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
lmer.mr= lme(Mac_Richness~Watershed+Gap_location, random=~1|Stream,
                  data=EAB_Paired_Macs_Stand,
                  method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.mr),col="darkgray")
shapiro.test(residuals(lmer.mr))
#W = 0.97617, p-value = 0.354, normal
summary(lmer.mr)
#intercept and downstream significant
#visualize
ggplot(EAB_Paired_Macs_Stand, aes(x=Watershed, y=Mac_Richness, color=Gap_location)) +
  geom_point(position=position_jitter(h=0.1, w=0.1),size=3)+
  ylab("Macroinvertebrate Family Richness")+
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

#Model diversity as mixed model like microbes
#Calculate diversity
EAB_Paired_Macs_Stand$Mac_Simp<-diversity(EAB_Paired_Macs_Stand_com, index="invsimpson")
#Change inf to 1
EAB_Paired_Macs_Stand$Mac_Simp[!is.finite(EAB_Paired_Macs_Stand$Mac_Simp)] <- 1
#build model
lmer.md= lme(Mac_Simp~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_Macs_Stand,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.md),col="darkgray")
shapiro.test(residuals(lmer.md))
#W = 0.98111, p-value = 0.5492, normal
summary(lmer.md)
#intercept and downstream significant
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

##Model populations as mixed model like microbes
#Calculate elmidae relative abundance
EAB_Paired_Macs_Stand$ElmidaeRA<-EAB_Paired_Macs_Stand$'InsectaColeopteraElmidae'/rowSums(EAB_Paired_Macs_Stand[167:213])
#build model
lmer.elm= lme(ElmidaeRA~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_Macs_Stand,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.elm),col="darkgray")
shapiro.test(residuals(lmer.elm))
#W = 0.97495, p-value = 0.3384, normal
summary(lmer.elm)
#nothing significant

#Calculate gammaridae relative abundance
EAB_Paired_Macs_Stand$GammaridRA<-EAB_Paired_Macs_Stand$'MalacostracaAmphipodaGammaridae'/rowSums(EAB_Paired_Macs_Stand[167:213])
#build model
lmer.gam= lme(GammaridRA~Watershed+Gap_location, random=~1|Stream,
              data=EAB_Paired_Macs_Stand,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.gam),col="darkgray")
shapiro.test(residuals(lmer.gam))
#W = 0.7636, p-value = 9.215e-08, not normal
EAB_Paired_Macs_Stand$log10GammaridRA<-log10(EAB_Paired_Macs_Stand$GammaridRA+1)
lmer.log10gam= lme(log10GammaridRA~Watershed+Gap_location, random=~1|Stream,
              data=EAB_Paired_Macs_Stand,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10gam),col="darkgray")
shapiro.test(residuals(lmer.log10gam))
#W = 0.81049, p-value = 1.058e-06 not normal try binomial
names(EAB_Paired_Macs_Stand)
EAB_Paired_Macs_Stand$Total_Macs<-rowSums(EAB_Paired_Macs_Stand[167:213])
glmer.gam<-glmer(GammaridRA~Watershed+Gap_location+1|Stream,family=binomial,
               weights=Total_Macs,
               data=EAB_Paired_Macs_Stand)
summary(glmer.gam)
#fail to converge, so see if Stream is adding to model
glm.gam<-glm(GammaridRA~Watershed+Gap_location,family=binomial,weights=Total_Macs,
           data=EAB_Paired_Macs_Stand)
#check assumptions
hist(residuals(glm.gam),col="darkgray")
shapiro.test(residuals(glm.gam))
#W = 0.97544, p-value = 0.3541, normal
summary(glm.gam)
#kalamazoo signifiant
EAB_Paired_Macs_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("GammaridRA", type = "mean_se")

#Calculate hydropsychidae relative abundance
EAB_Paired_Macs_Stand$HydroRA<-EAB_Paired_Macs_Stand$'InsectaTrichopteraHydropsychidae'/rowSums(EAB_Paired_Macs_Stand[167:213])
#build model
lmer.hy= lme(HydroRA~Watershed+Gap_location, random=~1|Stream,
              data=EAB_Paired_Macs_Stand,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.hy),col="darkgray")
shapiro.test(residuals(lmer.hy))
#W = 0.72775, p-value = 1.73e-08, not normal
EAB_Paired_Macs_Stand$log10HydroRA<-log10(EAB_Paired_Macs_Stand$HydroRA+1)
lmer.log10hy= lme(log10HydroRA~Watershed+Gap_location, random=~1|Stream,
                   data=EAB_Paired_Macs_Stand,
                   method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10hy),col="darkgray")
shapiro.test(residuals(lmer.log10hy))
#W = 0.75117, p-value = 5.073e-08, not significant, try binomial
glmer.hy<-glmer(HydroRA~Watershed+Gap_location+1|Stream,family=binomial,
                 weights=Total_Macs,
                 data=EAB_Paired_Macs_Stand)
summary(glmer.hy)
#fail to converge, so see if Stream is adding to model
glm.hy<-glm(HydroRA~Watershed+Gap_location,family=binomial,weights=Total_Macs,
             data=EAB_Paired_Macs_Stand)
#check assumptions
hist(residuals(glm.hy),col="darkgray")
shapiro.test(residuals(glm.hy))
#W = 0.96494, p-value = 0.1282, normal
summary(glm.hy)
#nothing signifiant

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
EAB_Paired_FFG_Stand_com<-EAB_Paired_FFG_Stand[,166:ncol(EAB_Paired_FFG_Stand)]
EABFFGComtotbray0Stand<-as.matrix(bray0(EAB_Paired_FFG_Stand_com))
EAB_Paired_FFG_Stand_env<-EAB_Paired_FFG_Stand[,1:165]
EAB_Paired_FFG_Stand_GA<-EAB_Paired_FFG_Stand_env$Gap_Area
EAB_Paired_FFG_Stand_GA_cat<-as.factor(EAB_Paired_FFG_Stand_GA)
EAB_Paired_FFG_Stand_YGF<-EAB_Paired_FFG_Stand_env$Year_Gapformation
EAB_Paired_FFG_Stand_YGF_cat<-as.factor(EAB_Paired_FFG_Stand_YGF)
EAB_Paired_FFG_Stand_GL<-EAB_Paired_FFG_Stand_env$Gap_location
EAB_Paired_FFG_Stand_W<-EAB_Paired_FFG_Stand_env$Watershed

#analyze community
adonis(as.dist(EABFFGComtotbray0Stand) ~ Watershed*Gap_location+Stream,
       data=EAB_Paired_FFG_Stand_env,permutations=9999)
#Watershed and stream significant

EAB_FFG_Com_W_indic<-signassoc(EAB_Paired_FFG_Stand_com, cluster=EAB_Paired_FFG_Stand_W,  mode=0, alternative = "two.sided",control = how(nperm=9999))
EAB_FFG_Com_W_indic_sig<-subset(EAB_FFG_Com_W_indic, psidak<=0.05)
#All 5 families indicate watershed
#Collector Filterer Grazer and Predator Kalamazoo
#Collector gatherer and shredder Clinton

#Gap Location
EAB_FFG_Com_GL_indic<-signassoc(EAB_Paired_FFG_Stand_com, cluster=EAB_Paired_FFG_Stand_GL,  mode=0, alternative = "two.sided",control = how(nperm=9999))
EAB_FFG_Com_GL_indic_sig<-subset(EAB_FFG_Com_GL_indic, psidak<=0.05)
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

#Model relative abundances of all functional feeding groups
#Combine FFGs with other metadata
EAB_Paired_FFG_Stand_RA<-merge(EAB_Ter_H2o_aqcwd_ALL,EAB_FFG_MatrixRA, by=c("Stream","Gap_location","Date"))
#collector gatherers
lmer.cg= lme(CollectorGatherer~Watershed+Gap_location, random=~1|Stream,
              data=EAB_Paired_FFG_Stand_RA,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cg),col="darkgray")
shapiro.test(residuals(lmer.cg))
#W = 0.9754, p-value = 0.3291, normal
summary(lmer.cg)
#nothing significant
EAB_Paired_FFG_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("CollectorGatherer", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Watershed) %>%
  get_summary_stats("CollectorGatherer", type = "mean_se")

#collector filterers
lmer.cf= lme(CollectorFilterer~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_FFG_Stand_RA,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.cf),col="darkgray")
shapiro.test(residuals(lmer.cf))
#W = 0.7749, p-value = 1.079e-07, not normal
EAB_Paired_FFG_Stand_RA$log10CF<-log10(EAB_Paired_FFG_Stand_RA$CollectorFilterer+1)
lmer.log10cf= lme(log10CF~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_FFG_Stand_RA,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10cf),col="darkgray")
shapiro.test(residuals(lmer.log10cf))
#W = 0.79245, p-value = 2.704e-07, not normal, try binomial
#Create weights by rowsums
EAB_Paired_FFG_Stand_RA$Total_Macs<-rowSums(EAB_Paired_FFG_Stand[166:170])
glmer.cf<-glmer(CollectorFilterer~Watershed+Gap_location+1|Stream,family=binomial,
                weights=Total_Macs,
                  data=EAB_Paired_FFG_Stand_RA)
summary(glmer.cf)
#fail to converge, so see if Stream is adding to model
glm.cf<-glm(CollectorFilterer~Watershed+Gap_location,family=binomial,weights=Total_Macs,
                data=EAB_Paired_FFG_Stand_RA)
#check assumptions
hist(residuals(glm.cf),col="darkgray")
shapiro.test(residuals(glm.cf))
#W = 0.95793, p-value = 0.05569, normal
summary(glm.cf)
#Kalamazoo significant
EAB_Paired_FFG_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("CollectorFilterer", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Watershed) %>%
  get_summary_stats("CollectorFilterer", type = "mean_se")

#shredder
lmer.sh= lme(Shredder~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_FFG_Stand_RA,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.sh),col="darkgray")
shapiro.test(residuals(lmer.sh))
#W = 0.89095, p-value = 0.0001403, not normal
EAB_Paired_FFG_Stand_RA$log10sh<-log10(EAB_Paired_FFG_Stand_RA$Shredder+1)
lmer.log10sh= lme(log10sh~Watershed+Gap_location, random=~1|Stream,
                  data=EAB_Paired_FFG_Stand_RA,
                  method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10sh),col="darkgray")
shapiro.test(residuals(lmer.log10sh))
#W = 0.93607, p-value = 0.006413, not normal, try binomial
glmer.sh<-glmer(Shredder~Watershed+Gap_location+1|Stream,family=binomial,
                weights=Total_Macs,
                data=EAB_Paired_FFG_Stand_RA)
summary(glmer.sh)
#fail to converge, so see if Stream is adding to model
glm.sh<-glm(Shredder~Watershed+Gap_location,family=binomial,weights=Total_Macs,
            data=EAB_Paired_FFG_Stand_RA)
#check assumptions
hist(residuals(glm.sh),col="darkgray")
shapiro.test(residuals(glm.sh))
#W = 0.97255, p-value = 0.249, normal
summary(glm.sh)
#Intercept and kalamazoo significant
EAB_Paired_FFG_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("Shredder", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Watershed) %>%
  get_summary_stats("Shredder", type = "mean_se")

#grazer
lmer.g= lme(Grazer~Watershed+Gap_location, random=~1|Stream,
             data=EAB_Paired_FFG_Stand_RA,
             method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.g),col="darkgray")
shapiro.test(residuals(lmer.g))
#W = 0.71551, p-value = 6.471e-09, not normal
EAB_Paired_FFG_Stand_RA$log10G<-log10(EAB_Paired_FFG_Stand_RA$Grazer+1)
lmer.log10g= lme(log10G~Watershed+Gap_location, random=~1|Stream,
            data=EAB_Paired_FFG_Stand_RA,
            method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10g),col="darkgray")
shapiro.test(residuals(lmer.log10g))
#W = 0.7705, p-value = 8.63e-08, not normal, try binomial
glmer.g<-glmer(Grazer~Watershed+Gap_location+1|Stream,family=binomial,
                weights=Total_Macs,
                data=EAB_Paired_FFG_Stand_RA)
summary(glmer.g)
#fail to converge, so see if Stream is adding to model
glm.g<-glm(Grazer~Watershed+Gap_location,family=binomial,weights=Total_Macs,
            data=EAB_Paired_FFG_Stand_RA)
#check assumptions
hist(residuals(glm.g),col="darkgray")
shapiro.test(residuals(glm.g))
#W = 0.96821, p-value = 0.1606, normal
summary(glm.g)
#Gap significant
EAB_Paired_FFG_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("Grazer", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Watershed) %>%
  get_summary_stats("Grazer", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Gap_location) %>%
  get_summary_stats("Grazer", type = "mean_se")

#predator
lmer.p= lme(Predator~Watershed+Gap_location, random=~1|Stream,
            data=EAB_Paired_FFG_Stand_RA,
            method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.p),col="darkgray")
shapiro.test(residuals(lmer.p))
#W = 0.69467, p-value = 2.638e-09, not normal
EAB_Paired_FFG_Stand_RA$log10P<-log10(EAB_Paired_FFG_Stand_RA$Predator+1)
lmer.log10p= lme(log10P~Watershed+Gap_location, random=~1|Stream,
                 data=EAB_Paired_FFG_Stand_RA,
                 method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10p),col="darkgray")
shapiro.test(residuals(lmer.log10p))
summary(lmer.log10p)
#W = 0.75305, p-value = 3.65e-08, not normal, try binomial
glmer.p<-glmer(Predator~Watershed+Gap_location+1|Stream,family=binomial,
               weights=Total_Macs,
               data=EAB_Paired_FFG_Stand_RA)
summary(glmer.p)
#fail to converge, so see if Stream is adding to model
glm.p<-glm(Predator~Watershed+Gap_location,family=binomial,weights=Total_Macs,
           data=EAB_Paired_FFG_Stand_RA)
#check assumptions
hist(residuals(glm.p),col="darkgray")
shapiro.test(residuals(glm.p))
#W = 0.94763, p-value = 0.01963, not normal
glm.log10p<-glm(log10P~Watershed+Gap_location,family=binomial,weights=Total_Macs,
                data=EAB_Paired_FFG_Stand_RA)
summary(glm.log10p)
#check assumptions
hist(residuals(glm.log10p),col="darkgray")
shapiro.test(residuals(glm.log10p))
#W = 0.9564, p-value = 0.04761, not normal
skewness(residuals(glm.log10p))
summary(glm.log10p)
#intercept signifiant
EAB_Paired_FFG_Stand %>%
  group_by(Watershed) %>%
  get_summary_stats("Predator", type = "mean_se")
EAB_Paired_FFG_Stand_RA %>%
  group_by(Watershed) %>%
  get_summary_stats("Predator", type = "mean_se")

#Model total abundance
lmer.tma= lme(Total_Macs~Watershed+Gap_location, random=~1|Stream,
            data=EAB_Paired_FFG_Stand_RA,
            method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.tma),col="darkgray")
shapiro.test(residuals(lmer.tma))
#W = 0.75244, p-value = 3.545e-08, not normal
EAB_Paired_FFG_Stand_RA$log10Total_Macs<-log10(EAB_Paired_FFG_Stand_RA$Total_Macs+1)
lmer.log10tma= lme(log10Total_Macs~Watershed+Gap_location, random=~1|Stream,
              data=EAB_Paired_FFG_Stand_RA,
              method="REML", na.action=na.omit)
#check assumptions
hist(residuals(lmer.log10tma),col="darkgray")
shapiro.test(residuals(lmer.log10tma))
#W = 0.96965, p-value = 0.1861, normal
summary(lmer.log10tma)
#intercept significant

