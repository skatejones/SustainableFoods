### Modelling effects using metafor ####
# Authors: SJ, AS

library(metafor) # for meta-analyses
library (meta) # for meta-analysis 
library(metagear) 
library(optimParallel) # for using parallel processing
library(dplyr)
library(ggplot2)
library(funModeling) 
library(tidyverse) 
library(Hmisc)
library(reshape2)
library(data.table)
library(gridExtra)
library(readxl)
library(openxlsx)
library(tidyr)
library(MuMIn)
eval(metafor:::.MuMIn)
#if (!require("devtools")) {
#  install.packages("devtools")
#}
#devtools::install_github("MathiasHarrer/dmetar")
library(dmetar) # https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/fitting-a-three-level-model.html
library(forestplot)
library(nnet)
library("DescTools")
library(forcats)
library(summarytools)
library(dunn.test)
library(car)
library("report")
library(broom) # for summary outputs
library(ggrepel)
library(cowplot) # arranging multiple plots
library(ggpubr) # ggarrange
library(ggExtra) # for marginal plots

wd <- readline() #at the prompt, copy and paste your filepath and press enter
D:\02_Bioversity\30_SustainableFoods\R
setwd(wd)

outpath <- "./Results_tradeoffs_wLER/"

d <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield")
d_yield <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_yield")
d_zeros <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_zeros")
d_zeros_yield <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_yield_zeros")
d_validity.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity")
d_yield_validity.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_yield_validity")
d_bioyield_LERplus <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_LERplus")
d_yield_LERplus <- d_bioyield_LERplus %>% select(colnames(d_yield)) %>% select(-c("Effect_ID","Validity_biodiversity_overall")) %>% unique()
#d_validity.location.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity.location")
#d_validity.N.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity.N")
#d_validity.time.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity.time")

# Add identifier and se columns
d <- d %>% mutate(Effect_ID = c(1:nrow(d))) %>% mutate(Effect_ID = as.character(Effect_ID))
d_yield_LERplus <- d_yield_LERplus %>% mutate(Effect_ID = c(1:nrow(d_yield_LERplus))) %>% mutate(Effect_ID = as.character(Effect_ID))

d <- d %>% mutate(se = sqrt(vi_B))
d_yield <- d_yield %>% mutate(se = sqrt(vi_Y))

#d <- d %>% mutate(Synergies_sig = as.character(Synergies_sig))
d <- d %>% mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ),"Win-Equiv",
                                              ifelse(Synergies_sig %in% c("Lose B-Equiv Y" , "Equiv B-Lose Y" ),"Lose-Equiv",Synergies_sig)))
d <- d %>% mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ,  "Win B-Win Y"),"Win B-Win Y",
                                              ifelse(Synergies_sig %in% c("Lose B-Lose Y" ,"Lose B-Equiv Y" , "Equiv B-Lose Y" ), "Lose B-Lose Y",Synergies_sig)))
d <- d %>%  mutate(Synergies_sig_binary = ifelse(((yi_B >0 & ci.lb_B>0) | (yi_B <0 & ci.ub_B<0)) & ((yi_Y >0 & ci.lb_Y>0)|(yi_Y<0 & ci.ub_Y<0)),1,0 ))

# Make functions used in code ####
transf = function(x){
  return((exp(x)-1)*100)
}

transf_paste <- function(model){
  paste0(round(transf(model$b),2),"% [",round(transf(model$ci.lb),2),",",round(transf(model$ci.ub),2),"]")
}

# I2 from Cheung (2014) formula 14
# typical within study variance using Q stat
sampling.variance <- function (data,vi) {  
  result<- ((nrow(data)-1) * sum(1/vi))/ (((sum(1/vi))^2)-(sum(1/(vi^2))))
  return(result)
}

I2_3level <- function(model,v.q){
  I2 <- (model$QE - model$QMdf[2])/model$QE*100 # from Abuzaid et al. 2020, I2
  I2.level1 <- v.q/(v.q + model$sigma2[1] + model$sigma2[2])*100
  I2.level2 <- model$sigma2[1]/(v.q + model$sigma2[1] + model$sigma2[2])*100
  I2.level3 <- model$sigma2[2]/(v.q + model$sigma2[1] + model$sigma2[2])*100
  return(paste0("I2 = ",round(I2,3),"% ",
                "Level 1 = ", round(I2.level1,5),"% ",
                "Level 2 = ", round(I2.level2,5),"% ",
                "Level 3 = ", round(I2.level3,5),"% "))}

# Model functions

# Using Study ID with nested Effect ID, as random effects
# Note, study ID as RE allows for variance between studies,
# while effect ID as RE allows for variance within studies.
# Level 1 RE refers to sampling variance between extracted effect sizes.
# Generally if Level 1 variance is <75%, moderator analyses should be done.

# MD = mean difference
# tau = sqrt of between-study variance
# I^2 = heterogeneity statistic

# removing the intercept when mods are categorical makes 
# estimators reflect response ratios (more useful for checking effects)
# otherwise estimates reflect difference from intercept
# see: http://www.metafor-project.org/doku.php/tips:models_with_or_without_intercept

# running rma.mv is a fixed effect model until moderators are specified (then it is a mixed effect model)
# or until random variables are specified (then it is a random effects model, unless there are moderators which makes it a mixed effects model)

# note REML should not be used if afterwards we want to compare models with anova - use ML instead

model1_function <- function(data,sigma1,sigma2,method){
  rma <- rma.mv(yi_Y, 
                vi_Y,
                ##variance between effect sizes within each study (LEVEL 2) and 
                # variance between studies (LEVEL 3) 
                random=list(~1|ID/Effect_ID), 
                intercept=TRUE,
                method=method,
                sigma2=c(sigma1,sigma2),
                test="t", # calculate CI using t-distribution
                tdist=TRUE, # slightly mimics Hartung-Knapp method to reduce bias in SE, CI, p values
                data=data)
  return(rma)}

model1_bio_function <- function(data,sigma1,sigma2,method){
  rma <- rma.mv(yi_B, 
                vi_B,
                random=list(~1|ID/Effect_ID), 
                intercept=TRUE,
                method=method,
                sigma2=c(sigma1,sigma2),
                test="t", 
                tdist=TRUE, 
                data=data)
  return(rma)}


model2_function <- function(data,mods,sigma1=NA,sigma2=NA,method="REML"){
  rma <- rma.mv(yi_Y, vi_Y, mods = mods,random=list(~1|ID/Effect_ID), 
                method=method,test="t",tdist=TRUE,
                data=data, sigma=c(sigma1,sigma2),
                verbose=FALSE,
                control=list(optimizer="optimParallel", ncpus=4)) 
  return(rma)}

model2_bio_function <- function(data,mods,sigma1=NA,sigma2=NA,method="REML"){
  rma <- rma.mv(yi_B, vi_B, mods = mods,random=list(~1|ID/Effect_ID), 
                method=method,test="t",tdist=TRUE,
                data=data, sigma=c(sigma1,sigma2),
                verbose=FALSE,
                control=list(optimizer="optimParallel", ncpus=4)) 
  return(rma)}

# Function to export model results ####
model = model2_bio_biome
model_omnibus = model2_bio_biome_omnibus
group = "Biome" 
mod = "Biome"
outcome = "Biodiversity"

model = model2_agrochem
model2_agrochem_omnibus
group="Agrochemicals"
mod="Agrochem_CT"
outcome = "Yield"
data = d

export_results <- function(data=d,vi,model,model_omnibus,group,mod="intrcpt",outcome,ci.ub.limit=250){
  model.predict <- predict.rma(model,transf=transf)   # compute prediction intervals
  data <- data.frame(data)
  
  if(mod=="intrcpt"){
    model.predict <- data.frame(model.predict) %>%
      mutate(mod = c(mod))  %>%
      mutate(group=group,outcome=outcome)}
  else{
    model.predict <- data.frame(model.predict) %>%
      cbind(mod= data[,c(mod)])  %>%
      mutate(group=group,outcome=outcome)
  }
  
  if(mod=="intrcpt"){
    results_N <- data %>%
      summarise(n_studies = n_distinct(ID),
                n_effectsizes = n_distinct(Effect_ID)) %>%
      mutate(mod=mod)
  } else{
    results_N <- data %>%
      group_by_at(mod) %>%
      summarise(n_studies = n_distinct(ID),
                n_effectsizes = n_distinct(Effect_ID))
    colnames(results_N)[1] <- "mod"
  }
  
  model.predict <- model.predict %>% rename("pred.ci.lb"="ci.lb") %>% rename("pred.ci.ub"="ci.ub")
  
  model.p <- data.frame(cbind(b=transf(model$b),ci.lb=transf(model$ci.lb),ci.ub=transf(model$ci.ub),pvalue = model$pval)) %>% 
    rename("estimate"="V1")
  model.p$term <- rownames(model.p)
  
  if(mod=="intrcpt"){
    model.predict$term <- paste0(mod)
  } else{
    model.predict$term <- paste0(mod,model.predict$mod)
  }
  model.predict <- model.predict %>% left_join(model.p,by="term") %>% select(-c(term))
  
  model.predict <- unique(model.predict) %>% 
    left_join(results_N,by="mod") %>%
    mutate(Label = paste0(mod," (",n_effectsizes,", ",n_studies,")")) %>%
    mutate(ci.ub.edit = ifelse(ci.ub>ci.ub.limit,pred,ci.ub)) %>%
    mutate(ci.ub.seg = ifelse(ci.ub>ci.ub.limit,ci.ub.limit,NA)) %>%
    mutate(ci.ub.seg = as.numeric(ci.ub.seg)) %>%
    mutate(ci.lb.seg = ifelse(ci.lb < -100 & pred<0,100-pred,
                              ifelse(ci.lb < -100 & pred>0,pred+100,NA))) %>%
    mutate(ci.lb.seg = as.numeric(ci.lb.seg)) %>%
    mutate(pi.ub.seg = ifelse(pi.ub>ci.ub.limit,ci.ub.limit,NA)) %>%
    mutate(pi.ub.seg = as.numeric(pi.ub.seg)) %>%
    mutate(pi.lb.seg = ifelse(pi.lb < -100,pred-100,NA)) %>%
    mutate(pi.lb.seg = as.numeric(pi.lb.seg))
  
  # get omnibus test results
  QM <- model_omnibus$QM
  df <- model_omnibus$QMdf 
  p <- model_omnibus$QMp
  omnibus <- paste0("F(",df[1],",",df[2], ") = ",round(QM,3),", p = ",round(p,5))[1]
  
  model.predict <- model.predict %>% 
    mutate(omnibus = c(omnibus),
           omnibus_p = c(round(p,5)))
  # Calculate I2
  v.q <- sampling.variance(data,vi)
  I2 <- I2_3level(model,v.q=v.q)
  model.predict <- model.predict %>% mutate(I2 = I2)
  print(paste0("For model including ",mod," to predict ",outcome," effect sizes"))
  print(omnibus)
  print(I2)
  
  return(model.predict)
}

# Quick check descriptive stats ####
#glimpse(d)
#freq(d[,c("Country","Crop_type_C","Crop_FAO_C","System_T","Agrochem_CT","Biome")])
#addmargins(table(d$System_T,d$Crop_FAO_C))
#unique(d$Country)

# Model 1: random effects model with no moderators #### 
model1 <- model1_function(d_yield,NA,NA,"REML") 
model1
fitstats(model1)
model1_predict <- predict.rma(model1,transf=transf)
#check <- model1_function(d,NA,NA,"REML") 
#predict.rma(check,transf=transf)
model1_bio <- model1_bio_function(d,NA,NA,"REML")
transf_paste(model1_bio)
fitstats(model1_bio)
model1_bio_predict <- predict.rma(model1_bio,transf=transf)

# Perform a likelihood-ratio-test to determine the
# significance of each level of variance
# the reported p values are for a two-sided test
# so need to be halved to confirm the variance level is significant using a one-sided test
# see Assink and Wibbelink 2016
#model1_novar2 <- model1_function(d,0,NA,"ML") # Set variance of level 2 to zero (so no variance between effect IDs is allowed)

model1.ML <- model1_function(d,NA,NA,"ML") 

model1_novar2 <- rma.mv(yi_Y, vi_Y, random=list(~1|ID/Effect_ID), 
                        intercept=TRUE,method = "ML",sigma2=c(0,NA),test="t", tdist=TRUE, data=d)
summary(model1_novar2, digits=3)
anova(model1.ML,model1_novar2) # if diff is sig, means level 2 variance is significant (because model fit is better when variance at this level is allowed to vary)

#model1_novar3 <- model1_function(d,NA,0,"ML") # Set variance of level 3 to zero
model1_novar3 <- rma.mv(yi_Y, vi_Y, random=list(~1|ID/Effect_ID), 
                        intercept=TRUE,method = "ML",sigma2=c(NA,0),test="t", tdist=TRUE, data=d)

summary(model1_novar3, digits=3)
anova(model1.ML,model1_novar3) # if difference is significant, means level 3 variance is significant

# Compute I2
mlm.variance.distribution(x = model1)

v.q <- sampling.variance(d_yield,d_yield$vi_Y)
v.q.bio <- sampling.variance(d,d$vi_B)
I2_3level(model1,v.q)
I2_3level(model1_bio,v.q.bio)

# Model 1 including intercropping treatments not using LER to measure yield ####
model1_LERplus <- model1_function(d_yield_LERplus,NA,NA,"REML")
transf_paste(model1_LERplus)

# Model 1_excZeros: RE model with no moderators ####
# run model
model1_zeros <- model1_function(d_zeros_yield,NA,NA,"REML")
model1_zeros_predict <- predict.rma(model1_zeros,transf=transf)
transf(model1_zeros$b)
summary(model1_zeros, digits=3)
summary(model1,digits=3)
plot(model1_predict$pred,model1_zeros_predict$pred)

model1_zeros_bio <- model1_bio_function(d_zeros,NA,NA,"REML")
model1_zeros_bio_predict <- predict.rma(model1_zeros_bio,transf=transf)
plot(model1_bio_predict$pred,model1_zeros_bio_predict$pred)

transf_paste(model1_zeros)
transf_paste(model1_zeros_bio)

# Model 1_validity.high: RE model with no moderators ####
# Run with only high quality effect sizes
# And for each criteria in turn

model1_validity.high <- model1_function(d_yield_validity.high,NA,NA,"REML")
transf_paste(model1_validity.high)
model1_bio_validity.high <- model1_bio_function(d_validity.high,NA,NA,"REML")
transf_paste(model1_bio_validity.high)

#model1_validity.time.high <- model1_function(d_validity.time.high,NA,NA,"REML")
#model1_validity.N.high <- model1_function(d_validity.N.high,NA,NA,"REML")
#model1_validity.location.high <- model1_function(d_validity.location.high,NA,NA,"REML")

# Model 1_noOutliers: RE model with no moderators #### 
hist(d$yi_Y)
boxplot(d$yi_Y)
hist(d$yi_B)
boxplot(d$yi_B)

#Identify possible effect size outliers
# Cook's distance
cooks.model1 <- cooks.distance(model1,parallel="multicore",reestimate=FALSE) # leverage and fit
cooks <- data.frame(cooks.model1)
cooks <- cooks %>% cbind(d_yield) 
cooks.outlier <- qchisq(0.5,df=1) # Chi sq = 0.455, df=1, alpha=0.5. 
cooks[cooks$cooks.model1 > cooks.outlier,] # no outliers for yield
cooks.outlier.visual <- 0.01

g <- ggplot(cooks,aes(x=as.numeric(Effect_ID),y=cooks.model1,label=paste(Effect_ID,sep="\n")))+
  geom_pointrange(aes(ymin=0,ymax=cooks.model1),shape=21,size=0.5)+
  #geom_pointrange(data=cooks_bio[cooks_bio$cooks.model1_bio>cooks.outlier,], 
  # aes(x=Effect_ID,ymin=0,y=cooks.model1_bio,ymax=cooks.model1_bio),colour="red")
  geom_text_repel(data=cooks[cooks$cooks.model1>0.01,],nudge_x=1,size=3)+
  geom_hline(yintercept=0.01,colour="red")+
  #geom_hline(yintercept=(mean(cooks$cooks.model1)*3),colour="blue")+
  #scale_y_continuous(limits=c(0,cooks.outlier))+
  #ggtitle("Biodiversity")+
  labs(x="Unique yield effect sizes",y="Cook's distance")+
  theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
g

png(paste0(outpath,"Cooks distance yield.png"),width=188,height=80,res=300,unit="mm")
g
dev.off()

tiff(paste0(outpath,"Cooks distance yield.tiff"),width=188,height=80,res=300,unit="mm",compression="lzw")
g
dev.off()

cooks.model1_bio <- cooks.distance(model1_bio,parallel="multicore",reestimate=FALSE) # leverage and fit
cooks_bio <- data.frame(cooks.model1_bio) %>% mutate(cooks.model1_bio = round(cooks.model1_bio,3))
cooks_bio <- cooks_bio %>% cbind(d) 
cooks.outlier <- qchisq(0.5,df=1) # Chi sq = 0.455, df=1, alpha=0.5. 
cooks_bio[cooks_bio$cooks.model1_bio > cooks.outlier,] 
cooks.outlier.visual.bio <- 0.015
#par(mfrow=c(1,1),oma = c(0, 0, 2, 0))#two plots in the same graph
#png(paste0(outpath,"Cooks distance biodiversity.png"),width=188,height=80,res=300,unit="mm")
plot(cooks.model1_bio, type="o", pch=19, xlab="Observed effect size", ylab="Cook's Distance")
abline(h = cooks.outlier, col="red")  # add cutoff line
#text(x=1:nrow(cooks.model1_bio)+1, y=x, labels=ifelse(cooks.model1_bio>cooks.outlier,names(cooks.model1),""), col="red")  
#dev.off()

g <- ggplot(cooks_bio,aes(x=Effect_ID,y=cooks.model1_bio,label=paste(Effect_ID,sep="\n")))+
  geom_pointrange(aes(ymin=0,ymax=cooks.model1_bio),shape=21,size=0.5)+
  geom_text(data=cooks_bio[cooks_bio$cooks.model1_bio>0.01,],nudge_x=20,size=3)+
  geom_hline(yintercept=0.015,colour="red")+
  #geom_hline(yintercept=(mean(cooks_bio$cooks.model1)*3),colour="blue")+
  #scale_y_continuous(limits=c(0,cooks.outlier))+
  #ggtitle("Biodiversity")+
  labs(x="Unique biodiversity effect sizes",y="Cook's distance")+
  theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid=element_blank())
g

png(paste0(outpath,"Cooks distance biodiversity.png"),width=188,height=80,res=300,unit="mm")
g
dev.off()

tiff(paste0(outpath,"Cooks distance biodiversity.tiff"),width=188,height=80,res=300,unit="mm",compression="lzw")
g
dev.off()

d_noOutlier_yield <- cooks %>% filter(cooks.model1 <= cooks.outlier.visual)
d_noOutlier_bio <- cooks_bio %>% filter(cooks.model1_bio <= cooks.outlier.visual.bio)

# Export data with outliers removed
write.xlsx(list("d_noOutlier_bio" =d_noOutlier_bio,"d_noOutlier_yield"=d_noOutlier_yield ),file=paste0(outpath,"bio_yield_outliers.xlsx"))

# Run models
model1_noOutlier_yield <- model1_function(d_noOutlier_yield,NA,NA,"REML")
model1_noOutlier_bio <- model1_bio_function(d_noOutlier_bio,NA,NA,"REML")

transf_paste(model1_noOutlier_yield)
transf_paste(model1_noOutlier_bio)

v.q.noOutlier <- sampling.variance(d_noOutlier_yield,d_noOutlier_yield$vi_Y)
I2_3level(model1_noOutlier_yield,v.q.noOutlier)

v.q.noOutlier <- sampling.variance(d_noOutlier_bio,d_noOutlier_bio$vi_Y)
I2_3level(model1_noOutlier_bio,v.q.noOutlier)

#Standardised residuals
#rs.model1 <- rstandard(model1)
#rs.model1<- as.data.frame(rs.model1)
#rs.model1$Effect_ID <- as.numeric(1:nrow(rs.model1)) #add unique effect size ID 

#rs.bio.model1 <- rstandard(model1_bio)
#rs.bio.model1<- as.data.frame(rs.bio.model1)
#rs.bio.model1$Effect_ID <- as.numeric(1:nrow(rs.bio.model1)) #add unique effect size ID 

#Hat values
#hat.model1 <- hatvalues.rma.mv(model1)
#hat.model1<- as.data.frame(hat.model1)
#hat.model1$Effect_ID <- as.numeric(1:nrow(hat.model1)) #add unique effect size ID 

#hat.outlier <- hat.model1 %>% summarise(threshold=mean(hat.model1)*2)
#hat.outlier <- 1/(nrow(hat.model1))*2
#hat.outlier

#hat.bio.model1 <- hatvalues.rma.mv(model1_bio)
#hat.bio.model1<- as.data.frame(hat.bio.model1)
#hat.bio.model1$Effect_ID <- as.numeric(1:nrow(hat.bio.model1)) #add unique effect size ID 

#hat.bio.outlier <- hat.bio.model1 %>% summarise(threshold=mean(hat.bio.model1)*2)
#hat.bio.outlier <- 1/(nrow(hat.bio.model1))*2
#hat.bio.outlier

# Using hat values and residuals
#d_noOutlier_hatres <- left_join(d,rs.model1, by="Effect_ID")%>%
#  left_join(hat.model1, by ="Effect_ID") %>%
#  rename(resid.Y= resid,se.Y=se) %>%
#  left_join(rs.bio.model1, by="Effect_ID")%>%
#  left_join(hat.bio.model1, by ="Effect_ID") %>%
#  rename(resid.bio= resid,se.bio=se) %>%  
#  filter(hat.model1 < hat.outlier) %>%
#  filter(hat.bio.model1 < hat.bio.outlier)%>% 
#  filter(resid.Y < 3) %>%
#  filter(resid.Y > -3)%>% 
#  filter(resid.bio < 3) %>%
#  filter(resid.bio > -3)

#Outliers_hatres <- left_join(d,rs.model1, by="Effect_ID")%>%
#  left_join(hat.model1, by ="Effect_ID")%>%
#  rename(resid.Y= resid,se.Y=se) %>%
#  left_join(rs.bio.model1, by="Effect_ID")%>%
#  left_join(hat.bio.model1, by ="Effect_ID")%>%
#  rename(resid.bio= resid,se.bio=se) %>% 
#  filter(hat.model1>=hat.outlier | resid.Y > 3 | resid.Y < -3) %>%
#  filter(hat.bio.model1>=hat.bio.outlier | resid.bio > 3 | resid.bio < -3)

#Plot hat values agains residual values
#par(mfrow=c(1,1),oma = c(0, 0, 2, 0))#two plots in the same graph
#tiff(paste0(outpath,"Hat values versus residuals 3 ",run,".tif"),width=10,height=4,res=120,unit="in")
#plot(x=hat.model1$hat.model1, y= rs.model1$resid, 
#     ylab= "Standardized residuals", xlab= "Hat values", 
#     main= "Database: Yield RE model without moderators", cex.main=0.6,col="red")
#points(x=d_noOutlier_hatres$hat.model1, y= d_noOutlier$resid.Y, cex.main=0.6,col="black")
#abline(v = hat.outlier, lty=2, lwd=2, col="grey50") 
#abline(h = 3, lty=2, lwd=2, col="grey50")
#abline(h = -3, lty=2, lwd=2, col="grey50")
#plot(x=hat.bio.model1$hat.bio.model1, y= rs.bio.model1$resid.bio, 
#     ylab= "Standardized residuals", xlab= "Hat values", 
#     main= "Database: Biodiversity RE model without moderators", cex.main=0.6,col="red")
#points(x=d_noOutlier_hatres$hat.bio.model1, y= d_noOutlier$resid.bio, cex.main=0.6,col="black")
#abline(v = hat.bio.outlier, lty=2, lwd=2, col="grey50") 
#abline(h = 3, lty=2, lwd=2, col="grey50")
#abline(h = -3, lty=2, lwd=2, col="grey50")
#dev.off()

# plot studentised residuals (residuals divided by their standard errors, calculating the residual for case i based on comparisons with model fit when the case i is removed)
rs.model1 <- residuals(model1,type="rstudent")
rs.model1<- as.data.frame(rs.model1)
rs.model1 <- rs.model1 %>% cbind(d_yield[,c("ID", "Effect_ID","yi_Y","vi_Y","se")])

g <- ggplot(data=rs.model1,aes(y=rs.model1,x=Effect_ID))+geom_point(alpha=0.4)+
  #geom_text(data=rs.model1, aes(y=rs.model1,x=Effect_ID,label=Effect_ID),nudge_y=0.2)+
  geom_hline(yintercept=0)+
  #geom_hline(yintercept=model1$beta,colour="grey50",linetype=2)+
  labs(y="Log RR studentized residuals","Yield effect sizes")+
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
g

# Make funnel plots ####
# check effect of removing one biodiversity effect size with extreme se
#d_check <- d %>% filter(se_B<100) # intercept only:  "9.63% [-18.4,47.29]"

g <- ggplot(data=d,aes(y=yi_B,x=se_B))+geom_point(alpha=0.4)+ #[which(d$se_B<60000),]
  #geom_text_repel(data=d[which(d$se<(1/600000)),], aes(y=yi_B,x=1/se,label=Effect_ID))+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=model1$beta,colour="grey50",linetype=2)+
  labs(y="Log RR for biodiversity")+
  theme_bw()
g
dev.off()
ggsave(paste0(outpath,"Funnel plot bio_ggplot.tiff"))
#ggsave(paste0(outpath,"Funnel plot bio_ggplot_one high se removed.tiff"))
check <- d %>% filter(Effect_ID == 70 | Effect_ID ==73)
#tiff(paste0(outpath,"Funnel plot bio_ggplot.tiff"),width=450,height=350,units="mm",res=300, compression="lzw")
#g
#dev.off
#while (!is.null(dev.list()))  dev.off()

g <- ggplot(data=d_yield,aes(y=yi_Y,x=se))+geom_point(alpha=0.4)+
  #geom_text_repel(data=d_yield[which(d_yield$se<(1/5000)),], aes(y=yi_Y,x=1/se,label=Effect_ID))+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=model1$beta,colour="grey50",linetype=2)+
  #coord_flip()+
  labs(y="Log RR for yield")+
  theme_bw()
g
ggsave(paste0(outpath,"Funnel plot yield_ggplot.tiff"))

#tiff(paste0(outpath,"Funnel plot biodiversity.tif"),width=6,height=4,units="in",res=300, compression="lzw")
metafor::funnel(x=d$yi_B,vi=d$vi_B,yaxis="sei",xlab="Log response ratio",
                steps=4,ylim=c(0,1.5),at=c(-10,-5,0,5,10),digits=1,pch=20,col="black",cex=0.8,level=c(0.95,0.99), 
                back="white",hlines="black",shade=c("grey70","grey50"))
#dev.off()

#tiff(paste0(outpath,"Funnel plot yields.tif"),width=6,height=4,units="in",res=300, compression="lzw")
metafor::funnel(x=d_yield$yi_Y,vi=d_yield$vi_Y,yaxis="sei",xlab="Log response ratio",
                steps=4,ylim=c(0,1.5),at=c(-10,-5,0,5,10),digits=1,pch=20,col="black",cex=0.8,level=c(0.95,0.99), 
                back="white",hlines="black",shade=c("grey70","grey50"))
#dev.off()


# Test for publication bias ####
# Use Egger's test by including SE as a moderator
# when the intercept is sig diff from zero, 
# it means the relationship between the mean and the sample size is asymmetrical
# and therefore considered biased. See Sterne and Egger 2005.
# Consider analyses to be biased when p<0.1 (Egger et al. 1997).
# Conduct Egger's regression test adapted by Nakagawa and Santos (2012) 
# for datasets with non-independent ESs. 
# reg test: https://wviechtb.github.io/metafor/reference/regtest.html
# https://www.bmj.com/content/315/7109/629
# funnel plots: https://rstudio-pubs-static.s3.amazonaws.com/28456_ea0b1faf0f4645cc8af81d81aaf0c1af.html 

# Identifying publication bias with a funnel plot (visual identification)
# Fairly symmetrical
# differences can also be caused by true heterogeneity, data irregularities, chance 
# (see Egger et al. 1997 and Nakagawa and Santos 2012)
# and problematic with random effects (Viectbauer and Cheng 2010)
#metafor::funnel(x=d$yi_Y,vi=d$vi_Y,yaxis="vi",atransf=exp,
#                steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=2,pch=20,col="black")
#metafor::funnel(x=d$yi_B,vi=d$vi_B,yaxis="vi",atransf=exp,
#                steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=2,pch=20,col="black")

# Identifying publication bias statistically with regression does not work when effect sizes are not independent 
# see Nakagawa and Santos 2012 equation 38 and p.1266
#regtest(d$yi_Y,d$vi_Y,model="lm",predictor="sei") 
# test for funnel plot asymmetry: t = -0.8431, df = 3379, p = 0.3992. Means not significant

model1_pub_bias <- rma.mv(yi_Y,vi_Y,mods= ~1/se, random=list(~1|ID/Effect_ID), data=d_yield, test="t",tdist=TRUE,control=list(optimizer="optimParallel", ncpus=4))
model1_pub_bias # not significant, Q(df = 224) = 203883.7200, p-val < .0001
#estimate      se    tval   df    pval    ci.lb   ci.ub 
#0.0211  0.0848  0.2491  224  0.8035  -0.1460  0.1883

model1_bio_pub_bias <- rma.mv(yi_B,vi_B,mods= ~1/se_B, random=list(~1|ID/Effect_ID), data=d, test="t",tdist=TRUE,control=list(optimizer="optimParallel", ncpus=4))
model1_bio_pub_bias # not significant
#estimate      se    tval   df    pval    ci.lb   ci.ub 
#0.0920  0.1504  0.6114  773  0.5411  -0.2033  0.3873 

model1_bio_pub_bias <- rma.mv(yi_B,vi_B,mods= ~Crop_type_C+Pest_group + 1/se_B, random=list(~1|Effect_ID,~1|ID), data=d, test="t",tdist=TRUE,control=list(optimizer="optimParallel", ncpus=4))
model1_bio_pub_bias  # not significant
# REDO !!! F(df1 = 14, df2 = 1199) = 1.4078, p-val = 0.1417

# Model 2: mixed effects model with moderators #### 
model2_treatment <- model2_function(data=d_yield,mods=~System_T-1,method="REML")
model2_treatment_omnibus <- model2_function(data=d_yield,mods=~System_T,method="REML")
fitstats(model2_treatment)

model2_treatment_action <- model2_function(data=d_yield,mods=~System_T_action-1,method="REML")
model2_treatment_action_omnibus <- model2_function(data=d_yield,mods=~System_T_action,method="REML")

model2_agrochem <- model2_function(data=d_yield,mods=~Agrochem_CT-1,method="REML")
model2_agrochem_omnibus <- model2_function(data=d_yield,mods=~Agrochem_CT,method="REML")
fitstats(model2_agrochem)

model2_agrochem_nond <- model2_function(data=d_yield[which(!(d_yield$Agrochem_CT %in% c("No data","Mixed"))),],mods=~Agrochem_CT-1,method="REML")
model2_agrochem_omnibus_nond <- model2_function(data=d_yield[which(!(d_yield$Agrochem_CT %in% c("No data","Mixed"))),],mods=~Agrochem_CT,method="REML")
fitstats(model2_agrochem_nond)

model2_crop <- model2_function(data=d_yield,mods=~Crop_type_C-1,method="REML")
model2_crop_omnibus <- model2_function(data=d_yield,mods=~Crop_type_C,method="REML")
fitstats(model2_crop)

model2_dev <- model2_function(data=d_yield,mods=~DevelopmentStatus-1,method="REML")
model2_dev_omnibus <- model2_function(data=d_yield,mods=~DevelopmentStatus,method="REML")
fitstats(model2_dev)

model2_reg <- model2_function(data=d_yield,mods=~Region.Name-1,method="REML")
model2_reg_omnibus <- model2_function(data=d_yield,mods=~Region.Name,method="REML")
fitstats(model2_dev)

model2_metric <- model2_function(data=d_yield,mods=~Yield_measure_group-1,method="REML")
model2_metric_omnibus <- model2_function(data=d_yield,mods=~Yield_measure_group,method="REML")
fitstats(model2_metric)

model2_biome <- model2_function(data=d_yield,mods=~Biome-1,method="REML")
model2_biome_omnibus <- model2_function(data=d_yield,mods=~Biome,method="REML")
fitstats(model2_biome)

model2_biome_simp <- model2_function(data=d_yield,mods=~Biome_simp-1,method="REML")
model2_biome_simp_omnibus <- model2_function(data=d_yield,mods=~Biome_simp,method="REML")

model2_crop_fao <- model2_function(data=d_yield,mods=~Crop_FAO_C-1,method="REML")
model2_crop_fao_omnibus <- model2_function(data=d_yield,mods=~Crop_FAO_C,method="REML")
fitstats(model2_crop_fao)

model2_crop_treatment <- model2_function(data=d_yield,mods=~Crop_type_C:System_T-1,method="REML")
model2_crop_treatment_omnibus <- model2_function(data=d_yield,mods=~Crop_type_C:System_T,method="REML")
model2_crop_treatment
fitstats(model2_crop_treatment)

model2_crop_fao_treatment <- model2_function(data=d_yield,mods=~Crop_FAO_C:System_T-1,method="REML")
model2_crop_fao_treatment_omnibus <- model2_function(data=d_yield,mods=~Crop_FAO_C:System_T,method="REML")
model2_crop_fao_treatment
fitstats(model2_crop_fao_treatment)

model2_crop_type_fao <- model2_function(data=d_yield,mods=~Crop_type_C:Crop_FAO_C-1,method="REML")
model2_crop_type_fao_omnibus <- model2_function(data=d_yield,mods=~Crop_type_C:Crop_FAO_C,method="REML")
model2_crop_type_fao
fitstats(model2_crop_type_fao)
#model2_reg_sub <- model2_function(data=d_yield,mods=~Sub.region.Name-1,method="REML")
#model2_reg_sub_omnibus <- model2_function(data=d_yield,mods=~Sub.region.Name,method="REML")

model2_treatment
model2_agrochem
model2_crop
model2_dev
model2_reg
model2_metric
model2_biome

model2_bio_treatment <- model2_bio_function(data=d,mods=~System_T-1,method="REML")
model2_bio_treatment_omnibus <- model2_bio_function(data=d,mods=~System_T,method="REML")

model2_bio_treatment_action <- model2_bio_function(data=d,mods=~System_T_action-1,method="REML")
model2_bio_treatment_action_omnibus <- model2_bio_function(data=d,mods=~System_T_action,method="REML")

model2_bio_agrochem <- model2_bio_function(data=d,mods=~Agrochem_CT-1,method="REML")
model2_bio_agrochem_omnibus <- model2_bio_function(data=d,mods=~Agrochem_CT,method="REML")

model2_bio_agrochem_nond <- model2_bio_function(data=d[which(!(d$Agrochem_CT %in% c("No data","Mixed"))),],mods=~Agrochem_CT-1,method="REML")
model2_bio_agrochem_omnibus_nond <- model2_bio_function(data=d[which(!(d$Agrochem_CT %in% c("No data","Mixed"))),],mods=~Agrochem_CT,method="REML")

model2_bio_crop <- model2_bio_function(data=d,mods=~Crop_type_C-1,method="REML")
model2_bio_crop_omnibus <- model2_bio_function(data=d,mods=~Crop_type_C,method="REML")

model2_bio_dev <- model2_bio_function(data=d,mods=~DevelopmentStatus-1,method="REML")
model2_bio_dev_omnibus <- model2_bio_function(data=d,mods=~DevelopmentStatus,method="REML")

model2_bio_pests <- model2_bio_function(data=d,mods=~Pest_group-1,method="REML")
model2_bio_pests_omnibus <- model2_bio_function(data=d,mods=~Pest_group,method="REML")

model2_bio_reg <- model2_bio_function(data=d,mods=~Region.Name-1,method="REML")
model2_bio_reg_omnibus <- model2_bio_function(data=d,mods=~Region.Name,method="REML")

model2_bio_taxa <- model2_bio_function(data=d,mods=~Taxa_group-1,method="REML")
model2_bio_taxa_omnibus <- model2_bio_function(data=d,mods=~Taxa_group,method="REML")

model2_bio_taxa_simp <- model2_bio_function(data=d,mods=~Taxa_group_simp -1,method="REML")
model2_bio_taxa_simp_omnibus <- model2_bio_function(data=d,mods=~Taxa_group_simp,method="REML")

model2_bio_taxa_ground <- model2_bio_function(data=d,mods=~B_ground-1,method="REML")
model2_bio_taxa_ground_omnibus <- model2_bio_function(data=d,mods=~B_ground,method="REML")

model2_bio_measure <- model2_bio_function(data=d,mods=~B_measure_group-1,method="REML")
model2_bio_measure_omnibus <- model2_bio_function(data=d,mods=~B_measure_group,method="REML")

model2_bio_biome <- model2_bio_function(data=d,mods=~Biome-1,method="REML")
model2_bio_biome_omnibus <- model2_bio_function(data=d,mods=~Biome,method="REML")

model2_bio_biome_simp <- model2_bio_function(data=d,mods=~Biome_simp-1,method="REML")
model2_bio_biome_simp_omnibus <- model2_bio_function(data=d,mods=~Biome_simp,method="REML")

model2_bio_crop_fao <- model2_bio_function(data=d,mods=~Crop_FAO_C-1,method="REML")
model2_bio_crop_fao_omnibus <- model2_bio_function(data=d,mods=~Crop_FAO_C,method="REML")

model2_bio_crop_fao_type <- model2_bio_function(data=d,mods=~Crop_type_C:Crop_FAO_C-1,method="REML")
model2_bio_crop_fao_type_omnibus <- model2_bio_function(data=d,mods=~Crop_type_C:Crop_FAO_C,method="REML")
fitstats(model2_bio_crop_fao_type)
model2_bio_crop_fao_type

model2_bio_crop_fao_treatment <- model2_bio_function(data=d,mods=~Crop_FAO_C:System_T-1,method="REML")
model2_bio_crop_fao_treatment_omnibus <- model2_bio_function(data=d,mods=~Crop_FAO_C:System_T,method="REML")
fitstats(model2_bio_crop_fao_treatment)
model2_bio_crop_fao_treatment

model2_bio_crop_type_treatment <- model2_bio_function(data=d,mods=~Crop_type_C:System_T-1,method="REML")
model2_bio_crop_type_treatment_omnibus <- model2_bio_function(data=d,mods=~Crop_type_C:System_T,method="REML")
model2_bio_crop_type_treatment
fitstats(model2_bio_crop_type_treatment)

#model2_bio_taxa_simp_biome <- model2_bio_function(data=d,mods=~Taxa_group_simp*Biome_simp-1,method="REML")
#model2_bio_taxa_simp_biome_omnibus <- model2_bio_function(data=d,mods=~Taxa_group_simp*Biome_simp,method="REML")

model2_bio_treatment
model2_bio_agrochem
model2_bio_crop 
model2_bio_dev
model2_bio_reg
model2_bio_pests 
model2_bio_taxa
model2_bio_taxa_ground
model2_bio_measure
model2_bio_biome

#### Profile plots ####

# profile.rma.mv(model1) # for post model fitting checks. Make sure surface is not flat.
# sig QE test suggests that the true effects or outcomes are heterogeneous (excluding variability accounted for by the moderators)

# identify best models ####
# http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti_and_mumin
# takes a very long time
# weights are the probability that this model is the best model
# importance is the sum of the weights for the models in which the variable appears

# Forest plots #####
d <- data.frame(d)

model1.results <- export_results(data=d_yield,vi=data$vi_Y,model= model1,model_omnibus = model1,group="Overall mean",mod="intrcpt",outcome="Yield")
model1_bio.results <- export_results(d,model1_bio,vi=data$vi_B,model_omnibus = model1_bio,"Overall mean","intrcpt","Biodiversity")
transf_paste(model1)

ggplot(model1.results,aes(y=pred))+
  geom_density(data=d_yield,aes(y=transf(yi_Y)),fill="lightgreen")+
  geom_point(aes(x=0.002,colour=mod),shape=18,size=8)+
  geom_errorbar(aes(x=0.002,ymin=ci.lb,ymax=ci.ub,colour=mod),size=1,width=0.001)+  geom_hline(yintercept=0)+
  labs(x="",y="Percentage difference (RR in %)")+coord_flip()+
  scale_y_continuous(breaks=seq(-150,250,50),limits=c(-150,250))+
  scale_colour_manual(values=terrain.colors(3))+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

g_intercept <- ggplot(model1.results,aes(y=pred))+
  geom_hline(yintercept=0)+
  #geom_histogram(data=d_yield,aes(y=yi_Y_pc,x=..density..,colour="Yield"),fill="white",bins=20,alpha=0.7,position="identity")+
  geom_density(data=d_yield,aes(y=transf(yi_Y),colour="Yield",fill="Yield"),alpha=0.4)+
  #geom_histogram(data=d,aes(y=yi_B_pc,x=..density..,colour="Biodiversity"),fill="white",bins=20,alpha=0.7,position="identity")+
  geom_density(data=d,aes(y=yi_B_pc,colour="Biodiversity",fill="Biodiversity"),alpha=0.4)+
  geom_point(aes(x=0.002,colour="Yield"),shape=18,size=4)+
  geom_errorbar(aes(x=0.002,ymin=ci.lb,ymax=ci.ub,colour="Yield"),size=0.5,width=0.001)+  
  #geom_errorbar(aes(x=0.002,ymin=pi.lb,ymax=pi.ub,colour="Yield"),size=0.5,linetype=2,width=0.001)+  
  geom_point(data=model1_bio.results, aes(x=0.004,colour="Biodiversity"),shape=18,size=4)+
  geom_errorbar(data=model1_bio.results, aes(x=0.004,ymin=ci.lb,ymax=ci.ub,colour="Biodiversity"),size=0.5,width=0.001)+  
  #geom_errorbar(data=model1_bio.results, aes(x=0.004,ymin=pi.lb,ymax=pi.ub,colour="Biodiversity"),size=1,linetype=2,width=0.001)+  
  #geom_segment(data=model1_bio.results,aes(x=0.004,xend=0.004,y=pred, yend=pi.lb.seg,colour="Biodiversity"),size=1,linetype=2,width=0.001)+
  #geom_segment(data=model1_bio.results,aes(x=0.004,xend=0.004,y=pred, yend=pi.ub.seg,colour="Biodiversity"),size=1,linetype=2,width=0.001,
  #             arrow=arrow(length=unit(5,"mm")))+
  labs(x="",y="Percentage difference (RR in %)")+coord_flip()+
  scale_y_continuous(breaks=seq(-150,250,50),limits=c(-150,250))+
  scale_fill_manual(values=c("lightgreen","grey70"),name="",labels=c("Biodiversity (n=744)","Yield (n=225)"))+
  scale_colour_manual(values=c("forestgreen","grey50"),name="",labels=c("Biodiversity (n=744)","Yield (n=225)"))+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.background=element_blank(),
        legend.position=c(0.8,y=0.8))
g_intercept

height=3
width=7.4
tiff(paste0(outpath,"Fig MA bio-yield forest ",nrow(d),"-",nrow(d_yield),".tif"),height=height,width=width,units="in",res=600,compression="lzw")
g_intercept
dev.off()

# Sensitivity forest plot ###
model1_validity.high.results <- export_results(data=d_yield_validity.high,vi=data$vi_Y,model= model1_validity.high,model_omnibus = model1_validity.high,group="Overall mean (low quality removed)",mod="intrcpt",outcome="Yield")
model1_zeros.results <- export_results(data=d_zeros,vi=data$vi_Y,model= model1_zeros,model_omnibus = model1_zeros,group="Overall mean (zeros removed)",mod="intrcpt",outcome="Yield")
model1_noOutlier_yield.results  <- export_results(data=d_noOutlier_yield,vi=data$vi_Y,model= model1_noOutlier_yield,model_omnibus = model1_noOutlier_yield,group="Overall mean (outliers removed)",mod="intrcpt",outcome="Yield")
model1_LERplus.results  <- export_results(data=d_yield_LERplus ,vi=data$vi_Y,model= model1_LERplus,model_omnibus = model1_LERplus,group="Overall mean (LER plus)",mod="intrcpt",outcome="Yield")

model1_bio.validity.high.results <- export_results(data=d_validity.high ,vi=data$vi_Y,model= model1_bio_validity.high ,model_omnibus = model1_bio_validity.high,group="Overall mean (low quality removed)",mod="intrcpt",outcome="Biodiversity")
model1_bio.zeros.results <- export_results(data=d_zeros ,vi=data$vi_Y,model= model1_zeros_bio ,model_omnibus = model1_zeros_bio,group="Overall mean (zeros removed)",mod="intrcpt",outcome="Biodiversity")
model1_bio.noOutlier_yield.results  <- export_results(data=d_noOutlier_bio,vi=data$vi_Y,model= model1_noOutlier_yield,model_omnibus = model1_noOutlier_yield,group="Overall mean (outliers removed)",mod="intrcpt",outcome="Biodiversity")

results_sensitivity <- rbind(model1.results,model1_bio.results, model1_validity.high.results,model1_zeros.results,
                             model1_noOutlier_yield.results,model1_LERplus.results, model1_bio.validity.high.results,model1_bio.zeros.results,model1_bio.noOutlier_yield.results)

write.xlsx(list(results_sensitivity),
           file=paste0(outpath,"MA results sensitivity_",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

g_sensitivity <- ggplot(results_sensitivity[results_sensitivity$outcome == "Biodiversity",],
                        aes(y=estimate,x=group))+
  geom_hline(yintercept=0)+
  geom_point(aes(x=group,colour=group),size=3)+
  geom_errorbar(aes(x=group,ymin=ci.lb,ymax=ci.ub,colour=group),size=0.5,width=0.05)+  
  geom_point(data=results_sensitivity[results_sensitivity$outcome == "Yield",],
             aes(x=group,colour=group),size=3)+
  geom_errorbar(data=results_sensitivity[results_sensitivity$outcome == "Yield",],
                aes(x=group,ymin=ci.lb,ymax=ci.ub,colour=group),size=0.5,width=0.05)+  
  
  labs(x="",y="Percentage difference (RR in %)")+
  #scale_y_continuous(breaks=seq(-150,250,50),limits=c(-150,250))+
  scale_fill_manual(values=hcl.colors(5, palette = "viridis"),name="")+
  scale_colour_manual(values=hcl.colors(5, palette = "viridis"),name="")+
  scale_x_discrete(labels=function(x)str_wrap(x,10))+
  theme_bw()+facet_wrap(~outcome,scales="free_x")+
  theme(#axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=8),
        strip.text=element_text(size=9),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.background=element_blank(),
        #legend.position=c(0.8,y=0.8),
        legend.position="none")
g_sensitivity
height=3.5
width=7.6
tiff(paste0(outpath,"Fig MA bio-yield forest sensitivity results.tif"),height=height,width=width,units="in",res=400,compression="lzw")
g_sensitivity
dev.off()

#freq(d,path=paste0(outpath))
unique(d$Biome)

model2_agrochem

model2_agrochem.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_agrochem,model2_agrochem_omnibus,group="Agrochemicals",mod="Agrochem_CT","Yield")
model2_biome.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_biome,model2_biome_omnibus,group="Biome",mod="Biome","Yield")
model2_biome_simp.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_biome_simp,model2_biome_simp_omnibus,group="Biome",mod="Biome_simp","Yield")
model2_crop.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_crop,model2_crop_omnibus,group="Crop",mod="Crop_type_C","Yield",ci.ub.limit=400)
model2_dev.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_dev,model2_dev_omnibus,group="Development",mod="DevelopmentStatus","Yield")
model2_reg.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_reg,model2_reg_omnibus,group="Region",mod="Region.Name","Yield")
model2_measure.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_metric,model2_metric_omnibus,group="Yield_metric",mod="Yield_measure_group","Yield")
model2_treatment.results <- export_results(d_yield,vi=d_yield$vi_Y,model2_treatment,model2_treatment_omnibus,group="Practice",mod="System_T","Yield")
model2_crop_fao.results  <- export_results(d_yield,vi=d_yield$vi_Y,model2_crop_fao,model2_crop_fao_omnibus,group="Crop_FAO",mod="Crop_FAO_C","Yield")
model2_crop_type_fao.results  <- export_results(d_yield,vi=d_yield$vi_Y,model2_crop_type_fao,model2_crop_type_fao_omnibus,group="Crop_FAO",mod="Crop_FAO_C","Yield")

model2_bio_agrochem.results <- export_results(d,vi=d$vi_B,model2_bio_agrochem,model2_bio_agrochem_omnibus,group="Agrochemicals",mod="Agrochem_CT","Biodiversity")
model2_bio_biome.results <- export_results(d,vi=d$vi_B,model2_bio_biome,model2_bio_biome_omnibus,group="Biome",mod="Biome","Biodiversity")
model2_bio_biome_simp.results <- export_results(d,vi=d$vi_B,model2_bio_biome_simp,model2_bio_biome_simp_omnibus,group="Biome",mod="Biome_simp","Biodiversity")
model2_bio_crop.results <- export_results(d,vi=d$vi_B,model2_bio_crop,model2_bio_crop_omnibus,group="Crop",mod="Crop_type_C","Biodiversity",ci.ub.limit=400)
model2_bio_dev.results <- export_results(d,vi=d$vi_B,model2_bio_dev,model2_bio_dev_omnibus,group="Development",mod="DevelopmentStatus","Biodiversity")
model2_bio_measure.results <- export_results(d,vi=d$vi_B,model2_bio_measure,model2_bio_measure_omnibus,group="Biodiversity_metric",mod="B_measure_group","Biodiversity")
model2_bio_pests.results <- export_results(d,vi=d$vi_B,model2_bio_pests,model2_bio_pests_omnibus,group="Pest_function",mod="Pest_group","Biodiversity")
model2_bio_reg.results <- export_results(d,vi=d$vi_B,model2_bio_reg,model2_bio_reg_omnibus,group="Region",mod="Region.Name","Biodiversity")
model2_bio_taxa.results <- export_results(d,vi=d$vi_B,model2_bio_taxa,model2_bio_taxa_omnibus,group="Taxa",mod="Taxa_group","Biodiversity",ci.ub.limit=400)
model2_bio_taxa_simp.results <- export_results(d,vi=d$vi_B,model2_bio_taxa_simp,model2_bio_taxa_simp_omnibus,group="Taxa_simp",mod="Taxa_group_simp","Biodiversity",ci.ub.limit=400)
model2_bio_taxa_ground.results <- export_results(d,vi=d$vi_B,model2_bio_taxa_ground,model2_bio_taxa_ground_omnibus,group="Ground_relation",mod="B_ground","Biodiversity")
model2_bio_treatment.results <- export_results(d,vi=d$vi_B,model2_bio_treatment,model2_bio_treatment_omnibus,group="Practice",mod="System_T","Biodiversity")
model2_bio_crop_fao.results  <- export_results(d,vi=d$vi_B,model2_bio_crop_fao,model2_bio_crop_fao_omnibus,group="Crop_FAO",mod="Crop_FAO_C","Biodiversity")

results_export_interaction <- function(data,model,model_omnibus,mod,group,outcome){
  results <- data.frame(estimate=model$b[,1],se=model$se,pvalue=model$pval, ci.lb=model$ci.lb,ci.ub=model$ci.ub) 
  results <- results %>%
    mutate(term=rownames(results))
  
  # get omnibus test results
  QM <- model_omnibus$QM
  df <- model_omnibus$QMdf 
  p <- model_omnibus$QMp
  omnibus <- paste0("F(",df[1],",",df[2], ") = ",round(QM,3),", p = ",round(p,5))[1]
  
  results <- results %>% 
    mutate(omnibus = c(omnibus),
           omnibus_p = c(round(p,5)))
  
  # add number of effect sizes and studies
  if(length(mod)==1){
    results_N <- data %>%
      group_by_at(mod) %>%
      summarise(n_studies = n_distinct(ID),
                n_effectsizes = n_distinct(Effect_ID))
    colnames(results_N)[1] <- "mod"
    } else{
    results_N <- data %>%
      group_by_at(mod) %>%
      summarise(n_studies = n_distinct(ID),
                n_effectsizes = n_distinct(Effect_ID))
    colnames(results_N)[1] <- "mod1"
    colnames(results_N)[2] <- "mod2"
    results_N <- results_N %>% mutate(mod = paste0(mod1,":",mod[2],mod2))
    }
  
  results <- results %>% 
    mutate(mod = ifelse(substr(term,1,7)=="Crop_ty", gsub("Crop_type_C*","",term),
                                               ifelse(substr(term,1,7)=="Crop_FA", gsub("Crop_FAO_C*","",term),term))) %>%
    mutate(group = group) %>%
    left_join(results_N,by="mod") %>%
    mutate(outcome=outcome)
  
  # Calculate I2
  v.q <- sampling.variance(data,data$vi)
  I2 <- I2_3level(model,v.q=v.q)
  results <- results %>% mutate(I2 = I2)
  print(omnibus)
  print(I2)
  
  return(results)
}

model2_crop_type_fao.results <- results_export_interaction(data=d_yield, model2_crop_type_fao,model2_crop_type_fao_omnibus,mod=c("Crop_type_C","Crop_FAO_C"),group="Crop type x taxon",outcome="Yield")
model2_crop_treatment.results <- results_export_interaction(data=d_yield, model2_crop_treatment,model2_crop_treatment_omnibus,mod=c("Crop_type_C","System_T"),group="Crop type x Practice",outcome="Yield")
model2_crop_fao_treatment.results <- results_export_interaction(data=d_yield, model2_crop_fao_treatment,model2_crop_fao_treatment_omnibus,mod=c("Crop_FAO_C","System_T"),group="Crop taxon x Practice",outcome="Yield")

model2_bio_crop_type_fao.results <- results_export_interaction(data=d, model2_bio_crop_fao ,model2_bio_crop_fao_omnibus,mod=c("Crop_type_C","Crop_FAO_C"),group="Crop type x taxon",outcome="Biodiversity")
model2_bio_crop_treatment.results <- results_export_interaction(data=d, model2_bio_crop_treatment,model2_bio_crop_treatment_omnibus,mod=c("Crop_type_C","System_T"),group="Crop type x Practice",outcome="Biodiversity")
model2_bio_crop_fao_treatment.results <- results_export_interaction(data=d, model2_bio_crop_fao_treatment,model2_bio_crop_fao_treatment_omnibus,mod=c("Crop_FAO_C","System_T"),group="Crop taxon x Practice",outcome="Biodiversity")

model2_crop_type_fao.results %>% select(c("group","omnibus","I2")) %>% unique()

results_interactions <- rbind(model2_crop_type_fao.results,model2_crop_treatment.results,model2_crop_fao_treatment.results,
                              model2_bio_crop_type_fao.results,model2_bio_crop_treatment.results,model2_bio_crop_fao_treatment.results)
write.xlsx(list(results_interactions),
           file=paste0(outpath,"MA results interactions_",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)


results_combined <- rbind(model1.results,model1_bio.results,
                          model2_agrochem.results, #model2_biome.results,
                          model2_biome_simp.results,
                          model2_crop.results, model2_crop_fao.results,
                          model2_dev.results, model2_reg.results,
                          model2_measure.results,model2_treatment.results,
                          model2_bio_agrochem.results,#model2_bio_biome.results,
                          model2_bio_biome_simp.results,
                          model2_bio_crop.results,model2_bio_crop_fao.results,
                          model2_bio_dev.results,model2_bio_reg.results,
                          model2_bio_measure.results, model2_bio_treatment.results,
                          model2_bio_pests.results,model2_bio_taxa_ground.results,
                          model2_bio_taxa.results,model2_bio_taxa_simp.results)

write.xlsx(list(results_combined),file=paste0(outpath,"MA results_",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

results_omnibus <- results_combined %>% select(outcome,group, omnibus,omnibus_p,I2) %>% unique() %>%
  mutate(omnibus_p = ifelse(omnibus_p < 0.0001,"<0.0001",round(omnibus_p,3)))

write.xlsx(list(results_omnibus),file=paste0(outpath,"MA omnibus tests_",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

# Model of all significant moderators ####
# keeping intercept then test for moderators is an omnibus test for significance of 
# the whole moderator (null hypothesis that the average Log RR is the same for all classes).
# removing intercept then test for moderator is NOT an omnibus test 
# but results are the same when using btt tests

# Yield
model2_all_sig <- model2_function(data=d_yield,mods=~Biome_simp+Crop_type_C-1,method="REML")
model2_all_sig_omnibus <- model2_function(data=d_yield,mods=~Biome_simp+Crop_type_C,method="REML")
model2_all_sig
model2_all_sig_omnibus
I2_3level(model2_all_sig,sampling.variance(d_yield,d_yield$vi_Y))

anova(model2_all_sig_omnibus,btt=2:6) # biome
anova(model2_all_sig_omnibus,btt=7:10) # crop type 

model2_all_sig_2 <- model2_function(data=d_yield,mods=~Biome_simp+Crop_FAO_C+System_T+Crop_FAO_C:System_T-1,method="REML")
model2_all_sig_2_omnibus <- model2_function(data=d_yield,mods=~Biome_simp+Crop_FAO_C+System_T+Crop_FAO_C:System_T,method="REML")
model2_all_sig_2_omnibus

anova(model2_all_sig_2,btt=2:6) # biome
anova(model2_all_sig_2,btt=7:13) # crop species
anova(model2_all_sig_2,btt=14:19) # practice
anova(model2_all_sig_2,btt=20:25) # interaction

I2_3level(model2_all_sig_2,sampling.variance(d_yield,d_yield$vi_Y))

# Biodiversity
model2_bio_all_sig <- model2_bio_function(data=d,mods=~Crop_type_C+Pest_group-1,method="REML")
model2_bio_all_sig_omnibus <- model2_bio_function(data=d,mods=~Crop_type_C+Pest_group ,method="REML")
model2_bio_all_sig
model2_bio_all_sig_omnibus

anova(model2_bio_all_sig,btt=1:5) # crop type 
anova(model2_bio_all_sig_omnibus,btt=2:5) # crop type 
anova(model2_bio_all_sig_omnibus,btt=6) # pests

anova(model2_bio_all_sig,btt=c(1,4)) 

model2_bio_test <- model2_bio_function(data=d,mods=~Crop_type_CT+Pest_group,method="REML")
I2_3level(model2_bio_all_sig_omnibus,v.q.bio)
I2_3level(model2_bio_test,v.q.bio)
model2_bio_all_sig_omnibus
model2_bio_test
anova(model2_bio_test,btt=2:5) # crop type 
anova(model2_bio_test,btt=6) # pests

model2_bio_test <- model2_bio_function(data=d,mods=~Pest_group+Crop_type_C:Crop_FAO_C,method="REML")
model2_bio_test
anova(model2_bio_test,btt=2) # crop pest 
anova(model2_bio_test,btt=3:11) # crop type x species

# check results ####
model1
model1.results
model1_bio
model1_bio.results

model2_crop
model2_crop_omnibus
model2_crop.results
model2_bio_crop
model2_bio_crop_omnibus
model2_bio_crop.results

model2_crop_fao
model2_crop_fao_omnibus
model2_bio_crop_fao
model2_bio_crop_fao_omnibus

model2_treatment
model2_treatment_omnibus
model2_treatment.results
model2_bio_treatment
model2_bio_treatment_omnibus
model2_bio_treatment.results

model2_biome
model2_biome.results
model2_bio_biome
transf(coef(model2_bio_biome))
model2_bio_biome.results

model2_agrochem
model2_agrochem_omnibus
model2_agrochem.results
model2_bio_agrochem
model2_bio_agrochem_omnibus

model2_bio_measure
model2_bio_measure_omnibus
transf(coef(model2_bio_measure))
anova(model2_bio_measure,btt = c(1,3))
anova(model2_bio_measure,btt = c(3,4))

model2_bio_taxa
model2_bio_taxa_omnibus
transf(coef(model2_bio_taxa))
table(d_full$Taxa_class,d_full$Taxa_group)
d %>% filter(Taxa_group =="Plants") %>% select(c(Taxa_order,Taxa_details,B_measure_group))

model2_bio_pests
model2_bio_pests_omnibus
transf(coef(model2_bio_pests))
d %>% filter(Pest_group =="Pests") %>% select(c(Taxa_order,Taxa_details,B_measure_group))

model2_bio_taxa_ground
model2_bio_taxa_ground_omnibus
transf(coef(model2_bio_taxa_ground))

model2_crop_treatment
model2_crop_treatment_omnibus
model2_bio_crop_treatment
model2_bio_crop_treatment_omnibus

model2_crop_fao_type
model2_crop_fao_type_omnibus
model2_bio_crop_fao_type
model2_bio_crop_fao_type_omnibus

model2_crop_fao_treatment_omnibus
model2_bio_crop_fao_treatment_omnibus

I2_3level(model2_treatment,v.q)
I2_3level(model2_agrochem,v.q)
I2_3level(model2_crop,v.q)
I2_3level(model2_crop_fao,v.q)
I2_3level(model2_dev,v.q)
I2_3level(model2_reg,v.q)
I2_3level(model2_metric,v.q)
I2_3level(model2_biome,v.q)
I2_3level(model2_crop_treatment,v.q)
I2_3level(model2_crop_fao_type,v.q)
I2_3level(model2_crop_fao_treatment ,v.q)
I2_3level(model2_all_sig_omnibus,v.q)

I2_3level(model2_bio_treatment,v.q.bio)
I2_3level(model2_bio_agrochem,v.q.bio)
I2_3level(model2_bio_crop,v.q.bio)
I2_3level(model2_bio_crop_fao,v.q.bio)
I2_3level(model2_bio_dev,v.q.bio)
I2_3level(model2_bio_reg,v.q.bio)
I2_3level(model2_bio_pests,v.q.bio)
I2_3level(model2_bio_taxa,v.q.bio)
I2_3level(model2_bio_taxa_ground,v.q.bio)
I2_3level(model2_bio_measure,v.q.bio)
I2_3level(model2_bio_biome,v.q.bio)
I2_3level(model2_bio_crop_treatment,v.q)
I2_3level(model2_bio_crop_fao_type,v.q)
I2_3level(model2_bio_crop_fao_treatment ,v.q)
I2_3level(model2_bio_all_sig_omnibus,v.q)

# Quick model fit check ####
fitstats(model1_bio)

## Sub-analyis per biome ####
#results_combined <- read.xlsx(paste0(outpath,"MA results_20211212.xlsx"))
#omnibus <- read.xlsx(paste0(outpath,"MA omnibus tests_20211205.xlsx"))

filter = "Mediterranean"
filter = "Temperate Forests"
filter= "Tropical & Subtropical Non-forest" 

#results_combined <- read.xlsx(paste0(outpath,"MA results_20211205.xlsx"))
model1.results <- results_combined %>% filter(group=="Overall mean" & outcome=="Yield") %>% mutate(I2=NA,Biome="All")
model1_bio.results <- results_combined %>% filter(group=="Overall mean" & outcome =="Biodiversity") %>% mutate(I2=NA,Biome="All")
model2_biome_simp.results <- results_combined %>% filter(group=="Biome_simp" & outcome=="Yield") %>% mutate(I2=NA,Biome=mod)
model2_bio_biome_simp.results <- results_combined %>% filter(group=="Biome_simp" & outcome =="Biodiversity") %>% mutate(I2=NA,Biome=mod)
#model2_treatment.results <- results_combined %>% filter(group=="Practice" & outcome=="Yield") %>% mutate(I2=NA,Biome="All")
#model2_bio_treatment.results <- results_combined %>% filter(group=="Practice" & outcome =="Biodiversity") %>% mutate(I2=NA,Biome="All")
#model2_crop_fao.results <- results_combined %>% filter(group=="Crop_FAO" & outcome =="Yield") %>% mutate(I2=NA,Biome="All")
#model2_bio_crop_fao.results <- results_combined %>% filter(group=="Crop_FAO" & outcome =="Biodiversity") %>% mutate(I2=NA,Biome="All")
#model2_crop_type.results <- results_combined %>% filter(group=="Crop" & outcome =="Yield") %>% mutate(I2=NA,Biome="All")
#model2_bio_crop_type.results <- results_combined %>% filter(group=="Crop" & outcome =="Biodiversity") %>% mutate(I2=NA,Biome="All")

# Biome x practice
results <- data.frame( "pred" =NA ,     "pred.ci.lb" =NA,   "pred.ci.ub" =NA,   "pi.lb" =NA,        "pi.ub"  =NA,
                       "cr.lb"   =NA,      "cr.ub" =NA,        "mod"  =NA,  "group"  =NA,       "outcome" =NA,      "estimate" =NA,     "ci.lb" =NA,        "ci.ub" =NA,
                       "pvalue"   =NA,     "n_studies"  =NA,   "n_effectsizes"=NA, "Label"  =NA,       "ci.ub.edit" =NA,   "ci.ub.seg"  =NA,   "ci.lb.seg" =NA,
                       "pi.ub.seg" =NA,    "pi.lb.seg"  =NA,   "omnibus" =NA,      "omnibus_p"  =NA,   "I2"=NA, "Biome" =NA)

for(filter in unique(d$Biome_simp)){
  print(paste0("For ",filter))
  d_sub <- d %>% filter(Biome_simp ==filter)
  d_yield_sub <- d_yield %>% filter(Biome_simp ==filter)

  model1_sub_bio <- model1_bio_function(data=d_sub,NA,NA,method="REML")
  model1_sub_yield <- model1_function(data=d_yield_sub,NA,NA,method="REML")
  model2_sub_bio <- model2_bio_function(data=d_sub,mods=~System_T-1,method="REML")
  model2_sub_yield <- model2_function(data=d_yield_sub,mods=~System_T-1,method="REML")
  model2_sub_bio_omnibus <- model2_bio_function(data=d_sub,mods=~System_T,method="REML")
  model2_sub_yield_omnibus <- model2_function(data=d_yield_sub,mods=~System_T,method="REML")
  
  results1_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model1_sub_bio,model_omnibus=model1_sub_bio,group="Overall",mod="intrcpt",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results1_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model1_sub_yield,model_omnibus=model1_sub_yield,group="Overall",mod="intrcpt",outcome="Yield") %>%
    mutate(Biome=filter)
  results2_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model2_sub_bio,model_omnibus=model2_sub_bio_omnibus,group="Practice",mod="System_T",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results2_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model2_sub_yield,model_omnibus=model2_sub_yield_omnibus,group="Practice",mod="System_T",outcome="Yield") %>%
    mutate(Biome=filter)
  results <- rbind(results,results1_bio,results1_yield,results2_bio,results2_yield)
  
  name = paste0("model.sub_biome_practice")
  assign(name,results,envir = .GlobalEnv)
}

model.sub_biome_practice <- model.sub_biome_practice  %>% filter(!(is.na(pred)))

temp <- model.sub_biome_practice %>%
  mutate(effect.size.ci = paste0(round(estimate,1),"% [",round(ci.lb,0),",",round(ci.ub,0),"]\n","(", n_effectsizes,",",n_studies,")"))

model.sub_biome_practice.table <- temp %>% setDT %>% dcast(group+mod~Biome+outcome,value.var="effect.size.ci") %>%
  select(-c("Other_Biodiversity","Other_Yield"),everything(),c("Other_Biodiversity","Other_Yield")) 

# Biome x Crop commodity
results <- data.frame( "pred" =NA ,     "pred.ci.lb" =NA,   "pred.ci.ub" =NA,   "pi.lb" =NA,        "pi.ub"  =NA,
                       "cr.lb"   =NA,      "cr.ub" =NA,        "mod"  =NA,  "group"  =NA,       "outcome" =NA,      "estimate" =NA,     "ci.lb" =NA,        "ci.ub" =NA,
                       "pvalue"   =NA,     "n_studies"  =NA,   "n_effectsizes"=NA, "Label"  =NA,       "ci.ub.edit" =NA,   "ci.ub.seg"  =NA,   "ci.lb.seg" =NA,
                       "pi.ub.seg" =NA,    "pi.lb.seg"  =NA,   "omnibus" =NA,      "omnibus_p"  =NA,   "I2"=NA, "Biome" =NA)

for(filter in unique(d$Biome_simp)){
  print(paste0("For ",filter))
  d_sub <- d %>% filter(Biome_simp ==filter)
  d_yield_sub <- d_yield %>% filter(Biome_simp ==filter)

  model1_sub_bio <- model1_bio_function(data=d_sub,NA,NA,method="REML")
  model1_sub_yield <- model1_function(data=d_yield_sub,NA,NA,method="REML")
  model2_sub_bio <- model2_bio_function(data=d_sub,mods=~Crop_FAO_C-1,method="REML")
  model2_sub_yield <- model2_function(data=d_yield_sub,mods=~Crop_FAO_C-1,method="REML")
  model2_sub_bio_omnibus <- model2_bio_function(data=d_sub,mods=~Crop_FAO_C,method="REML")
  model2_sub_yield_omnibus <- model2_function(data=d_yield_sub,mods=~Crop_FAO_C,method="REML")
  
  results1_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model1_sub_bio,model_omnibus=model1_sub_bio,group="Overall",mod="intrcpt",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results1_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model1_sub_yield,model_omnibus=model1_sub_yield,group="Overall",mod="intrcpt",outcome="Yield") %>%
    mutate(Biome=filter)
  results2_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model2_sub_bio,model_omnibus=model2_sub_bio_omnibus,group="Crop commodity",mod="Crop_FAO_C",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results2_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model2_sub_yield,model_omnibus=model2_sub_yield_omnibus,group="Crop commodity",mod="Crop_FAO_C",outcome="Yield") %>%
    mutate(Biome=filter)
  results <- rbind(results,results1_bio,results1_yield,results2_bio,results2_yield)
  
  name = paste0("model.sub_biome_crop_fao")
  assign(name,results,envir = .GlobalEnv)  
}

model.sub_biome_crop_fao <- model.sub_biome_crop_fao  %>% filter(!(is.na(pred)))

temp <- model.sub_biome_crop_fao %>%
  mutate(effect.size.ci = paste0(round(estimate,1),"% [",round(ci.lb,0),",",round(ci.ub,0),"]\n","(", n_effectsizes,",",n_studies,")"))

model.sub_biome_crop_fao.table <- temp %>% setDT %>% dcast(group+mod~Biome+outcome,value.var="effect.size.ci") %>%
  select(-c("Other_Biodiversity","Other_Yield"),everything(),c("Other_Biodiversity","Other_Yield")) 

# Biome x crop type (woodiness and duration)
results <- data.frame( "pred" =NA ,     "pred.ci.lb" =NA,   "pred.ci.ub" =NA,   "pi.lb" =NA,        "pi.ub"  =NA,
                       "cr.lb"   =NA,      "cr.ub" =NA,        "mod"  =NA,  "group"  =NA,       "outcome" =NA,      "estimate" =NA,     "ci.lb" =NA,        "ci.ub" =NA,
                       "pvalue"   =NA,     "n_studies"  =NA,   "n_effectsizes"=NA, "Label"  =NA,       "ci.ub.edit" =NA,   "ci.ub.seg"  =NA,   "ci.lb.seg" =NA,
                       "pi.ub.seg" =NA,    "pi.lb.seg"  =NA,   "omnibus" =NA,      "omnibus_p"  =NA,   "I2"=NA, "Biome" =NA)

filter=unique(d$Biome_simp)
filter <- filter[-2] # exclude 'Other' biome because all entries are from Annual Herb crop type so model cannot run with moderators
filter
for(filter in filter){
  print(paste0("For ",filter))
  d_sub <- d %>% filter(Biome_simp ==filter)
  d_yield_sub <- d_yield %>% filter(Biome_simp ==filter)
  
  model1_sub_bio <- model1_bio_function(data=d_sub,NA,NA,method="REML")
  model1_sub_yield <- model1_function(data=d_yield_sub,NA,NA,method="REML")
  model2_sub_bio <- model2_bio_function(data=d_sub,mods=~Crop_type_C-1,method="REML")
  model2_sub_yield <- model2_function(data=d_yield_sub,mods=~Crop_type_C-1,method="REML")
  model2_sub_bio_omnibus <- model2_bio_function(data=d_sub,mods=~Crop_type_C,method="REML")
  model2_sub_yield_omnibus <- model2_function(data=d_yield_sub,mods=~Crop_type_C,method="REML")
  
  results1_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model1_sub_bio,model_omnibus=model1_sub_bio,group="Overall",mod="intrcpt",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results1_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model1_sub_yield,model_omnibus=model1_sub_yield,group="Overall",mod="intrcpt",outcome="Yield") %>%
    mutate(Biome=filter)
  results2_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model2_sub_bio,model_omnibus=model2_sub_bio_omnibus,group="Crop",mod="Crop_type_C",outcome="Biodiversity") %>%
    mutate(Biome=filter)
  results2_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model2_sub_yield,model_omnibus=model2_sub_yield_omnibus,group="Crop",mod="Crop_type_C",outcome="Yield") %>%
    mutate(Biome=filter)
  results <- rbind(results,results1_bio,results1_yield,results2_bio,results2_yield)
  name = paste0("model.sub_biome_crop_type")
  assign(name,results,envir = .GlobalEnv)
}

filter="Other" #cannot be run as normal because all cases are for Annual Herbs, so no contrasts can be computed.
d_sub <- d %>% filter(Biome_simp ==filter)
d_yield_sub <- d_yield %>% filter(Biome_simp ==filter)

model1_sub_bio <- model1_bio_function(data=d_sub,NA,NA,method="REML")
model1_sub_yield <- model1_function(data=d_yield_sub,NA,NA,method="REML")

model1_sub_yield
transf(model1_sub_yield$b)
I2_3level(model1_sub_yield,sampling.variance(d_yield_sub,d_yield_sub$vi_Y))

model_check_main <- rma.mv(yi_Y, vi_Y, mods=~Biome_simp:Crop_type_C-1,random=list(~1|ID/Effect_ID), intercept=TRUE,method = "REML",test="t", tdist=TRUE, data=d_yield)
model_check_main
fitstats(model_check_main) #logLik = -169
transf(model_check_main$b)
I2_3level(model_check_main,sampling.variance(d_yield,d_yield$vi_Y))

ggplot(data=d_yield_sub,aes(y=transf(yi_Y),x=Effect_ID))+geom_point(colour="blue")+
  #geom_point(data=d_yield[which(d_yield$Biome_simp=="Other"),],colour="red")+
  geom_hline(yintercept=transf(model_check$b),colour="blue")+geom_hline(yintercept=transf(model_check_main$b[2]),colour="red")+
  facet_grid(~ID,scale="free") 
# shows effects are much higher for ID 638 than the other two ID's in d_yield_sub
# maybe random effects are influencing difference in overall mean for Other biome across the two models
ggplot(data=d_sub,aes(y=transf(yi_B),x=Effect_ID))+geom_point(colour="blue")+ facet_grid(~ID,scale="free") 

results1_bio <- export_results(data=d_sub,vi=d_sub$vi_B,model=model1_sub_bio,model_omnibus=model1_sub_bio,group="Overall",mod="intrcpt",outcome="Biodiversity") %>% mutate(Biome=filter)
results1_yield <- export_results(data=d_yield_sub,vi=d_yield_sub$vi_Y,model=model1_sub_yield,model_omnibus=model1_sub_yield,group="Overall",mod="intrcpt",outcome="Yield") %>%mutate(Biome=filter)
results <- rbind(results1_bio,results1_yield)

model.sub_biome_crop_type <- model.sub_biome_crop_type %>% rbind(results) %>% filter(!(is.na(pred)))  
model.sub_biome_crop_type <- model.sub_biome_crop_type %>% 
  mutate(Label = paste0(mod," (",n_effectsizes,",",n_studies,")"))
temp <- model.sub_biome_crop_type %>%
  mutate(effect.size.ci = paste0(round(estimate,1),"% [",round(ci.lb,0),",",round(ci.ub,0),"]\n","(", n_effectsizes,",",n_studies,")"))

model.sub_biome_crop_type.table <- temp %>% 
  setDT %>% dcast(group+mod~Biome+outcome,value.var="effect.size.ci") #%>%
  select(-c("Other_Biodiversity","Other_Yield"),everything(),c("Other_Biodiversity","Other_Yield"))
  
# export results
results_sub_biome <- rbind(model.sub_biome_crop_fao,
                           model.sub_biome_crop_type[which(model.sub_biome_crop_type$group != "Overall"),],
                           model.sub_biome_practice[which(model.sub_biome_practice$group != "Overall"),])

write.xlsx(x=list("biome_practice" = model.sub_biome_practice,
                  "biome_practice_table"= model.sub_biome_practice.table,
                  "biome_crop_fao"=model.sub_biome_crop_fao,
                  "biome_crop_fao_table" = model.sub_biome_crop_fao.table,
                  "biome_crop_type"=model.sub_biome_crop_type,
                  "biome_crop_type_table"=model.sub_biome_crop_type.table,
                  "results_sub_biome" =results_sub_biome),
           file=paste0(outpath,"Table coefficients per biome_MA.xlsx"),overwrite=TRUE)


# Forest plots ####

### Make forest plots for biome x other analysis (data subset) ####
results_sub_biome <- read.xlsx(paste0(outpath,"Table coefficients per biome_MA.xlsx"),sheet="results_sub_biome")

results_sub_biome <- results_sub_biome %>% unique() 
results_sub_biome <- results_sub_biome %>% 
  arrange(results_sub_biome,Biome,outcome,group,estimate) %>% 
  mutate(order = rep(1:nrow(results_sub_biome),1)) %>%
  mutate(order = ifelse(group=="Overall",-1,order)) %>%
  mutate(mod = ifelse(group =="Overall","Overall",mod))

#results_biome <- results_combined %>% filter(group=="Biome_simp") %>% mutate(I2=NA,Biome=mod,mod="OVERALL",order=rep(1:nrow(results_biome)),Label_bioyield=NA)
#results_sub_biome <- results_sub_biome %>% rbind(results_biome)

results_sub_biome <- results_sub_biome %>%mutate(Label = factor(Label,levels=unique(Label[order(order)])))

labels <- results_sub_biome %>% select(Biome,group,mod,outcome,n_studies,n_effectsizes) %>% unique() %>% setDT() %>%
  dcast.data.table(Biome+group+mod~outcome,value.var=c("n_studies","n_effectsizes")) %>%
  mutate(Label_bioyield = paste0(mod,"\n(",n_effectsizes_Biodiversity,",",n_effectsizes_Yield,",",n_studies_Yield, ")")) 

results_sub_biome <- results_sub_biome %>% left_join(labels[,c("Biome", "group","mod","Label_bioyield")])
results_sub_biome <- results_sub_biome %>% mutate(Label_bioyield = factor(Label_bioyield,levels=unique(Label_bioyield[order(order)])))

results_sub_biome <- results_sub_biome %>% mutate(Biome = fct_relevel(as.character(Biome),"Other",after=Inf))

ci.ub.limit = 400
results_sub_biome <- results_sub_biome %>%
  #mutate(estimate.edit = ifelse(estimate>ci.ub.limit,ci.ub.limit,estimate)) %>%
  mutate(ci.ub.edit = ifelse(ci.ub>ci.ub.limit,estimate,ci.ub)) %>%
  mutate(ci.ub.seg = ifelse(ci.ub>ci.ub.limit,ci.ub.limit,NA)) %>%
  mutate(ci.ub.seg = as.numeric(ci.ub.seg)) %>%
  mutate(ci.lb.seg = ifelse(ci.lb < -100 & estimate<0,100-estimate,
                            ifelse(ci.lb < -100 & estimate>0,estimate+100,NA))) %>%
  mutate(ci.lb.seg = as.numeric(ci.lb.seg)) %>%
  mutate(pi.ub.seg = ifelse(pi.ub>ci.ub.limit,ci.ub.limit,NA)) %>%
  mutate(pi.ub.seg = as.numeric(pi.ub.seg)) %>%
  mutate(pi.lb.seg = ifelse(pi.lb < -100,estimate-100,NA)) %>%
  mutate(pi.lb.seg = as.numeric(pi.lb.seg))

results_sub_biome <- results_sub_biome %>% mutate(Biome_mod = paste(Biome,Label_bioyield,sep="_")) %>%
  mutate(Biome_mod = factor(Biome_mod,levels=unique(Biome_mod[order(Biome,order)]))) 
#%>% mutate(Biome_mod = fct_relevel(as.character(Biome_mod),c("Other_Overall")))

Biome_mod_labels <- function(string) {
  str_wrap(sub(".+_", "", string),20)
}

label_wrap_gen <- function(width = 20) {
  function(variable, value) {
    lapply(strwrap(as.character(value), width=width, simplify=FALSE), 
           paste, collapse="\n")
  }
}

results_sub_biome <- results_sub_biome %>% mutate(size_group = ifelse(group=="Overall","Overall","Normal"))


g1 <- ggplot(results_sub_biome[which(results_sub_biome$group %in% c("Crop","Overall")),],
            aes(x=factor(outcome,levels=c("Yield","Biodiversity")),y=estimate))+
  #geom_histogram(data=d,aes(y=yi_Y_pc,x=..density..,colour="Yield"),fill="white",bins=20,alpha=0.7,position="identity")+
  #geom_density(data=d,aes(y=yi_Y_pc,colour=Agrochem_CT,fill=Agrochem_CT),alpha=0.4)+
  geom_point(aes(colour=Biome,shape=outcome,size=size_group))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=Biome),size=0.5,width=0.3,show_guide=TRUE)+  
    #geom_errorbar(aes(ymin=pi.lb,ymax=pi.ub,colour=Label),size=1,linetype=2,width=0.2)+  
  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")), 
                   y=estimate, yend=estimate-ci.lb.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")),
                   y=estimate, yend=ci.ub.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_hline(yintercept=0)+
  geom_text(aes(y=410,label=paste0(round(estimate,1)," [",round(ci.lb,0),",",round(ci.ub,0),"]")),hjust=0,size=3)+
  labs(x="",y="Effect size (%)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-100,470),expand = expansion(add = c(5, 5)))+
  #scale_colour_manual(values=hcl.colors(6, palette = "Earth"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_colour_manual(values=c("#A36B2B", "#BF9E66" ,"#FDE333", "#7CB0A1", "#2686A0","#4B0055"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_shape_manual(values=c("Biodiversity"=19,"Yield"=17),name="")+
  scale_size_manual(values=c("Overall"=5,"Normal"=2))+
  theme_minimal()+
  facet_grid(rows=vars(Biome_mod),scales="free",space="free",switch="y",labeller=labeller(Biome_mod=Biome_mod_labels))+
  #facet_grid(rows=vars(Biome,factor(Label_bioyield,levels=unique(Label_bioyield[order(order)]))),
  #           switch="y",   scales="free",space="free",labeller=labeller(Biome=biome.labels))+# labeller= label_wrap_gen(width=20))+ #multi_line=FALSE)
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.spacing = unit(-5,"mm"),
        text=element_text(size=10),
        axis.line.x=element_line(size=0.5,colour="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_blank(),
        #panel.border=element_rect(colour="black",fill=NA),
        panel.spacing = unit(0, "lines"),
        strip.text.y.left=element_text(angle=0,hjust=0),
        strip.switch.pad.grid=unit(0,"lines"),
        strip.background = element_blank(),
        strip.placement = "inside")+
  guides(colour=guide_legend(nrows=2,override.aes = list(linetype = 0)),
         shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0)),
         size="none")
g1


g2 <- ggplot(results_sub_biome[which(results_sub_biome$group %in% c("Overall","Practice")),],
            aes(x=factor(outcome,levels=c("Yield","Biodiversity")),y=estimate))+
  geom_point(aes(colour=Biome,shape=outcome,size=size_group))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=Biome),size=0.5,width=0.3,show_guide=TRUE)+  

  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")), 
                   y=estimate, yend=estimate-ci.lb.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")),
                   y=estimate, yend=ci.ub.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_hline(yintercept=0)+
  geom_text(aes(y=410,label=paste0(round(estimate,1)," [",round(ci.lb,0),",",round(ci.ub,0),"]")),hjust=0,size=3)+
  labs(x="",y="Effect size (%)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-100,500),expand = expansion(add = c(5, 5)))+
  scale_colour_manual(values=c("#A36B2B", "#BF9E66" ,"#FDE333", "#7CB0A1", "#2686A0","#4B0055"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_shape_manual(values=c("Biodiversity"=19,"Yield"=17),name="")+
  scale_size_manual(values=c("Overall"=5,"Normal"=2))+
  theme_minimal()+
  facet_grid(rows=vars(Biome_mod),scales="free",space="free",switch="y",labeller=labeller(Biome_mod=Biome_mod_labels))+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.spacing = unit(-5,"mm"),
        text=element_text(size=10),
        axis.line.x=element_line(size=0.5,colour="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_blank(),
        #panel.border=element_rect(colour="black",fill=NA),
        panel.spacing = unit(0, "lines"),
        strip.text.y.left=element_text(angle=0,hjust=0),
        strip.switch.pad.grid=unit(0,"lines"),
        #strip.background = element_blank(),
        strip.placement = "inside")+
  guides(colour=guide_legend(nrows=2,override.aes = list(linetype = 0)),
         shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0)),
         size="none")
g2

tiff(paste0(outpath,"Fig S10 forest plot BIOME x PRACTICE.tif"),height=9,width=7.4,units="in",res=600,compression="lzw")
g2
dev.off()

g3 <- ggplot(results_sub_biome[which(results_sub_biome$group %in% c("Crop commodity","Overall")),],
             aes(x=factor(outcome,levels=c("Yield","Biodiversity")),y=estimate))+
  geom_hline(yintercept=0)+
  geom_point(aes(colour=Biome,shape=outcome,size=size_group))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=Biome),size=0.5,width=0.3,show_guide=TRUE)+  
  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")), 
                   y=estimate, yend=estimate-ci.lb.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_segment(aes(x=factor(outcome,levels=c("Yield","Biodiversity")),xend=factor(outcome,levels=c("Yield","Biodiversity")),
                   y=estimate, yend=ci.ub.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_text(aes(y=410,label=paste0(round(estimate,1)," [",round(ci.lb,0),",",round(ci.ub,0),"]")),hjust=0,size=3)+
  
  labs(x="",y="Effect size (%)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-100,470),expand = expansion(add = c(5, 5)))+
  scale_colour_manual(values=c("#A36B2B", "#BF9E66" ,"#FDE333", "#7CB0A1", "#2686A0","#4B0055"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_shape_manual(values=c("Biodiversity"=19,"Yield"=17),name="")+
  scale_size_manual(values=c("Overall"=5,"Normal"=2))+
  theme_minimal()+
  facet_grid(rows=vars(Biome_mod),scales="free",space="free",switch="y",labeller=labeller(Biome_mod=Biome_mod_labels))+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.spacing = unit(-5,"mm"),
        text=element_text(size=10),
        axis.line.x=element_line(size=0.5,colour="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y.left=element_text(angle=0,hjust=0),
        strip.switch.pad.grid=unit(0,"lines"),
        strip.placement = "inside")+
  guides(colour=guide_legend(nrows=2,override.aes = list(linetype = 0)),
         shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0)),
         size="none")
g3

tiff(paste0(outpath,"Fig S11 forest plot BIOME x CROP.tif"),height=9,width=7.4,units="in",res=600,compression="lzw")
g3
dev.off()

results_sub_biome <- results_sub_biome %>% mutate(outcome_mod = paste0(outcome,"_",mod, " (",n_effectsizes,",",n_studies,")"))

g2 <- ggplot(results_sub_biome[which(results_sub_biome$group %in% c("Overall","Practice")),],
             aes(x=factor(outcome_mod,levels=unique(outcome_mod[order(mod,order)])),y=estimate))+
  geom_hline(yintercept=0)+
  geom_point(aes(colour=Biome,shape=outcome,size=size_group))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=Biome),size=0.5,width=0.3,show_guide=TRUE)+  
  geom_segment(aes(x=outcome_mod,xend=outcome_mod, 
                   y=estimate, yend=estimate-ci.lb.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_segment(aes(x=outcome_mod,xend=outcome_mod,
                   y=estimate, yend=ci.ub.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_text(aes(y=-400,label=str_wrap(paste0(mod," (",n_effectsizes,",",n_studies,")"),18)),hjust=0,size=2.7)+
  #geom_text(aes(y=-140,label=paste0(n_effectsizes,",",n_studies)),hjust=0,size=3)+
  #geom_text(aes(y=410,label=paste0(round(estimate,1),"\n[",round(ci.lb,0),",",round(ci.ub,0),"]")),hjust=0,size=2.7)+
  labs(x="",y="Effect size (%)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-400,420),expand = expansion(add = c(5, 5)))+
  #scale_x_discrete(labels=function(x)str_wrap(sub(".+_", "", x),20))+
  scale_colour_manual(values=c("#A36B2B", "#BF9E66" ,"#FDE333", "#7CB0A1", "#2686A0","#4B0055"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_shape_manual(values=c("Biodiversity"=19,"Yield"=17),name="")+
  scale_size_manual(values=c("Overall"=5,"Normal"=2))+
  theme_minimal()+
  facet_wrap(~Biome,scales="free",nrow=2)+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.spacing = unit(-5,"mm"),
        text=element_text(size=10),
        axis.line.x=element_line(size=0.5,colour="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_blank(),
        #panel.border=element_rect(colour="black",fill=NA),
        panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.placement = "inside")+
  guides(colour=guide_legend(nrows=2,override.aes = list(linetype = 0)),
         shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0)),
         size="none")
g2

tiff(paste0(outpath,"Fig forest plot BIOME x PRACTICE.tif"),height=8,width=7.4,units="in",res=600,compression="lzw")
g2
dev.off()

g3 <- ggplot(results_sub_biome[which(results_sub_biome$group %in% c("Overall","Crop")),],
             aes(x=factor(outcome_mod,levels=unique(outcome_mod[order(mod,order)])),y=estimate))+
  geom_hline(yintercept=0)+
  geom_point(aes(colour=Biome,shape=outcome,size=size_group))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=Biome),size=0.5,width=0.3,show_guide=TRUE)+  
  geom_segment(aes(x=outcome_mod,xend=outcome_mod, 
                   y=estimate, yend=estimate-ci.lb.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_segment(aes(x=outcome_mod,xend=outcome_mod,
                   y=estimate, yend=ci.ub.seg,colour=Biome),size=0.5,linetype=1,arrow=arrow(length=unit(2,"mm")),show_guide=TRUE)+
  geom_text(aes(y=-320,label=paste0(mod,"\n",n_effectsizes,",",n_studies)),hjust=0,size=3)+
  #geom_text(aes(y=-140,label=paste0(n_effectsizes,",",n_studies)),hjust=0,size=3)+
  geom_text(aes(y=410,label=paste0(round(estimate,1),"\n[",round(ci.lb,0),",",round(ci.ub,0),"]")),hjust=0,size=3)+
  labs(x="",y="Effect size (%)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-320,600),expand = expansion(add = c(5, 5)))+
  #scale_x_discrete(labels=function(x)str_wrap(sub(".+_", "", x),20))+
  scale_colour_manual(values=c("#A36B2B", "#BF9E66" ,"#FDE333", "#7CB0A1", "#2686A0","#4B0055"),name="",labels = function(x) str_wrap(x, width = 22))+#labels=c("Yield","Biodiversity"))+
  scale_shape_manual(values=c("Biodiversity"=19,"Yield"=17),name="")+
  scale_size_manual(values=c("Overall"=5,"Normal"=2))+
  theme_minimal()+
  facet_wrap(~Biome,scales="free")+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.spacing = unit(-5,"mm"),
        text=element_text(size=10),
        axis.line.x=element_line(size=0.5,colour="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_blank(),
        #panel.border=element_rect(colour="black",fill=NA),
        panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.placement = "inside")+
  guides(colour=guide_legend(nrows=2,override.aes = list(linetype = 0)),
         shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0)),
         size="none")
g3


legend_sub <- get_legend(
  g1 + theme(legend.position="bottom",legend.direction="horizontal",
                  #legend.box.margin = margin(0, 6, 0, 0),
                  legend.text=element_text(size=7))+
    guides(colour=guide_legend(nrow=2,override.aes = list(linetype = 0)),
           shape=guide_legend(nrow=2,order=1,override.aes = list(linetype = 0))))

plot_grid(legend_sub)


plot1 <- plot_grid(g3+theme(legend.position="none"),g2+theme(legend.position="none"),
                    align = "h", ncol=2, 
                    rel_widths = c(1/2, 1/2), rel_heights=c(1),
                    axis=c("bl"),labels=c("a","b"),label_size=10)
plot1


p <- plot_grid(plot1,legend_sub,ncol=1,rel_heights=c(0.5,0.06),axis=c("bl"))
p

tiff(paste0(outpath,"Fig forest plot BIOME x PRACTICE and CROP FAO.tif"),height=8,width=7.4,units="in",res=600,compression="lzw")
p
dev.off()


tiff(paste0(outpath,"Fig forest plot BIOME x PRACTICE.tif"),height=7.4,width=7.4,units="in",res=600,compression="lzw")
g2
dev.off()

### Make forest plots for each variable modelled in full dataset ####
results_combined <- read.xlsx(paste0(outpath,"MA results_20211212.xlsx"))
#omnibus <- read.xlsx(paste0(outpath,"MA omnibus tests_20211205.xlsx"))

results_combined <- arrange(results_combined,outcome,group,pred) %>% 
  #results_combined <- arrange(results_combined,mod) %>% 
  mutate(order = rep(1:nrow(results_combined),1)) %>%
  mutate(Label = factor(Label,levels=unique(Label[order(order)])))


labels <- results_combined %>% select(group,mod,outcome,n_studies,n_effectsizes) %>% setDT() %>%
  dcast.data.table(group+mod~outcome,value.var=c("n_studies","n_effectsizes")) %>%
  #  #mutate(n_effectsizes_Yield = ifelse(is.na(n_effectsizes_Yield),0,n_effectsizes_Yield)) %>%
  mutate(Label_bioyield = ifelse(is.na(n_effectsizes_Yield),paste0(mod,"\n(",n_effectsizes_Biodiversity,",",n_studies_Biodiversity, ")"),
                                 ifelse(is.na(n_effectsizes_Biodiversity),paste0(mod,"\n(",n_effectsizes_Yield,",",n_studies_Yield, ")"),
                                        paste0(mod,"\n(",n_effectsizes_Biodiversity,",",n_effectsizes_Yield,",",n_studies_Yield, ")"))))

results_combined <- results_combined %>% left_join(labels[,c("group","mod","Label_bioyield")])
#results_combined <- results_combined %>%   mutate(Label_bioyield_sep = paste0(n_effectsizes,",",n_studies))
#results_combined <- results_combined %>% mutate(outcome_test = factor(paste0(outcome,mod,Label_bioyield_sep))) %>%
#  mutate(outcome_test = factor(outcome_test,levels=outcome_test,labels=Label_bioyield_sep))#,labels=Label_bioyield_sep))
#results_combined <- results_combined %>% mutate(outcome=factor(outcome,levels=c("Yield","Biodiversity"))) %>%
#  mutate(Label_effectsize = paste0(round(pred,1)," [",round(ci.lb,1),",",round(ci.ub,1),"]")) 

g <- ggplot(results_combined[results_combined$group=="Biome",],aes(x=outcome,y=pred))+
  #geom_histogram(data=d,aes(y=yi_Y_pc,x=..density..,colour="Yield"),fill="white",bins=20,alpha=0.7,position="identity")+
  #geom_density(data=d,aes(y=yi_Y_pc,colour=Agrochem_CT,fill=Agrochem_CT),alpha=0.4)+
  geom_point(aes(colour=outcome),size=4)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=outcome),size=1,width=0.3)+  
  #geom_errorbar(aes(ymin=pi.lb,ymax=pi.ub,colour=Label),size=1,linetype=2,width=0.2)+  
  geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                   y=pred, yend=pred-ci.lb.seg,colour=outcome),size=1,linetype=1,arrow=arrow(length=unit(2.5,"mm")))+
  geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                   y=pred, yend=ci.ub.seg,colour=outcome),size=1,linetype=1,arrow=arrow(length=unit(2.5,"mm")))+
  geom_hline(yintercept=0)+
  #geom_text(aes(y=255,label=Label_bioyield,hjust=0),size=3,family="sans")+
  labs(x="",y="Percentage difference (Log-RR in %)")+
  coord_flip(clip="off")+
  scale_y_continuous(limits=c(-100,250),expand = expansion(add = c(5, 5)))+
  #scale_colour_manual(values=c("grey60","forestgreen"),name="",labels=c("Yield","Biodiversity"))+
  scale_colour_manual(values=c("grey60","black"),name="",labels=c("Yield","Biodiversity"))+
  #scale_colour_manual(values=c("gold","forestgreen"),name="",labels=c("Yield","Biodiversity"))+
  #scale_colour_manual( values=hcl.colors(2, palette = "cividis"),name="",labels=c("Yield","Biodiversity"))+
  theme_bw()+
  facet_grid(rows=vars(factor(str_wrap(Label_bioyield,20),levels=unique(str_wrap(Label_bioyield,20)[order(-order)]))),
             #cols=vars(Label_effectsize),switch="y",
             scales="free",space="free")+#,labeller= label_wrap_gen(multi_line=FALSE))+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y=element_text(angle=0,hjust=0),
        strip.background = element_blank(),
        strip.placement = "inside")
g

results_combined <- results_combined %>% mutate(outcome = factor(outcome))

variable="Biome"
group="Biome"
ci.ub.limit=250
wrap.limit=23
height=3
width=3.4


g_forest <- function(variable,group,ci.ub.limit=250,wrap.limit=20,height=3, width=7.4,
                     col.outcomes=c("Biodiversity"= "black","Yield"="grey60"),labels.outcomes=c("Biodiversity"="Biodiversity","Yield"="Yield")){
  data <- d
  results_N <- setDT(data) %>%
    group_by_at(variable) %>%
    summarise(n_studies = n_distinct(ID),
              n_effectsizes = n_distinct(Effect_ID))
  
  data <- data %>% left_join(results_N[,c(variable,"n_effectsizes")],by=c(variable))
  
  Synergies_count <- data %>% group_by_at(c("Synergies",variable)) %>%
    summarise(Synergies_count = n()) 
  
  data <- data %>% left_join(Synergies_count,by=c("Synergies",variable)) 
  data <- data %>%
    mutate(Synergies_prop = ifelse(Synergies=="Win B-Win Y",Synergies_count/n_effectsizes*100,
                                   ifelse(Synergies=="Lose B-Lose Y",Synergies_count/n_effectsizes*-1,0)))
  
  data <- data %>% select(c(variable, "Synergies_prop")) %>% unique() %>% 
    group_by_at(variable) %>% mutate(Synergies_prop = max(Synergies_prop)) %>% unique()
    
  data_results <- results_combined[results_combined$group==group,]
  data_results <- data_results %>% left_join(data,by=c("mod"=variable))
  
  g <- ggplot(data_results,aes(x=outcome,y=pred))+
    geom_hline(yintercept=0)+
  geom_point(aes(colour=outcome),size=2)+ #shape=15,size=2.5
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=outcome),size=0.5,width=0.1)+  
  geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                   y=pred, yend=pred-ci.lb.seg,colour=outcome),size=0.5,linetype=1,
               arrow=arrow(length=unit(2,"mm")),show.legend=FALSE)+
  geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                   y=pred, yend=ci.ub.seg,colour=outcome),size=0.5,linetype=1,
               arrow=arrow(length=unit(2,"mm")),show.legend=FALSE)+
    #geom_text(aes(y=ci.ub.limit+5,label=Label_bioyield,hjust=0),size=7*(5/14),family="sans")+
    labs(x="",y="Effect size (%)")+coord_flip()+
  scale_y_continuous(breaks=seq(-100,ci.ub.limit+100,100),limits=c(-100,ci.ub.limit),expand = expansion(add = c(5, 5)))+
  scale_colour_manual(values=col.outcomes,name="",labels=labels.outcomes)+
    #scale_colour_manual(values=c("forestgreen","grey50"),name="",labels=c("Biodiversity","Yield"))+
    theme_classic2()+
  #facet_grid(rows=vars(reorder(str_wrap(Label,wrap.limit),-order)),scales="free",space="free",margins=FALSE)+
  #facet_grid(rows=vars(factor(str_wrap(Label_bioyield,wrap.limit),levels=unique(str_wrap(Label_bioyield,wrap.limit)[order(-order)]))),
  #             scales="free",space="free")+
    facet_grid(rows=vars(factor(str_wrap(Label_bioyield,wrap.limit),levels=unique(str_wrap(Label_bioyield,wrap.limit)[order(-Synergies_prop)]))), #
               scales="free",space="free")+
    theme(legend.background=element_blank(),
          legend.position="none",
          line=element_line(size=0.3,colour="black"),
          axis.text.x=element_text(size=8,colour="black"),
          axis.title.x=element_text(size=8,colour="black"),
          text=element_text(size=8,colour="black"),
          legend.text=element_text(size=8,colour="black"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          #plot.margin=unit(c(1,1,1,1),"mm"),
          panel.spacing = unit(0, "lines"),
          strip.text.y=element_text(angle=0,hjust=0,size=8,colour="black"),
          strip.background = element_blank(),
          panel.background = element_blank(),
          strip.placement = "outside")
  print(g)
  assign(paste0("g_",group),g,envir = .GlobalEnv)
  
  tiff(paste0(outpath,"Fig forest plot each outcome by ",group,".tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g)
  dev.off()
}

results_combined <- results_combined %>% mutate(outcome = factor(outcome,levels=c("Yield","Biodiversity")))
  
g_forest(variable="Biome_simp", group="Biome",ci.ub.limit=250,wrap.limit=23,height=3, width=3.4)
g_forest(variable="Agrochem_CT",group="Agrochemicals",ci.ub.limit=250,wrap.limit=20,height=3, width=3.4)
g_forest(variable="System_T",group="Practice",ci.ub.limit=250,wrap.limit=23,height=3, width=3.4)
g_forest(variable="Crop_type_C",group="Crop",ci.ub.limit=400,wrap.limit=23,height=3, width=3.4)
g_forest(variable="Crop_FAO_C", group="Crop_FAO",ci.ub.limit=250,wrap.limit=20,height=3, width=3.4)

g_forest <- function(variable,group,ci.ub.limit=250,wrap.limit=20,height=3, width=7.4,
                     col.outcomes=c("Biodiversity"= "black","Yield"="grey60"),labels.outcomes=c("Biodiversity"="Biodiversity","Yield"="Yield")){
  data <- d
  results_N <- setDT(data) %>%
    group_by_at(variable) %>%
    summarise(n_studies = n_distinct(ID),
              n_effectsizes = n_distinct(Effect_ID))
  
  data <- data %>% left_join(results_N[,c(variable,"n_effectsizes")],by=c(variable))
  
  Synergies_count <- data %>% group_by_at(c("Synergies",variable)) %>%
    summarise(Synergies_count = n()) 
  
  data <- data %>% left_join(Synergies_count,by=c("Synergies",variable)) 
  data <- data %>%
    mutate(Synergies_prop = ifelse(Synergies=="Win B-Win Y",Synergies_count/n_effectsizes*100,
                                   ifelse(Synergies=="Lose B-Lose Y",Synergies_count/n_effectsizes*-1,0)))
  
  data <- data %>% select(c(variable, "Synergies_prop")) %>% unique() %>% 
    group_by_at(variable) %>% mutate(Synergies_prop = max(Synergies_prop)) %>% unique()
  
  data_results <- results_combined[results_combined$group==group,]
  data_results <- data_results %>% left_join(data,by=c("mod"=variable))
  
  g <- ggplot(data_results,aes(x=outcome,y=pred))+
    geom_hline(yintercept=0)+
    geom_point(aes(colour=outcome),size=2)+ #shape=15,size=2.5
    geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub.edit,colour=outcome),size=0.5,width=0.1)+  
    geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                     y=pred, yend=pred-ci.lb.seg,colour=outcome),size=0.5,linetype=1,
                 arrow=arrow(length=unit(2,"mm")),show.legend=FALSE)+
    geom_segment(aes(x=reorder(str_wrap(outcome,25),pred),xend=reorder(str_wrap(outcome,25),pred),
                     y=pred, yend=ci.ub.seg,colour=outcome),size=0.5,linetype=1,
                 arrow=arrow(length=unit(2,"mm")),show.legend=FALSE)+
    #geom_text(aes(y=ci.ub.limit+5,label=Label_bioyield,hjust=0),size=7*(5/14),family="sans")+
    labs(x="",y="Effect size (%)")+coord_flip()+
    scale_y_continuous(breaks=seq(-100,ci.ub.limit+100,100),limits=c(-100,ci.ub.limit),expand = expansion(add = c(5, 5)))+
    scale_colour_manual(values=col.outcomes,name="",labels=labels.outcomes)+
    theme_classic2()+
       facet_grid(rows=vars(factor(str_wrap(Label_bioyield,wrap.limit),levels=unique(str_wrap(Label_bioyield,wrap.limit)[order(-pred)]))), #
               scales="free",space="free",switch="y")+
    theme(legend.background=element_blank(),
          legend.position="none",
          text=element_text(size=8,colour="black"),
          line=element_line(size=0.3,colour="black"),
          axis.text.x=element_text(size=8,colour="black"),
          axis.title.x=element_text(size=8,colour="black"),
          legend.text=element_text(size=8,colour="black"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          #plot.margin=unit(c(1,1,1,1),"mm"),
          panel.spacing = unit(0, "lines"),
          strip.text.y.left=element_text(angle=0,hjust=1,size=8),
          strip.switch.pad.grid=unit(0,"lines"),
          #strip.text.y=element_text(angle=90,hjust=0,size=8,colour="black"),
          strip.background = element_blank(),
          panel.background = element_blank(),
          strip.placement = "outside")
  print(g)
  assign(paste0("g_",group),g,envir = .GlobalEnv)
  
  tiff(paste0(outpath,"Fig forest plot each outcome by ",group,".tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g)
  dev.off()
}

g_forest(variable="Taxa_group", group="Taxa",ci.ub.limit=400,wrap.limit=12,height=3, width=3.4,col.outcomes=c("black"),labels.outcomes=c("Biodiversity"))
g_forest(variable="Taxa_group_simp", group="Taxa_simp",ci.ub.limit=400,wrap.limit=15,height=3, width=3.4,col.outcomes=c("black"),labels.outcomes=c("Biodiversity"))
g_forest(variable="DevelopmentStatus", group="Development",ci.ub.limit=100,wrap.limit=15,height=1.5, width=3.4)
g_forest(variable="Region.Name",group="Region",ci.ub.limit=200,wrap.limit=15,height=3, width=3.4)
g_forest(variable="Yield_measure_group", group="Yield_metric",ci.ub.limit=200,wrap.limit=12,height=3, width=3.4,col.outcomes=c("grey50"),labels.outcomes=c("Yield"))
g_forest(variable="B_measure_group", group="Biodiversity_metric",ci.ub.limit=200,wrap.limit=10,height=3, width=3.4,col.outcomes=c("black"),labels.outcomes=c("Biodiversity"))
g_forest(variable="Pest_group", group="Pest_function",ci.ub.limit=100,wrap.limit=10,height=1.5, width=3.4,col.outcomes=c("black"),labels.outcomes=c("Biodiversity"))
g_forest(variable="B_ground",group="Ground_relation",ci.ub.limit=250,wrap.limit=10,height=1.5, width=3.4,col.outcomes=c("black"),labels.outcomes=c("Biodiversity"))

#write.xlsx(list("model1_bio.results" = model1_bio.results,"model1.results"=model1.results),
#           file=paste0(outpath,"MA results model1.xlsx"),overwrite=TRUE)

legend_meta <- get_legend(
  g_Biome + theme(legend.position="bottom",legend.direction="horizontal",
                   #legend.box.margin = margin(0, 6, 0, 0),
                   legend.text=element_text(size=9))+
    guides(fill=guide_legend(nrow=1,reverse=FALSE)))


### FIGURE: Probability plots, and proportions and forest plots side by side  ####

# Figure 2 for paper: combined plot with MA results ###
legend_multi <- get_legend(
  g_system + theme(legend.position="bottom",legend.direction="horizontal",
                   #legend.box.margin = margin(0, 6, 0, 0),
                   legend.text=element_text(size=8))+
    guides(fill=guide_legend(nrow=1)))

legend_meta <- get_legend(
  g_Biome + theme(legend.position="bottom",legend.direction="horizontal",
                  #legend.box.margin = margin(0, 6, 0, 0),
                  legend.text=element_text(size=8))+
    guides(fill=guide_legend(nrow=1)))

plot_legend <- plot_grid(NULL,legend_multi,NULL,legend_meta,
                         align = "vh", ncol=4, 
                         rel_widths = c(1/8,3/8, 1/8,3/8), rel_heights=c(1),
                         axis=c("b"),labels=c(""),label_size=10)
plot_legend

plot1 <- plot_grid(g2_biome +theme(legend.position="none"), g_Biome+theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("a","b"),label_size=10)
plot1

plot2 <- plot_grid(g2_agrochem +theme(legend.position="none"), g_Agrochemicals +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("c","d"),label_size=10)
plot2

plot2 <- plot_grid(g2_system +theme(legend.position="none"), g_Practice+theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("c","d"),label_size=10)
plot2

plot3 <- plot_grid(g2_crop +theme(legend.position="none"), g_Crop +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("e","f"),label_size=10)
plot3

p <- plot_grid(plot1,plot2,plot3,plot_legend,ncol=1,rel_heights=c(0.32,0.36,0.28,0.04),axis=c("lb"),align = "v")
p

tiff(paste0(outpath,"Fig 2 Proportions stacked_practice_crop_flip.tif"),height=18.5,width=18,units="cm",res=600,compression="lzw")
p
dev.off()

pdf(paste0(outpath,"Fig 2 Proportions stacked_practice_crop_flip.pdf"),height=18.5,width=18)
p
dev.off()

# Figure 3 for paper: biodiversity forest plots ###
library("jpeg")
library("patchwork")
icon_pests <- readJPEG("D:\\02_Bioversity\\30_SustainableFoods\\R\\pest_icon.jpg",native=TRUE)

plot1 <- plot_grid(g_Pest_function +theme(legend.position="none"), 
                   g_Ground_relation +theme(legend.position="none"),
                   g_Biodiversity_metric,
                   align = "vh", ncol=1, 
                   rel_widths = c(1), rel_heights=c(1/4,1/4,1/2),
                   axis=c("bt"),labels=c("a","b","c"),label_size=10)
plot1


plot2 <- plot_grid(plot1,g_Taxa,
                   align = "h", ncol=2, 
                   rel_widths = c(2/5,3/5), rel_heights=c(1),
                   axis=c("tlr"),labels=c("","d"),label_size=10)
plot2


tiff(paste0(outpath,"Fig 3 Biodiversity forest plots.tif"),height=115,width=115,units="mm",res=600,compression="lzw")
plot2
dev.off()

# Figure S4 ####
plot1 <- plot_grid(g2_crop_FAO +theme(legend.position="none"), g_Crop_FAO +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("a","b"),label_size=10)
plot1

plot2 <- plot_grid(g2_agrochem +theme(legend.position="none"), g_Agrochemicals+theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("c","d"),label_size=10)
plot2


p <- plot_grid(plot1,plot2,plot_legend,ncol=1,rel_heights=c(0.5,0.3,0.04),axis=c("lbtr"),align="v")
p

tiff(paste0(outpath,"Fig S4 Proportions stacked_flip.tif"),height=18,width=18,units="cm",res=600,compression="lzw")
p
dev.off()

# Fig S5 ####
plot1 <- plot_grid(g_Development,g_Yield_metric, 
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2,1/2), rel_heights=c(1),
                   axis=c("bl"),labels=c("a","b"),label_size=10)
plot1

p <- plot_grid(plot1,legend_meta,ncol=1,rel_heights=c(0.3,0.04),axis=c("lbtr"),align="v")
p

tiff(paste0(outpath,"Fig S5 MA forest plots.tif"),height=7,width=18,units="cm",res=600,compression="lzw")
p
dev.off()

# Figures - not used ####
plot1 <- plot_grid(g2_taxa +theme(legend.position="none"), g_Taxa +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("a","b"),label_size=10)
plot1

plot2 <- plot_grid(g2_pest +theme(legend.position="none"), g_Pest_function +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("c","d"),label_size=10)
plot2

plot3 <- plot_grid(g2_ground +theme(legend.position="none"), g_Ground_relation +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("e","f"),label_size=10)
plot3

plot23 <- plot_grid(plot2,plot3,ncol=2,rel_widths = c(0.5,0.5),axis=c("lb"),align="h")
plot23

plot4 <- plot_grid(g2_b_metric +theme(legend.position="none"), g_Biodiversity_metric +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("g","h"),label_size=10)
plot4

plot5 <- plot_grid(g2_y_metric +theme(legend.position="none"), g_Yield_metric +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("i","j"),label_size=10)
plot5

plot45 <- plot_grid(plot4,plot5,ncol=2,rel_widths = c(0.5,0.5),axis=c("lb"),align="h")
plot45


p <- plot_grid(plot1,plot23,plot45,plot_legend,ncol=1,rel_heights=c(0.4,0.15,0.25,0.04),axis=c("lb"),align = "v")
p



plot2 <- plot_grid(g2_crop_FAO +theme(legend.position="none"), g_Crop_FAO +theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("bt"),labels=c("c","d"),label_size=10)
plot2

plot1 <- plot_grid(g2_system_red +theme(legend.position="none"), g_Practice+theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(1/2, 1/2), rel_heights=c(1),
                   axis=c("b"),labels=c("a","b"),label_size=10)
plot1

#tiff(paste0(outpath,"Fig 1 Proportions stacked_practice_red_flip.tif"),height=4,width=7.4,units="in",res=600,compression="lzw")
#plot1
#dev.off()

plot1 <- plot_grid(g_system +theme(legend.position="none"), g_Practice+theme(legend.position="none"),
                   align = "vh", ncol=2, 
                   rel_widths = c(2/3, 1/3), rel_heights=c(1),
                   axis=c("b"),labels=c("a","b"),label_size=10)
plot1

#tiff(paste0(outpath,"Fig 1 Proportions stacked_practice.tif"),height=4,width=7.4,units="in",res=600,compression="lzw")
#plot1
#dev.off()

# Mega-regression ####

# Find best model
# https://www.statology.org/interpret-log-likelihood/
# https://www.statology.org/likelihood-ratio-test-in-r/

# Set reference levels to most frequently occurring 
d <- d %>% mutate(Biome = factor(Biome))
table(d$Biome)
d$Biome <- fct_relevel(d$Biome,"Tropical & Subtropical Non-forests")
d_yield$Biome <- fct_relevel(d_yield$Biome,"Tropical & Subtropical Non-forests")

table(d$Crop_FAO_C)
d$Crop_FAO_C <- fct_relevel(d$Crop_FAO_C,"Cereals")
d_yield$Crop_FAO_C <- fct_relevel(d_yield$Crop_FAO_C,"Cereals")

table(d$Crop_type_C)
d$Crop_type_C <- fct_relevel(d$Crop_type_C,"Annual Herb")
d_yield$Crop_type_C <- fct_relevel(d_yield$Crop_type_C,"Annual Herb")

table(d$System_T)
d$System_T <- fct_relevel(d$System_T,"Associated plant species")
d_yield$System_T <- fct_relevel(d_yield$System_T,"Associated plant species")

table(d$Agrochem_CT)
d$Agrochem_CT <- fct_relevel(d$Agrochem_CT,"Agrochemicals")
d_yield$Agrochem_CT <- fct_relevel(d_yield$Agrochem_CT,"Agrochemicals")

table(d$Taxa_group)
d$Taxa_group <- fct_relevel(d$Taxa_group,"Insect (other)")

table(d$Taxa_group_simp)
d$Taxa_group_simp <- fct_relevel(d$Taxa_group_simp,"Invertebrates")

table(d$B_measure_group)
d$B_measure_group <- fct_relevel(d$B_measure_group,"Abundance")

table(d$Pest_group)
d$Pest_group <- fct_relevel(d$Pest_group,"Other")

table(d$B_ground)
d$B_ground <- fct_relevel(d$B_ground,"Above")

table(d$Yield_measure_group)
d$Yield_measure_group <- fct_relevel(d$Yield_measure_group,"Mass per area")
d_yield$Yield_measure_group <- fct_relevel(d_yield$Yield_measure_group,"Mass per area")

# Yield #
# Using crop type not crop FAO
model_all <-  rma.mv(yi_Y, vi_Y, 
                     mods =~Crop_type_C+System_T + Biome_simp + Agrochem_CT + Yield_measure_group +
                       System_T:Crop_type_C,
                     random=list(~1|ID/Effect_ID), 
                     method="ML",test="t",tdist=TRUE,
                     data=d_yield, 
                     verbose=FALSE,control=list(optimizer="optimParallel", ncpus=4)) 
model_all
anova(model_all,btt=2:5) # crop type
anova(model_all,btt=6:10) # crop species
anova(model_all,btt=11:16) # practice
anova(model_all,btt=17:23) # biome
anova(model_all,btt=24:26) # agrochem
anova(model_all,btt=27:28) # yield metric
anova(model_all,btt=29:31) # crop x practice
anova(model_all,btt=32:34) # crop taxon x practice

res <- dredge(model_all, trace=2,rank="AICc",extra = alist(BIC,AIC))
plot(res)
head(res,10)
importance(res)
out1a <- data.frame(res)
out1b <- data.frame(importance(res))
out1b$covariate <- row.names(out1b)

# model average
model_all_average <- model.avg(res, revised.var=FALSE, subset= (delta <2)) # only one model
summary(model_all_average)
# use full estimates: https://uoftcoders.github.io/rcourse/lec09-model-selection.html 
res_avg_fit <- model_all_average$msTable[1,]
res_avg=summary(model_all_average)
res_avg=data.frame(cbind(summary(model_all_average)$coefmat.full,summary(model_all_average)$coefmat.nmod))
res_avg <- res_avg %>%  mutate(variable = rownames(res_avg)) %>%
  mutate(LogRR_pc = transf(Estimate),
         ci.lb = transf(Estimate-(1.96*Std..Error)),
         ci.ub = transf(Estimate+(1.96*Std..Error))) %>%
  rename("estimate"="Estimate", "se" = "Std..Error", "p.value"=   "Pr...z..")

# Using crop FAO not crop type
model_all_2 <-  rma.mv(yi_Y, vi_Y, 
                     mods =~Crop_FAO_C+System_T + Biome_simp + Agrochem_CT + Yield_measure_group +
                        System_T:Crop_FAO_C,
                     random=list(~1|ID/Effect_ID), 
                     method="ML",test="t",tdist=TRUE,
                     data=d_yield, 
                     verbose=FALSE,control=list(optimizer="optimParallel", ncpus=4)) 
model_all_2

res_2 <- dredge(model_all_2, trace=2,rank="AICc",extra = alist(BIC,AIC))
plot(res_2)
head(res_2,10)
head(res,10)
importance(res_2)
out1c <- data.frame(res_2)
out1d <- data.frame(importance(res_2))
out1d$covariate <- row.names(out1d)

# model average
model_all_2_average <- model.avg(res_2, revised.var=FALSE, subset= (delta <2)) #  'object' consists of only one model

# Biodiversity #
# Using crop type not crop FAO
model_bio_all <-  rma.mv(yi_B, vi_B, 
                         mods = ~Crop_type_C+System_T + Biome_simp + Agrochem_CT  +
                           B_measure_group + Taxa_group + Pest_group + B_ground+
                           System_T:Crop_type_C,
                         random=list(~1|ID/Effect_ID), 
                         method="ML",test="t",tdist=TRUE,
                         data=d, 
                         verbose=FALSE,control=list(optimizer="optimParallel", ncpus=4)) 
model_bio_all
res_bio <- dredge(model_bio_all, trace=2,rank="AICc",extra = alist(BIC,AIC)) # !!! This takes a long time to run
plot(res_bio)
head(res_bio,10)
importance(res_bio)
out2a_bio <- data.frame(res_bio)
out2b_bio <- data.frame(importance(res_bio))
out2b_bio$covariate <- row.names(out2b_bio)
model_bio_all_average <- model.avg(res_bio, revised.var=FALSE, subset= (delta <2))
summary(model_bio_all_average)
transf(summary(model_bio_all_average)$coefmat.full[,1:2])
res_bio_avg_fit <- model_bio_all_average$msTable[1,]
res_bio_avg_fit
res_bio_avg=summary(model_bio_all_average)
res_bio_avg=data.frame(cbind(summary(model_bio_all_average)$coefmat.full,summary(model_bio_all_average)$coefmat.nmod))
res_bio_avg <- res_bio_avg %>%  mutate(variable = rownames(res_bio_avg)) %>%
  mutate(LogRR_pc = transf(Estimate),
         ci.lb = transf(Estimate-(1.96*Std..Error)),
         ci.ub = transf(Estimate+(1.96*Std..Error))) %>%
  rename("estimate"="Estimate", "se" = "Std..Error", "p.value"=   "Pr...z..")

# Using crop FAO not crop type
model_bio_all_2 <-  rma.mv(yi_B, vi_B, 
                         mods = ~Crop_FAO_C+System_T + Biome_simp + Agrochem_CT  +
                           B_measure_group + Taxa_group + Pest_group + B_ground+
                           System_T:Crop_FAO_C,
                         random=list(~1|ID/Effect_ID), 
                         method="ML",test="t",tdist=TRUE,
                         data=d, 
                         verbose=FALSE,control=list(optimizer="optimParallel", ncpus=4)) 
model_bio_all_2

res_bio_2 <- dredge(model_bio_all_2, trace=2,rank="AICc",extra = alist(BIC,AIC)) # !!! This takes a long time to run
plot(res_bio_2)
head(res_bio_2,10)
importance(res_bio_2)
out2c_bio <- data.frame(res_bio_2)
out2d_bio <- data.frame(importance(res_bio_2))
out2d_bio$covariate <- row.names(out2d_bio)
model_bio_all_2_average <- model.avg(res_bio_2, revised.var=FALSE, subset= (delta <2))
summary(model_bio_all_2_average)
transf(summary(model_bio_all_2_average)$coefmat.full[,1:2])
res_bio_2_avg_fit <- model_bio_all_2_average$msTable[1,]
res_bio_2_avg_fit
res_bio_2_avg=summary(model_bio_all_2_average)
res_bio_2_avg=data.frame(cbind(summary(model_bio_all_2_average)$coefmat.full,summary(model_bio_all_2_average)$coefmat.nmod))
res_bio_2_avg <- res_bio_2_avg %>%  mutate(variable = rownames(res_bio_2_avg)) %>%
  mutate(LogRR_pc = transf(Estimate),
         ci.lb = transf(Estimate-(1.96*Std..Error)),
         ci.ub = transf(Estimate+(1.96*Std..Error))) %>%
  rename("estimate"="Estimate", "se" = "Std..Error", "p.value"=   "Pr...z..")


# export results
write.xlsx(list(res_yield_1=out1a, importance_yield_1=out1b,
                res_yield_2=out1c,importance_yield_2=out1d,
                res_bio_1=out2a_bio,importance_bio_1=out2b_bio,
                res_bio_2=out2c_bio,importance_bio_2=out2d_bio,
                #res_yield_1_avg= res_avg, res_yield_1_avg_fit = res_avg_fit,
                #res_yield_2_avg= res_2_avg, res_yield_2_avg_fit = res_2_avg_fit,
                res_bio_1_avg=res_bio_avg,res_bio_1_avg_fit = res_bio_avg_fit,
                res_bio_2_avg=res_bio_2_avg,res_bio_2_avg_fit = res_bio_2_avg_fit),
           paste0(outpath,"Multimodel comparison MA_",format(Sys.Date(),"%Y%m%d"), ".xlsx"),overwrite=TRUE)


# plot importance 
unique(out1b$covariate)
out1b <- out1b %>% mutate(labels = factor(covariate,levels=unique(covariate)),
                          labels=c("Biome","Crop type", " Yield metric","Agrochemical use","Practice","Crop type:Practice"))
g_imp_yields <- ggplot(data=out1b,aes(x=importance.res.,y=reorder(labels,importance.res.)))+
  geom_col(fill="grey70",width=0.8)+ggtitle("Yield (excluding crop commodity)")+
  geom_text(data=out1b[out1b$importance.res.>0.9,], aes(label=round(importance.res.,2)),nudge_x=-0.045,size=2.5) +
  geom_text(data=out1b[out1b$importance.res.<=0.9,], aes(label=round(importance.res.,2)),nudge_x=0.05,size=2.5) +
  coord_cartesian(expand=(0))+labs(x="Importance\n(sum of model weights)",y="")+ 
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +theme_bw()+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8), plot.title=element_text(hjust=0.5,size=8,face="bold"))
g_imp_yields

unique(out1d$covariate)
names(out1d)
out1d <- out1d %>% mutate(labels = factor(covariate,levels=unique(covariate)),
                          labels=c("Crop commodity","Practice", "Crop commodity:Practice","Biome","Yield metric","Agrochemical use"))
g_imp_yields_2 <- ggplot(data=out1d,aes(x=importance.res_2.,y=reorder(labels,importance.res_2.)))+
  geom_col(fill="grey70",width=0.8)+ggtitle("Yield (excluding crop type)")+
  geom_text(data=out1d[out1d$importance.res_2.>0.9,], aes(label=round(importance.res_2.,2)),nudge_x=-0.045,size=2.5) +
  geom_text(data=out1d[out1d$importance.res_2.<=0.9,], aes(label=round(importance.res_2.,2)),nudge_x=0.05,size=2.5) +
  coord_cartesian(expand=(0))+labs(x="Importance\n(sum of model weights)",y="")+ 
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +theme_bw()+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8), plot.title=element_text(hjust=0.5,size=8,face="bold"))
g_imp_yields_2

unique(out2b_bio$covariate)
out2b_bio <- out2b_bio %>% mutate(labels = factor(covariate,levels=unique(covariate)),
                          labels=c("Crop_type_C"= "Crop type","B_ground"= "Biodiversity ground relation","Pest_group"= "Biodiversity pest group", 
                                  "B_measure_group"= "Biodiversity metric", "Taxa_group" = "Biodiversity taxa",
                                  "System_T"= "Practice", "Agrochem_CT" = "Agrochemical use", 
                                  "Crop_type_C:System_T"= "Crop type:Practice",
                                  "Biome_simp"= "Biome"))
g_imp_bio <- ggplot(data=out2b_bio,aes(x=importance.res_bio.,y=reorder(labels,importance.res_bio.)))+
  geom_col(fill="grey70",width=0.8)+ggtitle("Biodiversity (excluding crop commodity)")+
  geom_text(data=out2b_bio[out2b_bio$importance.res_bio.>0.9,], aes(label=round(importance.res_bio.,2)),nudge_x=-0.03,size=2.5) +
  geom_text(data=out2b_bio[out2b_bio$importance.res_bio.<=0.9,], aes(label=round(importance.res_bio.,2)),nudge_x=+0.05,size=2.5) +
  coord_cartesian(expand=(0))+labs(x="Importance\n(sum of model weights)",y="")+ 
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +theme_bw()+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8), plot.title=element_text(hjust=0.5,size=8,face="bold"))
g_imp_bio

unique(out2d_bio$covariate)
out2d_bio <- out2d_bio %>% mutate(labels = factor(covariate,levels=unique(covariate)),
                                  labels=c("Crop_FAO_C"= "Crop commodity","B_ground"= "Biodiversity ground relation","Pest_group"= "Biodiversity pest group", 
                                           "B_measure_group"= "Biodiversity metric", "Taxa_group" = "Biodiversity taxa",
                                           "System_T"= "Practice", "Agrochem_CT" = "Agrochemical use", 
                                           "Crop_FAO_C:System_T"= "Crop commodity:Practice",
                                           "Biome_simp"= "Biome"))
g_imp_bio_2 <- ggplot(data=out2d_bio,aes(x=importance.res_bio_2.,y=reorder(labels,importance.res_bio_2.)))+
  geom_col(fill="grey70",width=0.8)+ggtitle("Biodiversity (excluding crop type)")+
  geom_text(data=out2d_bio[out2d_bio$importance.res_bio_2.>0.9,], aes(label=round(importance.res_bio_2.,2)),nudge_x=-0.03,size=2.5) +
  geom_text(data=out2d_bio[out2d_bio$importance.res_bio_2.<=0.9,], aes(label=round(importance.res_bio_2.,2)),nudge_x=+0.05,size=2.5) +
  coord_cartesian(expand=(0))+labs(x="Importance\n(sum of model weights)",y="")+ 
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +theme_bw()+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8), plot.title=element_text(hjust=0.5,size=8,face="bold"))
g_imp_bio_2

p1 <- plot_grid( g_imp_bio,g_imp_bio_2,
                    align = "v", ncol=1, nrow=2,
                    rel_widths = c(1), rel_heights=c(1/2,1/2),
                    axis=c("bl"),labels=c("a","c"),label_size=10)
p1

p2 <- plot_grid(g_imp_yields,g_imp_yields_2,
                 align = "v", ncol=1,nrow=2, 
                 rel_widths = c(1), rel_heights=c(1/2,1/2),
                 axis=c("bl"),labels=c("b", "d"),label_size=10)
p2

p <- plot_grid(p1,NULL,p2,ncol=3,rel_heights=c(1),rel_widths=c(1/2,0.03,1/2),axis=c("bl"),align="h")
p

tiff(paste0(outpath,"Fig multimodel importance.tif"),height=6,width=7.4,units="in",res=600,compression="lzw")
p
dev.off()

# BEST model ####
model1.ML <- rma.mv(yi_Y, vi_Y,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
model_yield_best <- get.models(res_2,1)[[1]]
model_yield_best
v.q = sampling.variance(d_yield,d_yield$vi_Y)
I2_3level(model_yield_best,v.q)
I2_3level(model_all,v.q)
I2_3level(model1,v.q)

model_bio_best <- get.models(res_bio,1)[[1]] #<- model2_function(data=d,mods=~System_T+Crop_FAO_extra+Pest_group-1,method="REML")
model_bio_best
I2_3level(model_bio_best,v.q.bio)
I2_3level(model1_bio,v.q.bio)

anova(model1,model_all)
logLik(model_all)

# get fitstats for the models included in supplementary and not generated by multimodel inference
library(lmtest)

lrtest(model_yield_best,model1.ML)
fitstats(model1.ML)-fitstats(model_yield_best)

model <- rma.mv(yi_Y, vi_Y, mods =~Crop_FAO_C,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
fitstats(model_yield_best)
fitstats(model)

model <- rma.mv(yi_Y, vi_Y, mods =~Crop_FAO_C+Yield_measure_group,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
model
lrtest(model_yield_best,model)

model <- rma.mv(yi_Y, vi_Y, mods =~Crop_type_C,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
model
lrtest(model_yield_best,model)

model <- rma.mv(yi_Y, vi_Y, mods =~Biome_simp+Yield_measure_group,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
lrtest(model_yield_best,model)

model <- rma.mv(yi_Y, vi_Y, mods =~Biome_simp+Crop_type_C+Yield_measure_group,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
lrtest(model_yield_best,model)

model <- rma.mv(yi_Y, vi_Y, mods =~Biome_simp,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
lrtest(model_yield_best,model)
fitstats(model)

model <- rma.mv(yi_Y, vi_Y, mods =~Biome_simp+Crop_FAO_C+System_T+Yield_measure_group,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
lrtest(model_yield_best,model)
fitstats(model)
fitstats(model)[5]-fitstats(model_yield_best)[5]

model <- rma.mv(yi_Y, vi_Y, mods =~Biome_simp+Crop_FAO_C+Yield_measure_group,random=list(~1|ID/Effect_ID), method="ML",test="t",tdist=TRUE, data=d_yield) 
lrtest(model_yield_best,model)
fitstats(model)
fitstats(model)[5]-fitstats(model_yield_best)[5]


### Sub-analysis...not used ####
table(d$Crop_type_C)
table(d$Crop_FAO_C) 
table(d$Biome) 

d_bio_crop_p_herb <- d %>% filter(Crop_type_C =="Perennial Herb")
d_bio_crop_p_shrub <- d %>% filter(Crop_type_C =="Perennial Shrub/Liana")
d_bio_crop_p_tree <- d %>% filter(Crop_type_C =="Perennial Tree")
d_bio_crop_a_herb <- d %>% filter(Crop_type_C =="Annual Herb")
d_bio_crop_a_shrub <- d %>% filter(Crop_type_C =="Annual Herb")

model2_bio_function(data=d_bio_crop_p_herb,mods=~Pest_group-1,sigma1=NA,sigma2=NA,method="REML")
model2_bio_function(data=d_bio_crop_p_shrub,mods=~Pest_group-1,sigma1=NA,sigma2=NA,method="REML")

submodel_crop_p_herb <- model1_bio_function(data=d_bio_crop_p_herb,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_p_shrub <- model1_bio_function(data=d_bio_crop_p_shrub,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_p_tree <- model1_bio_function(data=d_bio_crop_p_tree,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_a_herb <- model1_bio_function(data=d_bio_crop_a_herb,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_a_shrub <- model1_bio_function(data=d_bio_crop_a_shrub,sigma1=NA,sigma2=NA,method="REML")

submodel_crop_p_herb_biome <- model2_bio_function(data=d_bio_crop_p_herb,mods=~Biome-1,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_p_shrub_biome <- model2_bio_function(data=d_bio_crop_p_shrub,mods=~Biome-1,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_p_tree_biome <- model2_bio_function(data=d_bio_crop_p_tree,mods=~Biome-1,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_a_herb_biome <- model2_bio_function(data=d_bio_crop_a_herb,mods=~Biome-1,sigma1=NA,sigma2=NA,method="REML")
submodel_crop_a_shrub_biome <- model2_bio_function(data=d_bio_crop_a_shrub,mods=~Biome-1,sigma1=NA,sigma2=NA,method="REML")

# SENSITVITY ANALYSIS FOR MULTIVARIATE MODELS (NOT USED...) ####
model1.ML <- model1_function(data=d,sigma1=NA,sigma2=NA,method="ML")

model2_zeros_treatment <- model2_function(data=d_zeros,mods=~System_T-1,method="REML")
model2_zeros_agrochem <- model2_function(data=d_zeros,mods=~Agrochem_CT-1,method="REML")
model2_zeros_reg <- model2_function(data=d_zeros,mods=~Region.Name-1,method="REML")
model2_zeros_crop <- model2_function(data=d_zeros,mods=~Crop_type-1,method="REML")
model2_zeros_dev <- model2_function(data=d_zeros,mods=~DevelopmentStatus-1,method="REML")

model2_validity.high_treatment <- model2_function(data=d_validity.high,mods=~System_T-1,method="REML")
model2_validity.high_crop <- model2_function(data=d_validity.high,mods=~Crop_FAO_T-1,method="REML")
model2_validity.high_dev <- model2_function(data=d_validity.high,mods=~DevelopmentStatus-1,method="REML")

# note outliers were removed from d_noOutlier based on model 1 only
model2_noOutlier_treatment <- model2_function(data=d_noOutlier,mods=~System_T-1,method="REML")
model2_noOutlier_crop <- model2_function(data=d_noOutlier,mods=~Crop_FAO_T-1,method="REML")
model2_noOutlier_dev <- model2_function(data=d_noOutlier,mods=~DevelopmentStatus-1,method="REML")

