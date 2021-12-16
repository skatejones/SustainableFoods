# Probabilistic model of trade-offs ####
# multinomial logit to compute relative risk ratios
# https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/
# https://www.statology.org/interpret-relative-risk/
# Authors: SJ

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
library(lme4) # to apply GLMs
library(nnet) # to apply multinom
library("DescTools")
library(forcats)
library(MCMCglmm)
library(lmtest) # for likelihood ratio tests
library(foreign) # to apply logit models
library(stargazer) # for logit model outcome analysis
library(car)
library(report)
library(cowplot) # arranging multiple plots
library(ggpubr) # ggarrange
library(vcd)
library(MuMIn)

wd <- readline() #at the prompt, copy and paste your filepath and press enter
D:\02_Bioversity\30_SustainableFoods\R
setwd(wd)

outpath <- "./Results_tradeoffs_wLER/"

col.intervention = c("navy","forestgreen", "seagreen","gold","orange", "lightblue","purple")
col.synergies = c("#225ea8" ,"#41b6c4" ,"#a1dab4" , "#74c476") # blue to green
col.synergies = c("Lose B-Lose Y" = "darkred","Lose B-Win Y"= "lightblue","Win B-Lose Y" ="gold","Win B-Win Y"= "forestgreen")
col.synergies = c("Lose B-Lose Y" = "navy","Lose B-Win Y"= "lightblue","Win B-Lose Y" ="gold","Win B-Win Y"= "forestgreen")
col.synergies_sig_simp = c("navy", "lightblue" ,"grey70" ,"gold", "forestgreen") # blue to green with grey
col7 <- c("#00A600", "#63C600" ,"#E6E600", "#E9BD3A", "#ECB176", "#EFC2B3" ,"mediumpurple")
col2 <- c("#00A600","#E6E600")

d <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield")
d_zeros <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_zeros")
d_validity.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity")

d <- d %>% mutate(Synergies_sig = ifelse(yi_B >0 & ci.lb_B>0 & yi_Y >0 & ci.lb_Y>0,"Win B-Win Y",
                                         ifelse(yi_B >0 & ci.lb_B >0 & (yi_Y <0 & ci.ub_Y<0),"Win B-Lose Y",
                                                ifelse(yi_B >0 & ci.lb_B >0,"Win B-Equiv Y",
                                                       ifelse(yi_B <0 & ci.ub_B <0 & yi_Y>0 & ci.lb_Y>0,"Lose B-Win Y",
                                                              ifelse(yi_Y >0 & ci.lb_Y>0,"Equiv B-Win Y",
                                                                     ifelse((yi_B<0 & ci.ub_B<0) & (yi_Y<0 & ci.ub_Y<0),"Lose B-Lose Y",
                                                                            ifelse((yi_B<0 & ci.ub_B<0) & (yi_Y<=0 & ci.ub_Y>0),"Lose B-Equiv Y",
                                                                                   ifelse((yi_B<=0 & ci.ub_B>0) & (yi_Y<0 & ci.ub_Y<0),"Equiv B-Lose Y","Equiv B-Equiv Y")))))))))
  
  
d <- d %>% mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ),"Win-Equiv",
                                                ifelse(Synergies_sig %in% c("Lose B-Equiv Y" , "Equiv B-Lose Y" ),"Lose-Equiv",Synergies_sig)))
d <- d %>% mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ,  "Win B-Win Y"),"Win B-Win Y",
                                              ifelse(Synergies_sig %in% c("Lose B-Lose Y" ,"Lose B-Equiv Y" , "Equiv B-Lose Y" ), "Lose B-Lose Y",Synergies_sig)))
d <- d %>%  mutate(Synergies_sig_binary = ifelse(((yi_B >0 & ci.lb_B>0) | (yi_B <0 & ci.ub_B<0)) & ((yi_Y >0 & ci.lb_Y>0)|(yi_Y<0 & ci.ub_Y<0)),1,0 ))

unique(d$Synergies)
d$Synergies <- factor(d$Synergies,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))
unique(d_zeros$Synergies)
d_zeros$Synergies <- factor(d_zeros$Synergies,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

unique(d$Synergies_sig)
d$Synergies_sig <- factor(d$Synergies_sig,levels=c( "Lose B-Lose Y" ,"Lose B-Equiv Y" , "Equiv B-Lose Y" , "Equiv B-Equiv Y" , "Win B-Lose Y"  ,  "Lose B-Win Y" ,  "Win B-Equiv Y"  , "Equiv B-Win Y" ,  "Win B-Win Y"       ))

unique(d$Synergies_sig_simp)
d$Synergies_sig_simp <- factor(d$Synergies_sig_simp,levels=c( "Lose B-Lose Y" ,"Lose B-Win Y" , "Equiv B-Equiv Y" , "Win B-Lose Y" , "Win B-Win Y"))

table(d$Synergies_sig_simp)
table(d$Synergies_sig, d$Synergies_sig_simp)
d <- d %>% mutate(ID = as.numeric(as.character(ID)),Effect_ID = as.numeric(as.character(Effect_ID)))

# Set reference levels ####
# to most frequently occurring 
# and subset data to exclude empty cells
table(d$Biome_simp,d$Synergies) # complete
table(d$Biome_simp,d$Crop_type_C) # gaps in all except annual herb
table(d$Biome_simp,d$System_T) # complete for associated plant species only (nearly for embedded natural and intercropping)

d <- d %>% mutate(Biome = factor(Biome))
table(d$Biome)
table(d$Biome,d$Synergies)
d$Biome <- fct_relevel(d$Biome,"Tropical & Subtropical Non-forest")
#d$Biome <- fct_relevel(d$Biome,"Temperate Forests")
table(d_biome_simp$Biome,d_biome_simp$Synergies)

d <- d %>% mutate(Biome_simp = factor(Biome_simp))
table(d$Biome_simp)
table(d$Biome_simp,d$Synergies)
d$Biome_simp <- fct_relevel(d$Biome_simp,"Tropical & Subtropical Non-forest")
d$Biome_simp <- fct_relevel(d$Biome_simp,"Other",after=Inf)

table(d$Crop_FAO_C,d$Synergies)
table(d$Crop_FAO_C)
d$Crop_FAO_C <- fct_relevel(d$Crop_FAO_C,"Cereals")
d_crop_simp <- d %>% filter(Crop_FAO_C %in% c("Cereals","Fodder","Fruits","Oil-bearing crops","Vegetables")) %>%
  mutate(Crop_FAO_C = factor(Crop_FAO_C,levels=unique(Crop_FAO_C)))

table(d$Crop_type_C,d$Synergies)
table(d$Crop_type_C)
d$Crop_type_C <- fct_relevel(d$Crop_type_C,"Annual Herb")
d_crop_type_simp <- d %>% filter(!(Crop_type_C %in% c("Perennial Herb"))) %>% mutate(Crop_type_C = factor(Crop_type_C))
table(d_crop_type_simp$Crop_type_C,d_crop_type_simp$Synergies)
table(d_crop_type_simp$Biome_simp,d_crop_type_simp$Synergies)

table(d$System_T,d$Synergies)
table(d$System_T)
d$System_T <- fct_relevel(d$System_T,"Associated plant species")
d_system <- d %>% filter(!(System_T %in% c("Combined practices","Crop rotation","Cultivar mixture"))) %>%
  mutate(System_T = factor(System_T,levels=unique(System_T)))
table(d_system$System_T,d_system$Synergies)
table(d_system$Biome_simp,d_system$Synergies)

table(d$Agrochem_CT,d$Synergies)
d <- d  %>%  mutate(Agrochem_CT = factor(Agrochem_CT,levels=c("Agrochemicals","No agrochemicals","Mixed","No data")))

table(d$Taxa_group)
table(d$Taxa_group_simp)
d$Taxa_group <- fct_relevel(d$Taxa_group,"Insect (other)")

table(d$B_measure_group,d$Synergies)
d_bio_metric <- d %>% filter(B_measure_group != "Other")
d$B_measure_group <- fct_relevel(d$B_measure_group,"Abundance")
table(d_bio_metric$B_measure_group,d_bio_metric$Synergies)
table(d_bio_metric$Taxa_group_simp,d_bio_metric$Synergies) # gaps!

table(d$Pest_group,d$Synergies)
d$Pest_group <- fct_relevel(d$Pest_group,"Other")

table(d$B_ground,d$Synergies)
d$B_ground <- fct_relevel(d$B_ground,"Above")

table(d$Yield_measure_group,d$Synergies)
d$Yield_measure_group <- fct_relevel(d$Yield_measure_group,"Mass per area")

# Multinomial model is a type of GLM, 
# so the overall goodness-of-fit statistics and their interpretations and limitations apply.
# on multinomial logit models: https://www.princeton.edu/~otorres/LogitR101.pdf (I follwed this one mainly)
# on multinomial logit models in R: https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/
# on categorical predictors: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/
# example of paper using this approach: https://link.springer.com/article/10.1007/s10980-019-00775-1
# on how to interpret coefficients for categorical predictors: https://stats.stackexchange.com/questions/60817/significance-of-categorical-predictor-in-logistic-regression

# 4 categorical response outcomes, use a multinomial logit model
# reference category is LoseLose (first factor level)
# reset with: d$variable = relevel(d$variable, ref="optionA") or changing all factor levels

### Reasons and solutions for weird coef and p values ####
# Devika et al 2016: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944325/ 
# Separation problem is caused by zero or very small cell sizes
# Solution is Penalised Maximum Likelihood
# https://search.r-project.org/CRAN/refmans/MultBiplotR/html/RidgeMultinomialLogisticFit.html 
# https://rdrr.io/cran/brglm2/man/brmultinom.html 
# https://cran.r-project.org/web/packages/brglm2/vignettes/multinomial.html 
# the above good for multiple outcomes, but don't accept random effects
# Quick try here:
#install.packages("logistf")
#library(logistf) # only good for binomial outcome, and does not accept random effects 
#model <- logistf(Synergies_sig_simp~System_T, data=d[d$Synergies_sig_simp %in% c("Win B-Lose Y","Lose B-Lose Y"),],plcontrol=logistpl.control(maxit=100))
#summary(model)
#summary(multi1.system)
#exp(model$coef[1]);exp(model$ci.lower[1]);exp(model$ci.upper[1])
#exp(model$coef[2]);exp(model$ci.lower[2]);exp(model$ci.upper[2])

#install.packages("brglm2") # for multinomials
#library("brglm2")
#model <- brmultinom(Synergies_sig_simp ~ System_T ,  data = d,type = "AS_median")

#install.packages("detectseparation")
#library("detectseparation")
#update(model, method = "detect_separation") # check for separation 
#check_infinite_estimates(multi1.biome) # plot estimates against number of iterations https://cran.r-project.org/web/packages/brglm2/vignettes/multinomial.html 
#plot(y=multi1.biome$estimates,x=multi1.biome$residuals)
#plot(y=multi1.biome$estimates,x=multi1.biome$fitted.values)
#plot(x=multi1.biome$fitted.values,y=multi1.biome$residuals)

# Calculating the additional variance explained by a predictor can be done from the log-likelihood ratio values
# see: https://data.princeton.edu/wws509/notes/c6s2 

# Sort model datasets to have adequate count in each cell formed by the factors in your analysis 
# All cell frequencies should be greater than 1 and 80% or more of cells are should be greater than 5 in count.
# The presence of small or empty cells may cause the logistic model to become unstable, 
# reporting implausibly large b coefficients and odds ratios for independent variables. 
table(d$Synergies) # fine
table(d$Synergies,d$System_T) # 3 out of 7 have some zeros
table(d$Synergies,d$Crop_type_C) # complete, except two zeros for perennial herb
table(d$Synergies,d$Crop_FAO_C) # complete, except nuts/stimulants, roots tubers, and fibres
table(d$Synergies,d$Agrochem_CT) # complete!
table(d$Synergies,d$Biome) # Temperate forests, tropical grasslands, tropical forests, temperate grasslands, mediterranean ok, the rest have at least one zero cell. 
table(d$Synergies,d$Biome_simp) 

### Run models ####

# multi0: intercept only
# multi1: all predictors
# multi2: single predictors

multi0 <- multinom(Synergies ~ 1+  (1+Effect_ID|ID),data=d) 
#multi0 <- MCMCglmm(Synergies ~ 1,random=~(1|Effect_ID)+(1|ID),data=d,family="multinomial")
summary(multi0)
multi0_result <- data.frame(broom::tidy(multi0),broom::glance(multi0)) # provides Wald statistic and p-value
multi0_result
multi0_result$rrr <- round(exp(multi0_result$estimate),3)

multi0.rrr = exp(coef(multi0))
multi0.rrr # e.g. rrr = 1.176 for win B Lose Y which means this outcome is 1.17 times more likely than a lose-lose outcome (the reference category)
stargazer(multi0, type="text",out=paste(outpath,"multi0.txt"),coef=list(multi0.rrr),p.auto=FALSE) # p-values calculated with Wald tests

multi0_pp <- data.frame(Effect_ID=1,ID=1)
multi0_pp <- cbind(multi0_pp,predict(multi0,newdata=multi0_pp,type="probs"))
multi0_pp
summary(predict(multi0,newdata=multi0_pp,type="probs"))
multi0_pp$Outcome <- row.names(multi0_pp)
row.names(multi0_pp) <- NULL
colnames(multi0_pp) <- c("Effect_ID","ID","Probability","Outcome")
multi0_pp

# biodiversity-yield outcomes separating significant effect sizes
multi0_sig <- multinom(Synergies_sig ~ 1+  (1+Effect_ID|ID),data=d) 
summary(multi0_sig)
multi0_sig_result <- data.frame(broom::tidy(multi0_sig),broom::glance(multi0_sig)) # provides Wald statistic and p-value
multi0_sig_result
multi0_sig_result$rrr <- round(exp(multi0_sig_result$estimate),3)

multi0_sig.rrr = exp(coef(multi0_sig))
multi0_sig.rrr # e.g. rrr = 1.176 for win B Lose Y which means this outcome is 1.17 times more likely than a lose-lose outcome (the reference category)
stargazer(multi0_sig, type="text",out=paste(outpath,"multi0_sig.txt"),coef=list(multi0_sig.rrr),p.auto=FALSE) # p-values calculated with Wald tests

multi0_sig_pp <- data.frame(Effect_ID=1,ID=1)
multi0_sig_pp <- cbind(multi0_sig_pp,predict(multi0_sig,newdata=multi0_sig_pp,type="probs"))
multi0_sig_pp
summary(predict(multi0_sig,newdata=multi0_sig_pp,type="probs"))
multi0_sig_pp$Outcome <- row.names(multi0_sig_pp)
row.names(multi0_sig_pp) <- NULL
colnames(multi0_sig_pp) <- c("Effect_ID","ID","Probability","Outcome")
multi0_sig_pp


# Run intercept only models ####
model= nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_biome_simp)
data=d_biome_simp
run_name="d_sub_biome"

pp_multi0_function <- function(model=multi0.var,data,run_name){
  
  name <-paste0("multi0.", run_name) 
  assign(name,model,envir = .GlobalEnv)

  multi0.var.result <- data.frame(broom::tidy(model),broom::glance(model)) 
  multi0.var.result <- multi0.var.result %>% rename(Outcome="y.level")
  multi0.var.result$rrr <- exp(multi0.var.result$estimate)
  multi0.var.result <- multi0.var.result %>%
    mutate(estimate = round(estimate,3),
           std.error = round(std.error,3),
           statistic  = round(statistic,3),
           p.value = round(p.value,4),
           rrr = round(rrr,3))

  data <- data.frame(data)
  
  multi0.var.pp <- data.frame(Effect_ID=1,ID=1)
  multi0.var.pp <- cbind(multi0.var.pp,predict(model,newdata=multi0.var.pp,type="probs"))
  colnames(multi0.var.pp)[length(multi0.var.pp)] <- "Probability"
  multi0.var.pp$Outcome <- rownames(multi0.var.pp)
  multi0.var.pp <- cbind(multi0.var.pp,pred.class =predict(model,newdata=multi0.var.pp,type="class"))
  
  #results_N <- setDT(data) %>%summarise(n_studies = n_distinct(ID), n_effectsizes = n_distinct(Effect_ID))
  #multi0.var.pp <- multi0.var.pp %>% left_join(results_N,by=variable)
  #multi1.var_pp$Label <- paste0(multi1.var_pp[,c(variable)]," (",multi1.var_pp[,c("n_effectsizes")],", ",multi1.var_pp[,c("n_studies")],")")
  # setDT(multi0.var.pp)
  # multi0.var.pp <- melt(multi0.var.pp,id.vars=c("Effect_ID","ID","pred.class"),variable.name="Outcome",value.name="Probability")

  name <-paste0("multi0.", run_name,"_pp") 
  assign(name,multi0.var.pp,envir = .GlobalEnv)
  
  multi0.var.result <- multi0.var.result %>% full_join(multi0.var.pp,by=c("Outcome")) %>%
    mutate(Probability = round(Probability,3))
  multi0.var.result <- multi0.var.result[,c(1,2,ncol(multi0.var.result),3:(ncol(multi0.var.result)-1))]
  
  name <-paste0("multi0.", run_name,".result") 
  assign(name,multi0.var.result,envir = .GlobalEnv)
}

pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d), data=d, run_name="d")

pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_system), data=d_system, run_name="d_sub_system")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_crop_simp), data=d_crop_simp, run_name="d_sub_crop_fao")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_crop_type_simp), data=d_crop_type_simp, run_name="d_sub_crop_type")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_bio_metric), data=d_bio_metric, run_name="d_bio_metric")

# Model on 2 variables ####
model= nnet::multinom(Synergies ~ Biome + Taxa_group_simp + (1 + Effect_ID|ID),data=d_biome_simp)

table(d_system$Biome_simp,d_system$System_T,d_system$Synergies) # many empty cells at this level...

data=d_biome_system
run_name="d_biome1_system"
filter="Temperate Forests"
variable="System_T"

pp_multi1_sub_function <- function(data,run_name, filter="Temperate Forests",variable="System_T"){
  
  data <- data %>% filter(Biome_simp ==filter)
  
  model = nnet::multinom(Synergies ~ System_T + (1 + Effect_ID|ID),data=data)
  print(summary(model))
  
  multi.var.result <- data.frame(broom::tidy(model),broom::glance(model)) 
  multi.var.result <- multi.var.result %>% rename(Outcome="y.level") %>%
    mutate(rrr = exp(estimate),
           Biome_simp = filter) %>%
    mutate(estimate = round(estimate,3),
           std.error = round(std.error,3),
           statistic  = round(statistic,3),
           p.value = round(p.value,4),
           rrr = round(rrr,3))
  
  data <- data.frame(data)
  multi.var.pp <- data.frame(unique(data[,c(variable)])) %>% mutate(Effect_ID=1,ID=1)
  colnames(multi.var.pp)[1] <- variable
  multi.var.pp <- cbind(multi.var.pp,predict(model,newdata=multi.var.pp,type="probs"))
  multi.var.pp <- cbind(multi.var.pp,pred.class =predict(model,newdata=multi.var.pp,type="class")) %>%
    mutate(Biome_simp = filter)
  #multi.var.table <- multi.var.pp %>% select(c("Biome", variable,"pred.class"))
  
  multi.var.result <- multi.var.result %>%
    mutate(level = term) %>%mutate(level = ifelse(level %in% c("(Intercept)","1 + Effect_ID | IDTRUE"),level,
                                                  ifelse(substr(term,1,3)=="Bio",gsub("Biome*","",term),
                                                         ifelse(substr(term,1,3)=="Sys",gsub("System_T*","",term),
                                                                ifelse(substr(term,1,6)=="Crop_t",gsub("Crop_type_C*","",term),
                                                                       ifelse(ifelse(substr(term,1,6)=="Crop_F",gsub("Crop_FAO_C*","",term),level)))))))
  
  
  name <-paste0("multi1.", run_name) 
  assign(name,model,envir = .GlobalEnv)
  
  name <-paste0("multi1.", run_name,".pp") 
  assign(name,multi.var.pp,envir = .GlobalEnv)
  
  #name <-paste0("multi1.", run_name,".result") 
  #assign(name,multi.var.result,envir = .GlobalEnv)
}

unique(d$Biome_simp)
#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_biome_system), data=d_biome_system, run_name="d_sub_biome_practice")
#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_biome_crop_simp), data=d_biome_crop_simp, run_name="d_sub_biome_crop_fao")
#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_biome_crop_type_simp ), data=d_biome_crop_type_simp, run_name="d_sub_biome_crop_type")

pp_multi1_sub_function(data=d,run_name="d_biome1_system",filter="Temperate Forests",variable="System_T")
pp_multi1_sub_function(data=d,run_name="d_biome2_system",filter="Temperate Non-forest",variable="System_T")
pp_multi1_sub_function(data=d,run_name="d_biome3_system",filter="Tropical & Subtropical Forests",variable="System_T")
pp_multi1_sub_function(data=d,run_name="d_biome4_system",filter="Tropical & Subtropical Non-forest",variable="System_T")
pp_multi1_sub_function(data=d,run_name="d_biome5_system",filter="Mediterranean",variable="System_T")
pp_multi1_sub_function(data=d,run_name="d_biome6_system",filter="Other",variable="System_T") # warning empty groups
#pp_multi1_sub_function(data=d,run_name="d_biome6_system",filter="Montane",variable="System_T") # warning empty groups
#pp_multi1_sub_function(data=d,run_name="d_biome7_system",filter="Deserts",variable="System_T") # warning empty groups
#pp_multi1_sub_function(data=d,run_name="d_biome8_system",filter="Boreal",variable="System_T") # warning empty groups

sub_biome_practice_combined <- rbind(multi1.d_biome1_system.pp,multi1.d_biome2_system.pp,multi1.d_biome3_system.pp,multi1.d_biome4_system.pp,multi1.d_biome5_system.pp,multi1.d_biome6_system.pp)
sub_biome_practice_combined <- sub_biome_practice_combined %>% select(c("Biome_simp", "System_T","pred.class"))

#https://www.statmethods.net/stats/frequencies.html
sub_biome_practice_combined <- d %>% group_by(Biome_simp,System_T) %>% mutate(count_group=n()) %>% ungroup() %>%
  group_by(Biome_simp,System_T,Synergies,count_group) %>% summarise(count=n(),freq=round(n()/count_group,3)) %>% unique() %>%
  group_by(Biome_simp,System_T) %>% mutate(pred.class = ifelse(freq == max(freq),as.character(Synergies),NA)) %>% filter(!is.na(pred.class)) %>%
  mutate(pred.label = paste0(freq*100,"%, n=",count_group)) %>% 
  mutate(pred.class = ifelse((System_T == "Intercropping" & Biome_simp == "Temperate Forests"),"Lose B-Win Y,Win B-Win Y",pred.class)) %>%
  select(Biome_simp,System_T,pred.class,pred.label,freq,count,count_group) %>%  unique() %>%
  mutate(pred.class = ifelse(pred.class =="Win B-Win Y","WW",ifelse(pred.class=="Win B-Lose Y","WL",ifelse(pred.class=="Lose B-Win Y","LW",ifelse(pred.class=="Lose B-Lose Y","LL",pred.class)))))
temp1 <- sub_biome_practice_combined %>% setDT() %>% dcast(System_T~Biome_simp,value.var="pred.class") %>% select(-Other,everything(),Other)
temp2 <- sub_biome_practice_combined %>% setDT() %>% dcast(System_T~Biome_simp,value.var="pred.label") %>% select(-Other,everything(),Other)
sub_biome_practice_combined <- cbind(temp1,temp2)

sub_biome_crop_fao_combined <- d %>% group_by(Biome_simp,Crop_FAO_C) %>% mutate(count_group=n()) %>% ungroup() %>%
  group_by(Biome_simp,Crop_FAO_C,Synergies,count_group) %>% summarise(count=n(),freq=round(n()/count_group,3)) %>% unique() %>%
  group_by(Biome_simp,Crop_FAO_C) %>% mutate(pred.class = ifelse(freq == max(freq),as.character(Synergies),NA)) %>% filter(!is.na(pred.class)) %>%
  mutate(pred.label = paste0(freq*100,"%, n=",count_group)) %>% 
  select(Biome_simp,Crop_FAO_C,pred.class,pred.label,freq,count,count_group) %>%  unique() %>%
  mutate(pred.class = ifelse(pred.class =="Win B-Win Y","WW",ifelse(pred.class=="Win B-Lose Y","WL",ifelse(pred.class=="Lose B-Win Y","LW",ifelse(pred.class=="Lose B-Lose Y","LL",pred.class)))))
temp1 <- sub_biome_crop_fao_combined %>% setDT() %>% dcast(Crop_FAO_C~Biome_simp,value.var="pred.class") %>% select(-Other,everything(),Other)
temp2 <- sub_biome_crop_fao_combined %>% setDT() %>% dcast(Crop_FAO_C~Biome_simp,value.var="pred.label") %>% select(-Other,everything(),Other)
sub_biome_crop_fao_combined <- cbind(temp1,temp2)

sub_biome_crop_type_combined <- d %>% group_by(Biome_simp,Crop_type_C) %>% mutate(count_group=n()) %>% ungroup() %>%
  group_by(Biome_simp,Crop_type_C,Synergies,count_group) %>% summarise(count=n(),freq=round(n()/count_group,3)) %>% unique() %>%
  group_by(Biome_simp,Crop_type_C) %>% mutate(pred.class = ifelse(freq == max(freq),as.character(Synergies),NA)) %>% filter(!is.na(pred.class)) %>%
  mutate(pred.label = paste0(freq*100,"%, n=",count_group)) %>% 
  select(Biome_simp,Crop_type_C,pred.class,pred.label,freq,count,count_group) %>%  unique() %>%
  mutate(pred.class = ifelse(pred.class =="Win B-Win Y","WW",ifelse(pred.class=="Win B-Lose Y","WL",ifelse(pred.class=="Lose B-Win Y","LW",ifelse(pred.class=="Lose B-Lose Y","LL",pred.class)))))
temp1 <- sub_biome_crop_type_combined %>% setDT() %>% dcast(Crop_type_C~Biome_simp,value.var="pred.class") %>% select(-Other,everything(),Other)
temp2 <- sub_biome_crop_type_combined %>% setDT() %>% dcast(Crop_type_C~Biome_simp,value.var="pred.label") %>% select(-Other,everything(),Other)
sub_biome_crop_type_combined <- cbind(temp1,temp2)


write.xlsx(x=list("biome_practice"= sub_biome_practice_combined,
                  "biome_crop_fao"=sub_biome_crop_fao_combined,
                  "biome_crop_type"=sub_biome_crop_type_combined),
           file=paste0(outpath,"Table coefficients per biome.xlsx"),overwrite=TRUE)

## Double check associations (should be moved to code 1_...)
test <- prop.table(table(as.character(d$Biome_simp),as.character(d$System_T),as.character(d$Synergies)))
test <- xtabs(~Biome_simp+System_T+Synergies, data=d)
ftable(test)
summary(test)
assocstats(xtabs(~Biome_simp+System_T,data=d))
assocstats(xtabs(~Crop_FAO_C+Crop_type_C,data=d))
assocstats(xtabs(~Crop_FAO_C+System_T,data=d))
assocstats(xtabs(~System_T+Crop_type_C,data=d))
mosaic(xtabs(~Biome_simp+System_T,data=d))
assocplot(xtabs(~Biome_simp+System_T,data=d))
assocplot(xtabs(~Crop_FAO_C+Crop_type_C,data=d))

# Find best model ####
# Use AIC if N/K > 40. 
nrow(d)/6 # 129
nrow(d)/40 #19. So ok to use AIC with up to 10 variables, BUT AICc converges to AIC at large N so better to use AICc all the time # https://stats.stackexchange.com/questions/319769/is-aicc-ever-worse-than-aic
# Use AICc: https://towardsdatascience.com/introduction-to-aic-akaike-information-criterion-9c9ba1c96ced
# Use BIC to find model with fewest variables: https://www.researchgate.net/publication/8588301_AIC_model_selection_using_Akaike_weights

# 1- all variables we could consider

model.all <- nnet::multinom(Synergies ~ System_T + Crop_type_C + Biome_simp + Agrochem_CT + (1+Effect_ID|ID),
                            data=d_crop_type_simp,na.action="na.fail") # check out pdredge (parallel processing)
res <- dredge(model.all, trace=2,rank="AICc",extra=c("BIC","AIC")) 

head(res,10)
print(res)
importance(res)
out1 <- data.frame(res)
out1b <- data.frame(importance(res))
out1b$variable <- row.names(out1b)
model.average <- model.avg(res, subset= (BIC <1760), revised.var=FALSE)
model.average <- model.avg(res,cumsum(weight) <= .95)
summary(model.average)
model.average.summary <- summary(model.average)
model.average.summary$coefmat.subset
model.average.summary$coef.nmod
model.average[1]
model.average[2]
model.average_result <- data.frame(model.average.summary$coefmat.subset)
model.average_result$rrr <- round(exp(model.average_result$Estimate),3) 
model.average_result$variable <- row.names(model.average_result)
model.average$msTable

# 2- like 1 but now with Crop_FAO_C instead of Crop_type_C
model.all2 <- nnet::multinom(Synergies ~ System_T + Crop_FAO_C + Biome_simp + Agrochem_CT + (1+Effect_ID|ID),
                            data=d_crop_simp,na.action="na.fail")
res2 <- dredge(model.all2, trace=2,rank="AICc",extra=c("BIC","AIC")) 
head(res2,10)
importance(res2)
print(res2)
out2 <- data.frame(res2)
out2b <- data.frame(importance(res2))

# 3- like 1 but include biodiversity variables instead
model.all3 <- nnet::multinom(Synergies ~ System_T + Crop_type_C + Biome_simp + Agrochem_CT + 
                               B_ground + Pest_group +  (1+Effect_ID|ID),
                             data=d,na.action="na.fail")
res3 <- dredge(model.all3, trace=2,rank="AICc",extra=c("BIC","AIC")) 
head(res3,10)
importance(res3)
print(res3)
out3 <- data.frame(res3)
out3b <- data.frame(importance(res3))

write.xlsx(list(res1=out1,res1_importance=out1b,
                res2=out2,res2_importance=out2b),
                #res3=out3,res3_importance=out3b),
                #average1=model.average_result),
           paste0(outpath,"Multimodel comparison Multinom ",format(Sys.Date(),"%Y%m%d"), ".xlsx"),overwrite=TRUE)

### FIGURE: Relative risk ratios from multivariate models ####
table(d_crop_type_simp$System_T,d_crop_type_simp$Synergies)
d_crop_type_simp <- d_crop_type_simp %>% filter(!(System_T %in% c("Combined practices","Crop rotation","Cultivar mixture")))
table(d_crop_type_simp$Biome_simp,d_crop_type_simp$Synergies)
table(d_crop_type_simp$System_T)
table(d_crop_type_simp$System_T,d_crop_type_simp$Synergies)
table(d_crop_type_simp$Crop_type_C,d_crop_type_simp$Synergies)
table(d_crop_type_simp$Crop_FAO_C,d_crop_type_simp$Synergies) # gaps
table(d_crop_type_simp$Agrochem_CT,d_crop_type_simp$Synergies)
table(d_crop_type_simp$Pest_group,d_crop_type_simp$Synergies)
table(d_crop_type_simp$B_ground,d_crop_type_simp$Synergies) # gaps
d_crop_type_simp <- d_crop_type_simp %>% mutate(System_T = fct_relevel(System_T,"Associated plant species"))

table(d_crop_type_simp$Pest_group)
table(d_crop_type_simp$Pest_group,d_crop_type_simp$Synergies)

ggplot(d,aes(y=System_T,fill=System_T))+geom_bar()+ scale_fill_manual(values=hcl.colors(7, palette = "Viridis"),name="")+
  #stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5))+
  geom_text(colour = "black", size = 3.5,aes(label = paste0(round(..count../sum(..count..)*100,1),"%")),stat='count', nudge_x=30)+
  labs(x="Number of cases",y="")+scale_x_continuous(expand = expansion(add = c(0, 30)))+
  theme_pubr()+
  theme(legend.position="none")

ggplot(d,aes(y=Taxa_class,fill=Taxa_class))+geom_bar()+ scale_fill_manual(values=hcl.colors(11, palette = "Viridis"),name="")+
  #stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5))+
  geom_text(colour = "black", size = 3.5,aes(label = paste0(round(..count../sum(..count..)*100,1),"%")),stat='count', nudge_x=20)+
  labs(x="Number of cases",y="")+scale_x_continuous(expand = expansion(add = c(0, 20)))+
  theme_pubr()+
  theme(legend.position="none")

#multi1 <- nnet::multinom(Synergies ~ Biome_simp + Crop_type_C + System_T + (1+Effect_ID|ID),data=d_crop_type_simp)
multi1 <- nnet::multinom(Synergies ~ Biome_simp + Crop_type_C + System_T + Agrochem_CT + (1+Effect_ID|ID),data=d_crop_type_simp)
multi1_result <- data.frame(broom::tidy(multi1),broom::glance(multi1)) 
multi1_result$rrr <- round(exp(multi1_result$estimate),3)

model0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_crop_type_simp)
lrtest(model0,multi1)
paste0("Variance explained (%) = ",round(((model0$deviance-multi1$deviance)/model0$deviance*100),2))

multi1_result <- multi1_result %>%  mutate(ci.lb = estimate - 1.96*std.error,ci.ub = estimate + 1.96*std.error) %>%
  mutate(ci.lb.exp = exp(ci.lb),ci.ub.exp = exp(ci.ub)) %>%
  mutate(mod = ifelse(substr(term,1,3)=="Sys", gsub("System_T*","",term),
                        ifelse(substr(term,1,3)=="Bio",gsub("Biome_simp*","",term),
                               ifelse(substr(term,1,6)=="Crop_F",gsub("Crop_FAO_C*","",term),
                                      ifelse(substr(term,1,6)=="Crop_t",gsub("Crop_type_C*","",term),
                                             ifelse(substr(term,1,6)=="Agroch",gsub("Agrochem_CT*","",term),
                                                    ifelse(substr(term,1,4)=="Pest",gsub("Pest_group*","",term),term))))))) %>%
  mutate(group = ifelse(substr(term,1,3)=="Sys", "Practice",
                        ifelse(substr(term,1,3)=="Bio","Biome",
                               ifelse(substr(term,1,6)=="Crop_F","Crop commodity",
                                             ifelse(substr(term,1,6)=="Crop_t","Crop type",
                                                    ifelse(substr(term,1,6)=="Agroch","Agrochemical use",
                                                           ifelse(substr(term,1,4)=="Pest","Pest group",
                                                                  ifelse(term=="(Intercept)","Reference",term)))))))) %>%
  mutate(Label = ifelse(rrr>100 & p.value<0.001,paste0(">100,p<0.001"),
                        ifelse(rrr>100,paste0(">100,p=",round(p.value,3)),
                               ifelse(p.value<0.001,paste0(round(rrr-1,2),",p<0.001"),paste0(round(rrr-1,2),",p=",round(p.value,3)))))) %>%
  mutate(Label_sig = ifelse(is.nan(p.value),"*",ifelse(p.value<0.05,"*","")))

check <- d %>% select(Synergies, System_T,Biome_simp,Crop_type_C,Crop_FAO_C,Agrochem_CT)

round(addmargins(table(d$System_T,d$Synergies)),2)
round(addmargins(table(d$Crop_type_C,d$Synergies)),2)

write.xlsx(list(multi1=multi1_result),paste0(outpath,"Multimodel multinom results.xlsx"),overwrite=TRUE)


multi1_result <- multi1_result %>% mutate(mod = factor(mod,levels=unique(mod[order(desc(y.level),estimate)]))) # fct_relevel(mod,c("No data","Other")))
multi1_result <- multi1_result %>% mutate(group = fct_relevel(group,"Biome"))

g <- ggplot(multi1_result[ !(multi1_result$term %in% c("1 + Effect_ID | IDTRUE","Agrochem_CTNo data","Agrochem_CTMixed")),],
            aes(x=mod,y=estimate))+
  #geom_histogram(data=d,aes(y=yi_Y_pc,x=..density..,colour="Yield"),fill="white",bins=20,alpha=0.7,position="identity")+
  #geom_density(data=d,aes(y=yi_Y_pc,colour=Agrochem_CT,fill=Agrochem_CT),alpha=0.4)+
  geom_segment(aes(x=mod,xend=mod,y=0, yend=estimate))+
  geom_point(aes(colour=y.level),shape=16,size=3)+
  #geom_col(aes(fill=y.level),position=position_dodge())+
  #geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),colour="black",size=0.7,width=0.2)+  
  geom_hline(yintercept=0,linetype=2)+
  #geom_text(aes(y=ci.ub,label=str_wrap(mod,30),hjust=0),nudge_x=0,nudge_y=1, size=2.5,family="sans")+
  geom_text(aes(y=estimate,label=Label_sig,hjust=0),nudge_x=0.1,nudge_y=0.2, size=5,family="sans")+
  labs(x="",# x="Practice, crop type and biome,\n(reference associated plant species,\ncereal crops, tropical grasslands",
       y="Log odds relative to lose-lose")+
  coord_flip(clip="off")+
  #scale_y_continuous(limits=c(-44,56),expand = expansion(add = c(5, 5)))+
  scale_y_continuous(expand = expansion(add = c(2, 2)))+
  #scale_colour_manual(values=c("black","grey50","grey65","grey80"),name="")+
  scale_colour_manual(values=col.synergies[2:4],name="")+
  scale_fill_manual(values=col.synergies[2:4],name="")+
  scale_x_discrete(labels=function(x)str_wrap(x,28))+
  theme_bw()+
  facet_grid(cols=vars(factor(y.level,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y"))),
             rows=vars(group),scales="free",space="free_y")+
  #facet_grid(rows=vars(group),scales="free",space="free")+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=9,colour="black"),
        #axis.text.y=element_blank(),
        text=element_text(size=9,colour="black"),
        axis.title=element_text(size=9,colour="black"),
        axis.text.y=element_text(size=9,colour="black"),
        axis.text.x=element_text(size=9,colour="black"),
        panel.grid.minor.y = element_line(size=0.5,colour="black"),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x=element_text(size=9,face="bold",colour="black"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
        strip.background = element_blank(),
        strip.placement = "outside")
g

height=7
width=9
tiff(paste0(outpath,"Fig probability each outcome best model_wCropType.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig probability each outcome best model_wCropType.pdf"),height=height,width=width)
g
dev.off()

table(d_crop_simp$Crop_FAO_C)
table(d_crop_type_simp$Crop_FAO_C)
table(d$Crop_FAO_C)
table(d_crop_type_simp$Crop_FAO_C,d_crop_type_simp$Synergies) # gaps

d_crop_simp_2 <- d_crop_simp %>% filter(!(System_T %in% c("Combined practices","Cultivar mixture","Crop rotation")))
table(d_crop_simp_2$Crop_FAO_C,d_crop_simp_2$Synergies) # fine
table(d_crop_simp_2$Agrochem_CT,d_crop_simp_2$Synergies)
table(d_crop_simp_2$Biome_simp,d_crop_simp_2$Synergies)

d_crop_simp_2 <- d_crop_simp_2 %>% mutate(System_T = fct_relevel(System_T,"Associated plant species"))

multi1_fao <- nnet::multinom(Synergies ~ Biome_simp + Crop_FAO_C + System_T + Agrochem_CT + (1+Effect_ID|ID),data=d_crop_simp_2)
multi1_fao_result <- data.frame(broom::tidy(multi1_fao),broom::glance(multi1_fao)) 
multi1_fao_result$rrr <- round(exp(multi1_fao_result$estimate),3)

model0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_crop_simp_2)
lrtest(model0,multi1_fao)
paste0("Variance explained (%) = ",round(((model0$deviance-multi1_fao$deviance)/model0$deviance*100),2))

multi1_fao_result <- multi1_fao_result %>%  mutate(ci.lb = estimate - 1.96*std.error,ci.ub = estimate + 1.96*std.error) %>%
  mutate(ci.lb.exp = exp(ci.lb),ci.ub.exp = exp(ci.ub)) %>%
  mutate(mod = ifelse(substr(term,1,3)=="Sys", gsub("System_T*","",term),
                      ifelse(substr(term,1,3)=="Bio",gsub("Biome_simp*","",term),
                             ifelse(substr(term,1,6)=="Crop_F",gsub("Crop_FAO_C*","",term),
                                    ifelse(substr(term,1,6)=="Crop_t",gsub("Crop_type_C*","",term),
                                           ifelse(substr(term,1,6)=="Agroch",gsub("Agrochem_CT*","",term),term)))))) %>%
  mutate(group = ifelse(substr(term,1,3)=="Sys", "Practice",
                        ifelse(substr(term,1,3)=="Bio","Biome",
                               ifelse(substr(term,1,6)=="Crop_F","Crop commodity",
                                      ifelse(substr(term,1,6)=="Crop_t","Crop type",
                                             ifelse(substr(term,1,6)=="Agroch","Agrochemical use",
                                                    ifelse(term=="(Intercept)","Reference",term))))))) %>%
  mutate(Label = ifelse(rrr>100 & p.value<0.001,paste0(">100,p<0.001"),
                        ifelse(rrr>100,paste0(">100,p=",round(p.value,3)),
                               ifelse(p.value<0.001,paste0(round(rrr-1,2),",p<0.001"),paste0(round(rrr-1,2),",p=",round(p.value,3)))))) %>%
  mutate(Label_sig = ifelse(is.nan(p.value),"*",ifelse(p.value<0.05,"*","")))

multi1_fao_result <- multi1_fao_result %>% mutate(mod = factor(mod,levels=unique(mod[order(desc(y.level),estimate)]))) # fct_relevel(mod,c("No data","Other")))
multi1_fao_result <- multi1_fao_result %>% mutate(group = fct_relevel(group,"Biome"))

g <- ggplot(multi1_fao_result[ !(multi1_fao_result$term %in% c("1 + Effect_ID | IDTRUE","Agrochem_CTNo data","Agrochem_CTMixed")),],
            aes(x=mod,y=estimate))+
  #geom_histogram(data=d,aes(y=yi_Y_pc,x=..density..,colour="Yield"),fill="white",bins=20,alpha=0.7,position="identity")+
  #geom_density(data=d,aes(y=yi_Y_pc,colour=Agrochem_CT,fill=Agrochem_CT),alpha=0.4)+
  geom_segment(aes(x=mod,xend=mod,y=0, yend=estimate))+
  geom_point(aes(colour=y.level),shape=16,size=3)+
  #geom_col(aes(fill=y.level),position=position_dodge())+
  #geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),colour="black",size=0.7,width=0.2)+  
  geom_hline(yintercept=0,linetype=2)+
  #geom_text(aes(y=ci.ub,label=str_wrap(mod,30),hjust=0),nudge_x=0,nudge_y=1, size=2.5,family="sans")+
  geom_text(aes(y=estimate,label=Label_sig,hjust=0),nudge_x=0.1,nudge_y=0.2, size=5,family="sans")+
  labs(x="",# x="Practice, crop type and biome,\n(reference associated plant species,\ncereal crops, tropical grasslands",
       y="Log odds relative to lose-lose")+
  coord_flip(clip="off")+
  #scale_y_continuous(limits=c(-44,56),expand = expansion(add = c(5, 5)))+
  scale_y_continuous(expand = expansion(add = c(2, 2)))+
  #scale_colour_manual(values=c("black","grey50","grey65","grey80"),name="")+
  scale_colour_manual(values=col.synergies[2:4],name="")+
  scale_fill_manual(values=col.synergies[2:4],name="")+
  scale_x_discrete(labels=function(x)str_wrap(x,28))+
  theme_bw()+
  facet_grid(cols=vars(factor(y.level,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y"))),
             rows=vars(group),scales="free",space="free_y")+
  #facet_grid(rows=vars(group),scales="free",space="free")+
  theme(legend.background=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=9,colour="black"),
        #axis.text.y=element_blank(),
        text=element_text(size=9,colour="black"),
        axis.title=element_text(size=9,colour="black"),
        axis.text.y=element_text(size=9,colour="black"),
        axis.text.x=element_text(size=9,colour="black"),
        panel.grid.minor.y = element_line(size=0.5,colour="black"),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x=element_text(size=9,face="bold",colour="black"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
        strip.background = element_blank(),
        strip.placement = "outside")
g

height=7
width=9
tiff(paste0(outpath,"Fig probability each outcome best model_wCropFAO.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig probability each outcome best model_wCropFAO.pdf"),height=height,width=width)
g
dev.off()

# run likelihood ratio test
# log-likelihood is a measure of how much unexplained variability there is in the data
# so higher values are worse, lower values are better.
# The difference or change in log-likelihood (Chi sq value) indicates how much new variance has been explained by the model
# https://www.statology.org/likelihood-ratio-test-in-r/
# https://www.bookdown.org/chua/ber642_advanced_regression/multinomial-logistic-regression.html - following this one
# if result is significant, we should reject the null hypothesis that model 1 and 2 are equally good
# and use model 2 which is a better fit
# if not significant, use model that has fewest variables
# Then check if individual variables explain a significant amount of the variance

# Goodness of fit by comparing the residual deviance statistic 
# for full model against null (intercept only) model , e.g. using ANOVA.
# a model with lower deviance is better
# See https://www.research-training.net/addedfiles/2013aManchester/POmnlNotes.pdf p.9
# https://online.stat.psu.edu/stat504/book/export/html/788 

anova(multi0,multi1) # LR 800, p=0, so first model is NOT at least as good as the second model (null hypothesis rejected)
lrtest(multi0,multi1) # significant, and large reduction in unexplained variability under multi1

# Classification table to see how well the predicted values match the actual values
# Can also check the percentage of variation in the response variable that is explained by the model
# using Cox and Snell R square - but this is less reliable
# PseudoR2(multi1, which = c("CoxSnell")) # doesn't work for nnet models
chisq.test(d$Synergies,predict(multi1))
ctable <- table(d$Synergies,predict(multi1))
ctable
table(d$Synergies)

## All sig variables - NO LONGER USED ####
unique(d$Crop_FAO_C)
d$Crop_FAO_C <- fct_relevel(d$Crop_FAO_C,"Fibres")
multi1_sig <- nnet::multinom(Synergies ~ Biome + Crop_type_C + System_T + Agrochem_CT -1 + (1+Effect_ID|ID),data=d)
summary(multi1_sig)
multi1_sig_result <- data.frame(broom::tidy(multi1_sig),broom::glance(multi1_sig)) 
multi1_sig_result
multi1_sig_result$rrr <- round(exp(multi1_sig_result$estimate),3)
multi1_sig_result$odds <- round(multi1_sig_result$rrr - 1,3)

d$Crop_FAO_C <- fct_relevel(d$Crop_FAO_C,"Cereals")
multi1_sig_1 <- nnet::multinom(Synergies ~ Biome + Crop_FAO_C + System_T + Agrochem_CT -1 + (1+Effect_ID|ID),data=d)
multi1_sig_1_result <- data.frame(broom::tidy(multi1_sig_1),broom::glance(multi1_sig_1)) 
multi1_sig_1_result <- multi1_sig_1_result %>%
  mutate(rrr = round(exp(estimate),3),
         odds = round(rrr - 1,3))


multi1_sig.pp <- predict(multi1_sig,type="probs")
summary(multi1_sig.pp)
unique(d$Biome)
unique(d$Crop_FAO_C)
unique(d$System_T)
unique(d$Agrochem_CT)
newdata <- data.frame(Biome=rep("Tropical & Subtropical Forests",56),
                      Crop_FAO_C=rep(unique(d$Crop_FAO_C),7),
                      System_T = rep(unique(d$System_T),8),
                      Agrochem_CT = rep("No pesticides or fertilisers",56),
                      Effect_ID = rep(1,56),ID = rep(1,56))
newdata[,c("pred.prob")] <- predict(multi1_sig,newdata=newdata,type="probs")
newdata[,c("pred.class")] <- predict(multi1_sig,newdata=newdata,type="class")


write.xlsx(list(multi1_sig=multi1_sig_result,multi1_sig_1=multi1_sig_1_result,multi1_BIC_result=multi1_result),
           paste0(outpath,"Multimodel sig and best models multinom.xlsx"),overwrite=TRUE)

#### run univariate models ####

# compare against intercept only model (same as in the dredge above)
# If excluding intercept, we are testing if the coefficients of each level are equal to zero 
# and if significant, can conclude the level has an effect on outcome.
# If including intercept, we are testing pairwise differences between the coefficients of each level and 
# the coefficient of the reference class, so can only 
# conclude that effects on outcome are not equal across classes (called  DUMMY CODING OR TREATMENT CONTRASTS). 
# See here for dummy coding options: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
# Change reference level with relevel(variable, ref="reference level name").
# Wald tests are used to test for differences between coefficients in both cases.
# To test if whole variable improves model fit, need to run a likelihood ratio test using anova or lrtest.

table(d$Synergies, predict(model))
table(d$Synergies, predict(multi0))
table(d$Synergies, predict(multi1))
table(d$Synergies, predict(multi1_sig))
freq(d$Synergies)
multi0_pp
fitstats(multi0_pp)

addmargins(table(d$Biome,d$Synergies))
multi1.biome_pp[,c("Biome","Outcome","Probability")] # gives exactly the same probabilities as proportions in original data.

model= nnet::multinom(Synergies ~ Biome - 1 + (1 + Effect_ID|ID),data=d)
variable="Biome"
variable_short="biome"

data=d_temp
model= nnet::multinom(Synergies ~ Taxa_group - 1+ (1+Effect_ID|ID),data=d_temp)
model0=nnet::multinom(Synergies ~ 1+ (1+Effect_ID|ID),data=d_temp)
variable="Taxa_group"
variable_short="Biome_tf_taxa_simp"


pp_function <- function(data,model=multi1.var,model0=multi0,variable,variable_short){
 
  print(paste0("Multinomial model for ",variable," with cases n=",nrow(data)))
  #print(summary(model))
  print(anova(model0,model))
  print(lrtest(model0,model))
  print(paste0("Variance explained (%) = ",round(((model0$deviance-model$deviance)/model0$deviance*100),2)))
  name <-paste0("multi1.", variable_short) 
  assign(name,model,envir = .GlobalEnv)
  
  #multi1.var.lr <- data.frame(broom::tidy(lrtest(model0,model)),model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d)",model2=paste0("multinom(Synergies ~ )",variable," + (1+Effect_ID|ID),data=d")) 
  #name <-paste0("multi1.", variable_short,".lr") 
  #assign(name,multi1.var.lr,envir = .GlobalEnv)
  
  multi1.var.result <- data.frame(broom::tidy(model),broom::glance(model)) 
  multi1.var.result <- multi1.var.result %>% rename(Outcome="y.level")
  multi1.var.result$rrr <- exp(multi1.var.result$estimate)
  multi1.var.result <- multi1.var.result %>%
    mutate(estimate = round(estimate,3),
           std.error = round(std.error,3),
           statistic  = round(statistic,3),
           p.value = round(p.value,4),
           rrr = round(rrr,3))
  multi1.var.result <- multi1.var.result %>%
    mutate(level = term) %>%mutate(level = ifelse(level %in% c("(Intercept)","1 + Effect_ID | IDTRUE"),
                                                  level,substr(level,nchar(variable)+1,100))) %>%
    mutate(level = ifelse(level=="(Intercept)",levels(d[[variable]])[1],level))
  
  data <- data.frame(data)
  multi1.var_pp <- data.frame(unique(data[,c(variable)]), Effect_ID=1,ID=1) # if code blocks here, type in variable name
  colnames(multi1.var_pp)[1] <- variable
  multi1.var_pp <- cbind(multi1.var_pp,predict(model,newdata=multi1.var_pp,type="probs"))
  multi1.var_pp <- cbind(multi1.var_pp,pred.class =predict(model,newdata=multi1.var_pp,type="class"))
  #multi1.var_pp
  
  #results_N <- setDT(data) %>%group_by_at(variable) %>%summarise(n_studies = n_distinct(ID),n_effectsizes = n_distinct(Effect_ID))
  
  #multi1.var_pp <- multi1.var_pp %>% left_join(results_N,by=variable)
  #multi1.var_pp$Label <- paste0(multi1.var_pp[,c(variable)]," (",multi1.var_pp[,c("n_effectsizes")],", ",multi1.var_pp[,c("n_studies")],")")
  multi1.var_pp$Label_simp <- paste0(multi1.var_pp[,c(variable)])
  setDT(multi1.var_pp)
  #multi1.var_pp <- melt(multi1.var_pp,id.vars=c(variable, "Effect_ID","ID","pred.class","n_effectsizes","n_studies", "Label","Label_simp"),variable.name="Outcome",value.name="Probability")
  multi1.var_pp <- melt(multi1.var_pp,id.vars=c(variable, "Effect_ID","ID","pred.class","Label_simp"),variable.name="Outcome",value.name="Probability")
  multi1.var_pp$level <- multi1.var_pp[[variable]]
  
  name <-paste0("multi1.", variable_short,"_pp") 
  assign(name,multi1.var_pp,envir = .GlobalEnv)
  
  multi1.var.result <- multi1.var.result %>% full_join(multi1.var_pp,by=c("level","Outcome")) %>%
    mutate(Probability = round(Probability,3))
  multi1.var.result <- multi1.var.result[,c(1,2,ncol(multi1.var.result),3:(ncol(multi1.var.result)-1))]
  
  name <-paste0("multi1.", variable_short,".result") 
  assign(name,multi1.var.result,envir = .GlobalEnv)
}

unique(d$Biome)
#d$Biome <- fct_relevel(d$Biome,"Temperate Grasslands, Savannas & Shrublands")
table(d$Synergies,d$System_T)

pp_function(data=d,model= nnet::multinom(Synergies ~ Biome_simp -1 + (1+Effect_ID|ID),data=d),model0=multi0,variable="Biome_simp",variable_short="biome_simp")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ Biome -1 + (1+Effect_ID|ID),data=d),model0=multi0_sig, variable="Biome",variable_short="biome_sig")

#pp_function(data=d,model= nnet::multinom(Synergies ~ System_T -1 + (1+Effect_ID|ID),data=d),model0=multi0,variable="System_T",variable_short="system")
pp_function(data=d_system,model= nnet::multinom(Synergies ~ System_T -1 + (1+Effect_ID|ID),data=d_system),model0=nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_system),variable="System_T",variable_short="system")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ System_T -1 + (1+Effect_ID|ID),data=d),model0=multi0_sig, variable="System_T",variable_short="system_sig")

#pp_function(data=d,model= nnet::multinom(Synergies ~ Crop_FAO_C - 1+ (1+Effect_ID|ID),data=d),variable="Crop_FAO_C",variable_short="crop_taxon")
pp_function(data=d_crop_simp,model= nnet::multinom(Synergies ~ Crop_FAO_C - 1+ (1+Effect_ID|ID),data=d_crop_simp),model0=nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_crop_simp),variable="Crop_FAO_C",variable_short="crop_taxon")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ Crop_FAO_C - 1+ (1+Effect_ID|ID),data=d),model0=multi0_sig, variable="Crop_FAO_C",variable_short="crop_taxon_sig")

#pp_function(data=d,model= nnet::multinom(Synergies ~ Crop_type_C -1 + (1+Effect_ID|ID),data=d),variable="Crop_type_C",variable_short="crop")
pp_function(data=d_crop_type_simp,model= nnet::multinom(Synergies ~ Crop_type_C -1 + (1+Effect_ID|ID),data=d_crop_type_simp),model0=nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_crop_type_simp),variable="Crop_type_C",variable_short="crop")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ Crop_type_C -1 + (1+Effect_ID|ID),data=d),model0=multi0_sig,variable="Crop_type_C",variable_short="crop_sig")

pp_function(data=d,model= nnet::multinom(Synergies ~ Agrochem_CT -1 + (1+Effect_ID|ID),data=d),variable="Agrochem_CT",variable_short="agrochem")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ Agrochem_CT -1 + (1+Effect_ID|ID),data=d),model0=multi0_sig,variable="Agrochem_CT",variable_short="agrochem_sig")

pp_function(data=d,model= nnet::multinom(Synergies ~ System_T_action - 1+ (1+Effect_ID|ID),data=d),variable="System_T_action",variable_short="system_action")

pp_function(data=d,model= nnet::multinom(Synergies ~ Region.Name -1+ (1+Effect_ID|ID),data=d),variable="Region.Name",variable_short="Region")

pp_function(data=d,model= nnet::multinom(Synergies ~ DevelopmentStatus -1 + (1+Effect_ID|ID),data=d),variable="DevelopmentStatus",variable_short="development")

pp_function(data=d,model= nnet::multinom(Synergies ~ B_ground -1 + (1+Effect_ID|ID),data=d),variable="B_ground",variable_short="b_ground")

pp_function(data=d,model= nnet::multinom(Synergies ~ Pest_group - 1+ (1+Effect_ID|ID),data=d),variable="Pest_group",variable_short="pest_group")

pp_function(data=d,model= nnet::multinom(Synergies ~ Yield_measure_group - 1 + (1+Effect_ID|ID),data=d),variable="Yield_measure_group",variable_short="yield_metric")

#pp_function(model= nnet::multinom(Synergies ~ B_measure_group - 1+ (1+Effect_ID|ID),data=d),variable="B_measure_group",variable_short="b_metric")
pp_function(data=d_bio_metric, model= nnet::multinom(Synergies ~ B_measure_group - 1+ (1+Effect_ID|ID),data=d_bio_metric),model0=nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_bio_metric),variable="B_measure_group",variable_short="b_metric")
#pp_function(data=d,model= nnet::multinom(Synergies_sig ~ B_measure_group - 1+ (1+Effect_ID|ID),data=d),model0=multi0_sig,variable="B_measure_group",variable_short="b_metric_sig")

pp_function(data=d,model= nnet::multinom(Synergies ~ Taxa_group_simp - 1+ (1+Effect_ID|ID),data=d),model0=multi0,variable="Taxa_group_simp",variable_short="taxa_simp")

      
anova(multi0,multi1.system)
anova(multi0,multi1.biome)
anova(multi0,multi1.crop)
anova(multi0,multi1.crop_taxon)
anova(multi0,multi1.agrochem)


write.xlsx(x=list("multi0" = multi0_result,
                  "biome"=multi1.biome_simp.result,
                  "system" = multi1.system.result,
                  "crop_type" = multi1.crop.result,
                  "crop_taxon" = multi1.crop_taxon.result,
                  "agrochem" = multi1.agrochem.result,
                  "region"=multi1.Region.result,
                  "development"=multi1.development.result,
                  "yield_metric" = multi1.yield_metric.result,
                  "bio_metric" = multi1.b_metric.result,
                  "pests" = multi1.pest_group.result,
                  "ground_relation"=multi1.b_ground.result,
                  "taxa_simp"=multi1.taxa_simp.result),
           file=paste0(outpath,"Multinom results_synergies",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

model= nnet::multinom(Synergies_sig_simp ~ Crop_type_C:Taxa_group_simp + (1+Effect_ID|ID),data=d)
model
summary(model)
anova(multi0,model)
lrtest(multi0,model)
result <- data.frame(broom::tidy(model),broom::glance(model)) 
result <- result %>% rename(Outcome="y.level")
result$rrr <- exp(result$estimate)
table(d$Crop_type_C,d$Taxa_group_simp)
table(d$Crop_type_C,d$Taxa_group_simp,d$Synergies_sig_simp)
table(d$Crop_type_C,d$Synergies_sig_simp)

model= nnet::multinom(Synergies ~ Crop_type_C:System_T + (1+Effect_ID|ID),data=d)
model
summary(model)
anova(multi0,model)
lrtest(multi0,model)

model= nnet::multinom(Synergies ~ Crop_type_C:Crop_FAO_C + (1+Effect_ID|ID),data=d)
model
summary(model)
anova(multi0,model)
lrtest(multi0,model)

multi1.lat <- nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d)
summary(multi1.lat) 
lrtest(multi0,multi1.lat) # significant, Chi sq  186.9
multi1.lat.lr <- data.frame(broom::tidy(lrtest(multi0,multi1.system)),
                            model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d)",
                            model2="multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d)") 
multi1.lat.result <- data.frame(broom::tidy(multi1.lat),broom::glance(multi1.lat)) 
multi1.lat.result$rrr <- round(exp(multi1.lat.result$estimate),3)
multi1.lat_pp <- data.frame( Lat_T=rep(c(-30:70)), Effect_ID=1,ID=1)
multi1.lat_pp <- cbind(multi1.lat_pp,predict(multi1.lat,newdata=multi1.lat_pp,type="probs",se=TRUE))
setDT(multi1.lat_pp)
multi1.lat_pp <- melt(multi1.lat_pp,id.vars=c("Lat_T","Effect_ID","ID"),variable.name="Outcome",value.name="Probability")

g <- ggplot(multi1.lat_pp,aes(x=Lat_T,y=Probability,colour=Outcome))+geom_line(size=1)+
  scale_colour_manual(values=col.synergies,name="")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  labs(x="Latitude")+theme_minimal()+theme_bw()+
  theme(line=element_line(size=0.3),
        #axis.text.x=element_text(angle=90,hjust=1,vjust=0.2),
        panel.grid.major.y=element_line(size=0.3,colour="grey50",linetype="dashed"),
        strip.background=element_rect(fill="white"),
        text=element_text(size=10),
        legend.position="bottom",legend.title=element_blank())+
  guides(colour=guide_legend(nrow=2))
g
height=3
width=3.6
tiff(paste0(outpath,"Fig probability each outcome by lat.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig probability each outcome by lat.pdf"),height=height,width=width)
g
dev.off()

#### Figure: proportion of each outcome with significance as per univariate multinomial models ####
data=d
synergies_variable="Synergies"
multi_data=multi1.system.result
variable="System_T"
name="system"
pos.legend="bottom"
nrow.legend=1
x_name="Diversification practice"

data=d
multi_data=multi1.b_ground.result
variable="B_ground"
name="ground"
height=3
width=3.6
x_name="Biodiversity relation to ground"

multi1.stacked <- function(data, synergies_variable="Synergies", multi_data,variable,name, x_name,height=3, width=7.4,pos.legend="none",nrow.legend=1){
  
  results_N <- setDT(data) %>%
    group_by_at(variable) %>%
    summarise(n_studies = n_distinct(ID),
              n_effectsizes = n_distinct(Effect_ID))
  results_N$Label <- paste0(results_N[[variable]]," (",results_N$n_effectsizes, ";",results_N$n_studies,")")
  results_N$Label_simp <- paste0(results_N[[variable]])
  
  data <- data %>% left_join(results_N[,c(variable,"Label","Label_simp", "n_effectsizes")],by=c(variable))
  
  #data$Synergies <- factor(data$Synergies,levels = c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y"))
  #data$Synergies <- factor(data$Synergies,levels = c("Lose B-Lose Y" ,"Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))
  
  Synergies_count <- data %>% group_by_at(c(synergies_variable,variable)) %>%
    summarise(Synergies_count = n()) 
  
  data <- data %>% left_join(Synergies_count,by=c(synergies_variable,variable)) 
  data <- data %>%
    mutate(Synergies_prop = ifelse(Synergies=="Win B-Win Y",Synergies_count/n_effectsizes*100,
                                   ifelse(Synergies=="Lose B-Lose Y",Synergies_count/n_effectsizes*-1,0)))
  #data$Label <- factor(data$Label,levels=data$Label[order(data$Synergies_prop,decreasing=TRUE)])
  
  # add pvalue info
  multi_data <- multi_data %>% filter(term !="1 + Effect_ID | IDTRUE" & estimate>0  & p.value<0.05) %>% #& Outcome != "Lose B-Lose Y"
    mutate(Label_rrr_p = ifelse(Outcome =="Win B-Win Y","WW*",
                                ifelse(Outcome == "Win B-Lose Y","WL*",
                                       ifelse(Outcome =="Equiv B-Equiv Y","EE*",
                                              ifelse(Outcome=="Lose B-Win Y","LW*",
                                                     ifelse(Outcome =="Lose B-Lose Y","LL*",NA))))))%>%
    mutate(Outcome = ifelse(Outcome =="Win B-Win Y","WW",
                            ifelse(Outcome == "Win B-Lose Y","WL",
                                   ifelse(Outcome=="Lose B-Win Y","LW",
                                          ifelse(Outcome =="Equiv B-Equiv Y","EE",
                                                 ifelse(Outcome =="Lose B-Lose Y","LL","check"))))))
  multi_data <- setDT(multi_data) %>% dcast(level~Outcome,value.var=c("Label_rrr_p")) 
  multi_data[is.na(multi_data)] <- ""
  colnames(multi_data)[1] <- variable
  cols_to_paste <- colnames(multi_data[,2:length(multi_data)])
  multi_data <- data.frame(multi_data)
  dim(multi_data)
  if(dim(multi_data)[2]<3){
    multi_data$Label_rrr_p <-multi_data[,c(cols_to_paste)]
  } else{
    multi_data$Label_rrr_p <- apply(multi_data[,c(cols_to_paste)],1, paste,collapse=" ")
    }

  data <- data %>% left_join(multi_data[,c(variable,"Label_rrr_p")])
  
  
  if(synergies_variable=="Synergies"){
    g <- ggplot(data,aes(x=reorder(str_wrap(Label_simp,12),-Synergies_prop),fill=Synergies),colour=NA)+
      geom_bar(position="fill",width=0.8)+
      geom_text(aes(y=1.03,label=Label_rrr_p, hjust=0.5),size=2,fontface="plain",colour="black")+
      #geom_text(aes(y=..count..,label=Label_rrr_p, hjust=0.5),position=position_stack(), size=2,fontface="plain",colour="black")+
      labs(y="Proportion of cases",x="")+#x=str_wrap(x_name,33),
      scale_y_continuous(expand=c(0,0))+
      scale_fill_manual(values=col.synergies)+#scale_colour_manual(values=col.synergies)+
      coord_cartesian(clip="off")+
      theme_classic2()+
      theme(line=element_line(size=0.3,colour="black"),
            axis.text=element_text(size=8,colour="black"),
            axis.title=element_text(size=8,colour="black"),
            text=element_text(size=8,colour="black"),
            legend.position=pos.legend,legend.title=element_blank(),
            plot.margin = unit(c(6,0,0,0), "pt"))+
      guides(fill=guide_legend(nrow=nrow.legend))
    print(g)
    
    g2 <- ggplot(data,aes(x=reorder(str_wrap(Label_simp,12),Synergies_prop),fill=Synergies),colour=NA)+
      geom_bar(position="fill",width=0.8)+
      geom_text(aes(y=1.03,label=Label_rrr_p, hjust=0.5,angle=270),size=2,fontface="plain",colour="black")+
      labs(y="Proportion of cases",x="")+#x=str_wrap(x_name,33),
      scale_y_continuous(expand=c(0,0))+
      scale_fill_manual(values=col.synergies)+#scale_colour_manual(values=col.synergies)+
      coord_flip(clip="off")+
      theme_classic2()+
      theme(line=element_line(size=0.3,colour="black"),
            axis.text=element_text(size=8,colour="black"),
            axis.title=element_text(size=8,colour="black"),
            text=element_text(size=8,colour="black"),
            legend.position=pos.legend,legend.title=element_blank(),
            plot.margin = unit(c(0,6,0,0), "pt"))+
      guides(fill=guide_legend(nrow=nrow.legend))
    print(g2)
    
  }
  
  if(synergies_variable=="Synergies_sig_simp"){
    g <- ggplot(data,aes(x=reorder(str_wrap(Label_simp,12),-Synergies_prop),fill=Synergies_sig_simp),colour=NA)+
      geom_bar(position="fill",width=0.8)+
      geom_text(aes(y=1.03,label=Label_rrr_p, hjust=0.5),size=2,fontface="plain",colour="black")+
      labs(y="Proportion of cases",x="")+#x=str_wrap(x_name,33),
      scale_y_continuous(expand=c(0,0))+
      scale_fill_manual(values=col.synergies_sig_simp)+#scale_colour_manual(values=col.synergies_sig_simp)+
      coord_cartesian(clip="off")+
      theme_classic2()+
      theme(line=element_line(size=0.3,colour="black"),
            axis.text=element_text(size=8,colour="black"),
            axis.title=element_text(size=8,colour="black"),
            text=element_text(size=8,colour="black"),
            legend.position=pos.legend,legend.title=element_blank(),
            plot.margin = unit(c(6,0,0,0), "pt")
      )+
      guides(fill=guide_legend(nrow=nrow.legend))
    print(g)
  }
  assign(paste0("g_",name),g,envir = .GlobalEnv)
  assign(paste0("g2_",name),g2,envir = .GlobalEnv)
  
  tiff(paste0(outpath,"Fig proportion each outcome by ",name,"_",synergies_variable,".tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g)
  dev.off()
  
  tiff(paste0(outpath,"Fig proportion each outcome by ",name,"_",synergies_variable,"_flip.tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g2)
  dev.off()
  #pdf(paste0(outpath,"Fig proportion each outcome by ",name,"_stacked_test.pdf"),height=height,width=width)
  #print(g)
  #dev.off()
}


# https://stats.idre.ucla.edu/stata/output/multinomial-logistic-regression/
# https://stats.stackexchange.com/questions/63222/getting-p-values-for-multinom-in-r-nnet-package 
# https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/
# The test statistic z is the ratio of the Coef. to the Std. Err. of the respective predictor, 
# and the p-value P>|z| is the probability the z test statistic (or a more extreme test statistic) 
# would be observed under the null hypothesis. Tested using a 2-tailed Wald test.

#z <- summary(multi1.agrochem)$coefficients/summary(multi1.agrochem)$standard.errors
#p <- (1 - pnorm(abs(z), 0, 1)) * 2

multi1.stacked(data=d,multi_data=multi1.system.result, variable="System_T",name="system",pos.legend="bottom",nrow.legend=1,height=4, width=5,x_name="Diversification practice")
#multi1.stacked(data=d_system,multi_data=multi1.system.result, variable="System_T",name="system_red",pos.legend="bottom",nrow.legend=1,height=4, width=5,x_name="Diversification practice")

multi1.stacked(data=d,multi_data=multi1.biome_simp.result, variable = "Biome_simp", name="biome",x_name="Biome",height=4, width=5)

multi1.stacked(data=d,multi_data=multi1.crop.result, variable = "Crop_type_C",name="crop",x_name="Cropping system",pos.legend="bottom",height=4, width=5)
#multi1.stacked(data=d_crop_type_simp,multi_data=multi1.crop.result, variable = "Crop_type_C",name="crop_red",x_name="Cropping system",pos.legend="bottom",height=4, width=5)

multi1.stacked(data=d,multi_data=multi1.crop_taxon.result, variable = "Crop_FAO_C",name="crop_FAO",x_name="Crop taxonomic group",pos.legend="bottom",height=4, width=5)
#multi1.stacked(data=d_crop_simp,multi_data=multi1.crop_taxon.result, variable = "Crop_FAO_C",name="crop_FAO_red",x_name="Crop taxonomic group",pos.legend="bottom",height=4, width=5)


multi1.stacked(data=d,multi_data=multi1.agrochem.result, variable = "Agrochem_CT",name="agrochem",pos.legend="bottom",height=3, width=3.6,x_name="Agrochemical inputs")
multi1.stacked(data=d,multi_data=multi1.system_action.result, variable = "System_T_action",name="system_action", x_name="Diversification strategy",height=4, width=5)
multi1.stacked(data=d,multi_data=multi1.taxa_simp.result, variable="Taxa_group_simp",  name="taxa_simp",x_name="Biodiversity taxanomic group",height=4, width=5)
multi1.stacked(data=d,multi_data=multi1.taxa_simp.result, variable="Taxa_group",  name="taxa",x_name="Biodiversity taxanomic group",height=4, width=5)
multi1.stacked(data=d,multi_data=multi1.b_ground.result , variable="B_ground",  name="ground",height=3, width=3.6,x_name="Biodiversity relation to ground")
multi1.stacked(data=d,multi_data=multi1.development.result, variable = "DevelopmentStatus",name="dev",height=3, width=3.6,x_name="Country development status")
multi1.stacked(data=d,multi_data=multi1.pest_group.result, variable = "Pest_group",name="pest",height=3, width=3.6,x_name="Non-domesticated biodiversity crop pest status")
multi1.stacked(data=d,multi_data=multi1.Region.result, variable = "Region.Name",name="reg",height=3, width=3.6,x_name="Region")

multi1.stacked(data=d,multi_data=multi1.b_metric.result, variable = "B_measure_group",name="b_metric",height=4, width=5,x_name="Biodiversity metric")
multi1.stacked(data=d_bio_metric,multi_data=multi1.b_metric.result, variable = "B_measure_group",name="b_metric_red",x_name="Biodiversity metric",height=4, width=5)

multi1.stacked(data=d,multi_data=multi1.yield_metric.result, variable="Yield_measure_group",name="y_metric",x_name="Yield metric",height=4, width=5)


# Figure: Probability of each outcome as calculated by univariate multinomial models (stacked and dodged) ####
multi1.stacked_prob <- function(data, name, x_name,height=3, width=7.4,pos.legend="none",nrow.legend=1){
  g <- ggplot(data,aes(x=str_wrap(Label_simp,9),y=Probability,fill=Outcome))+
    geom_col(position="stack",width=0.8,colour="grey50")+
    labs(x=str_wrap(x_name,33),y="Probability")+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(values=col.synergies)+
    theme_bw()+
    theme(line=element_line(size=0.3),
          axis.text=element_text(size=8),
          axis.title=element_text(size=9),
          text=element_text(size=9),
          legend.position=pos.legend,legend.title=element_blank(),
          plot.margin = unit(c(12,6,12,6), "pt")
    )+
    guides(fill=guide_legend(nrow=nrow.legend))
  print(g)
  
  assign(paste0("g_",name),g,envir = .GlobalEnv)
  
  tiff(paste0(outpath,"Fig probability each outcome by ",name,"_stacked.tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g)
  dev.off()
  
  pdf(paste0(outpath,"Fig probability each outcome by ",name,"_stacked.pdf"),height=height,width=width)
  print(g)
  dev.off()
}



multi1.stacked_prob(data=multi1.system_pp,name="system",pos.legend="bottom",nrow.legend=1,x_name="Diversification practice")
multi1.stacked_prob(data=multi1.system_action_pp,name="System_T_action", x_name="Diversification strategy")
multi1.stacked_prob(data=multi1.taxa_simp_pp,name="taxa",x_name="Non-domesticated biodiversity taxonomic group")
multi1.stacked_prob(data=multi1.b_ground_pp,name="taxa_ground",height=3, width=3.6,x_name="Non-domesticated biodiversity relation to ground")
multi1.stacked_prob(data=multi1.crop_pp,name="crop",x_name="Cropping system",pos.legend="bottom")
multi1.stacked_prob(data=multi1.crop_taxon_pp,name="crop_taxon",x_name="Crop taxonomic group",pos.legend="bottom")
multi1.stacked_prob(data=multi1.agrochem_pp,name="agrochem",pos.legend="bottom",height=3, width=3.6,x_name="Agrochemical inputs")
multi1.stacked_prob(data=multi1.development_pp,name="dev",height=3, width=3.6,x_name="Country development status")
multi1.stacked_prob(data=multi1.pest_group_pp,name="pest",height=3, width=3.6,x_name="Non-domesticated biodiversity crop pest status")
multi1.stacked_prob(data=multi1.Region_pp,name="reg",height=3, width=3.6,x_name="Region")
multi1.stacked_prob(data=multi1.biome_pp, name="biome",x_name="Biome")
multi1.stacked_prob(data=multi1.b_metric_pp, name="b_metric",x_name="Biodiversity metric")
multi1.stacked_prob(data=multi1.yield_metric_pp, name="y_metric",x_name="Yield metric")


multi1.dodge <- function(data, x, name, col=col7,height=3, width=7.4){
  
  g <- ggplot(data,aes(x=x,y=Probability,fill=str_wrap(Label,15)))+
    geom_col(position="dodge",width=0.8,colour="grey50")+labs(x="")+
    scale_y_continuous(expand=c(0,0),limits=c(0,1))+
    scale_fill_manual(values=col)+
    facet_grid(~Outcome)+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          line=element_line(size=0.3),
          panel.grid.major.y=element_line(size=0.3,colour="grey50",linetype="dashed"),
          strip.background=element_rect(fill="white"),
          text=element_text(size=10),
          legend.position="bottom",legend.title=element_blank())+
    guides(fill=guide_legend(nrow=1))
  print(g)
  
  tiff(paste0(outpath,"Fig probability each outcome by ",name,"_dodge.tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g)
  dev.off()
  
  pdf(paste0(outpath,"Fig probability each outcome by ",name,"_dodge.pdf"),height=height,width=width)
  print(g)
  dev.off()
}

col8 <- c(col7,"thistle")
col12 <- c(col7,"thistle","lightblue","lightblue3","steelblue","skyblue4")
multi1.dodge(data=multi1.system_pp,x="System_T", name="system",col=col7)
multi1.dodge(data=multi1.system_action_pp,x="System_T_action", name="system",col=col7)
multi1.dodge(data=multi1.taxa_pp,x="Taxa", name="taxa",col=col12)
multi1.dodge(data=multi1.taxa_ground_pp,x="B_ground", name="taxa_ground",col=col2)
multi1.dodge(data=multi1.crop_pp,x="Crop_type_C", name="crop",col=col7)
multi1.dodge(data=multi1.crop_taxon_pp,x="Crop_FAO_C", name="crop",col=col8)
multi1.dodge(data=multi1.agrochem_pp,x="Agrochem_CT", name="agrochem",col=col7)
multi1.dodge(data=multi1.dev_pp,x="DevelopmentStatus", name="dev",col=col2)
multi1.dodge(data=multi1.pest_pp,x="Pest_group", name="pest",col=col2)
multi1.dodge(data=multi1.reg_pp,x="Region.Name", name="reg",col=col7)
multi1.dodge(data=multi1.biome_pp,x="Biome", name="biome",col=col8)
multi1.dodge(data=multi1.b_metric_pp,x="B_measure_group", name="b_metric",col=col7)
multi1.dodge(data=multi1.yield_metric_pp,x="Yield_measure_group", name="y_metric",col=col7)

### Biome sub-analysis ####

# RELATIVE RISK (or LOG ODDS) is probability of choosing one outcome category 
# over the probability of choosing the baseline category.
# The relative risk ratios are exponentiated coefficients.
# Predicted probabilities are predictions of which outcome will occur, 
# together with the probability of each outcome class. 

# multi1.lat shows that a 1 unit increase in latitude is associated with 
# a small increase (0.054) in the log odds of a win-win versus a lose-lose outcome
# relative risk (or odds) of having a win-win versus lose lose outcome increases 
# by 1.05 for a 1 unit increase in latitude, and the coefficient is significant
# (1.05 times more likely to have win-win outcome if there is a 1 unit increase in lat)

# multi1.system has coefficient 1.064 (for win-win in intercropping. 
# This shows that log odds of a win win versus lose lose outcome are exp(1.064) = 2.89 higher 
# in agroforestry (the first factor level) compared to intercropping 
# (Can say if you move from agroforestry to intercropping, you are 
# 1.89 times more likely to have win-win outcome compared to a lose lose outcome,
# or the odds of being in the win-win class as opposed to the lose-lose class are 
# 2.89 - 1 = 1.89% higher in intercropping systems compared to agroforestry.)
# In agroforestry systems, the most likely outcome is win B lose Y.
# In intercropping the most likely outcome is lose B win Y.

# Examine changes in predicted probabilities for two variables at a time ####
# BUT need to be careful, as class sizes become too small and unreliable 

### Biome x crop commodity ####
d_biome_crop_simp <- d_biome_crop_simp %>% mutate(Biome = factor(Biome,unique(Biome)),Crop_FAO_C=factor(Crop_FAO_C,unique(Crop_FAO_C)))
multi0_biome_crop <- nnet::multinom(Synergies ~ 1 +  (1+Effect_ID|ID),data=d_biome_crop_simp) 
multi1_biome_crop <- nnet::multinom(Synergies ~ Biome + Crop_FAO_C -1 +  (1+Effect_ID|ID),data=d_biome_crop_simp) 
multi1_biome_crop.lr <- data.frame(broom::tidy(lrtest(multi0_biome_crop,multi1_biome_crop)),
                            model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d_biome_crop_simp)",
                            model2="multinom(Synergies ~ Biome + Crop_FAO_C + (1+Effect_ID|ID),data=d_biome_crop_simp)") 
multi1_biome_crop_result <- data.frame(broom::tidy(multi1_biome_crop),broom::glance(multi1_biome_crop)) # provides Wald statistic and p-value
multi1_biome_crop_result$rrr <- round(exp(multi1_biome_crop_result$estimate),3)
multi1_biome_crop_result <- multi1_biome_crop_result %>%
  mutate(p.value = round(p.value,4),
         estimate = round(estimate,4),
         rrr= round(rrr,4))

multi1_biome_crop_pp <- d_biome_crop_simp %>% expand(Biome,Crop_FAO_C) %>% mutate(Effect_ID=1,ID=1)
multi1_biome_crop_pp <- cbind(multi1_biome_crop_pp,predict(multi1_biome_crop,newdata=multi1_biome_crop_pp,type="probs",se.fit=TRUE))
multi1_biome_crop_pp
summary(predict(multi1_biome_crop,newdata=multi1_biome_crop_pp,type="probs"))

setDT(multi1_biome_crop_pp)
multi1_biome_crop_pp <- melt(multi1_biome_crop_pp,id.vars=c("Biome","Crop_FAO_C", "Effect_ID","ID"),variable.name="Outcome",value.name="Probability")
multi1_biome_crop_pp$Outcome <- factor(multi1_biome_crop_pp$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

multi1_biome_crop_pp <- multi1_biome_crop_pp %>% mutate(Probability = round(Probability,4))
multi1_biome_crop_pp$Outcome <- factor(multi1_biome_crop_pp$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

table(d_biome_crop_simp$Synergies,predict(multi1_biome_crop))

multi1_biome_crop_pp <- multi1_biome_crop_pp %>% group_by(Outcome,Biome) %>%
  mutate(Outcome_max = max(Probability))

multi1_biome_crop_pp <- multi1_biome_crop_pp %>% 
  mutate(Order = ifelse(Outcome =="Win B-Win Y",(Probability/Outcome_max)*10000,
                        ifelse(Outcome =="Win B-Lose Y",(Probability/Outcome_max)*1000,
                               ifelse(Outcome=="Lose B-Lose Y",(Probability/Outcome_max)*100,(Probability/Outcome_max))))) %>%
  mutate(Order = round(Order,2))

g <- ggplot(multi1_biome_crop_pp,aes(y=Probability,x=reorder(as.character(Crop_FAO_C),-Order) ,fill=Outcome))+
  geom_col(position="stack",width=0.5)+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  scale_x_discrete(labels=function(x)str_wrap(x,10))+
  labs(x="")+theme_minimal()+
  facet_wrap(~str_wrap(Biome,22),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g

Synergy_count <- d_biome_crop_simp %>% select(Biome,Synergies) %>% group_by(Biome,Synergies) %>% mutate(Synergies_count=n()) %>% unique()

d_biome_crop_simp <- d_biome_crop_simp %>% left_join(Synergy_count,by=c("Synergies","Biome")) 

d_biome_crop_simp <- d_biome_crop_simp %>%
  group_by(Biome,Crop_FAO_C,Synergies) %>% mutate(Synergies_sub_count = n()) %>% ungroup() %>%
  mutate(Order = ifelse(Synergies =="Win B-Win Y",(Synergies_sub_count/Synergies_count)*10000,
                        ifelse(Synergies =="Win B-Lose Y",(Synergies_sub_count/Synergies_count)*1000,
                               ifelse(Synergies=="Lose B-Lose Y",(Synergies_sub_count/Synergies_count)*100,(Synergies_sub_count/Synergies_count))))) %>%
  mutate(Order = round(Order,2))

g <- ggplot(d_biome_crop_simp,aes(x=reorder(as.character(Crop_FAO_C),-Order) ,fill=Synergies))+
  geom_bar(position="stack",width=0.5)+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  scale_x_discrete(labels=function(x)str_wrap(x,10))+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,22),nrow=5,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g


### Biome x System ####
unique(d$Biome)
multi.test <- nnet::multinom(formula = Synergies ~ Biome + System_T -1+ (1 + Effect_ID | ID), data = d_biome_system)
multi.test.lr <- data.frame(broom::tidy(lrtest(multi0,multi.test)),
                            model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d_biome_system)",
                            model2="multinom(Synergies ~ Biome + System_T + (1+Effect_ID|ID),data=d_biome_system)") 
multi.test.result <- data.frame(broom::tidy(multi.test),broom::glance(multi.test)) 
multi.test.result$rrr <- round(exp(multi.test.result$estimate),3)
multi.test.result <- multi.test.result %>%
  mutate(p.value = round(p.value,4),
         estimate = round(estimate,4),
         rrr= round(rrr,4))

pp_a <- data.frame(System_T = unique(d$System_T),
                   Biome = c(rep("Tropical & Subtropical Grasslands, Savannas & Shrublands",7),
                             rep("Tropical & Subtropical Forests",7),
                             rep("Temperate Grasslands, Savannas & Shrublands",7),
                             rep("Temperate Forests",7),
                             rep("Mediterranean Forests, Woodlands & Scrub",7)), 
                    Effect_ID=1,ID=1)
pp_a <- cbind(pp_a,predict(multi.test,newdata=pp_a,type="probs"))
pp_a
class_a <- cbind(pp_a,pred.class =predict(multi.test,newdata=pp_a,type="class"))
class_a

setDT(class_a)
class_a <- melt(class_a,id.vars=c("System_T","Biome", "Effect_ID","ID","pred.class"),variable.name="Outcome",value.name="Probability")
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"), labels=c("Lose B-Lose Y","Win B-Lose Y","Lose B-Win Y","Win B-Win Y"))
results_N <- setDT(d) %>%
  group_by(System_T,Biome) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  filter(Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                      "Tropical & Subtropical Forests",
                      "Temperate Grasslands, Savannas & Shrublands",
                      "Temperate Forests","Mediterranean Forests, Woodlands & Scrub"))
class_a <- class_a %>% left_join(results_N) %>%
  mutate(Label = paste0(System_T," (",n_effectsizes, ",",n_studies,")"))
class_a <- class_a %>% mutate(Probability = round(Probability,4))
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))
table(d$Synergies,predict(multi.test))
check <- data.frame(predict(multi.test))
ggplot(data=cbind(d,data.frame(pred.class=predict(multi.test))),aes(x=Synergies))+geom_bar()+
  geom_bar(aes(x=pred.class),width=0.5,fill="orange")
 
g_tropical <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                               "Tropical & Subtropical Forests"),],
            aes(x=str_wrap(System_T,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~Biome,nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_tropical

height=3.4
width=7.4

tiff(paste0(outpath,"Fig probability by system and biome tropical.tif"),height=height,width=width,units="in",res=300,compression="lzw")
g_tropical
dev.off()

pdf(paste0(outpath,"Fig probability by system and biome tropical.pdf"),height=height,width=width)
g_tropical
dev.off()

results_N <- setDT(d) %>%
  group_by(System_T,Biome) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  mutate(Label_biome_system = paste0(System_T," (",n_effectsizes, ",",n_studies,")"))

d <- d %>% left_join(results_N[,c("System_T","Biome","Label_biome_system")],by=c("System_T","Biome"))

g_tropical <- ggplot(d[ Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                     "Tropical & Subtropical Forests"),],
                     aes(x=str_wrap(Label_biome_system,14) ,fill=Synergies))+geom_bar(position="fill",width=0.8)+
  labs(y="Proportion")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_tropical

height=3.4
width=7.4

tiff(paste0(outpath,"Fig proportion by system and biome tropical.tif"),height=height,width=width,units="in",res=300,compression="lzw")
g_tropical
dev.off()

pdf(paste0(outpath,"Fig proportion by system and biome tropical.pdf"),height=height,width=width)
g_tropical
dev.off()


g_temperate <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Temperate Grasslands, Savannas & Shrublands",
                                                                         "Mediterranean Forests, Woodlands & Scrub"),],
                     aes(x=str_wrap(Label,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~Biome,nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="none",
        legend.text=element_text(size=8))
g_temperate

g_temperate2 <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Temperate Forests"),],
                      aes(x=str_wrap(Label,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~Biome,nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="none",
        legend.text=element_text(size=8))
g_temperate2

height=5
width=7.4

tiff(paste0(outpath,"Fig probability by system and biome temperate.tif"),height=height,width=width,units="in",res=300,compression="lzw")
ggarrange(g_temperate,g_temperate2,nrow=2,ncol=1, 
          #heights=c(12,12,1),
          common.legend=TRUE,legend="bottom", 
          labels=c("",""),font.label=list(size=10))
dev.off()

pdf(paste0(outpath,"Fig probability by system and biome temperate.pdf"),height=height,width=width)
ggarrange(g_temperate,g_temperate2,nrow=2,ncol=1, 
          #heights=c(12,12,1),
          common.legend=TRUE,legend="bottom", 
          labels=c("",""),font.label=list(size=10))
dev.off()

p <- plot_grid(g_tropical,g_temperate,g_temperate2,legend,ncol=1,rel_heights=c(0.30,0.30,0.36,0.06),axis=c("bl"))
p

height=7
width=7.4

tiff(paste0(outpath,"Fig probability by system and biome.tif"),height=height,width=width,units="in",res=300,compression="lzw")
p
dev.off()

pdf(paste0(outpath,"Fig probability by system and biome.pdf"),height=height,width=width)
p
dev.off()

## Biome x crop type ###
unique(d$Crop_type_C)
unique(d$Biome)
d$Biome <- fct_relevel(d$Biome,"Tropical & Subtropical Grasslands, Savannas & Shrublands")

multi.test <- nnet::multinom(formula = Synergies ~ Biome + Crop_type_C + (1 + Effect_ID | ID), data = d)
multi.test.lr <- data.frame(broom::tidy(lrtest(multi0,multi.test)),
                            model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d)",
                            model2="multinom(Synergies ~ Biome + Crop_type_C + (1+Effect_ID|ID),data=d)") 
multi.test.result <- data.frame(broom::tidy(multi.test),broom::glance(multi.test)) 
multi.test.result$rrr <- round(exp(multi.test.result$estimate),3)
multi.test.result <- multi.test.result %>%
  mutate(p.value = round(p.value,4),
         estimate = round(estimate,4),
         rrr= round(rrr,4))


pp_a <- data.frame(Crop_type_C = unique(d$Crop_type_C),
                   Biome = c(rep("Tropical & Subtropical Grasslands, Savannas & Shrublands",5),
                             rep("Tropical & Subtropical Forests",5),
                             rep("Temperate Grasslands, Savannas & Shrublands",5),
                             rep("Temperate Forests",5),
                             rep("Mediterranean Forests, Woodlands & Scrub",5)), 
                   Effect_ID=1,ID=1)
pp_a <- cbind(pp_a,predict(multi.test,newdata=pp_a,type="probs"))
pp_a
class_a <- cbind(pp_a,pred.class =predict(multi.test,newdata=pp_a,type="class"))
class_a
setDT(class_a)
class_a <- melt(class_a,id.vars=c("Crop_type_C","Biome", "Effect_ID","ID","pred.class"),variable.name="Outcome",value.name="Probability")
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

g_biome_crop <- ggplot(class_a,
                     aes(x=str_wrap(Crop_type_C,10) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_biome_crop


results_N <- setDT(d) %>%
  group_by(Crop_type_C,Biome) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  filter(Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                      "Tropical & Subtropical Forests",
                      "Temperate Grasslands, Savannas & Shrublands",
                      "Temperate Forests","Mediterranean Forests, Woodlands & Scrub"))
class_a <- class_a %>% left_join(results_N) %>%
  mutate(Label = paste0(Crop_type_C," (",n_effectsizes, ",",n_studies,")"))

unique(class_a$Outcome)
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

g_tropical <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                                        "Tropical & Subtropical Forests"),],
                     aes(x=str_wrap(Label,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_tropical

results_N <- setDT(d) %>%
  group_by(Crop_type_C,Biome) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  mutate(Label_biome_crop_type = paste0(Crop_type_C," (",n_effectsizes, ",",n_studies,")"))

d <- d %>% left_join(results_N[,c("Crop_type_C","Biome","Label_biome_crop_type")],by=c("Crop_type_C","Biome"))

g_tropical_crop <- ggplot(d[ Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                     "Tropical & Subtropical Forests"),],
                      aes(x=str_wrap(Label_biome_crop_type,14) ,fill=Synergies))+geom_bar(position="fill",width=0.8)+
  labs(y="Proportion")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_tropical_crop

legend <- get_legend(
  g_tropical_crop + theme(legend.position="bottom",legend.direction="horizontal",
                   #legend.box.margin = margin(0, 6, 0, 0),
                   legend.text=element_text(size=9))+
    guides(fill=guide_legend(nrow=1)))


plot1 <- plot_grid(g_tropical_crop+theme(legend.position="none"), 
                   align = "vh", ncol=1, 
                   #rel_widths = c(1/2,1/2), rel_heights=c(1),
                   axis=c("bl"),labels=c("a"),label_size=10)
plot1

plot2 <- plot_grid(g_tropical+theme(legend.position="none"), 
                   align = "vh", ncol=1, 
                   #rel_widths = c(1/2,1/2), rel_heights=c(1),
                   axis=c("bl"),labels=c("b"),label_size=10)
plot2



p <- plot_grid(plot1,plot2,legend,ncol=1,rel_heights=c(0.30,0.30,0.06),axis=c("bl"))
p

height=7
width=7.4

tiff(paste0(outpath,"Fig S5 proportions by biome system crop tropical.tif"),height=height,width=width,units="in",res=300,compression="lzw")
p
dev.off()

pdf(paste0(outpath,"Fig S5 proportions by biome system crop tropical.pdf"),height=height,width=width)
p
dev.off()

g_tropical <- g_tropical_crop

height=3.5
width=7.4

tiff(paste0(outpath,"Fig probability by crop type and biome.tif"),height=height,width=width,units="in",res=300,compression="lzw")
g_tropical
dev.off()

pdf(paste0(outpath,"Fig probability by crop type and biome.pdf"),height=height,width=width)
g_tropical
dev.off()

g_temperate <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Temperate Grasslands, Savannas & Shrublands",
                                                                         "Mediterranean Forests, Woodlands & Scrub"),],
                      aes(x=str_wrap(Label,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="none",
        legend.text=element_text(size=8))
g_temperate


g_temperate <- ggplot(d[ Biome %in% c("Temperate Grasslands, Savannas & Shrublands",
                                       "Mediterranean Forests, Woodlands & Scrub"),],
                       aes(x=str_wrap(Label_biome_crop_type,14) ,fill=Synergies))+geom_bar(position="fill",width=0.8)+
  labs(y="Proportion")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_temperate


g_temperate2 <- ggplot(class_a[(class_a$n_effectsizes !=0) & Biome %in% c("Temperate Forests"),],
                       aes(x=str_wrap(Label,14) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_temperate2

g_temperate2 <- ggplot(d[ Biome %in% c("Temperate Forests"),],
                       aes(x=str_wrap(Label_biome_crop_type,14) ,fill=Synergies))+geom_bar(position="fill",width=0.8)+
  labs(y="Proportion")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g_temperate2

height=5
width=7.4

tiff(paste0(outpath,"Fig probability by crop and biome temperate.tif"),height=height,width=width,units="in",res=300,compression="lzw")
ggarrange(g_temperate,g_temperate2,nrow=2,ncol=1, 
          #heights=c(12,12,1),
          common.legend=TRUE,legend="bottom", 
          labels=c("",""),font.label=list(size=10))
dev.off()

pdf(paste0(outpath,"Fig probability by crop and biome temperate.pdf"),height=height,width=width)
ggarrange(g_temperate,g_temperate2,nrow=2,ncol=1, 
          #heights=c(12,12,1),
          common.legend=TRUE,legend="bottom", 
          labels=c("",""),font.label=list(size=10))
dev.off()

p <- plot_grid(g_tropical+theme(legend.position="none"),g_temperate,g_temperate2,legend,ncol=1,rel_heights=c(0.30,0.30,0.36,0.06),axis=c("bl"))
p

height=7
width=7.4

tiff(paste0(outpath,"Fig probability by crop and biome.tif"),height=height,width=width,units="in",res=300,compression="lzw")
p
dev.off()

pdf(paste0(outpath,"Fig probability by crop and biome.pdf"),height=height,width=width)
p
dev.off()


#### System action x crop species ###
unique(d$Crop_FAO_C)
unique(d$Crop_ann_pen_C)
unique(d$System_T_action)
table(d$Crop_FAO_C)

multi.test <- nnet::multinom(formula = Synergies ~ Crop_FAO_C + System_T_action + (1 + Effect_ID | ID), data = d)
# get prediction accuracy https://www.guru99.com/r-generalized-linear-model.html 
addmargins(table(d$Synergies,predict(multi.test)))
multi.accuracy <- sum(diag(table(d$Synergies,predict(multi.test))))/sum(table(d$Synergies,predict(multi.test))) #54%

table(d$Synergies,predict(multi1))
multi.accuracy <- sum(diag(table(d$Synergies,predict(multi1))))/sum(table(d$Synergies,predict(multi1))) # 62%

table(d$Synergies,predict(multi1_sig))
multi.accuracy <- sum(diag(table(d$Synergies,predict(multi1_sig))))/sum(table(d$Synergies,predict(multi1_sig))) # 64%

pp_a <- data.frame(System_T_action = unique(d$System_T_action),
                   Crop_FAO_C = c(rep("Cereals",5),
                             rep("Fruits",5),
                             rep("Oil-bearing crops",5),
                             rep("Vegetables",5),
                             rep("Fibres",5)), 
                   Effect_ID=1,ID=1)
pp_a <- cbind(pp_a,predict(multi.test,newdata=pp_a,type="probs"))
pp_a
class_a <- cbind(pp_a,pred.class =predict(multi.test,newdata=pp_a,type="class"))
class_a

setDT(class_a)
class_a <- melt(class_a,id.vars=c("System_T_action","Crop_FAO_C", "Effect_ID","ID","pred.class"),variable.name="Outcome",value.name="Probability")
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

results_N <- setDT(d) %>%
  group_by(Crop_FAO_C,System_T_action) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  filter(Crop_FAO_C %in% c("Cereals","Fruits","Oil-bearing crops","Vegetables","Fibres"))

class_a <- class_a %>% left_join(results_N) %>%
  mutate(Label = paste0(System_T_action," (",n_effectsizes, ",",n_studies,")"))


g <- ggplot(class_a[(class_a$n_effectsizes !=0),],
            aes(x=str_wrap(Label,10) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~Crop_FAO_C,nrow=2,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g

results_N <- setDT(d) %>%
  group_by(System_T_action,Biome) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  mutate(Label_biome_system_action = paste0(System_T_action," (",n_effectsizes, ",",n_studies,")"))

d <- d %>% left_join(results_N[,c("System_T_action","Biome","Label_biome_system_action")],by=c("System_T_action","Biome"))

g <- ggplot(d[ Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                     "Tropical & Subtropical Forests"),],
                     aes(x=str_wrap(Label_biome_system_action,14) ,fill=Synergies))+geom_bar(position="fill",width=0.8)+
  labs(y="Proportion")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~str_wrap(Biome,25),nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g



height=7.4
width=7.4

tiff(paste0(outpath,"Fig probability by crop and system action.tif"),height=height,width=width,units="in",res=300,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig probability by crop and system action.pdf"),height=height,width=width)
g
dev.off()

# Agrochemical x diversification practice ####
unique(d$Agrochem_CT)
multi.test <- nnet::multinom(formula = Synergies ~ Agrochem_CT + System_T + (1 + Effect_ID | ID), data = d)
pp_a <- data.frame(System_T = unique(d$System_T),
                   Agrochem_CT = c(rep("No pesticides or fertilisers",7),
                                  rep("No fertiliser"  ,7),
                                  rep("No pesticide" ,7),
                                  rep("Pestcides & fertilisers",7)), 
                   Effect_ID=1,ID=1)
pp_a <- cbind(pp_a,predict(multi.test,newdata=pp_a,type="probs"))
pp_a
class_a <- cbind(pp_a,pred.class =predict(multi.test,newdata=pp_a,type="class"))
class_a

setDT(class_a)
class_a <- melt(class_a,id.vars=c("System_T","Agrochem_CT", "Effect_ID","ID","pred.class"),variable.name="Outcome",value.name="Probability")
class_a$Outcome <- factor(class_a$Outcome,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

results_N <- setDT(d) %>%
  group_by(Agrochem_CT,System_T) %>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(Effect_ID)) %>%
  filter(Agrochem_CT != "No data")

class_a <- class_a %>% left_join(results_N) %>%
  mutate(Label = paste0(Agrochem_CT," (",n_effectsizes, ",",n_studies,")"))


g <- ggplot(class_a[(class_a$n_effectsizes !=0),],
            aes(x=str_wrap(Label,10) ,y=Probability,fill=Outcome))+geom_col(position="stack",width=0.8)+
  #coord_flip()+
  #ggtitle(")")+
  scale_fill_manual(values=col.synergies,name="")+
  labs(x="")+theme_minimal()+facet_wrap(~System_T,nrow=1,scales="free_x",drop=TRUE)+
  theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        legend.position="bottom",
        legend.text=element_text(size=8))
g

height=3.4
width=7.4

tiff(paste0(outpath,"Fig probability by system and agrochem.tif"),height=height,width=width,units="in",res=300,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig probability by system and agrochem.pdf"),height=height,width=width)
g
dev.off()

# Look at the averaged predicted probabilities for 
# different values of continuous predictor variables within each level of categorical ones.
test <- data.frame(System_T = rep(unique(d$System_T), each=101), Agrochem_CT = c("No"),Lat_T=rep(c(-30:70),7), Effect_ID=1,ID=1)
pp.test <- cbind(test,predict(multi1,newdata=test,type="probs",se=TRUE))
by(pp.test[,5:9],pp.test$System_T,colMeans)
setDT(pp.test)
lpp <- melt(pp.test,id.vars=c("System_T","Agrochem_T","Lat_T","Effect_ID","ID"),variable.name="Outcome",value.name="Probability")
ggplot(lpp,aes(x=Lat_T,y=Probability,colour=Outcome))+geom_line()+
  ggtitle("")+
  scale_colour_manual(values=col.synergies,name="")+
  labs(x="Latitude")+theme_minimal()+ facet_grid(~System_T)+
  theme(plot.title=element_text(hjust=0.5))
