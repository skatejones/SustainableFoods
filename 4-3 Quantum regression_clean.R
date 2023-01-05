### Run quantile regression: biodiversity on yield effect sizes ####

# Authors: SJ
# Last updated 12/12/2022

library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(readxl)
library(openxlsx)
library(quantreg)
library(qrLMM)

#### Set file paths ####

wd <- readline() #at the prompt, copy and paste your filepath and press enter

setwd(wd)

outpath <- "./Results_tradeoffs_wLER/"

#### Import data ####
d <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield")
d_validity.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity")

#### Format data for analysis ####
# from multinomial code 
d <- reclassify(d)
d_validity.high <- reclassify(d_validity.high)

# split by biodiversity metric type
d_abun <- d %>% filter(B_measure_group =="Abundance")
d_rich <- d %>% filter(B_measure_group =="Richness")
d_even <- d %>% filter(B_measure_group =="Richness-Evenness")


#### Check relationship between yield and biodiversity ####

# Quantile regression is a type of linear regression 
# to understand relationship between one or more predictors and a response variable.
# Linear regression is normally used to estimate the mean value of a response variable, but in 
# quantile regression we estimate any percentile of the response value, e.g. 70th, 90th.
# https://www.statology.org/quantile-regression-in-r/ 
# 50th percentile (median) is robust to outliers so can be better than linear regression (based on mean)

# interpreting results when dependent and independent variables are logarithms:
# we interpret the coefficient as the expected change in percent in the dependent variable 
# when the independent variable is increased by 1%
# e.g. an increase in yield of 1% is associated with an increase in biodiversity of ##%.
# https://www.stathelp.se/en/regression_logarithm_en.html

# QRLMM package allows for mixed effects modelling with quantile regression
# https://rdrr.io/cran/qrLMM/man/QRLMM.html
# paper with results from this package: https://www.researchgate.net/publication/303811593_Quantile_regression_for_mixed-effects_models/figures

# Quick plots
ggplot(d,aes(yi_Y,yi_B))+geom_point()+
  geom_abline(intercept=exp(coef(qmodel)[1]),slope=coef(qmodel)[2])+
  geom_smooth(method="lm",se=F)+theme_bw()

# set data to each of these in turn
data <- d
b_measure <- "Biodiversity"
b_measure_short <- "Bio"

data <- d_abun
b_measure <- "Abundance"
b_measure_short <- "Abun"

data <- d_rich
b_measure <- "Richness"
b_measure_short <- "Rich"

data <- d_even
b_measure <- "Richness-Evenness"
b_measure_short <- "RichEven"

data <- d_validity.high
b_measure <- "Biodiversity (HQ)"
b_measure_short <- "Bio_HQ"

data <- data %>% mutate(ID = as.factor(ID),Effect_ID=as.factor(Effect_ID))

attach(data)

y=yi_B # response
x=cbind(1,yi_Y) # design matrix for fixed effects
z=cbind(1,Effect_ID) # design matrix for random effects
#groups <- rep(1,nrow(data)) 
groups=ID
qr_values <- seq(0.1,0.9,0.1)

qmodel <- QRLMM(y=y,x=x,z=z,groups=groups,p=qr_values,show.convergence=TRUE,MaxIter=200,M=10) 
# for Bio, says convergence not reached but graphically it is, processing time 8-12 mins. 
# for Rich, says convergence not reached but graphically it is at about 50 iterations, processing time 2 mins.

# save convergence plot for quantile 0.5 for reference !!!

summary(qmodel)
check <- data.frame(t(qmodel[5][[1]]$res$beta))
qmodel[5][[1]]$res$AIC
qmodel[5][[1]]$res$BIC
qmodel[5][[1]]$res$loglik

qr_results <- data.frame(qmodel[1][[1]]$res$table) %>% mutate(quantile = 0.1)
qr_results <- qr_results %>%
  mutate(variable = rownames(qr_results)) %>%
  mutate(processing_time = qmodel[1][[1]]$res$time) %>%
  mutate(AIC = qmodel[1][[1]]$res$AIC, BIC = qmodel[1][[1]]$res$BIC,
         loglik = qmodel[1][[1]]$res$loglik,sigma = qmodel[1][[1]]$res$sigma)

for(i in c(2:9)){
  qr_results_0 <- data.frame(qmodel[i][[1]]$res$table) %>% mutate(quantile = paste0("0.",i))
  qr_results_0 <- qr_results_0 %>%
    mutate(quantile = as.numeric(quantile)) %>%
    mutate(variable = rownames(qr_results_0)) %>%
    mutate(processing_time = qmodel[i][[1]]$res$time) %>%
    mutate(AIC = qmodel[i][[1]]$res$AIC, BIC = qmodel[i][[1]]$res$BIC,
           loglik = qmodel[i][[1]]$res$loglik,sigma = qmodel[i][[1]]$res$sigma) 
  qr_results <- rbind(qr_results,qr_results_0)}

assign(paste0("qr_results_",b_measure_short),qr_results)

# for a list of quantiles, the result is a list of the same dimension where 
# each element corresponds to each quantile as detailed above.  

# Model results: 

# On figures, beta 0 is the intercept, 
# beta 1 (b1) represents the coefficients for the fixed effects (yield) 
# so is the most important graph
# Beta 2 (b2) represents coefficients for the random effects (effect ID)
# sigma represents standard errors for all parameters
# [export figures manually with height=500, width=800, name Fig X quan reg...]

beta <- qmodel[5][[1]]$res$beta  #fixed effects
weights = qmodel[5][1]$res$weights  #random weights
nj = c(as.data.frame(table(groups))[,2]) #obs per subject
fixed = tcrossprod(x,t(beta))
random = rep(0,dim(x)[1])  #initializing random shift
group.plot(x=fixed,y=yi_B,groups=groups,type="l")

for (j in 1:length(nj)){ 
  z1=matrix(z[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=dim(z)[2])
  random[(sum(nj[1:j-1])+1):(sum(nj[1:j]))] = tcrossprod(z1,t(weights[j,]))
}

pred = fixed + random  #predictions
group.plot(yi_Y,pred,groups,type = "l",xlab="Yields",ylab="Biodiversity")
group.points(yi_Y,yi_B,groups)
#colnames(pred) = c("yi_B_pred")
#qmodel.pred <- data.frame(pred) %>% cbind(yi_Y)

col.quantiles <- c("0.1"= "grey85","0.2" = "grey80","0.3" =  "grey75","0.4"= "grey70","0.5"= "orange",
                   "0.6" =  "grey65","0.7"= "grey60","0.8" = "grey55","0.9" = "grey50","Mean"= "blue")

g <- ggplot(data,aes(x=yi_Y,y=yi_B))+geom_point(pch=1)+
  labs(x="Yield RR",y=paste0(b_measure," RR"))+
  geom_smooth(formula=y~x,method="lm",se=F,aes(colour="Mean"),size=0.8)+
  geom_text(aes(x=1.5,y=6, label=paste0("y = ",round(qmodel[[5]]$res$beta[1,1],2),round((slope=qmodel[[5]]$res$beta[2,1]),2),"x, p = ",round(qmodel[[5]]$res$table$`Pr(>|z|)`[2],4)),colour="0.5"),size=3.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[1]]$res$beta[1,1], slope=qmodel[[1]]$res$beta[2,1],colour="0.1"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[2]]$res$beta[1,1], slope=qmodel[[2]]$res$beta[2,1],colour="0.2"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[3]]$res$beta[1,1], slope=qmodel[[3]]$res$beta[2,1],colour="0.3"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[4]]$res$beta[1,1], slope=qmodel[[4]]$res$beta[2,1],colour="0.4"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[5]]$res$beta[1,1], slope=qmodel[[5]]$res$beta[2,1],colour="0.5"),size=0.8,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[6]]$res$beta[1,1], slope=qmodel[[6]]$res$beta[2,1],colour="0.6"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[7]]$res$beta[1,1], slope=qmodel[[7]]$res$beta[2,1],colour="0.7"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[8]]$res$beta[1,1], slope=qmodel[[8]]$res$beta[2,1],colour="0.8"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qmodel[[9]]$res$beta[1,1], slope=qmodel[[9]]$res$beta[2,1],colour="0.9"),size=0.5,show.legend=FALSE)+
  scale_colour_manual(name="Quantiles",values=col.quantiles)+
  theme_classic()+
  theme(text=element_text(size=10))
g


tiff(paste0(outpath,"Quantile regression curves ",b_measure_short,".tif"),width=6,height=4,units="in",res=300,compression="lzw")
g
dev.off()

g <- ggplot(qr_results,aes(x=quantile,y=Estimate,group=variable,ymin=`Inf.CI95.`,ymax=`Sup.CI95.`,se=`Std..Error`))+
  geom_hline(yintercept=0)+
  geom_point(pch=1)+
  #geom_path()+
  geom_smooth(colour="black")+
  labs(y="Estimate",x="Yield quantiles")+
  scale_x_continuous(breaks=seq(0.1,0.9,0.1))+
  facet_grid(cols=vars(variable))+
  theme_bw()+
  theme(panel.grid.minor=element_blank())
g

tiff(paste0(outpath,"Quantile regression point estimates ",b_measure_short,".tif"),width=6,height=4,units="in",res=300,compression="lzw")
g
dev.off()

detach(data)

write.xlsx(list(qr_results_Bio,
                  qr_results_Abun,
                  qr_results_Rich,
                  qr_results_RichEven),
           file=paste0(outpath,"qr results.xlsx"),overwrite=TRUE)
write.xlsx(list(qr_results_Bio_HQ),
           file=paste0(outpath,"qr results_HQ.xlsx"),overwrite=TRUE)

### Make one figure with regression curves for all metrics ####

qr_results_Bio <- read.xlsx(paste0(outpath,"qr results.xlsx"),sheet=1)
qr_results_Abun <- read.xlsx(paste0(outpath,"qr results.xlsx"),sheet=2)
qr_results_Rich <- read.xlsx(paste0(outpath,"qr results.xlsx"),sheet=3)
qr_results_RichEven <- read.xlsx(paste0(outpath,"qr results.xlsx"),sheet=4)
qr_results_Bio_HQ <- read.xlsx(paste0(outpath,"qr results_HQ.xlsx"),sheet=4)

qr_results_all <- qr_results_Bio %>% mutate(B_measure = "Biodiversity") %>%
  rbind(
    qr_results_Abun %>% mutate(B_measure = "Abundance")) %>%
  rbind(
    qr_results_Rich %>% mutate(B_measure = "Richness")) %>%
  rbind(
    qr_results_RichEven %>% mutate(B_measure = "Richness-Evenness")) 

g <- ggplot(qr_results_all,aes(x=quantile,y=Estimate,group=variable,ymin=`Inf.CI95.`,ymax=`Sup.CI95.`,se=`Std..Error`))+
  geom_hline(yintercept=0)+
  geom_point(pch=1)+
  #geom_path()+
  geom_smooth(colour="black",size=0.5)+
  labs(y="Estimate",x="Yield quantiles")+
  scale_x_continuous(breaks=seq(0.1,0.9,0.4))+
  facet_grid(rows=vars(factor(B_measure,levels=c("Biodiversity","Abundance","Richness","Richness-Evenness","Biodiversity (HQ)"))),
             cols=vars(variable))+
  theme_bw()+
  theme(text=element_text(size=10,colour="black"),
        axis.title=element_text(size=10,colour="black"),
        axis.text=element_text(size=10,colour="black"),
        panel.grid.minor=element_blank(),
        panel.spacing=unit(0.7,"lines"),
        strip.background=element_blank(),
        strip.text=element_text(face="bold",size=10))
g

gb <- ggplot_build(g)
lay <- gb$layout$layout %>% mutate(B_measure = `factor(...)`)
tags <- cbind(lay, label = paste0(LETTERS[lay$PANEL]), x = -Inf, y = Inf)
tags <- cbind(lay,label=c("Ai","Aii", "Bi","Bii","Ci","Cii","Di","Dii"),x=-Inf,y=Inf)
g <- g + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), hjust = -0.5, 
              vjust = 1.5, fontface = "bold",size=3.5, inherit.aes = FALSE) 
g

tiff(paste0(outpath,"Quantile regression point estimates all.tif"),width=7,height=6.5,units="in",res=500,compression="lzw")
g
dev.off()

# Plot quantile regression slopes using exported data ####

qr_results <- qr_results_Bio %>% mutate(B_measure = "Biodiversity") %>%
  mutate(p.value = ifelse(`Pr...z..`<0.001,"< 0.001",paste0(" = ",round(`Pr...z..`,3))))

qr_results <- qr_results_Abun %>% mutate(B_measure = "Abundance") %>%
  mutate(p.value = ifelse(`Pr...z..`<0.001,"< 0.001",paste0(" = ",round(`Pr...z..`,3))))

qr_results <- qr_results_Rich %>% mutate(B_measure = "Richness") %>%
  mutate(p.value = ifelse(`Pr...z..`<0.001,"< 0.001",paste0(" = ",round(`Pr...z..`,3))))

qr_results <- qr_results_RichEven %>% mutate(B_measure = "Richness-Evenness") %>%
  mutate(p.value = ifelse(`Pr...z..`<0.001,"< 0.001",paste0(" = ",round(`Pr...z..`,3))))

qr_results <- qr_results_Bio_HQ %>% mutate(B_measure = "Biodiversity (HQ)") %>%
  mutate(p.value = ifelse(`Pr...z..`<0.001,"< 0.001",paste0(" = ",round(`Pr...z..`,3))))

g <- ggplot(d,aes(x=yi_Y,y=yi_B))+geom_point(pch=1)+
  labs(x="Yield RR",y=paste0(qr_results$B_measure," RR"))+
  geom_smooth(formula=y~x,method="lm",se=F,aes(colour="Mean"),size=0.8)+
  
  #geom_text(aes(x=1.3,y=8, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.1 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.1 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.1 & qr_results$variable =="beta 2"),]$p.value),colour="0.1"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=7, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.2 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.2 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.2 & qr_results$variable =="beta 2"),]$p.value),colour="0.2"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=6, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.3 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.3 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.3 & qr_results$variable =="beta 2"),]$p.value),colour="0.3"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=5, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.4 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.4 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.4 & qr_results$variable =="beta 2"),]$p.value),colour="0.4"),size=3.5,show.legend=FALSE)+
  geom_text(aes(x=1.1,y=8, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.5 & qr_results$variable =="beta 1"),]$Estimate,3)," + ",round((slope=qr_results[which(qr_results$quantile==0.5 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p ",qr_results[which(qr_results$quantile==0.5 & qr_results$variable =="beta 2"),]$p.value),colour="0.5"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=-6, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.6 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.6 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.6 & qr_results$variable =="beta 2"),]$p.value),colour="0.6"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=-7, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.7 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.7 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.7 & qr_results$variable =="beta 2"),]$p.value),colour="0.7"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=-8, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.8 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.8 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.8 & qr_results$variable =="beta 2"),]$p.value),colour="0.8"),size=3.5,show.legend=FALSE)+
  #geom_text(aes(x=1.3,y=-9, label=paste0("y = ",round(qr_results[which(qr_results$quantile==0.9 & qr_results$variable =="beta 1"),]$Estimate,2)," + ",round((slope=qr_results[which(qr_results$quantile==0.9 & qr_results$variable =="beta 2"),]$Estimate),2),"x, p = ",qr_results[which(qr_results$quantile==0.9 & qr_results$variable =="beta 2"),]$p.value),colour="0.9"),size=3.5,show.legend=FALSE)+
  
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.1 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.1 & qr_results$variable =="beta 2"),]$Estimate,colour="0.1"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.2 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.2 & qr_results$variable =="beta 2"),]$Estimate,colour="0.2"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.3 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.3 & qr_results$variable =="beta 2"),]$Estimate,colour="0.3"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.4 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.4 & qr_results$variable =="beta 2"),]$Estimate,colour="0.4"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.5 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.5 & qr_results$variable =="beta 2"),]$Estimate,colour="0.5"),size=0.8,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.6 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.6 & qr_results$variable =="beta 2"),]$Estimate,colour="0.6"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.7 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.7 & qr_results$variable =="beta 2"),]$Estimate,colour="0.7"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.8 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.8 & qr_results$variable =="beta 2"),]$Estimate,colour="0.8"),size=0.5,show.legend=FALSE)+
  geom_abline(aes(intercept=qr_results[which(qr_results$quantile==0.9 & qr_results$variable =="beta 1"),]$Estimate, slope=qr_results[which(qr_results$quantile==0.9 & qr_results$variable =="beta 2"),]$Estimate,colour="0.9"),size=0.5,show.legend=FALSE)+
  
  scale_colour_manual(name="Quantiles",values=col.quantiles)+
  theme_classic()+
  theme(text=element_text(size=10))
g

tiff(paste0(outpath,"Quantile regression point estimates ",qr_results$B_measure,".tif"),width=6,height=3.5,units="in",res=300,compression="lzw")
g
dev.off()

