### Run quantile regression ####
# Authors: SJ

library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(readxl)
library(openxlsx)
library(quantreg)
library(qrLMM)

wd <- readline() #at the prompt, copy and paste your filepath and press enter
D:\02_Bioversity\30_SustainableFoods\R
setwd(wd)

outpath <- "./Results_tradeoffs_wLER/"

d <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield")

# Relationship between yield and biodiversity ###

# Quantile regression is a type of linear regression 
# to understand relationship between one or more predictors and a response variable.
# Linear regression is normally used to estimate the mean value of a response variable, but in 
# quantile regression we estimate any percentile of the response value, e.g. 70th, 90th.
# https://www.statology.org/quantile-regression-in-r/ 
# 50th percentile (median) is robust to outliers so can be better than linear regression (based on mean)

ggplot(d,aes(y=yi_Y,x=System_T))+geom_boxplot()
ggplot(d,aes(y=yi_Y,x=Lat_T))+geom_point()+geom_smooth()
ggplot(d,aes(y=yi_Y,x=Long_T))+geom_point()+geom_smooth()
ggplot(d,aes(y=yi_Y,x=Agrochem_CT))+geom_boxplot()

ggplot(d,aes(y=yi_Y,x=yi_B))+geom_point()+geom_smooth()

# run quantile regression without multi-level model (so don't use this, just for checking as is fast)
qmodel <- rq(yi_Y~yi_B,data=d,tau=0.9)
summary(qmodel)
qmodel <- rq(yi_Y~yi_B,data=d,tau=0.1)
summary(qmodel)
qmodel <- rq(yi_Y~yi_B,data=d,tau=0.5)
summary(qmodel)

# interpreting results when dependent and independent variables are logarithms:
# we interpret the coefficient as the expected change in percent in the dependent variable 
# when the independent variable is increased by 1%
# e.g. an increase in biodiversity of 1% is associated with an increase in yield by ##%.
# https://www.stathelp.se/en/regression_logarithm_en.html

# Including all intercropping data points:
# 90th percentile of yield outcome  = exp(0.834)  -0.0309*biodiversity outcome - WEAK NEGATIVE
# 50th percentile of yield outcome  = exp(0.095)  -0.038*biodiversity outcome -  WEAK NEGATIVE 
# 10th percentile of yield outcome  = exp(-0.044)  -0.0086*biodiversity outcome -  WEAK NEGATIVE 
# Means using the mean we find no relationship
# using the median we find a very weak positive relationship

# Including only LER intercropping data points:
# 90th percentile of yield outcome  = exp(0.581)  -0.06725*biodiversity outcome -  WEAK NEGATIVE 
# 50th percentile of yield outcome  = exp(0.465)  -0.0114*biodiversity outcome -  WEAK NEGATIVE
# 10th percentile of yield outcome  = exp(-0.053)  -0.027*biodiversity outcome -  WEAK NEGATIVE 
# Means using the mean we find no relationship
# using the median we find a very weak positive relationship

ggplot(d,aes(yi_Y,yi_B))+geom_point()+
  geom_abline(intercept=exp(coef(qmodel)[1]),slope=coef(qmodel)[2])+
  geom_smooth(method="lm",se=F)+theme_bw()

# QRLMM package allows for mixed effects modelling with quantile regression
# https://rdrr.io/cran/qrLMM/man/QRLMM.html
# paper with results from this package: https://www.researchgate.net/publication/303811593_Quantile_regression_for_mixed-effects_models/figures

d <- d %>% mutate(ID = as.factor(ID),Effect_ID=as.factor(Effect_ID))
attach(d)

y=yi_B # response
x=cbind(1,yi_Y) # design matrix for fixed effects
z=cbind(1,Effect_ID) # design matrix for random effects
#groups <- rep(1,nrow(d)) 
groups=ID
qr_values <- seq(0.1,0.9,0.1)

qmodel <- QRLMM(y=y,x=x,z=z,groups=groups,p=qr_values,show.convergence=TRUE,MaxIter=200,M=10)
qmodel
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
  qr_results_0 <- data.frame(qmodel[i][[1]]$res$table) %>% mutate(quantile = paste0("0.",i)) %>%
    mutate(quantile = as.numeric(quantile)) %>%
    mutate(variable = rownames(qr_results_0)) %>%
    mutate(processing_time = qmodel[i][[1]]$res$time) %>%
    mutate(AIC = qmodel[i][[1]]$res$AIC, BIC = qmodel[i][[1]]$res$BIC,
           loglik = qmodel[i][[1]]$res$loglik,sigma = qmodel[i][[1]]$res$sigma) 
  qr_results <- rbind(qr_results,qr_results_0)}

write.xlsx(qr_results,file=paste0(outpath,"qr results.xlsx"),overwrite=TRUE)

# for a list of quantiles, the result is a list of the same dimension where 
# each element corresponds to each quantile as detailed above.  

# Model results: 

# On figures, beta 1 (b1) represents the coefficients for the fixed effects (biodiversity response) 
# so is the most important graph
# Beta 2 (b2) represents coefficients for the random effects (effect ID)
# sigma represents standard errors for all parameters
# [export figures manually with height=500, width=800, name Fig X quan reg...]

beta <- qmodel[5][1]$res$beta  #fixed effects
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
g <- ggplot(d,aes(x=yi_Y,y=yi_B))+geom_point(pch=1)+
  labs(x="Yield RR",y="Biodiversity RR")+
  geom_smooth(formula=y~x,method="lm",se=F,aes(colour="Mean"),size=0.8)+
  geom_abline(aes(intercept=qmodel[[1]]$res$beta[1,1], slope=qmodel[[1]]$res$beta[2,1],colour="0.1"),size=0.5,show.legend=FALSE)+
  #geom_text(aes(x=-3,y=qmodel[[1]]$res$beta[1,1], label=c("0.1")))+
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

tiff(paste0(outpath,"Quantile regression curves.tif"),width=6,height=4,units="in",res=300,compression="lzw")
g
dev.off()

g <- ggplot(qr_results,aes(x=quantile,y=Estimate,group=variable,ymin=`Inf.CI95.`,ymax=`Sup.CI95.`,se=`Std..Error`))+
  geom_hline(yintercept=0)+
  geom_point(pch=1)+
  #geom_path()+
  geom_smooth(colour="black")+
  labs(y="Estimate",x="Quantiles")+
  scale_x_continuous(breaks=seq(0.1,0.9,0.1))+
  facet_grid(~variable)+
  theme_bw()+
  theme(panel.grid.minor=element_blank())
g

tiff(paste0(outpath,"Quantile regression point estimates.tif"),width=6,height=4,units="in",res=300,compression="lzw")
g
dev.off()

detach(d)
