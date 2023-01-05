#### Multinomial logit models to compute relative risk ratios ####
#### associated with biodiversity and yield outcomes and      ####
#### random forest to assess predictor variable importance    ####

# Author: SJ
# Last updated 12/12/2022

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
library(AICcmodavg)
library(ggallin)
library(metR)
library(ggnewscale)  
library(jmv)# Descriptive stats
library(metaforest)
library(randomForest)
library(randomForestExplainer)
library(ggrepel)

#### Set file paths ####

wd <- readline() #at the prompt, copy and paste your filepath and press enter

setwd(wd)

outpath <- "./Results_tradeoffs_wLER/"

#### Set parameters ####
col.intervention = c("navy","forestgreen", "seagreen","gold","orange", "lightblue","purple")
col.synergies = c("Lose B-Lose Y" = "navy","Other"="grey50", "Win B-Win Y"= "forestgreen", "Win B-Lose Y" ="gold","Lose B-Win Y"= "lightblue")
col.synergies_sig_simp = c("navy", "lightblue" ,"grey70" ,"gold", "forestgreen") # blue to green with grey
col7 <- c("#00A600", "#63C600" ,"#E6E600", "#E9BD3A", "#ECB176", "#EFC2B3" ,"mediumpurple")
col2 <- c("#00A600","#E6E600")

#### Import data ####
d <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield") # full dataset
d_validity.high <- read.xlsx(paste0(outpath,"Bio_yields_effects.xlsx"),sheet="d_bioyield_validity") # reduced dataset with only high quality data points

#### Format data for analysis ####
reclassify <- function(data){
  
  # exclude simplified controls that are not either monocultures or have an absence of embedded natural vegetation
  data <- data %>% mutate(System_C = ifelse(System_C =="Simplified other" & System_T %in% c("Associated plant species","Agroforestry"),"Exclude",System_C)) %>%
    filter(System_C !="Exclude")
  
  # categorize biodiversity metrics
  table(data$B_measure)
  data <- data %>% mutate(B_measure_group = ifelse(B_measure %in% c("Abundance", "Activity-density", "Area under disease progress curve (AUDPC)",
                                                                    "Colonization Percent","Colonization Percent ", "Vistitation frequency"),"Abundance",
                                                   ifelse(B_measure %in% c("Chao1 Index","Species Richness"),"Richness",
                                                          ifelse(B_measure %in% c("Pielou Index","Shannon-Wiener Index","Shannon Index","Shannon Index ","Species Eveness","Species Eveness ","Simpson Index","Simpson Index "),"Richness-Evenness",B_measure))))
  
  # categorize biodiversity taxa
  data <- data %>% mutate(Taxa_group_simp = ifelse(Taxa_group %in% c("Annelids & millipedes","Arachnids","Arthropods","Coleoptera","Hymenoptera","Insect (other)","Lepidoptera"),"Invertebrates",ifelse(Taxa_group %in% c("Birds","Mammals"),"Vertebrates",as.character(Taxa_group))))
  
  # categorize biomes
  data <- data %>% mutate(Biome = as.character(Biome)) %>% 
    mutate(Biome = ifelse(Biome=="Tropical & Subtropical Non-forest", "Tropical & Subtropical Grasslands",
                          ifelse(Biome == "Temperate Non-forest", "Temperate Grasslands",Biome)))
  data <- data %>% mutate(Biome_simp = as.character(Biome_simp)) %>% 
    mutate(Biome_simp = ifelse(Biome_simp=="Tropical & Subtropical Non-forest", "Tropical & Subtropical Grasslands",
                               ifelse(Biome_simp == "Temperate Non-forest", "Temperate Grasslands",
                                      ifelse(Biome_simp=="Other","Boreal, Montane & Deserts",Biome_simp))))
  }

d <- reclassify(d)
d_validity.high <- reclassify(d_validity.high)

set_synergies <- function(data){
  # identify and order pairwise biodiversity-yield outcome categories
  # based on direction 
  data <- data %>%
    mutate(Synergies_withEquiv = ifelse(yi_B >0 & yi_Y >0,"Win B-Win Y",
                                      ifelse(yi_B>0 & yi_Y<0,"Win B-Lose Y",
                                             ifelse(yi_B>0 & yi_Y==0,"Win B-Equiv Y",
                                                    ifelse(yi_B<0 & yi_Y >0,"Lose B-Win Y",
                                                           ifelse(yi_B<0 & yi_Y == 0,"Lose B-Equiv Y",
                                                                  ifelse(yi_B==0 & yi_Y>0,"Equiv B-Win Y",
                                                                         ifelse(yi_B==0 & yi_Y<0,"Equiv B-Lose Y",
                                                                                ifelse(yi_B==0 & yi_Y==0,"Equiv B-Equiv Y",
                                                                                       ifelse(yi_B<0 & yi_Y<0,"Lose B-Lose Y","check"))))))))))

  data <- data %>%
    mutate(Synergies = ifelse(as.character(Synergies_withEquiv) %in% c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"),Synergies_withEquiv,"Other"))
  data$Synergies <- factor(data$Synergies,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y","Other"))
  
  data$Synergies_withEquiv <- fct_relevel(data$Synergies_withEquiv,"Lose B-Lose Y")

  # based on direction and significance
  data <- data %>% mutate(Synergies_sig = ifelse(as.character(Synergies_sig) %in% c("Win B-Equiv Y"  , "Equiv B-Win Y","Lose B-Equiv Y" , "Equiv B-Lose Y","Equiv B-Equiv Y" ),"Other",as.character(Synergies_sig)))
  data$Synergies_sig <- factor(data$Synergies_sig,levels=c("Lose B-Lose Y","Lose B-Win Y","Win B-Lose Y","Win B-Win Y","Other"))
  
  data <- data %>% mutate(ID = as.numeric(as.character(ID)),Effect_ID = as.numeric(as.character(Effect_ID)))
  return(data)
}

d <- set_synergies(d)
d_validity.high <- set_synergies(d_validity.high)

table(d$Synergies)
table(d$Synergies_withEquiv)
table(d$Synergies_sig)

d_synergies_sig <- d %>% filter(Synergies_sig %in% c( "Lose B-Lose Y",  "Lose B-Win Y","Win B-Lose Y" ,"Win B-Win Y"))
d_synergies_sig$Synergies_sig <- factor(d_synergies_sig$Synergies_sig,levels=c( "Lose B-Lose Y" ,"Lose B-Win Y" ,  "Win B-Lose Y"  ,  "Win B-Win Y"))

#### Set reference levels ####
# to most frequently occurring 

set_ref <- function(data){
  data <- data %>% mutate(Taxa_group_simp = ifelse(Taxa_group %in% c("Annelids & millipedes","Arachnids","Arthropods","Coleoptera","Hymenoptera","Insect (other)","Lepidoptera"),"Invertebrates",ifelse(Taxa_group %in% c("Birds","Mammals"),"Vertebrates",as.character(Taxa_group))))
  table(data$Taxa_group_simp)
  data$Taxa_group_simp <- fct_relevel(data$Taxa_group_simp,"Invertebrates")
  
  table(data$Biome_simp,data$Synergies) # complete
  table(data$Biome_simp,data$Crop_type_C) # gaps in all except annual herb
  table(data$Biome_simp,data$System_T) # complete for associated plant species only (nearly for embedded natural and intercropping)
  
  data <- data %>% mutate(Biome = factor(Biome))
  table(data$Biome)
  table(data$Biome,data$Synergies)
  data$Biome <- fct_relevel(data$Biome,"Tropical & Subtropical Grasslands")
  #data$Biome <- fct_relevel(data$Biome,"Temperate Forests")

  data <- data %>% mutate(Biome_simp = factor(Biome_simp))
  table(data$Biome_simp)
  table(data$Biome_simp,data$Synergies)
  data$Biome_simp <- fct_relevel(data$Biome_simp,"Tropical & Subtropical Grasslands")
  data$Biome_simp <- fct_relevel(data$Biome_simp,"Boreal, Montane & Deserts",after=Inf)

  table(data$Crop_FAO_C,data$Synergies)
  table(data$Crop_FAO_C)
  data$Crop_FAO_C <- fct_relevel(data$Crop_FAO_C,"Cereals")

  table(data$Crop_type_C,data$Synergies)
  table(data$Crop_type_C)
  data$Crop_type_C <- fct_relevel(data$Crop_type_C,"Annual Herb")

  table(data$System_T,data$Synergies)
  table(data$System_T)
  data$System_T <- fct_relevel(data$System_T,"Associated plant species")

  table(data$Agrochem_CT,data$Synergies)
  data <- data  %>%  mutate(Agrochem_CT = factor(Agrochem_CT,levels=c("Agrochemicals","No agrochemicals","Mixed","No data")))

  table(data$Taxa_group)
  table(data$Taxa_group_simp)
  data$Taxa_group <- fct_relevel(data$Taxa_group,"Insect (other)")

  data$Taxa_group_simp <- fct_relevel(data$Taxa_group_simp,"Invertebrates")
  
  table(data$B_measure_group,data$Synergies)
  data$B_measure_group <- fct_relevel(data$B_measure_group,"Abundance")

  table(data$Pest_group,data$Synergies)
  data$Pest_group <- fct_relevel(data$Pest_group,"Other")

  table(data$B_ground,data$Synergies)
  data$B_ground <- fct_relevel(data$B_ground,"Above")

  table(data$Yield_measure_group,data$Synergies)
  data$Yield_measure_group <- fct_relevel(data$Yield_measure_group,"Mass per area")
  
  return(data)
}

d <- set_ref(d)
d_validity.high <- set_ref(d_validity.high)
set_ref(d_synergies_sig)

#### Split dataset by biodiversity metric ####

d_abun <- d %>% filter(B_measure_group =="Abundance")
d_rich <- d %>% filter(B_measure_group =="Richness")
d_even <- d %>% filter(B_measure_group =="Richness-Evenness")

#d_sig_abun <- d_synergies_sig %>% filter(B_measure_group =="Abundance")
#d_sig_rich <- d_synergies_sig %>% filter(B_measure_group =="Richness")
#d_sig_even <- d_synergies_sig %>% filter(B_measure_group =="Richness-Evenness")

#### Random forest to identify important moderators ####

data <- d %>% select(Synergies,System_T,Crop_type_C,Crop_FAO_C,Agrochem_CT, Biome_simp,Region.Name,Lat_T, Yield_measure_group,B_ground, B_measure_group,Taxa_group_simp,Pest_group)
set.seed(10)

rf <- randomForest(Synergies~.,data=data,
                   type="classification",importance=TRUE,importanceSD=TRUE, localImp=TRUE,ntree=10000)

rf.imp <- data.frame(importance(rf)) %>% setDT(keep.rownames=TRUE) %>% rename(Variable=rn) 

rf.imp <- rf.imp %>% 
  mutate(Overall_accuracy = (MeanDecreaseAccuracy-min(MeanDecreaseAccuracy))/(max(MeanDecreaseAccuracy-min(MeanDecreaseAccuracy)))) %>%
  mutate(Variable_labels = ifelse(Variable =="Crop_FAO_C","Crop commodity",
                                  ifelse(Variable == "System_T","Practice",
                                         ifelse(Variable=="Lat_T","Latitude",
                                                ifelse(Variable=="Biome_simp","Biome",
                                                       ifelse(Variable=="Pest_group","Pest group",
                                                              ifelse(Variable == "Agrochem_CT","Agrochemical use",
                                                                     ifelse(Variable =="Region.Name","Continent",
                                                                            ifelse(Variable == "Crop_type_C","Crop type",
                                                                                   ifelse(Variable =="Yield_measure_group","Yield metric",
                                                                                          ifelse(Variable=="Taxa_group_simp","Taxonomic group",
                                                                                                 ifelse(Variable=="B_measure_group","Biodiversity metric",
                                                                                                        ifelse(Variable=="B_ground","Ground relation",Variable)))))))))))))

ggplot(rf.imp,aes(x=reorder(Variable_labels,Overall_accuracy),y=Overall_accuracy))+geom_col(fill="grey50")+
  coord_flip(expand=c(0))+
  labs(x="",y="Relative importance")+
  theme_pubr()+theme(text= element_text(size=10),plot.margin = unit(c(0,1,0,0),"lines"))

ggsave(paste0(outpath,"Fig random forest importance.tif"),device="tiff",width=6,height=4.5,units="in")

explain_forest(rf,path=paste0(outpath,"randomforest_explained.html"))#,interactions=TRUE,data=data)
#forest_interactions <- min_depth_interactions(rf) # long to run
#plot_min_depth_interactions(forest_interactions)
rf_imp_explainer <- measure_importance(rf)


#### Descriptive stats ####

# check associations
test <- prop.table(table(as.character(d$Biome_simp),as.character(d$System_T),as.character(d$Synergies)))
test <- xtabs(~Biome_simp+System_T+Synergies, data=d)
ftable(test)
summary(test)
assocstats(xtabs(~Biome_simp+Region.Name,data=d))
assocstats(xtabs(~Biome_simp+System_T,data=d))
assocstats(xtabs(~Crop_FAO_C+Crop_type_C,data=d))
assocstats(xtabs(~System_T+Yield_measure_group,data=d))
assocstats(xtabs(~Crop_FAO_C+Yield_measure_group,data=d))
assocstats(xtabs(~Yield_measure_group+Crop_type_C,data=d))
assocstats(xtabs(~Crop_FAO_C+System_T,data=d))
assocstats(xtabs(~System_T+Crop_type_C,data=d))
assocstats(xtabs(~Biome_simp+B_measure_group,data=d))
assocstats(xtabs(~Taxa_group_simp+B_measure_group,data=d))
mosaic(xtabs(~Biome_simp+System_T,data=d))
assocplot(xtabs(~Biome_simp+System_T,data=d))
assocplot(xtabs(~Crop_FAO_C+Crop_type_C,data=d))

# tabulate frequencies
freq(d$Biome_simp)
freq(d$Region.Name)
table(d$Region.Name,d$Country)
table(d$Region.Name,d$Biome_simp)

ggplot(d,aes(y=System_T))+geom_bar(aes(fill=Biome_simp))+ scale_fill_manual(values=hcl.colors(7, palette = "Viridis"),name="")+
  #stat_count(geom = "text", colour = "grey70", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5))+
  stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),vjust=0.5,hjust=0)+#,position=position_stack(vjust=0.5))+
  #geom_text(colour = "black", size = 3.5,aes(label = paste0(round(..count../sum(..count..)*100,1),"%")),stat='count', nudge_x=30)+
  labs(x="Number of cases",y="")+scale_x_continuous(expand = expansion(add = c(0, 0)))+
  coord_cartesian(clip="off")+
  theme_pubr()+
  facet_wrap(~Region.Name,scales="free_x",nrow=1)+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        axis.text=element_text(size=10),axis.title=element_text(size=10),
        panel.spacing.x=unit(1.5,units="lines"),#strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggplot(d,aes(y=System_T))+geom_bar(aes(fill=Crop_FAO_C))+ scale_fill_manual(values=hcl.colors(8, palette = "Viridis"),name="")+
  stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),vjust=0.5,hjust=0)+#,position=position_stack(vjust=0.5))+
  labs(x="Number of cases",y="")+scale_x_continuous(expand = expansion(add = c(0, 0)))+
  coord_cartesian(clip="off")+
  theme_pubr()+
  facet_wrap(~Agrochem_CT,scales="free_x",nrow=1)+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        axis.text=element_text(size=10),axis.title=element_text(size=10),
        panel.spacing.x=unit(1.5,units="lines"),#strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))


ggplot(d,aes(y=Biome_simp))+geom_bar(aes(fill=System_T))+ scale_fill_manual(values=hcl.colors(8, palette = "Viridis"),name="")+
  stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),vjust=0.5,hjust=0)+#,position=position_stack(vjust=0.5))+
  labs(x="Number of cases",y="")+scale_x_continuous(expand = expansion(add = c(0, 0)))+
  coord_cartesian(clip="off")+
  theme_pubr()+
  facet_wrap(~Synergies_withEquiv,scales="free_x",nrow=1)+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        axis.text=element_text(size=10),axis.title=element_text(size=10),
        panel.spacing.x=unit(1.5,units="lines"),#strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

Synergies_n_cases <- d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n()) %>% select(c("Synergies","n_cases")) %>% unique()

ggplot(d %>%setDT() %>% group_by_at(c("System_T","Biome_simp","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("System_T","Biome_simp","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("System_T","Biome_simp","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("System_T","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Biome_simp = "Sum") %>% select(c("System_T","Biome_simp","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Biome_simp","Synergies")) %>% mutate(n_cases = n()) %>% mutate(System_T = "Sum") %>% select(c("System_T","Biome_simp","Synergies","n_cases"))) %>% unique() %>% rbind(
           d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),System_T="Sum",Biome_simp="Sum") %>% select(c("System_T","Biome_simp","Synergies","n_cases")) %>% unique()) %>%
         mutate(System_T = factor(System_T),
                System_T = fct_relevel(System_T,"Sum"),
                Biome_simp = factor(Biome_simp),
                Biome_simp = fct_relevel(Biome_simp,"Sum",after=Inf)) %>%
  mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(System_T=="Sum" | Biome_simp=="Sum",1,0)),
aes(y=System_T,x=Biome_simp))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 biome x practice.tif"),device="tiff",width=10,height=6,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("Crop_type_C","Biome_simp","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("Crop_type_C","Biome_simp","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("Crop_type_C","Biome_simp","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("Crop_type_C","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Biome_simp = "Sum") %>% select(c("Crop_type_C","Biome_simp","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Biome_simp","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Crop_type_C = "Sum") %>% select(c("Crop_type_C","Biome_simp","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),Crop_type_C="Sum",Biome_simp="Sum") %>% select(c("Crop_type_C","Biome_simp","Synergies","n_cases")) %>% unique()) %>%
         mutate(Crop_type_C = factor(Crop_type_C),
                Crop_type_C = fct_relevel(Crop_type_C,"Sum"),
                Biome_simp = factor(Biome_simp),
                Biome_simp = fct_relevel(Biome_simp,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(Crop_type_C=="Sum" | Biome_simp=="Sum",1,0)),
       aes(y=Crop_type_C,x=Biome_simp))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 biome x crop type.tif"),device="tiff",width=10,height=6,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("System_T","Region.Name","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("System_T","Region.Name","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("System_T","Region.Name","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("System_T","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Region.Name = "Sum") %>% select(c("System_T","Region.Name","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Region.Name","Synergies")) %>% mutate(n_cases = n()) %>% mutate(System_T = "Sum") %>% select(c("System_T","Region.Name","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),System_T="Sum",Region.Name="Sum") %>% select(c("System_T","Region.Name","Synergies","n_cases")) %>% unique()) %>%
         mutate(System_T = factor(System_T),
                System_T = fct_relevel(System_T,"Sum"),
                Region.Name = factor(Region.Name),
                Region.Name = fct_relevel(Region.Name,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(System_T=="Sum" | Region.Name=="Sum",1,0)),
       aes(y=System_T,x=Region.Name))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 practice x region.tif"),device="tiff",width=10,height=5,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("Crop_FAO_C","Agrochem_CT","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("Crop_FAO_C","Agrochem_CT","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("Crop_FAO_C","Agrochem_CT","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("Crop_FAO_C","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Agrochem_CT = "Sum") %>% select(c("Crop_FAO_C","Agrochem_CT","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Agrochem_CT","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Crop_FAO_C = "Sum") %>% select(c("Crop_FAO_C","Agrochem_CT","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),Crop_FAO_C="Sum",Agrochem_CT="Sum") %>% select(c("Crop_FAO_C","Agrochem_CT","Synergies","n_cases")) %>% unique()) %>%
         mutate(Crop_FAO_C = factor(Crop_FAO_C),
                Crop_FAO_C = fct_relevel(Crop_FAO_C,"Sum"),
                Agrochem_CT = factor(Agrochem_CT),
                Agrochem_CT = fct_relevel(Agrochem_CT,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(Crop_FAO_C=="Sum" | Agrochem_CT=="Sum",1,0)),
       aes(y=Crop_FAO_C,x=Agrochem_CT))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 crop x agrochem.tif"),device="tiff",width=10,height=5,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("Taxa_group_simp","Biome_simp","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("Taxa_group_simp","Biome_simp","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("Taxa_group_simp","Biome_simp","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("Taxa_group_simp","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Biome_simp = "Sum") %>% select(c("Taxa_group_simp","Biome_simp","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Biome_simp","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Taxa_group_simp = "Sum") %>% select(c("Taxa_group_simp","Biome_simp","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),Taxa_group_simp="Sum",Biome_simp="Sum") %>% select(c("Taxa_group_simp","Biome_simp","Synergies","n_cases")) %>% unique()) %>%
         mutate(Taxa_group_simp = factor(Taxa_group_simp),
                Taxa_group_simp = fct_relevel(Taxa_group_simp,"Sum"),
                Biome_simp = factor(Biome_simp),
                Biome_simp = fct_relevel(Biome_simp,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(Taxa_group_simp=="Sum" | Biome_simp=="Sum",1,0)),
       aes(y=Taxa_group_simp,x=Biome_simp))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 biome x taxa.tif"),device="tiff",width=10,height=6,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("Taxa_group_simp","B_measure_group","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("Taxa_group_simp","B_measure_group","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("Taxa_group_simp","B_measure_group","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("Taxa_group_simp","Synergies")) %>% mutate(n_cases = n()) %>% mutate(B_measure_group = "Sum") %>% select(c("Taxa_group_simp","B_measure_group","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("B_measure_group","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Taxa_group_simp = "Sum") %>% select(c("Taxa_group_simp","B_measure_group","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),Taxa_group_simp="Sum",B_measure_group="Sum") %>% select(c("Taxa_group_simp","B_measure_group","Synergies","n_cases")) %>% unique()) %>%
         mutate(Taxa_group_simp = factor(Taxa_group_simp),
                Taxa_group_simp = fct_relevel(Taxa_group_simp,"Sum"),
                B_measure_group = factor(B_measure_group),
                B_measure_group = fct_relevel(B_measure_group,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(Taxa_group_simp=="Sum" | B_measure_group=="Sum",1,0)),
       aes(y=Taxa_group_simp,x=B_measure_group))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 metric x taxa.tif"),device="tiff",width=10,height=5,units="in")

ggplot(d %>%setDT() %>% group_by_at(c("B_measure_group","Pest_group","Synergies","ID")) %>% mutate(n_studies = length(unique(ID))) %>% 
         group_by_at(c("B_measure_group","Pest_group","Synergies"))  %>% mutate(n_cases=n()) %>% select(c("B_measure_group","Pest_group","Synergies","n_cases")) %>% rbind(
           d %>% group_by_at(c("B_measure_group","Synergies")) %>% mutate(n_cases = n()) %>% mutate(Pest_group = "Sum") %>% select(c("B_measure_group","Pest_group","Synergies","n_cases"))) %>% rbind(
             d %>% group_by_at(c("Pest_group","Synergies")) %>% mutate(n_cases = n()) %>% mutate(B_measure_group = "Sum") %>% select(c("B_measure_group","Pest_group","Synergies","n_cases"))) %>% unique() %>% rbind(
               d %>% group_by_at(c("Synergies")) %>% mutate(n_cases = n(),B_measure_group="Sum",Pest_group="Sum") %>% select(c("B_measure_group","Pest_group","Synergies","n_cases")) %>% unique()) %>%
         mutate(B_measure_group = factor(B_measure_group),
                B_measure_group = fct_relevel(B_measure_group,"Sum"),
                Pest_group = factor(Pest_group),
                Pest_group = fct_relevel(Pest_group,"Sum",after=Inf)) %>%
         mutate(n_cases_cat = ifelse(n_cases>10,"1","0"),n_cases_type = ifelse(B_measure_group=="Sum" | Pest_group=="Sum",1,0)),
       aes(y=B_measure_group,x=Pest_group))+geom_tile(colour=NA,fill=NA,lwd=1.5)+
  geom_point(colour="white",fill="white",size=10)+
  geom_text(aes(label=n_cases,colour=n_cases_cat,fontface=ifelse(n_cases_type==1,"bold","plain")),size=3.5,show.legend=FALSE)+
  scale_colour_manual(values=c("grey50","black"),name="")+
  labs(x="",y="")+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~factor(Synergies,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other","Equiv B-Win Y","Equiv B-Lose Y")),
             nrow=1,scales="fixed")+
  theme(legend.position="bottom",legend.text=element_text(size=10),legend.margin=margin(c(0,0,0,-3),unit="lines"), 
        panel.grid.major = element_line(linetype=1,colour="grey"),
        axis.text=element_text(size=10,colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_text(size=10),axis.ticks=element_line(colour="black"),
        panel.spacing.x=unit(0,units="lines"),
        strip.text=element_text(colour="black"), strip.background=element_blank(),
        plot.margin=margin(c(0,1,0),unit="lines"))

ggsave(paste0(outpath,"Fig SX1 metric x pest group.tif"),device="tiff",width=11,height=5,units="in")

stats <- descriptives(d, vars=vars(Synergies, Biome_simp, System_T, Crop_type_C,Region.Name), freq = TRUE)
stats <- as.data.frame(do.call(rbind,stats$frequencies[[x]]))
stats
data.frame(stats$frequencies[[1]])

stats <- d %>% select(Synergies,System_T,Biome_simp,Agrochem_CT,Crop_type_C,Crop_FAO_C,Taxa_group_simp,B_measure_group,Yield_measure_group,Pest_group,B_ground) %>%
  setDT() %>%
  dcast.data.table(...~Synergies)
stats

stats <- d %>% select(Synergies,System_T,Biome_simp,Agrochem_CT,Crop_type_C,Crop_FAO_C,Taxa_group_simp,B_measure_group,Yield_measure_group,Pest_group,B_ground) %>%
  setDT() %>%
  melt.data.table(id.vars="Synergies")

data <- setDT(d)
n_group <- c("Synergies")

freq_tables <- function(data,n_group){
  n_data <- data.table(table(setDT(data)[,..n_group]))
  n_data <- n_data %>% dcast.data.table(...~Synergies,value.var="N")
  return(n_data)
}

for(x in c("System_T","Biome_simp","Agrochem_CT", "Crop_type_C","Crop_FAO_C","B_measure_group","Yield_measure_group",
           "Taxa_group_simp","Pest_group","B_ground")){
  x <- c("Synergies",x)
  output <- data.frame(freq_tables(d,x)) %>% 
    rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
           `Lose B-Win Y`=Lose.B.Win.Y,
           `Win B-Lose Y`=Win.B.Lose.Y,
           `Win B-Win Y`=Win.B.Win.Y) %>%
    mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
    select(-c(Other,N),c(Other,N))
  print(output)
  assign(paste0("freq_",x[2]),output,envir = .GlobalEnv)
}

freq_biome_crop_type_system <- d %>% select(Synergies,Biome_simp,Crop_type_C,System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_biome_crop_type_system

freq_reg_crop_type_system <- d %>% select(Synergies,Region.Name,Crop_type_C,System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_reg_crop_type_system

freq_biome_crop_type <- d %>% select(Synergies,Biome_simp,Crop_type_C,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_biome_crop_type

freq_reg_crop_type <- d %>% select(Synergies,Region.Name,Crop_type_C,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_reg_crop_type

freq_biome_system <- d %>% select(Synergies,Biome_simp, System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_biome_system

freq_reg_system <- d %>% select(Synergies,Region.Name,System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_reg_system

freq_crop_type_system <- d %>% select(Synergies,Crop_type_C,System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_crop_type_system

freq_crop_type_y_measure <- d %>% select(Synergies,Crop_type_C,Yield_measure_group,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_crop_type_y_measure

freq_biome_crop_fao_system <- d %>% select(Synergies,Biome_simp,Crop_FAO_C,System_T,Agrochem_CT) %>%
  setDT() %>%
  dcast.data.table(...~Synergies,length,drop=FALSE) %>%
  mutate(N = `Lose B-Lose Y` + `Lose B-Win Y` + `Win B-Lose Y` + `Win B-Win Y` +   `Other`) %>%
  data.frame() %>% 
  rename(`Lose B-Lose Y`=Lose.B.Lose.Y,
         `Lose B-Win Y`=Lose.B.Win.Y,
         `Win B-Lose Y`=Win.B.Lose.Y,
         `Win B-Win Y`=Win.B.Win.Y) %>%
  select(-c(Other,N),c(Other,N))
freq_biome_crop_fao_system

write.xlsx(list(freq_System_T,freq_Biome_simp,freq_Agrochem_CT,freq_Crop_type_C,freq_Crop_FAO_C,
                freq_B_measure_group,freq_Yield_measure_group,freq_Taxa_group_simp,freq_Pest_group,freq_B_ground,
                freq_biome_crop_type_system,freq_reg_crop_type_system, 
                freq_biome_crop_type,freq_reg_crop_type,freq_biome_system,freq_reg_system,
                freq_crop_type_system,freq_crop_type_y_measure,
                freq_biome_crop_fao_system),
           paste0(outpath,"Synergies freq tables.xlsx"),overwrite = TRUE)

#### Run multimomial models ####

# Notes on assumptions, use and limitations of multinomial models:
# Multinomial model is a type of GLM, 
# so the overall goodness-of-fit statistics and their interpretations and limitations apply.
# on multinomial logit models: https://www.princeton.edu/~otorres/LogitR101.pdf 
# on multinomial logit models in R: https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/
# on categorical predictors: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/
# example of paper using this approach: https://link.springer.com/article/10.1007/s10980-019-00775-1
# on how to interpret coefficients for categorical predictors: https://stats.stackexchange.com/questions/60817/significance-of-categorical-predictor-in-logistic-regression
# on multinomial models: Hilbe (2009)

# 4 categorical response outcomes, use a multinomial logit model
# reference category is LoseLose (first factor level)
# reset with: d$variable = relevel(d$variable, ref="optionA") or changing all factor levels

# Reasons and solutions for weird coef and p values:
# Devika et al 2016: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944325/ 
# Separation problem is caused by zero or very small cell sizes
# Solution is Penalised Maximum Likelihood
# https://search.r-project.org/CRAN/refmans/MultBiplotR/html/RidgeMultinomialLogisticFit.html 
# https://rdrr.io/cran/brglm2/man/brmultinom.html 
# https://cran.r-project.org/web/packages/brglm2/vignettes/multinomial.html 
# the above good for multiple outcomes, but don't accept random effects

# Sort model datasets to have adequate count in each cell formed by the factors in your analysis 
# All cell frequencies should be greater than 1 and 80% or more of cells should be greater than 5.
# The presence of small or empty cells may cause the logistic model to become unstable, 
# reporting implausibly large b coefficients and odds ratios for independent variables. 

# Calculating the additional variance explained by a predictor can be done from the log-likelihood ratio values
# see: https://data.princeton.edu/wws509/notes/c6s2 

table(d$Synergies) # fine
table(d$Synergies_sig) # fine
table(d$Synergies_withEquiv) # fine
table(d$Synergies,d$System_T) # 3 out of 7 have some zeros
table(d$Synergies,d$Crop_type_C) # complete, except two zeros for perennial herb
table(d$Synergies,d$Crop_FAO_C) # complete, except nuts/stimulants, roots tubers, and fibres
table(d$Synergies,d$Agrochem_CT) # complete
table(d$Synergies,d$Biome) # Temperate forests, tropical grasslands, tropical forests, temperate grasslands, mediterranean ok, the rest have at least one zero cell. 
table(d$Synergies,d$Biome_simp) # fine

# multi0: intercept only
# multi1: all predictors
# multi2: single predictors

multi0 <- multinom(Synergies ~ 1+  (1+Effect_ID|ID),data=d) 
summary(multi0)
multi0_result <- data.frame(broom::tidy(multi0),broom::glance(multi0)) # provides Wald statistic and p-value
multi0_result
multi0_result$rrr <- round(exp(multi0_result$estimate),3) 
multi0_result <- multi0_result %>% mutate(ci.lb = estimate - 1.96*std.error,ci.ub = estimate + 1.96*std.error)

multi0.rrr = exp(coef(multi0))
multi0.rrr # e.g. rrr = 1.28 for win B Lose Y which means this outcome is 1.28 times more likely than a lose-lose outcome (the reference category)
stargazer(multi0, type="text",out=paste(outpath,"multi0.txt"),coef=list(multi0.rrr),p.auto=FALSE) # p-values calculated with Wald tests

multi0_pp <- data.frame(Effect_ID=1,ID=1)
multi0_pp <- cbind(multi0_pp,predict(multi0,newdata=multi0_pp,type="probs"))
multi0_pp
summary(predict(multi0,newdata=multi0_pp,type="probs"))
multi0_pp$Outcome <- row.names(multi0_pp)
row.names(multi0_pp) <- NULL
colnames(multi0_pp) <- c("Effect_ID","ID","Probability","Outcome")
multi0_pp

# Intercept only models ####

model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=d_sig_other)
data=d_sig_other
run_name = "d_sig_other"
run_label="Only sig. - Other"
n_group="Synergies_sig"

model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=d)
data=d
run_name="d_synergies_sig"
run_label="Only significant"
n_group="Synergies_sig"
type="Z_Overall"

data = d
n_group = c("Synergies","Crop_type_C")
sig="Synergies"
cutoff_cases = 1
cutoff_studies = 1

# filter out levels with no lose-lose cases
n_studies_filter <- function(data,n_group,sig="Synergies",cutoff_cases=1,cutoff_studies=1){
 
  n_data <- data.table(table(setDT(data)[,..n_group]))

  colnames(n_data)[1] <- n_group[1]
  if(sig=="Synergies"){
    exclude <- n_data %>% filter(Synergies=="Lose B-Lose Y" & N ==0) %>% 
      dplyr::select(n_group[length(n_group)]) %>% 
      mutate(exclude="exclude")
    #exclude <- exclude[,-1]
  }
  if(sig=="Synergies_sig"){
    exclude <- n_data %>% filter(Synergies_sig=="Lose B-Lose Y" & N ==0) %>% 
      dplyr::select(n_group[length(n_group)]) %>% 
      mutate(exclude="exclude")
    #exclude <- exclude[,-1]
  }
  data <- data %>% left_join(exclude) %>% filter(is.na(exclude))
  
  # get count of unique studies and cases
  n_data <- data %>% 
    group_by_at(n_group) %>% summarise(n_cases=n(),n_studies=length(unique(ID)))

  data <- data %>% left_join(n_data,by=c(n_group)) 
  data <- data %>% filter(!(n_studies<cutoff_studies)) %>% filter(!(n_cases<cutoff_cases))
}

pp_multi0_function <- function(model,data,run_name,run_label,n_group,type){
  
  name <-paste0("multi0.", run_name) 
  assign(name,model,envir = .GlobalEnv)
  
  n_data <- data %>% 
    group_by_at(n_group) %>% summarise(n_cases=n(),n_studies=length(unique(ID)))

  multi0.var.result <- data.frame(broom::tidy(model),broom::glance(model)) 
  multi0.var.result <- multi0.var.result %>% rename(Outcome="y.level")
  multi0.var.result$rrr <- exp(multi0.var.result$estimate)
  multi0.var.result <- multi0.var.result %>%
    mutate(estimate = round(estimate,3),
           std.error = round(std.error,3),
           statistic  = round(statistic,3),
           p.value = round(p.value,4),
           rrr = round(rrr,4)) %>% 
    mutate(rrr_pc = 100*(rrr-1)) %>%
    mutate(ci.lb = estimate - 1.96*std.error,ci.ub = estimate + 1.96*std.error) %>%
    left_join(n_data,by=c("Outcome"=n_group))
  
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
  #assign(name,multi0.var.pp,envir = .GlobalEnv)
  
  multi0.var.result <- multi0.var.result %>% full_join(multi0.var.pp,by=c("Outcome")) %>%
    mutate(Probability = round(Probability,3))
  multi0.var.result <- multi0.var.result[,c(1,2,ncol(multi0.var.result),3:(ncol(multi0.var.result)-1))] %>%
    mutate(run = run_name) %>%
    mutate(run_label = paste0(run_label),
           run_label_n = paste0(run_label," (",n_cases,")")) %>%
    mutate(type = type) %>%
    mutate(mod = ifelse(substr(term,1,3)=="Sys", gsub("System_T*","",term),
                      ifelse(term=="(Intercept)","Reference",term))) %>%
    mutate(Label = ifelse(rrr>100 & p.value<0.001,paste0(">100,p<0.001"),
                          ifelse(rrr>100,paste0(">100,p=",round(p.value,3)),
                                 ifelse(p.value<0.001,paste0(round(rrr-1,2),",p<0.001"),paste0(round(rrr-1,2),",p=",round(p.value,3)))))) %>%
    mutate(Label_sig = ifelse(is.nan(p.value),"*",ifelse(p.value<0.05,"*",""))) %>%
    mutate(run_label_n = factor(run_label_n,levels=unique(run_label_n[order(type)]))) %>%
    filter(!is.na(term))
  
  name <-paste0("multi0.", run_name,".result") 
  assign(name,multi0.var.result,envir = .GlobalEnv)
}

#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d), data=d, run_name="d",run_label="All d",n_group="Synergies")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d), data=d,run_name="d_synergies", run_label="All biodiversity metrics",n_group="Synergies",type="Z_Overall")
pp_multi0_function(model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=d), data=d,"d_synergies_sig", run_label="All biodiversity metrics (Sig)",n_group="Synergies_sig",type="Z_Overall")

pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group="Synergies")), data=n_studies_filter(d_abun,n_group="Synergies"),run_name="d_abun", run_label="Abundance",n_group="Synergies",type="Abundance")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group="Synergies")), data=n_studies_filter(d_rich,n_group="Synergies"),run_name="d_rich", run_label="Richness",n_group="Synergies",type="B_Richness")
pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group="Synergies")), data=n_studies_filter(d_even,n_group="Synergies"),run_name="d_even", run_label="Richness-Evenness",n_group="Synergies",type="Richness-Evenness")

pp_multi0_function(model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group="Synergies_sig",sig="Synergies_sig")), data=n_studies_filter(d_abun,n_group="Synergies_sig",sig="Synergies_sig"),run_name="d_sig_abun", run_label="Abundance (Sig)",n_group="Synergies_sig",type="Abundance")
pp_multi0_function(model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group="Synergies_sig",sig="Synergies_sig")), data=n_studies_filter(d_rich,n_group="Synergies_sig",sig="Synergies_sig"),run_name="d_sig_rich", run_label="Richness (Sig)",n_group="Synergies_sig",type="B_Richness")
pp_multi0_function(model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group="Synergies_sig",sig="Synergies_sig")), data=n_studies_filter(d_even,n_group="Synergies_sig",sig="Synergies_sig"),run_name="d_sig_even", run_label="Richness-Evenness (Sig)",n_group="Synergies_sig",type="Richness-Evenness")

pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d_validity.high), data=d_validity.high,run_name="d_validity.high_synergies", run_label="All biodiversity metrics (HQ)",n_group="Synergies",type="Z_Overall")
pp_multi0_function(model=nnet::multinom(Synergies_sig ~ 1 + (1 + Effect_ID|ID),data=d_validity.high), data=d_validity.high,"d_validity.high_synergies_sig", run_label="All biodiversity metrics (HQ,Sig)",n_group="Synergies_sig",type="ZZ_Overall")

#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group="Synergies",cutoff_cases=5,cutoff_studies=1)), data=n_studies_filter(d,n_group="Synergies",cutoff_cases=5,cutoff_studies=1),run_name="d_synergies_cases5", run_label="All biodiversity metrics (>5 cases)",n_group="Synergies",type="ZZ_Overall")
#pp_multi0_function(model=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group="Synergies",cutoff_cases=5,cutoff_studies=5)), data=n_studies_filter(d,n_group="Synergies",cutoff_cases=5,cutoff_studies=5),run_name="d_synergies_cases5_studies5", run_label="All biodiversity metrics (>5 cases,>5 studies)",n_group="Synergies",type="ZZ_Overall")

multi0.results <- rbind(multi0.d_synergies.result,
                        multi0.d_abun.result,multi0.d_rich.result,multi0.d_even.result,
                        multi0.d_synergies_sig.result,
                        multi0.d_sig_abun.result,multi0.d_sig_rich.result,multi0.d_sig_even.result,
                        multi0.d_validity.high_synergies.result,multi0.d_validity.high_synergies_sig.result)

write.xlsx(list(multi0=multi0.results),paste0(outpath,"Multimodel multi0 results.xlsx"),overwrite=FALSE)

# wrap strip labels
#label_wrap_gen <- function(width = 20) {
#  function(variable, value) {
#    lapply(strwrap(as.character(value), width=width, simplify=FALSE), 
#           paste, collapse="\n")
#  }
#}

#labels <- multi0.results %>% select(c("run_label","run_label_n"))
#run_labeller <- function(variable,value){
#  return(labels[value])
#}

#### Plot intercept only results ####
multi0.results <- read.xlsx(paste0(outpath,"Multimodel multi0 results.xlsx"),sheet=2)

unique(multi0.results$Outcome)
unique(multi0.results$run_label)
unique(multi0.results$type)

multi0.results <- multi0.results %>%
  mutate(run_label = factor(run_label,levels=unique(run_label[order(desc(type))]))) %>%
  mutate(Outcome = ifelse(Outcome =="1","Win B-Win Y",Outcome)) %>%
  mutate(estimate = ifelse(is.na(n_cases),NA,estimate),
         ci.lb = ifelse(is.na(n_cases),NA,ci.lb),
         ci.ub = ifelse(is.na(n_cases),NA,ci.ub),
         Label_sig = ifelse(is.na(n_cases),NA,Label_sig)) %>%
  mutate(sensitivity_n = ifelse(n_studies<3 | n_cases < 10,"k<10","k>=10"))

multi0.results$Outcome <- fct_relevel(multi0.results$Outcome,c("Other","Lose B-Win Y","Win B-Lose Y","Win B-Win Y"))

transf <- function(x){
  x <- exp(x)
  return(x)
}

col.synergies = c( "Win B-Win Y"= "forestgreen", "Win B-Lose Y" ="gold","Lose B-Win Y"= "lightblue","Lose B-Lose Y" = "navy","Other"="grey50")

g <- ggplot(multi0.results[ !(multi0.results$term %in% c("1 + Effect_ID | IDTRUE")) & 
                              multi0.results$run %in% c("d_synergies","d_abun","d_rich","d_even") &
                              grepl("*Sig*",multi0.results$run_label)==FALSE  & #!(multi0.results$n_cases<5) &
                              grepl("*HQ*",multi0.results$run_label)==FALSE  &
                              grepl("*HQ*",multi0.results$run_label)==FALSE  &
                              !is.na(multi0.results$n_cases),],
            aes(y=run_label,x=transf(estimate)))+
  geom_vline(xintercept=1,linetype=2)+
  geom_linerange(aes(xmin=transf(ci.lb), xmax=transf(ci.ub),group=Outcome,colour=Outcome),position=position_dodge2(width=0.8),show.legend=FALSE)+
  geom_point(aes(colour=Outcome,shape=type,group=Outcome),size=3,position=position_dodge2(width=0.8))+
  scale_colour_manual(values=col.synergies[1:3,5],name="")+
  new_scale_color()+
  geom_text(aes(x=transf(ci.ub),label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=transf(ci.ub),label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+
  scale_colour_manual(values=c("k<10"="grey50","k>=10"="black"),name="")+
  labs(y="",x="Likelihood relative to lose-lose")+
  scale_x_continuous( expand = expansion(add = c(0.1,0.5)))+
  scale_shape_manual(values=c("Z_Overall"=15,"Abundance" = 16,"B_Richness" = 16,"Richness-Evenness" = 16,"Other" = 16))+
  scale_fill_manual(values=col.synergies[1:3,5],name="")+
  scale_y_discrete(labels=function(x)gsub("All cases - ","",x))+ #labels=function(x)str_wrap(x,28))+
  coord_cartesian(clip="off")+
  theme_classic()+
  #facet_wrap(~factor(Outcome,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")),
  #           scales="free_x",ncol=2)+
  theme(legend.background=element_blank(),
        #legend.position=c(0.4,-0.15),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=9,colour="black"),
        text=element_text(size=9,colour="black"),
        axis.title=element_text(size=9,colour="black"),
        axis.text.y=element_text(size=9,colour="black"),
        axis.text.x=element_text(size=9,colour="black"),
        #panel.grid.major.x = element_line(size=0.5,colour="grey60"),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,3,0.5,0.5),"lines"),
        panel.spacing = unit(1, "lines"),
        strip.text.x=element_text(size=9,face="bold",colour="black"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
        strip.background = element_blank(),
        strip.placement = "outside")+
  #guides(shape="none",colour=guide_legend(override.aes=(colour=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")))
  guides(shape="none")#colour=guide_legend(override.aes=(colour=col.synergies)))
g

#jpeg(paste0(outpath,"Fig X probability each outcome.jpg"),height=height,width=width,units="in",res=600)
#height=7
#width=8
#tiff(paste0(outpath,"Fig X probability each outcome.tif"),height=height,width=width,units="in",res=600,compression="lzw")

height=5
width=7.5
tiff(paste0(outpath,"Fig X probability each outcome_dodge.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig X probability each outcome_dodge.pdf"),height=height,width=width)
g
dev.off()

# truncate large confidence intervals to increase readability
ci.edit <- function(data,cutoff=7){
  data<- data %>% mutate(estimate.transf = transf(estimate),
                                            ci.lb.edit = transf(ci.lb),
                                            ci.ub.edit = transf(ci.ub)) %>% 
  mutate(ci.ub.edit = ifelse(estimate>0 & (transf(estimate)+transf(ci.ub)>cutoff),NA,ci.ub.edit)) %>%
  mutate(ci.ub.edit.arrow =  ifelse(transf(estimate)+transf(ci.ub)>cutoff,cutoff,NA))
}

multi0.results <- ci.edit(multi0.results)

g <- ggplot(multi0.results[ !(multi0.results$term %in% c("1 + Effect_ID | IDTRUE")) & #grepl("*(Sig)",multi0.results$run_label)==FALSE  & #!(multi0.results$n_cases<5) &
                              grepl("*(HQ)",multi0.results$run_label)==FALSE  &
                              multi0.results$run %in% c("d_synergies","d_abun","d_rich","d_even",
                                                        "d_synergies_sig","d_sig_abun","d_sig_rich","d_sig_even") &
                               !is.na(multi0.results$n_cases),],
            aes(y=run_label,x=transf(estimate)))+
  geom_vline(xintercept=1,linetype=2)+
  geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.ub.edit))+
  geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.lb.edit))+
  geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.ub.edit.arrow),arrow=arrow(length=unit(0.25,"cm"),angle=25))+
  #geom_linerange(aes(xmin=transf(ci.lb), xmax=transf(ci.ub),group=Outcome,colour=Outcome),position=position_dodge2(width=0.8),show.legend=FALSE)+
  geom_point(aes(colour=Outcome,shape=type,group=Outcome),size=3,position=position_dodge2(width=0.8))+
  scale_colour_manual(values=col.synergies[1:3,5],name="")+
  new_scale_color()+
  geom_text(aes(x=ci.ub.edit,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit.arrow,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit.arrow,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  scale_colour_manual(values=c("k<10"="grey50","k>=10"="black"),name="")+
  labs(y="",x="Likelihood relative to lose-lose")+
  scale_x_continuous( expand = expansion(add = c(0.2, 0.5)),labels=function(x)round(x,1))+
  scale_shape_manual(values=c("Z_Overall"=15,"Abundance" = 16,"B_Richness" = 16,"Richness-Evenness" = 16,"Other" = 16))+
  scale_fill_manual(values=col.synergies[1:3,5],name="")+
  scale_y_discrete(labels=function(x)gsub("All cases - ","",x))+ #labels=function(x)str_wrap(x,28))+
  coord_cartesian(clip="off")+
  theme_classic()+
  facet_wrap(~factor(Outcome,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")),
             scales="free_x",ncol=2)+
  theme(legend.background=element_blank(),
        #legend.position=c(0.4,-0.15),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=9,colour="black"),
        text=element_text(size=9,colour="black"),
        axis.title=element_text(size=9,colour="black"),
        axis.text.y=element_text(size=9,colour="black"),
        axis.text.x=element_text(size=9,colour="black"),
        #panel.grid.major.x = element_line(size=0.5,colour="grey60"),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,3,0.5,0.5),"lines"),
        panel.spacing = unit(1, "lines"),
        strip.text.x=element_text(size=9,face="bold",colour="black"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
        strip.background = element_blank(),
        strip.placement = "outside")+
  #guides(shape="none",colour=guide_legend(override.aes=(colour=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")))
  guides(shape="none")#colour=guide_legend(override.aes=(colour=col.synergies)))
g
height=7
width=8
tiff(paste0(outpath,"Fig X probability each outcome_with sig.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig X probability each outcome_with sig.pdf"),height=height,width=width)
g
dev.off()

# sensitivity analysis

g <- ggplot(multi0.results[ !(multi0.results$term %in% c("1 + Effect_ID | IDTRUE")) & #grepl("*(Sig)",multi0.results$run_label)==FALSE  & #!(multi0.results$n_cases<5) &
                              #grepl("*All biodiversity*",multi0.results$run_label)==TRUE  &
                              multi0.results$run %in% c("d_synergies","d_validity.high_synergies") &
                              !is.na(multi0.results$n_cases),],
            aes(y=factor(run_label,labels=c("All cases","Only high quality")),x=transf(estimate)))+
  geom_vline(xintercept=1,linetype=2)+
  #geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.ub.edit))+
  #geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.lb.edit))+
  #geom_segment(aes(y=run_label, yend=run_label, x=transf(estimate), xend=ci.ub.edit.arrow),arrow=arrow(length=unit(0.25,"cm"),angle=25))+
  geom_linerange(aes(xmin=transf(ci.lb), xmax=transf(ci.ub),group=Outcome,colour=Outcome),position=position_dodge2(width=0.8),show.legend=FALSE)+
  geom_point(aes(colour=Outcome,group=Outcome),shape=16,size=3,position=position_dodge2(width=0.8))+
  scale_colour_manual(values=col.synergies[1:3,5],name="")+
  new_scale_color()+
  geom_text(aes(x=ci.ub.edit,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit.arrow,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  geom_text(aes(x=ci.ub.edit.arrow,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+#nudge_x=0.1,nudge_y=0,
  scale_colour_manual(values=c("k<10"="grey50","k>=10"="black"),name="")+
  labs(y="",x="Likelihood relative to lose-lose")+
  scale_x_continuous( expand = expansion(add = c(0.2, 0.5)),labels=function(x)round(x,1))+
  #scale_shape_manual(16)+#values=c("Z_Overall"=15,"Abundance" = 16,"B_Richness" = 16,"Richness-Evenness" = 16,"Other" = 16))+
  scale_fill_manual(values=col.synergies[1:3,5],name="")+
  scale_y_discrete(labels=function(x)gsub("All cases - ","",x))+ #labels=function(x)str_wrap(x,28))+
  coord_cartesian(clip="off")+
  theme_classic()+
  #facet_wrap(~factor(Outcome,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")),scales="free_x",ncol=2)+
  theme(legend.background=element_blank(),
        #legend.position=c(0.4,-0.15),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=9,colour="black"),
        text=element_text(size=9,colour="black"),
        axis.title=element_text(size=9,colour="black"),
        axis.text.y=element_text(size=9,colour="black"),
        axis.text.x=element_text(size=9,colour="black"),
        #panel.grid.major.x = element_line(size=0.5,colour="grey60"),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,3,0.5,0.5),"lines"),
        panel.spacing = unit(1, "lines"),
        strip.text.x=element_text(size=9,face="bold",colour="black"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
        strip.background = element_blank(),
        strip.placement = "outside")+
  #guides(shape="none",colour=guide_legend(override.aes=(colour=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y","Lose B-Lose Y","Other")))
  guides(shape="none")#colour=guide_legend(override.aes=(colour=col.synergies)))
g
height=3.5
width=7
tiff(paste0(outpath,"Fig X probability each outcome_validity.tif"),height=height,width=width,units="in",res=600,compression="lzw")
g
dev.off()

pdf(paste0(outpath,"Fig X probability each outcome_validity.pdf"),height=height,width=width)
g
dev.off()

#### Model 1: univariate models ####

# If excluding intercept, we are testing if the coefficients of each level are equal to zero 
# and if significant, can conclude the level has an effect on outcome.
# If including intercept, we are testing pairwise differences between the coefficients of each level and 
# the coefficient of the reference class, so can only 
# conclude that effects on outcome are not equal across classes (called  DUMMY CODING OR TREATMENT CONTRASTS). 
# See here for dummy coding options: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
# Change reference level with relevel(variable, ref="reference level name").
# Wald tests are used to test for differences between coefficients in both cases.
# To test if whole variable improves model fit, need to run a likelihood ratio test using anova or lrtest.
# If one variable in model, pp gives exactly the same probabilities as proportions in original data.

# Useful refs:
# https://www.statology.org/interpret-log-likelihood/ 
# https://www.rdocumentation.org/packages/AICcmodavg/versions/2.3-1/topics/AICc

table(d$Synergies, predict(multi0))
freq(d$Synergies)

model=nnet::multinom(Synergies ~ Biome_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Biome_simp")))
model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Biome_simp")))
data=n_studies_filter(d,n_group=c("Synergies","Biome_simp"))
variable="Biome_simp"
variable_short="biome"
run_label="Biome - All cases"
run_name=paste0("d_","biome")
type=" All biodiversity metrics"
n_group="Synergies"

pp_function <- function(data,model=multi1.var,model0=multi0,variable,variable_short,run_name,run_label,n_group="Synergies",type){
 
  print(paste0("Multinomial model for ",variable," with cases #=",nrow(data)," and studies #=",length(unique(data$ID))))
  #print(summary(model))
  print(anova(model0,model))
  print(lrtest(model0,model))
  lr.model <- paste0(model$call$formula[3])
  lr <- data.frame(lrtest(model0,model))[2,]
  lr.LogLik <- lr$LogLik
  lr.Df <- lr$Df
  lr.Chisq <- lr$Chisq
  lr.Pvalue <- lr$Pr..Chisq.
  lr.Pvalue.bon <- lr.Pvalue*13
  
  #aicc <-  -2*(log(lr.LogLik,base=exp(1))) + (3*(nrow(data)/(nrow(data)-2)))
  aicc <- round(AICc(model),2)
  var.explained <- round(((model0$deviance-model$deviance)/model0$deviance*100),3)
  print(paste0("Variance explained (%) = ",round(((model0$deviance-model$deviance)/model0$deviance*100),2)))
  name <-paste0("multi1.", variable_short) 
  #assign(name,model,envir = .GlobalEnv)
  
  n_data <- data %>% 
    group_by_at(c(n_group,variable)) %>% summarise(n_cases=n(),n_studies=length(unique(ID)))
  
  #multi1.var.lr <- data.frame(broom::tidy(lrtest(model0,model)),model1="multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = d)",model2=paste0("multinom(Synergies ~ )",variable," + (1+Effect_ID|ID),data=d")) 
  #name <-paste0("multi1.", variable_short,".lr") 
  #assign(name,multi1.var.lr,envir = .GlobalEnv)
  
  multi1.var.result <- data.frame(broom::tidy(model),broom::glance(model)) 
  multi1.var.result <- multi1.var.result %>% rename(Outcome="y.level") %>%
    mutate(var.explained = var.explained,
           lr.model = lr.model,
           lr.LogLik = lr.LogLik,
           lr.Df = lr.Df,
           lr.Chisq = lr.Chisq,
           lr.Pvalue = lr.Pvalue,
           lr.Pvalue.bon = lr.Pvalue.bon,
           aicc = aicc)
  multi1.var.result$rrr <- exp(multi1.var.result$estimate)
  multi1.var.result <- multi1.var.result %>%
    mutate(estimate = round(estimate,3),
           std.error = round(std.error,3),
           statistic  = round(statistic,3),
           p.value = round(p.value,4),
           rrr = round(rrr,4)) %>% 
    mutate(rrr_pc = 100*(rrr-1)) %>% 
    mutate(ci.lb = estimate - 1.96*std.error,ci.ub = estimate + 1.96*std.error) 
  
  multi1.var.result <- multi1.var.result %>%
    mutate(level = term) %>%mutate(level = ifelse(level %in% c("(Intercept)","1 + Effect_ID | IDTRUE"),
                                                  level,substr(level,nchar(variable)+1,100))) %>%
    mutate(level = ifelse(level=="(Intercept)",levels(d[[variable]])[1],level))
  
  multi1.var.result <- multi1.var.result %>%
    left_join(n_data,by=c("Outcome"=n_group,"level"=variable))
  
  data <- data.frame(data)
  multi1.var_pp <- data.frame(unique(data[,c(variable)]), Effect_ID=1,ID=1) # if code blocks here, type in variable name
  multi1.var_pp
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
  #assign(name,multi1.var_pp,envir = .GlobalEnv)
  
  multi1.var.result <- multi1.var.result %>% full_join(multi1.var_pp,by=c("level","Outcome")) %>%
    mutate(Probability = round(Probability,3))
  
  multi1.var.result <- multi1.var.result[,c(1,2,ncol(multi1.var.result),3:(ncol(multi1.var.result)-1))] %>%
    mutate(run = run_name) %>%
    mutate(run_label = paste0(run_label)) %>%
    mutate(type = type) %>%
    mutate(mod = ifelse(substr(term,1,3)=="Sys", gsub("System_T*","",term),
                        ifelse(term=="(Intercept)","Reference",term))) %>%
    mutate(Label = ifelse(rrr>100 & p.value<0.001,paste0(">100,p<0.001"),
                          ifelse(rrr>100,paste0(">100,p=",round(p.value,3)),
                                 ifelse(p.value<0.001,paste0(round(rrr-1,2),",p<0.001"),paste0(round(rrr-1,2),",p=",round(p.value,3)))))) %>%
    mutate(Label_sig = ifelse(is.nan(p.value),"*",ifelse(p.value<0.05,"*",""))) %>%
    mutate(Label_simp = factor(Label_simp,levels=unique(Label_simp[order(desc(Outcome),estimate)]))) %>%
    mutate(run_label_n = paste0(Label_simp," (",n_cases,")")) %>% 
    dplyr::select(-c(variable)) %>% 
    mutate(run_label_n = factor(run_label_n,levels=unique(run_label_n[order(type)]))) %>%
    filter(!is.na(n_cases))
  
  name <-paste0("multi1.", run_name,".result") 
  assign(name,multi1.var.result,envir = .GlobalEnv)
}

#d$Biome <- fct_relevel(d$Biome,"Temperate Grasslands, Savannas & Shrublands")
#table(d$Synergies,d$System_T)

# Biome_simp ####
pp_function(model=nnet::multinom(Synergies ~ Biome_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Biome_simp"))),data=n_studies_filter(d,n_group=c("Synergies","Biome_simp")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Biome_simp"))), 
            variable="Biome_simp",  variable_short="biome", run_label="Biome - All cases",
            run_name=paste0("d_","biome"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Biome_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Biome_simp"))), data=n_studies_filter(d_abun,n_group=c("Synergies","Biome_simp")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Biome_simp"))), 
            variable="Biome_simp",  variable_short="biome", run_label="Biome - All cases",
            run_name=paste0("d_abun_","biome"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Biome_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Biome_simp"))), data=n_studies_filter(d_rich,n_group=c("Synergies","Biome_simp")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Biome_simp"))),
            variable="Biome_simp",  variable_short="biome", run_label="Biome - All cases",
            run_name=paste0("d_rich_","biome"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Biome_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Biome_simp"))), data=n_studies_filter(d_even,n_group=c("Synergies","Biome_simp")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Biome_simp"))),
            variable="Biome_simp",  variable_short="biome", run_label="Biome - All cases",
            run_name=paste0("d_even_","biome"),type="Richness-Evenness")

# System_T ####
pp_function(model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","System_T"))),data=n_studies_filter(d,n_group=c("Synergies","System_T")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","System_T"))),
            variable="System_T",  variable_short="practice", run_label="Practice - All cases",
            run_name=paste0("d_","practice"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","System_T"))), data=n_studies_filter(d_abun,n_group=c("Synergies","System_T")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","System_T"))),
            variable="System_T",  variable_short="practice", run_label="Practice - All cases",
            run_name=paste0("d_abun_","practice"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","System_T"))), data=n_studies_filter(d_rich,n_group=c("Synergies","System_T")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","System_T"))),
            variable="System_T",  variable_short="practice", run_label="Practice - All cases",
            run_name=paste0("d_rich_","practice"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","System_T"))), data=n_studies_filter(d_even,n_group=c("Synergies","System_T")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","System_T"))),
            variable="System_T",  variable_short="practice", run_label="Practice - All cases",
            run_name=paste0("d_even_","practice"),type="Richness-Evenness")

# Agrochemicals ####
pp_function(model=nnet::multinom(Synergies ~ Agrochem_CT -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Agrochem_CT"))),data=n_studies_filter(d,n_group=c("Synergies","Agrochem_CT")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Agrochem_CT"))),
            variable="Agrochem_CT",  variable_short="agrochem", run_label="Agrochemicals - All cases",
            run_name=paste0("d_","agrochem"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Agrochem_CT -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Agrochem_CT"))), data=n_studies_filter(d_abun,n_group=c("Synergies","Agrochem_CT")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Agrochem_CT"))),
            variable="Agrochem_CT",  variable_short="agrochem", run_label="Agrochemicals - All cases",
            run_name=paste0("d_abun_","agrochem"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Agrochem_CT -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Agrochem_CT"))), data=n_studies_filter(d_rich,n_group=c("Synergies","Agrochem_CT")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Agrochem_CT"))),
            variable="Agrochem_CT",  variable_short="agrochem", run_label="Agrochemicals - All cases",
            run_name=paste0("d_rich_","agrochem"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Agrochem_CT -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Agrochem_CT"))), data=n_studies_filter(d_even,n_group=c("Synergies","Agrochem_CT")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Agrochem_CT"))),
            variable="Agrochem_CT",  variable_short="agrochem", run_label="Agrochemicals - All cases",
            run_name=paste0("d_even_","agrochem"),type="Richness-Evenness")

# Crop FAO ####

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=d),
            data=d,  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=d),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_","crops"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Crop_FAO_C"))),
            data=n_studies_filter(d,n_group=c("Synergies","Crop_FAO_C")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Crop_FAO_C"))),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_","crops"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_FAO_C"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_FAO_C")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_FAO_C"))),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_abun_","crops"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_FAO_C"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_FAO_C")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_FAO_C"))),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_rich_","crops"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Crop_FAO_C"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","Crop_FAO_C")),             
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Crop_FAO_C"))),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_even_","crops"),type="Richness-Evenness")

# Crop type ####
data <- n_studies_filter(d,n_group=c("Synergies","Crop_type_C"))
data <- n_studies_filter(d,n_group=c("Synergies","Crop_FAO_C"))

pp_function(model=nnet::multinom(Synergies ~ Crop_type_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Crop_type_C"))),
            data=n_studies_filter(d,n_group=c("Synergies","Crop_type_C")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Crop_type_C"))),
            variable="Crop_type_C",  variable_short="crop_type", run_label="Crop type - All cases",
            run_name=paste0("d_","crop_type"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Crop_type_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_type_C"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_type_C")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Crop_type_C"))),
            variable="Crop_type_C",  variable_short="crop_type", run_label="Crop type - All cases",
            run_name=paste0("d_abun_","crop_type"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Crop_type_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_type_C"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_type_C")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Crop_type_C"))),
            variable="Crop_type_C",  variable_short="crop_type", run_label="Crop type - All cases",
            run_name=paste0("d_rich_","crop_type"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Crop_type_C -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Crop_type_C"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","Crop_type_C")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Crop_type_C"))),
            variable="Crop_type_C",  variable_short="crop_type", run_label="Crop type - All cases",
            run_name=paste0("d_even_","crop_type"),type="Richness-Evenness")

# Region ####

pp_function(model=nnet::multinom(Synergies ~ Region.Name -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Region.Name"))),
            data=n_studies_filter(d,n_group=c("Synergies","Region.Name")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Region.Name"))),
            variable="Region.Name",  variable_short="region", run_label="Region - All cases",
            run_name=paste0("d_","region"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Region.Name -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Region.Name"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","Region.Name")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Region.Name"))),
            variable="Region.Name",  variable_short="region", run_label="Region - All cases",
            run_name=paste0("d_abun_","region"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Region.Name -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Region.Name"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","Region.Name")),
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Region.Name"))),
            variable="Region.Name",  variable_short="region", run_label="Region - All cases",
            run_name=paste0("d_rich_","region"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Region.Name -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Region.Name"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","Region.Name")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Region.Name"))),
            variable="Region.Name",  variable_short="region", run_label="Region - All cases",
            run_name=paste0("d_even_","region"),type="Richness-Evenness")

# Taxa group ####

pp_function(model=nnet::multinom(Synergies ~ Taxa_group_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Taxa_group_simp"))),
            data=n_studies_filter(d,n_group=c("Synergies","Taxa_group_simp")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Taxa_group_simp"))),
            variable="Taxa_group_simp",  variable_short="taxa", run_label="Taxa - All cases",
            run_name=paste0("d_","taxa"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Taxa_group_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Taxa_group_simp"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","Taxa_group_simp")),
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Taxa_group_simp"))),
            variable="Taxa_group_simp",  variable_short="taxa", run_label="Taxa - All cases",
            run_name=paste0("d_abun_","taxa"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Taxa_group_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Taxa_group_simp"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","Taxa_group_simp")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Taxa_group_simp"))),
            variable="Taxa_group_simp",  variable_short="taxa", run_label="Taxa - All cases",
            run_name=paste0("d_rich_","taxa"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Taxa_group_simp -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Taxa_group_simp"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","Taxa_group_simp")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Taxa_group_simp"))),
            variable="Taxa_group_simp",  variable_short="taxa", run_label="Taxa - All cases",
            run_name=paste0("d_even_","taxa"),type="Richness-Evenness")

# Pest group ####

pp_function(model=nnet::multinom(Synergies ~ Pest_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Pest_group"))),
            data=n_studies_filter(d,n_group=c("Synergies","Pest_group")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Pest_group"))),
            variable="Pest_group",  variable_short="pests", run_label="Pest group - All cases",
            run_name=paste0("d_","pests"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ Pest_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Pest_group"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","Pest_group")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","Pest_group"))),
            variable="Pest_group",  variable_short="pests", run_label="Pest group - All cases",
            run_name=paste0("d_abun_","pests"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ Pest_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Pest_group"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","Pest_group")),
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","Pest_group"))),
            variable="Pest_group",  variable_short="pests", run_label="Pest group- All cases",
            run_name=paste0("d_rich_","pests"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ Pest_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Pest_group"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","Pest_group")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","Pest_group"))),
            variable="Pest_group",  variable_short="pests", run_label="Pest group - All cases",
            run_name=paste0("d_even_","pests"),type="Richness-Evenness")

# Biodiversity metric ####
pp_function(model=nnet::multinom(Synergies ~ B_measure_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","B_measure_group"))),
            data=n_studies_filter(d,n_group=c("Synergies","B_measure_group")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","B_measure_group"))),
            variable="B_measure_group",  variable_short="b_measure", run_label="Biodiversity metric - All cases",
            run_name=paste0("d_","b_measure"),type=" All biodiversity metrics")

# Yield metric ####
pp_function(model=nnet::multinom(Synergies ~ Yield_measure_group -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Yield_measure_group"))),
            data=n_studies_filter(d,n_group=c("Synergies","Yield_measure_group")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","Yield_measure_group"))),
            variable="Yield_measure_group",  variable_short="yield_measure", run_label="Yield metric - All cases",
            run_name=paste0("d_","y_measure"),type=" All biodiversity metrics")

# Ground group ####

pp_function(model=nnet::multinom(Synergies ~ B_ground -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","B_ground"))),
            data=n_studies_filter(d,n_group=c("Synergies","B_ground")),  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","B_ground"))),
            variable="B_ground",  variable_short="b_ground", run_label="Ground relation - All cases",
            run_name=paste0("d_","b_ground"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ B_ground -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","B_ground"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","B_ground")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","B_ground"))),
            variable="B_ground",  variable_short="b_ground", run_label="Ground relation - All cases",
            run_name=paste0("d_abun_","b_ground"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ B_ground -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","B_ground"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","B_ground")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","B_ground"))),
            variable="B_ground",  variable_short="b_ground", run_label="Ground relation - All cases",
            run_name=paste0("d_rich_","b_ground"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ B_ground -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","B_ground"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","B_ground")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","B_ground"))),
            variable="B_ground",  variable_short="b_ground", run_label="Ground relation - All cases",
            run_name=paste0("d_even_","b_ground"),type="Richness-Evenness")

# Development group ####

pp_function(model=nnet::multinom(Synergies ~ DevelopmentStatus -1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","DevelopmentStatus"))),
            data=n_studies_filter(d,n_group=c("Synergies","DevelopmentStatus")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d,n_group=c("Synergies","DevelopmentStatus"))),
            variable="DevelopmentStatus",  variable_short="dev", run_label="Development status - All cases",
            run_name=paste0("d_","dev"),type=" All biodiversity metrics")

pp_function(model=nnet::multinom(Synergies ~ DevelopmentStatus -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","DevelopmentStatus"))), 
            data=n_studies_filter(d_abun,n_group=c("Synergies","DevelopmentStatus")),
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_abun,n_group=c("Synergies","DevelopmentStatus"))),
            variable="DevelopmentStatus",  variable_short="dev", run_label="Development status- All cases",
            run_name=paste0("d_abun_","dev"),type="Abundance")

pp_function(model=nnet::multinom(Synergies ~ DevelopmentStatus -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","DevelopmentStatus"))), 
            data=n_studies_filter(d_rich,n_group=c("Synergies","DevelopmentStatus")), 
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_rich,n_group=c("Synergies","DevelopmentStatus"))),
            variable="DevelopmentStatus",  variable_short="dev", run_label="Development status - All cases",
            run_name=paste0("d_rich_","dev"),type="Richness")

pp_function(model=nnet::multinom(Synergies ~ DevelopmentStatus -1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","DevelopmentStatus"))), 
            data=n_studies_filter(d_even,n_group=c("Synergies","DevelopmentStatus")),
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=n_studies_filter(d_even,n_group=c("Synergies","DevelopmentStatus"))),
            variable="DevelopmentStatus",  variable_short="dev", run_label="Development status - All cases",
            run_name=paste0("d_even_","dev"),type="Richness-Evenness")

# Latitude ####

multi0.d <- nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_rich)
multi1.lat <- nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d_rich)
variable_short="d_rich"
type="Richness"

multi0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d)
multi1 = nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d)
variable_short="d_"
type=""

lat_loop <- function(multi0,multi1,variable_short,type){
  summary(multi1) 
  print(lrtest(multi0,multi1)) # significant, Chi sq  71.37
  
  lr.model <- paste0(multi1$call$formula[3])
  lr <- data.frame(lrtest(multi0,multi1))[2,]
  lr.LogLik <- lr$LogLik
  lr.Df <- lr$Df
  lr.Chisq <- lr$Chisq
  lr.Pvalue <- lr$Pr..Chisq.
  lr.Pvalue.bon <- lr.Pvalue*13
  aicc <- round(AICc(multi1),2)
  var.explained <- round(((multi0$deviance-multi1$deviance)/multi0$deviance*100),3)
  print(paste0("Variance explained (%) = ",round(((multi0$deviance-multi1$deviance)/multi0$deviance*100),2)))
  
  multi1.lat.lr <- data.frame(broom::tidy(lrtest(multi0,multi1)),
                              model1=paste0("multinom(formula = Synergies ~ 1 + (1 + Effect_ID | ID), data = ",variable_short),
                              model2=paste0("multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=",variable_short)) 
  multi1.lat.result <- data.frame(broom::tidy(multi1),broom::glance(multi1)) 
  
  multi1.lat.result <- multi1.lat.result %>% #rename(Outcome="y.level") %>%
    mutate(var.explained = var.explained,
           lr.model = lr.model,
           lr.LogLik = lr.LogLik,
           lr.Df = lr.Df,
           lr.Chisq = lr.Chisq,
           lr.Pvalue = lr.Pvalue,
           lr.Pvalue.bon = lr.Pvalue.bon,
           aicc = aicc)
  
  multi1.lat.result <- multi1.lat.result %>% 
    mutate(rrr = round(exp(estimate),3)) %>%
    mutate(Label_sig = ifelse(is.nan(p.value),"*",ifelse(p.value<0.05,"*","")))
  
  assign("multi1.results.lat",multi1.lat.result,envir = .GlobalEnv )
  
  multi1.lat_pp <- data.frame( Lat_T=rep(c(-30:70)), Effect_ID=1,ID=1)
  multi1.lat_pp <- cbind(multi1.lat_pp,predict(multi1,newdata=multi1.lat_pp,type="probs",se=TRUE))
  
  setDT(multi1.lat_pp)
  multi1.lat_pp <- melt(multi1.lat_pp,id.vars=c("Lat_T","Effect_ID","ID"),variable.name="Outcome",value.name="Probability")
  multi1.lat_pp <- multi1.lat_pp %>% left_join(multi1.lat.result[!(multi1.lat.result$term %in% c("1 + Effect_ID | IDTRUE","(Intercept)")), 
                                                                 c("y.level","estimate","p.value","Label_sig")],by=c("Outcome"="y.level")) 
  labels <- multi1.lat_pp %>% filter(Lat_T==70) %>% mutate(p.value = ifelse(p.value<0.001,"<0.001",as.character(round(p.value,3))))
  print(labels)
  
  col.synergies = c( "Win B-Win Y"= "forestgreen", "Win B-Lose Y" ="gold","Lose B-Win Y"= "lightblue","Lose B-Lose Y" = "navy","Other"="grey50")
  
  g <- ggplot(multi1.lat_pp,aes(x=Lat_T,y=Probability,colour=Outcome))+geom_line(size=1)+
    scale_colour_manual(values=col.synergies,name="")+
    scale_y_continuous(expand=c(0,0),limits=c(0,1))+
    geom_text(aes(colour="Win B-Win Y", x=50,y=0.95,label=paste0("RRR=",round(exp(labels[labels$Outcome=="Win B-Win Y",]$estimate),3),",p=",labels[labels$Outcome=="Win B-Win Y",]$p.value,labels[labels$Outcome=="Win B-Win Y",]$Label_sig)),size=3.2,family="sans",show.legend=FALSE)+
    geom_text(aes(colour="Win B-Lose Y", x=50,y=0.90,label=paste0("RRR=",round(exp(labels[labels$Outcome=="Win B-Lose Y",]$estimate),3),",p=",labels[labels$Outcome=="Win B-Lose Y",]$p.value,labels[labels$Outcome=="Win B-Lose Y",]$Label_sig)),size=3.2,family="sans",show.legend=FALSE)+
    geom_text(aes(colour="Lose B-Win Y", x=50,y=0.85,label=paste0("RRR=",round(exp(labels[labels$Outcome=="Lose B-Win Y",]$estimate),3),",p=",labels[labels$Outcome=="Lose B-Win Y",]$p.value,labels[labels$Outcome=="Lose B-Win Y",]$Label_sig)),size=3.2,family="sans",show.legend=FALSE)+
    geom_text(aes(colour="Other", x=50,y=0.8,label=paste0("RRR=",round(exp(labels[labels$Outcome=="Other",]$estimate),3),",p=",labels[labels$Outcome=="Other",]$p.value,labels[labels$Outcome=="Other",]$Label_sig)),size=3.2,family="sans",show.legend=FALSE)+
    geom_text(data=labels,aes(x=Lat_T,label=Label_sig,hjust=0),nudge_x=0,nudge_y=0, size=6,family="sans",show.legend=FALSE)+
    ggtitle(type)+
    labs(x="Latitude",y="Probability")+
    theme_classic()+
    theme(line=element_line(size=0.3,colour="black"),
          #axis.text.x=element_text(angle=90,hjust=1,vjust=0.2),
          panel.grid.major.y=element_line(size=0.3,colour="grey50",linetype="dashed"),
          strip.background=element_rect(fill="white"),
          plot.title=element_text(hjust=0.5,size=10,colour="black"),
          axis.text.y=element_text(size=9,colour="black"),
          axis.text.x=element_text(size=9,colour="black"),
          text=element_text(size=10,colour="black"),
          legend.position="bottom",legend.title=element_blank())+
    guides(colour=guide_legend(nrow=1))
  
  print(g)
  
  #temp <- d %>% group_by(Lat_T,Long_T,Synergies) %>% summarise(n_cases = n())
  #g2 <- ggplot(temp,aes(x=Lat_T,y=Long_T,colour=Synergies,size=n_cases))+geom_point(alpha=0.5)
  #g2
  
  assign(paste0("fig_lat_",variable_short),g,envir = .GlobalEnv)
  
  height=4
  width=5.5
  tiff(paste0(outpath,"Fig probability each outcome by lat_",variable_short,".tif"),height=height,width=width,units="in",res=600,compression="lzw")
  print(g+ggtitle(""))
  dev.off()
  
  pdf(paste0(outpath,"Fig probability each outcome by lat_",variable_short,".pdf"),height=height,width=width)
  print(g+ggtitle(""))
  dev.off()
}

lat_loop(multi0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d),
         multi1 = nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d),variable_short="d",type="All biodiversity metrics")
lat_loop(multi0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_abun),
         multi1 = nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d_abun),variable_short="d_abun",type="Abundance")
lat_loop(multi0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_rich),
         multi1 = nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d_rich),variable_short="d_rich",type="Richness")
lat_loop(multi0 = nnet::multinom(Synergies ~ 1 + (1+Effect_ID|ID),data=d_even),
         multi1 = nnet::multinom(Synergies ~ Lat_T + (1+Effect_ID|ID),data=d_even),variable_short="d_even",type="Richness-Evenness")

plot.legend <- get_legend(fig_lat_d+theme(legend.position="bottom")+guides(colour=guide_legend(nrow=1)))

p <- plot_grid(fig_lat_d+theme(legend.position="none"),
               fig_lat_d_abun+theme(legend.position="none"),
               fig_lat_d_rich+theme(legend.position="none"),
               fig_lat_d_even+theme(legend.position="none"),
               ncol=2,
               align= "hv",axis="l",
               labels = c("A)","B)","C)","D)"),label_size=10,rel_heights=c(1,1))
p

p <- plot_grid(p,plot.legend,
               ncol=1,
               align= "hv",axis="l",
               labels = c(""),rel_heights=c(1,0.05))
p  


tiff(paste0(outpath,"Fig X probability each outcome lat.tif"), width=7.5,height=8,units="in",res=600,compression="lzw")
p
dev.off()

pdf(paste0(outpath,"Fig X probability each outcome lat.pdf"), width=7.5,height=8)
p
dev.off()


# Export multi1 results ####

multi1.results.biome <- rbind(multi1.d_biome.result, multi1.d_abun_biome.result,
                              multi1.d_rich_biome.result,multi1.d_even_biome.result)

multi1.results.practice <- rbind(multi1.d_practice.result, multi1.d_abun_practice.result,
                              multi1.d_rich_practice.result,multi1.d_even_practice.result)

multi1.results.agrochem <- rbind(multi1.d_agrochem.result, multi1.d_abun_agrochem.result,
                              multi1.d_rich_agrochem.result,multi1.d_even_agrochem.result)

multi1.results.crops <- rbind(multi1.d_crops.result, multi1.d_abun_crops.result,
                              multi1.d_rich_crops.result)

multi1.results.crop_type <- rbind(multi1.d_crop_type.result, multi1.d_abun_crop_type.result,
                              multi1.d_rich_crop_type.result,multi1.d_even_crop_type.result)

multi1.results.region <- rbind(multi1.d_region.result, multi1.d_abun_region.result,
                              multi1.d_rich_region.result,multi1.d_even_region.result)

multi1.results.taxa <- rbind(multi1.d_taxa.result, multi1.d_abun_taxa.result,
                              multi1.d_rich_taxa.result,multi1.d_even_taxa.result)

multi1.results.pests <- rbind(multi1.d_pests.result, multi1.d_abun_pests.result,
                              multi1.d_rich_pests.result,multi1.d_even_pests.result)

multi1.results.b_ground <- rbind(multi1.d_b_ground.result,multi1.d_even_b_ground.result)

multi1.results.dev<- rbind(multi1.d_dev.result,multi1.d_abun_dev.result,multi1.d_rich_dev.result,multi1.d_even_dev.result)

multi1.results.b_measure <- multi1.d_b_measure.result

multi1.results.y_measure <- multi1.d_y_measure.result


multi1.results.all <- rbind(multi1.d_biome.result,multi1.d_practice.result,multi1.d_agrochem.result,multi1.d_crops.result,
                        multi1.d_crop_type.result,multi1.d_region.result,multi1.d_taxa.result,multi1.d_pests.result,
                        multi1.d_b_ground.result,multi1.d_b_measure.result,multi1.d_y_measure.result) %>%
  mutate(lr.Pvalue.bon =lr.Pvalue*12) # multiply by 12 because this is the number of univariate tests (excluding development)

var.explained <- multi1.results.all %>% dplyr::select(lr.model,nobs,lr.Df,lr.Chisq,lr.Pvalue,lr.Pvalue.bon, lr.LogLik,var.explained,aicc) %>% 
  unique() %>% arrange(desc(var.explained)) %>%
  mutate(lr.Chisq = round(lr.Chisq,1),
         lr.Pvalue = ifelse(lr.Pvalue<0.001,"<0.001",round(lr.Pvalue,3)),
         lr.Pvalue.bon = ifelse(lr.Pvalue.bon<0.001,"<0.001",round(lr.Pvalue.bon,3)),
         lr.LogLik = round(lr.LogLik,0),
         var.explained = round(var.explained,1),
    aicc = round(aicc,1))

write.xlsx(x=list("all" = multi1.results.all,
                  "var explained" = var.explained,
                  "biome" = multi1.results.biome,"practice"=multi1.results.practice,"agrochem"=multi1.results.agrochem,
                  "crops"=multi1.results.crops,"crop_type"=multi1.results.crop_type,
                  "region" = multi1.results.region, "taxa"=multi1.results.taxa,"pests"=multi1.results.pests,
                  "b_ground" = multi1.results.b_ground,"dev"=multi1.results.dev,
                  "b_measure"=multi1.results.b_measure,"y_measure"= multi1.results.y_measure,"lat"=multi1.results.lat),
           #file=paste0(outpath,"Multinom multi1 results nlt5 excluded",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)
           file=paste0(outpath,"Multinom multi1 results",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

# Repeat using a different outcome as reference, to include variable levels with no lose-lose ####

# System_T ####
data <- d %>% filter(System_T %in% c("Combined practices","Crop rotation","Cultivar mixture"))
#data <- d %>% filter(!(System_T %in% c("Cultivar mixture")))
#data <- d
table(data$System_T,data$Synergies)
data <- data %>% mutate(Synergies = fct_relevel(Synergies,"Lose B-Win Y")) %>% droplevels()
table(data$System_T,data$Synergies)
data <- droplevels(data)

pp_function(model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=data),
            data=data,  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=data),
            variable="System_T",  variable_short="practice", run_label="Practice - All cases",
            run_name=paste0("d_refLW_","practice"),type=" All biodiversity metrics")

model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=data)
model=nnet::multinom(Synergies ~ System_T -1 + (1 + Effect_ID|ID),data=data)
anova(model,model0)
lrtest(model,model0)
AICc(model0)

# Taxa group ####
data <- d 
table(data$Taxa_group_simp,data$Synergies)
data <- data %>% mutate(Synergies = fct_relevel(Synergies,"Win B-Lose Y")) %>% droplevels()

pp_function(model=nnet::multinom(Synergies ~ Taxa_group_simp  - 1 + (1 + Effect_ID|ID),data=data),
            data=data,  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=data),
            variable="Taxa_group_simp",  variable_short="practice", run_label="Taxa - All cases",
            run_name=paste0("d_refWL_","taxa"),type=" All biodiversity metrics")

model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=data)
model=nnet::multinom(Synergies ~ Taxa_group_simp  - 1 + (1 + Effect_ID|ID),data=data)
AICc(model0)
AICc(model)
anova(model,model0)
print(paste0("Variance explained (%) = ",round(((model0$deviance-model$deviance)/model0$deviance*100),2)))

# Crop commodity ####
data <- d
data <- d %>% filter(!(Crop_FAO_C %in% c("Roots & tubers","Nuts/Stimulants","Fibres")))
table(data$Crop_FAO_C,data$Synergies)
#data <- data %>% mutate(Synergies = fct_relevel(Synergies,"Win B-Lose Y")) 
data <- data %>% droplevels()

pp_function(model=nnet::multinom(Synergies ~ Crop_FAO_C -1 + (1 + Effect_ID|ID),data=data),
            data=data,  
            model0=nnet::multinom(Synergies ~ 1 + (1 + Effect_ID|ID),data=data),
            variable="Crop_FAO_C",  variable_short="crops", run_label="Crops - All cases",
            run_name=paste0("d_refLL_","crops"),type=" All biodiversity metrics")

#check quality %
freq(d$Validity_biodiversity_overall)
freq(d$Validity_yield_overall)
freq(d$Validity_location)
freq(d$Validity_biodiversity_overall)
freq(d$Validity_biodiversity_overall)
check <- d %>% filter(Validity_location =="nd") # these ones had no precise location info. Counted as nd.
freq(d$Validity_time_C)
freq(d$Validity_time_T)
check <- d %>% filter(Validity_time_C !=Validity_time_T) # these represent sites with nd for control, >1 year for treatment. Counted as nd.
freq(d$Validity_N_biodiversity_C) # use this one, but remove one case from 1
freq(d$Validity_N_biodiversity_T)
check <- d %>% filter(Validity_N_biodiversity_C !=Validity_N_biodiversity_T) # these represent sites where either C (10 cases) or T (1 case) have <5. Counted as 0.
freq(d$Validity_N_yield_C) # use this one, but remove one case from 1
freq(d$Validity_N_yield_T)
check <- d %>% filter(Validity_N_yield_C !=Validity_N_yield_T) # same as bio, these represent sites where either C (10 cases) or T (1 case) have <5. Counted as 0.

# Combine, export ####

multi1.results.all.refLW <- rbind(multi1.d_refLW_practice.result,multi1.d_refWL_taxa.result) %>%
  mutate(lr.Pvalue.bon =lr.Pvalue*1) # multiply by 1 because this is the number of univariate tests on the data subset

var.explained.refLW <- multi1.results.all.refLW %>% dplyr::select(lr.model,nobs,lr.Df,lr.Chisq,lr.Pvalue,lr.Pvalue.bon, lr.LogLik,var.explained,aicc) %>% 
  unique() %>% arrange(desc(var.explained)) %>%
  mutate(lr.Chisq = round(lr.Chisq,1),
         lr.Pvalue = ifelse(lr.Pvalue<0.001,"<0.001",round(lr.Pvalue,3)),
         lr.Pvalue.bon = ifelse(lr.Pvalue.bon<0.001,"<0.001",round(lr.Pvalue.bon,3)),
         lr.LogLik = round(lr.LogLik,0),
         var.explained = round(var.explained,1),
         aicc = round(aicc,1))

write.xlsx(x=list("all" = multi1.results.all.refLW,
                  "var explained" = var.explained.refLW),
           file=paste0(outpath,"Multinom multi1 results_refLW_",format(Sys.Date(),"%Y%m%d"),".xlsx"),overwrite=TRUE)

# Make plots ####

data=multi1.results.crops
run_label=c("Crops - All cases")
variable_short="crops"
type=c(" All biodiversity metrics")
height=5
width=7.5
graph_limits=c(10, 10)
ncol=4
cutoff=2500

data <- multi1.results.biome
run_label=c("Biome - All cases")
variable_short="biome"
type=c(" All biodiversity metrics")
ncol=2
graph_limits=c(0.1, 0.2)
height=6
width=7
cutoff=7

data = multi1.d_refLW_practice.result
run_label=c("Practice - All cases")
variable_short="practice_refLW"
type=c(" All biodiversity metrics")
graph_limits=c(0.1, 0.2)
height=4.5
col.points=col.synergies[c(1,2,3,5)]
ref="Lose B-Win Y"
height=6
width=7
cutoff=7

multi1_plot <- function(data=multi1.results.all,run_label,variable_short,type,height=7,width=7,graph_limits=c(0.2, 0.2),ncol=3,cutoff=7,col.points=col.synergies[c(1:3,5)],ref="lose-lose"){
  
  data <- data %>% 
    mutate(estimate.transf = transf(estimate),
                          ci.lb.edit = transf(ci.lb),
                          ci.ub.edit = transf(ci.ub)) %>%
    mutate(estimate.edit = ifelse(estimate.transf > cutoff,NA,round(estimate.transf,3))) %>%
    mutate(ci.lb.edit = ifelse((transf(ci.lb)>cutoff),cutoff-0.1,ci.lb.edit)) %>%
    mutate(ci.ub.edit = ifelse((transf(estimate)+transf(ci.ub))>cutoff,as.numeric(NA),ci.ub.edit)) %>%
    mutate(ci.ub.edit.arrow =  ifelse((transf(estimate)+transf(ci.ub))>cutoff,cutoff,NA)) %>%
    #ci.ub.edit = ifelse(ci.ub>1.791759,1.791759,ci.ub),ifelse(ci.ub.limit = ifelse(ci.ub>1.791759,1.791759,NA))
    mutate(sensitivity_n = ifelse(n_studies<3 | n_cases < 10,"k<10","k>=10"))
  
  data <- data[which(data$type %in% c(type)),]
  
  
  g <- ggplot(data[ !(data$term %in% c("1 + Effect_ID | IDTRUE")) & 
                      #!(data$Outcome =="Other" & data$ci.ub>50) &
                      !is.na(data$term) & (data$run_label %in% run_label) &
                      (data$type %in% type),],
            aes(y=Label_simp,x=estimate.edit),position=position_dodge2(width=0.8))+
    geom_vline(xintercept=1,linetype=2)+
    geom_linerange(aes(y=Label_simp,xmin=ci.lb.edit, xmax=ci.ub.edit,colour=Outcome),position=position_dodge2(width=0.8),show.legend=FALSE)+
    geom_linerange(aes(y=Label_simp,xmin=ci.lb.edit, xmax=ci.ub.edit.arrow,colour=Outcome),position=position_dodge2(width=0.8),show.legend=FALSE)+
    geom_point(aes(colour=Outcome),shape=16,size=3,position=position_dodge2(width=0.8))+
    scale_colour_manual(values=col.points,name="")+
    new_scale_color() +
    geom_text(aes(x=ci.ub.edit,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+
    geom_text(aes(x=ci.ub.edit,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+
    geom_text(aes(x=ci.ub.edit.arrow,label=paste0(Label_sig),hjust=0,group=Outcome,colour=sensitivity_n), size=6,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+
    geom_text(aes(x=ci.ub.edit.arrow,label=paste0("   (",n_cases,",",n_studies, ")"),hjust=0,group=Outcome,colour=sensitivity_n), size=3.1,family="sans",position=position_dodge2(width=0.8),show.legend=FALSE)+
    scale_colour_manual(values=c("k<10"="grey50","k>=10"="black"),name="")+
    labs(y="",# x="Practice, crop type and biome,\n(reference associated plant species,\ncereal crops, tropical grasslands",
       #x="Log odds relative to lose-lose")+
       x=paste0("Likelihood relative to ",ref))+
    coord_cartesian(clip="off")+
    #scale_x_continuous(labels=function(x)round(transf(x),0),expand = expansion(add = graph_limits))+
    scale_x_continuous(labels=function(x)round(x,1),breaks=seq(0,cutoff,1), expand = expansion(add = graph_limits))+
    scale_fill_manual(values=col.points,name="")+
    scale_y_discrete(labels=function(x)str_wrap(x,15))+
    theme_classic()+
    #facet_wrap(~type+factor(Outcome,levels=c("Win B-Win Y","Win B-Lose Y","Lose B-Win Y")),scales="free_x",ncol=ncol,drop=T)+
    #facet_wrap(~type,scales="free_x",ncol=ncol, drop=T)+
    theme(legend.background=element_blank(),
          legend.position="bottom",
          legend.text=element_text(size=9,colour="black"),
          text=element_text(size=9,colour="black"),
          axis.title=element_text(size=9,colour="black"),
          axis.text.y=element_text(size=9,colour="black"),
          axis.text.x=element_text(size=9,colour="black"),
          #panel.grid.minor.y = element_line(size=0.5,colour="black"),
          axis.ticks.y=element_blank(),
          plot.margin=unit(c(0,3,0,0),"lines"),
          panel.spacing = unit(1,"lines"),
          strip.text.x=element_text(size=9,face="bold",colour="black"),
          strip.text.y=element_text(angle=0,hjust=0,face="bold",size=9,colour="black"),
          strip.background = element_blank(),
          strip.placement = "outside")
print(g)
  
gb <- ggplot_build(g)
lay <- data.frame(gb$data[3]) # from the third geom, so the one using ci.ub.edit.arrow
g <- g + geom_segment(data=lay,aes(y=y,yend=y,x=xmax-1,xend=xmax),colour=lay$colour,
             arrow=arrow(length=unit(0.25,"cm"),angle=25,type="open"),show.legend=FALSE)
print(g)

tiff(paste0(outpath,"Fig X probability each outcome ",variable_short,".tif"),height=height,width=width,units="in",res=600,compression="lzw")
print(g)
dev.off()

assign(paste0("g_",variable_short),g,envir = .GlobalEnv)
}
col.synergies

multi1_plot(data=multi1.results.biome, run_label=c("Biome - All cases"),variable_short="biome",type=c(" All biodiversity metrics"),ncol=2,graph_limit=c(0.1, 0.2),height=6) # type=c(" All biodiversity metrics", "Abundance","Richness","Richness-Evenness")
multi1_plot(data = multi1.results.practice, run_label=c("Practice - All cases"),variable_short="practice",type=c(" All biodiversity metrics"),height=4.5)
multi1_plot(multi1.results.agrochem, run_label=c("Agrochemicals - All cases"),variable_short="agrochem",type=c(" All biodiversity metrics"),  height=4.5)
multi1_plot(multi1.results.crops, run_label=c("Crops - All cases"),variable_short="crops",type=c(" All biodiversity metrics"),height=6)
multi1_plot(multi1.results.crop_type, run_label=c("Crop type - All cases"),variable_short="crop_type",type=c(" All biodiversity metrics"),height=6)
multi1_plot(multi1.results.region, run_label=c("Region - All cases"),variable_short="region",type=c(" All biodiversity metrics"),height=4.5)
multi1_plot(multi1.results.pests, run_label=c("Pest group - All cases"),variable_short="pests",type=c(" All biodiversity metrics"),height=3)
multi1_plot(multi1.results.b_ground, run_label=c("Ground relation - All cases"),variable_short="b_ground",type=c(" All biodiversity metrics"),height=3)
multi1_plot(multi1.results.taxa, run_label=c("Taxa - All cases"),variable_short="taxa_simp",type=c(" All biodiversity metrics"),height=4.5)
multi1_plot(multi1.results.dev, run_label=c("Development status - All cases"),variable_short="development",type=c(" All biodiversity metrics"),height=3)
multi1_plot(multi1.results.b_measure, run_label=c("Biodiversity metric - All cases"),variable_short="b_measure",type=c(" All biodiversity metrics"),height=3)
multi1_plot(multi1.results.y_measure, run_label=c("Yield metric - All cases"),variable_short="y_measure",type=c(" All biodiversity metrics"),height=4.5)

col.synergies
multi1_plot(data = multi1.d_refLW_practice.result, run_label=c("Practice - All cases"),variable_short="practice_refLW",type=c(" All biodiversity metrics"),height=4.5,col.points=col.synergies[c(1,2,3,5)],ref="Lose B-Win Y")
multi1_plot(data = multi1.d_refWL_taxa.result, run_label=c("Taxa - All cases"),variable_short="taxa_refWL",type=c(" All biodiversity metrics"),height=4.5,col.points=col.synergies[c(1,3:5)],ref="Win B-Lose Y")
multi1_plot(data = multi1.d_refWL_crops.result, run_label=c("Crops - All cases"),variable_short="crops_refWL",type=c(" All biodiversity metrics"),height=4.5,col.points=col.synergies[c(1,3:5)],ref="Win B-Lose Y")
multi1_plot(data = multi1.d_refLL_crops.result, run_label=c("Crops - All cases"),variable_short="crops_refLL",type=c(" All biodiversity metrics"),height=4.5,col.points=col.synergies[c(1:3,5)],ref="Lose B-Lose Y")

height=5.5
width=7.5
  
plot.legend <- get_legend(g_biome+theme(legend.position="bottom")+guides(colour=guide_legend(nrow=1)))

p <- plot_grid(g_biome+theme(legend.position="none"),
               g_practice+theme(legend.position="none"),
               ncol=1,
               align= "hv",axis="lb",
               labels = c("A)","B)"),label_size=10,rel_heights=c(1,1),rel_widths=c(1,1))
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.9,0.1))
p


p <- plot_grid(g_practice+theme(legend.position="none"),
               g_practice_refLW+theme(legend.position="none"),
               ncol=2,
               align= "hv",axis="bl",
               labels = c("A)","B)"),label_size=10,rel_heights=c(1,1),rel_widths=c(1,1))
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.9,0.1))
p

tiff(paste0(outpath,"Fig X probability each outcome practice_cowplot.tif"),height=height,width=width,units="in",res=600,compression="lzw")
print(p)
dev.off()

pdf(paste0(outpath,"Fig X probability each outcome practice_cowplot.pdf"),height=height,width=width)
print(p)
dev.off()

p <- plot_grid(g_crop_type+theme(legend.position="none"),
               g_crops+theme(legend.position="none"),
               ncol=2,
               align= "hv",axis="b",
               labels = c("A)","B)"),label_size=10,rel_heights=c(1,1),rel_widths=c(1,1))
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.92,0.08))
p
tiff(paste0(outpath,"Fig X2 probability each outcome cowplot.tif"),height=height,width=width,units="in",res=600,compression="lzw")
print(p)
dev.off()

p <- plot_grid(g_region +theme(legend.position="none"),
               g_agrochem+theme(legend.position="none"),
               ncol=2,
               align= "hv",axis="b",
               labels = c("A)","B)"),label_size=10,rel_heights=c(1,1),rel_widths=c(1,1))
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.92,0.08))
p
tiff(paste0(outpath,"Fig X3 probability each outcome cowplot.tif"),height=5,width=width,units="in",res=600,compression="lzw")
print(p)
dev.off()

p <- plot_grid(g_region+theme(legend.position="none"),
               g_development+theme(legend.position="none"),
               ncol=2,
               align= "hv",axis="bl",
               labels = c("A)","B)"),label_size=10,rel_heights=c(1,1),rel_widths=c(1,1))
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.92,0.08))
p
tiff(paste0(outpath,"Fig X4 probability each outcome cowplot.tif"),height=height,width=width,units="in",res=600,compression="lzw")
print(p)
dev.off()

p <- plot_grid(g_b_ground +theme(legend.position="none"),
              g_pests+theme(legend.position="none"),
               ncol=1,
               align= "hv",axis="bl",
               labels = c("G)","H)"),label_size=12, rel_heights=c(1,1),rel_widths=c(1,1))
p
p <- plot_grid(g_taxa_simp+theme(legend.position="none"),p,ncol=2,axis="lb",labels=c("F)",""),label_size=12,rel_widths=c(1,1),rel_heights=c(1,1))
p
p <- plot_grid(p,plot.legend,ncol=1,rel_heights=c(0.92,0.08))
p
tiff(paste0(outpath,"Fig X5 probability each outcome cowplot.tif"),height=height,width=width,units="in",res=600,compression="lzw")
print(p)
dev.off()