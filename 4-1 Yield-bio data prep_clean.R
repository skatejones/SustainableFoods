### Explore trade-offs between biodiversity and yield in diversified versus simplified farming systems ####
### Data prep ####

# Authors: SJ
# Last updated 12/12/2022

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
#library(dmetar) 
library(forestplot)
library(nnet)
library("DescTools")
library(forcats)
library(MCMCglmm)
library(lmtest) # for likelihood ratio tests
library(summarytools)
library(dunn.test)
library(car)
library("report") # nice model summaries
library(foreign) # for reading dbase files

#### Set file paths ####

wd <- readline() #at the prompt, copy and paste your filepath and press enter

setwd(wd)

#### Set parameters ####
outpath <- "./Results_tradeoffs_wLER/"

#### Instructions to access data needed for this code ####

# Dataset 1_sources.xlsx and Dataset_2_outcomes.xlsx, download from Harvard DataVerse at https://doi.org/10.7910/DVN/XIDI1X 
# UNSD Methodology.csv, download from UN Stats at https://unstats.un.org/unsd/methodology/m49/overview/
# Biodiversity_yield_data_all_T_ecoreg_v2.dbf, available from corresponding author, or reproducible using the identity tool in ArcGIS to join the metadataset (using Lat_T and Long_T) to Dinerstein et al. (2017) ecoregion dataset (at: https://ecoregions.appspot.com/) and save as a shapefile

#### Import data ####
# Info on ID columns:
# C_ID represents rows using the same control across treatments (crop diversity practices), e.g. monoculture versus agroforestry versus crop rotation, where insect richness is the outcome measure
# C_T_ID represents rows with same control and treatment, so rows sharing this ID represent repeated measures (e.g. multiple timesteps in sampling period) 
# Effect_ID = unique effects 

bioclim <- read.dbf("./Biodiversity_yield_data_all_T_ecoreg_v2.dbf") %>% 
  mutate(ID = as.character(ID),Location=as.character(Location)) %>%
  mutate(Merge_ID = substr(Merge_ID,1,50))

unsd <- read.csv("./UNSD Methodology.csv") %>%
  mutate(Country.or.Area = ifelse(M49.Code %in% c(" Hong Kong Special Administrative Region"," Macao Special Administrative Region"),
                                  M49.Code,Country.or.Area)) # Do this so china is not repeated on three rows

articles <- read.xlsx("Dataset 1_sources.xlsx",sheet="Literature_screened") 

data_meta <- read.xlsx("Dataset_2_outcomes.xlsx",sheet="Data")
data_meta <- data_meta %>%
  mutate(Taxa_group = Taxa_order) %>%
  mutate(Taxa_group = ifelse(Taxa_group == "Aves", "Birds",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Arachnida","Araneae"), "Arachnids",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group == "Chiroptera", "Bats",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Nematoda","Annelida"), "Nematodes and annelids",Taxa_group))%>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Lagomorpha"), "Mammals",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Arthropods", "Insect n.s","Insects ns/oth",  "Collembola",  "Diptera",   "Dermaptera","Homoptera",       "Isoptera",   "Isopoda",      "Orthoptera", "Neuroptera",       "Neoptera",  "Thysanoptera"),  "Arthropods other",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Plant vascular non-woody" ,"Plant nonvascular" , "Plants n.s.","Plants n.s","Plant vascular n.s."), "Plants other",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Plant vascular woody"), "Plants woody",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Lepidoptera"), "Butterflies, moths",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Hymenoptera"), "Sawflies, wasps, bees, ants",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Hemiptera"), "True bugs",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Coleoptera"), "Beetles",Taxa_group)) %>%
  mutate(Taxa_group = ifelse(Taxa_group %in% c("Myriapoda"), "Millipedes, centipedes, symphylans",Taxa_group)) 
data_meta <- data_meta %>%
  mutate(Merge_ID = paste(ID,Experiment_stage,Taxa_group,Taxa_details,Functional_group, B_measure,Country, sep="_")) %>% unique() %>%
  mutate(Yield_data = ifelse(!is.na(Yield_value_T)&!is.na(Yield_SD_T) & !is.na(Yield_N_T) & 
                               Yield_value_T !="nd"& Yield_SD_T !="nd" & Yield_N_T !="nd",1,0)) %>%
  mutate(Yield_data = ifelse(Yield_measure!="LER" & (is.na(Yield_value_C) | Yield_value_C=="nd" | is.na(Yield_SD_C)|Yield_SD_C=="nd"),0,Yield_data)) %>%
  mutate(Yield_data = ifelse(is.na(Yield_data),0,Yield_data))

d_full <- data_meta[which(data_meta$System_C %in% c("Monoculture","Simplified other")),] %>%
  mutate(Merge_ID = substr(Merge_ID,1,50)) 
d_full$Effect_ID <- factor(1:nrow(d_full)) 
d_full <- d_full %>% filter(Yield_data==1)

# correct error in 1036 comparison ID_T column 
d_full <- d_full %>%
  mutate(Comparison_ID_T = ifelse(ID == "1036" & System_details_T == "S18: 18 months tree growth with Sesbania fallow","T3",
                                  ifelse(ID == "1036" & System_details_T == "N18: 18 months natural regrowth of vegetation fallow without cultivation","T4",Comparison_ID_T)))
d_full <- d_full %>% 
  mutate(C_T_ID = paste0(Comparison_ID_C,"_",Comparison_ID_T),
         Experiment_stage = ifelse(is.na(Experiment_stage),99,as.numeric(Experiment_stage))) %>%
  select(ID,C_T_ID,Experiment_stage,System_C,System_details_C,System_T,System_details_T, everything())

# Make functions used in code ####
transf = function(x){
  return((exp(x)-1)*100)
}

# Format data for analysis ####

data_formating = function(data,name,run=0,zeros=0,outcome="bio_yield"){

  print(paste0("# of unique biodiversity comparisons at start = ", 
               data %>% select(ID,C_T_ID,System_C,System_T,B_value_C,B_SD_C,B_N_C, B_value_T,B_SD_T,B_N_T,
                  Taxa_group,Functional_group,B_measure,Experiment_stage) %>% unique() %>% nrow())) 
  
  print(paste0("# of unique yield comparisons at start = ", 
               data %>% select(ID,C_T_ID,System_C,System_T,Yield_value_C,Yield_SD_C,Yield_N_C, Yield_value_T,Yield_SD_T,Yield_N_T,
                  Crop_C,Crop_T,Yield_measure) %>% unique() %>% nrow()))

  if(outcome %in% c("yield","bio_yield")){
    # remove rows that have no yield data
    data <- data %>% filter(Yield_data==1)
  }
  
  # Without adjusting zero mean or SD to allow their inclusion (so these are excluded).
  # This is to check how diff the results are if keeping zeros
  # instead of adding 0.01 to all zero values
  if(zeros=="keep_zeros"){
    data = data
  }
  if(zeros!="keep_zeros"){
    if(outcome %in% c("bio","bio_yield")){
      # Add small increment to zeros to allow their inclusion
      data$B_value_C[data$B_value_C == 0] <- 0.01
      data$B_SD_C[data$B_SD_C == 0] <- 0.01
      data$B_value_T[data$B_value_T == 0] <- 0.01
      data$B_SD_T[data$B_SD_T == 0] <- 0.01
    }
    if(outcome %in% c("yield","bio_yield")){
      data$Yield_value_C[data$Yield_value_C == 0] <- 0.01
      data$Yield_SD_C[data$Yield_SD_C == 0] <- 0.01
      data$Yield_value_T[data$Yield_value_T == 0] <- 0.01
      data$Yield_SD_T[data$Yield_SD_T == 0] <- 0.01
    }
  }
  
  data <- data %>% left_join(bioclim[,c("Merge_ID", "BIOME_NUM", "BIOME_NAME","REALM")],by=c("Merge_ID")) %>% 
    unique()

  if(outcome %in% c("bio_yield","bio")){
    data <- data %>% mutate(B_measure_group = ifelse(B_measure %in% c("Activity-density","Abundance","Vistitation frequency"),"Abundance",
                                                     ifelse(B_measure %in% c("Species Richness"),"Richness",
                                                            ifelse(B_measure %in%  c("Species Eveness ","Shannon Eveness Index", "Shannon Index" , "Shannon Index " ,   "Shannon-Wiener Index"),"Shannon diversity","Other")))) 
    
        data <- data %>% mutate(B_measure_group = ifelse(B_measure %in% c("Abundance", "Activity-density", "Area under disease progress curve (AUDPC)",
                                                                      "Colonization Percent","Colonization Percent ", "Vistitation frequency"),"Abundance",
                                                     ifelse(B_measure %in% c("Chao1 Index","Species Richness"),"Richness",
                                                            ifelse(B_measure %in% c("Pielou Index","Shannon-Wiener Index","Shannon Index","Shannon Index ","Species Eveness","Species Eveness ","Simpson Index","Simpson Index "),"Richness-Evenness",B_measure))))
    
    #freq(data$B_measure_group)
    #table(data$B_measure,data$B_measure_group)
    #check <- data %>% filter(B_measure =="Area under disease progress curve (AUDPC)")
    
    # Compute effect sizes ####
    data <- escalc(measure="ROM", # ROM = log ratio, SMD = standardised mean difference (Hedge's g)
                   data=data,
                   m2i=B_value_C,m1i=B_value_T, # mean
                   sd2i=B_SD_C,sd1i=B_SD_T, # SD
                   n2i=B_N_C,n1i=B_N_T, # sample size
                   var.names=c("yi_B","vi_B"), 
                   vtype="LS", # assume sampling variances are not homoscedastic (not the same across groups)
                   append=T) # append results to input data
    
    # remove inf and NA values
    # meta-analysis doesn't work if var=0, so remove these too
    data <- data %>% 
      filter(is.finite(yi_B)) %>%
      filter(is.finite(vi_B)) %>%
      filter(vi_B != 0) %>%
      filter(!(is.na(yi_B)))
    
    #tiff(paste0(outpath,"Funnel plot biodiversity ",name," n=",nrow(data),".tif"),width=6,height=4,units="in",res=120)
    metafor::funnel(x=data$yi_B,vi=data$vi_B,yaxis="sei",atransf=exp,
                    steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=1,pch=20,col="black")
    #dev.off()
    
  }
  
  # Classify response metrics into a few groups
  # ID 1647 kg yield represents average mean dry weight per 1m2.
  # ID 744 kg yield represents average above ground collards biomass in 120m2 plots
  # each of the others checked and classified below
  if(outcome %in% c("bio_yield","yield")){
    data <- data %>% mutate(Yield_measure_group = ifelse(Yield_measure %in% c("kg/ha", "kg/ha (cattle)","Mg/ha dry matter",
                                                                              "kg per plot","kg/plot","kg per row m", "kg per row",
                                                                              "kg per 3.6m of row" ,"kg" ),"Mass per area",
                                                         ifelse(Yield_measure %in% c("g/plant","fruits per branch" , 
                                                                                     "fruit set (proportion of flowers on branches that produce fruit)",
                                                                                     "kg per plant" , "g per plant","g per cluster",
                                                                                     "kg per plants","kg nut DM per tree" ,
                                                                                     "ear weight(g)","nuts/tree","fruit set (%)"),"Mass or count per unit",
                                                                ifelse(Yield_measure %in% c("plant weight, g","plant weight in g (dry)","whole plant dry weight (g)"),"Biomass",Yield_measure))))
    
    #check <- data %>% filter(!Yield_measure_group %in% c("LER","Mass per area","Mass or count per unit")) %>%
    #  select(ID,System_T,Crop_C,Crop_T,Yield_measure,Sampling_unit_C, Sampling_unit_T,Notes)
    
    #check <- data %>% filter(Yield_measure %in% c("nuts/tree")) %>%
    #  select(ID,System_T,Crop_C,Crop_T,Yield_measure,Sampling_unit_C, Sampling_unit_T,Notes)
    
    #check <- data %>% filter(!Yield_measure_group %in% c("LER") & System_T == "Intercropping") %>%
    #  select(ID,System_T,Crop_C,Crop_T,Yield_measure,Sampling_unit_C, Sampling_unit_T,Notes)
    
    # Compute effect sizes ####
    data <- escalc(measure="ROM", # ROM = log ratio, SMD = standardised mean difference (Hedge's g)
                   data=data,
                   m2i=Yield_value_C,m1i=Yield_value_T, # mean
                   sd2i=Yield_SD_C,sd1i=Yield_SD_T, # SD
                   n2i=Yield_N_C,n1i=Yield_N_T, # sample size
                   var.names=c("yi_Y","vi_Y"), 
                   vtype="LS", # assume sampling variances are not homoscedastic (not the same across groups)
                   append=T) # append results to input data
    
    data <- data %>%
      mutate(yi_Y = ifelse(Yield_measure =="LER",log(Yield_value_T),yi_Y),
             vi_Y = ifelse(Yield_measure =="LER",Yield_SD_T,vi_Y)) 
    #check_vi_Y = data$Yield_SD_C^2/(data$Yield_value_C^2*data$Yield_N_C) + data$Yield_SD_T^2/(data$Yield_value_T^2*data$Yield_N_T))
    
    # remove inf and NA values
    # meta-analysis doesn't work if var=0, so remove these too
    data <- data %>% 
      filter(is.finite(yi_Y)) %>%
      filter(is.finite(vi_Y)) %>%
      filter(vi_Y != 0) %>% 
      filter(!(is.na(yi_Y)))
    
    # Remove intercropping yield data points that don't use LER #
    if(run=="LER"){
      data <- data %>% filter(!(System_T == "Intercropping" & Yield_measure != "LER")) # 774 from 44 studies compared to 1214 from 57 studies
    }
    
    #tiff(paste0(outpath,"Funnel plot yields ",name," n=",nrow(data),".tif"),width=6,height=4,units="in",res=120)
    metafor::funnel(x=data$yi_Y,vi=data$vi_Y,yaxis="sei",atransf=exp,steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=1,pch=20,col="black")
    #dev.off()
    
  }
  
  print(paste0("# of unique biodiversity comparisons after removing no LER = ", 
               data %>% select(ID,C_T_ID,System_C,System_T,B_value_C,B_SD_C,B_N_C, B_value_T,B_SD_T,B_N_T,
                               Taxa_group,Functional_group,B_measure,Experiment_stage) %>% unique() %>% nrow())) 
  
  print(paste0("# of unique yield comparisons after removing no LER = ", 
               data %>% select(ID,C_T_ID,System_C,System_T,Yield_value_C,Yield_SD_C,Yield_N_C, Yield_value_T,Yield_SD_T,Yield_N_T,
                               Crop_C,Crop_T,Yield_measure) %>% unique() %>% nrow()))

  
  if(outcome %in% c("bio_yield")){
    data <- data %>%
      mutate(yi_B_pc = transf(yi_B),
             yi_Y_pc = transf(yi_Y))
    
    data <- data %>%
      mutate(Synergies = ifelse(yi_B >=0 & yi_Y >=0,"Win B-Win Y",
                                ifelse(yi_B>=0,"Win B-Lose Y",
                                       ifelse(yi_Y>=0,"Lose B-Win Y","Lose B-Lose Y"))),
             Synergies = factor(Synergies,levels=c("Lose B-Lose Y","Win B-Lose Y","Lose B-Win Y","Win B-Win Y"))) 
    
    
    labels <- data.frame(table(data$Synergies)) %>% rename(Synergies = Var1) %>%
      mutate(Freq_pc = paste0(Synergies," ",round(Freq/sum(Freq)*100,1),"%"))
    assign(paste0("labels_",name),labels,envir=.GlobalEnv)
  }
  
  
    # Make moderator variables #### 
  data <- data %>% mutate(Fertiliser_chem_C = ifelse(Fertiliser_chem_C %in% c("no"),"No",
                                                     ifelse(Fertiliser_chem_C %in% c("yes"),"Yes", Fertiliser_chem_C))) %>%
    mutate(Fertiliser_chem_T = ifelse(Fertiliser_chem_T %in% c("no"),"No",
                                      ifelse(Fertiliser_chem_T %in% c("yes"),"Yes",Fertiliser_chem_T))) %>%
    mutate(Fertiliser_CT = paste0(as.character(Fertiliser_chem_C),":",as.character(Fertiliser_chem_T))) %>%
    mutate(Fertiliser_CT = ifelse(Fertiliser_CT %in% c("nd:nd"),"nd",
                                  ifelse(Fertiliser_CT %in% c("nd:No","No:nd","No:No","No:NA"),"No",
                                         ifelse(Fertiliser_CT %in% c("nd:Yes","Yes:nd","Yes:Yes","Yes:NA"),"Yes",
                                                ifelse(Fertiliser_CT %in% c("No:Yes","Yes:No"),"Mixed",Fertiliser_CT)))))
  
  
  data <- data %>% mutate(Pesticide_CT = paste0(as.character(Pesticide_C),":",as.character(Pesticide_T)))
  data$Pesticide_CT[data$Pesticide_CT %in% c("nd:nd","nd:NA")] <- "nd"
  data$Pesticide_CT[data$Pesticide_CT %in% c("nd:No","No:nd","No:No","No:NA")] <- "No"
  data$Pesticide_CT[data$Pesticide_CT %in% c("nd:Yes","Yes:nd","Yes:Yes")] <- "Yes"
  data$Pesticide_CT[data$Pesticide_CT %in% c("No:Yes","Yes:No")] <- "Mixed"
  data$Pesticide_CT <- factor(data$Pesticide_CT,levels=c("No","Mixed","Yes", "nd"))
  
  data <- data  %>% mutate(Agrochem_CT = ifelse(Pesticide_CT =="Yes" & Fertiliser_CT=="Yes","Pesticides and fertilizers",
                                                ifelse(Fertiliser_CT=="Mixed"|Pesticide_CT=="Mixed","Mixed",
                                                       ifelse(Pesticide_CT %in% c("Yes"), "Pesticides only",
                                                              ifelse(Fertiliser_CT %in% c("Yes"),"Fertilizers only",
                                                                     ifelse(Pesticide_CT == "No" | Pesticide_CT == "No","No pesticides and/or no fertilizers","No data")))))) %>%
    mutate(Agrochem_CT = factor(Agrochem_CT)) %>%
    mutate(Agrochem_CT = fct_relevel(Agrochem_CT,"No data",after=Inf))
  
  data <- data  %>% mutate(Agrochem_CT = ifelse(Pesticide_CT =="Yes" | Fertiliser_CT=="Yes","Agrochemicals",
                                                ifelse(Fertiliser_CT=="Mixed"|Pesticide_CT=="Mixed","Mixed",
                                                       ifelse(Pesticide_CT == "No" | Fertiliser_CT =="No","No agrochemicals","No data")))) %>%
    mutate(Agrochem_CT = factor(Agrochem_CT)) %>%
    mutate(Agrochem_CT = fct_relevel(Agrochem_CT,"No data",after=Inf))
  
  data <- data %>% mutate(Tillage_CT = paste0(as.character(Soil_management_C),":",as.character(Soil_management_T)))
  data$Tillage_CT[data$Tillage_CT %in% c("nd:nd","Burning:Burning","nd:NA","NA:NA","NA:nd")] <- "nd"
  data$Tillage_CT[data$Tillage_CT %in% c("nd:No Tillage","No Tillage:nd","No Tillage:No Tillage","NA:No Tillage")] <- "No"
  data$Tillage_CT[data$Tillage_CT %in% c("nd:Tillage","Tillage:nd","Tillage:Tillage","NA:Tillage")] <- "Yes"
  data$Tillage_CT[data$Tillage_CT %in% c("No Tillage:Tillage","Tillage:No Tillage")] <- "Mixed"
  data$Tillage_CT <- factor(data$Tillage_CT,levels=c("No","Mixed","Yes", "nd"))
  
  data <- data %>% mutate(Pest_group = as.character(Functional_group)) %>%
    mutate(Pest_group = ifelse(Pest_group %in% c("Pests"),Pest_group, "Other"))
  
  data <- data %>% mutate(Crop_type = paste0(Crop_ann_pen_T," ",Crop_woodiness_T))
  data <- data %>% mutate(Crop_type = ifelse(Crop_type %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type))
  
  if(outcome %in% c("bio_yield","bio")){
    data <- data %>% mutate(Taxa_group = ifelse(Taxa_class%in% c("Diplopoda","Annelida"),"Annelids & millipedes",
                                                ifelse(Taxa_class == "Arachnida","Arachnids",
                                                       ifelse(Taxa_class %in% c("Tree","Herb"),"Plants",
                                                              ifelse(Taxa_order %in% c("Coleoptera","Hymenoptera","Lepidoptera"),Taxa_order,Taxa_class)))))%>%
      mutate(Taxa_group = ifelse(Taxa_group =="Insect","Insect (other)",Taxa_group))
    
  }
  
  # Add development status
  data <- data %>% mutate(Country = ifelse(Country =="United Kingdom","United Kingdom of Great Britain and Northern Ireland",Country)) %>%
    left_join(unsd,by=c("Country"="Country.or.Area"))
  data <- data %>% mutate(DevelopmentStatus = Developed...Developing.Countries)
  
  # Add continent
  data <- data %>%
    mutate(Country = as.character(Country)) %>%
    mutate(Country = ifelse(Country %in% c("United States of America"),"USA",
                            ifelse(Country %in% c("United Kingdom","United Kingdom of Great Britain and Northern Ireland"),"UK",Country))) %>%
    mutate(Continent= if_else(Country == "Argentina"| Country =="Brazil" | Country =="Ecuador" | Country =="Uruguay" |  Country == "Peru"| Country =="Colombia", "South America",
                              if_else(Country == "Canada"|Country =="USA" |Country =="Mexico", "North America", 
                                      if_else(Country =="Germany"| Country =="Belgium"|Country == "France"|Country == "Swiss"| Country == "Sweden"| Country == "Poland"| Country =="Spain"| Country =="Italy"| Country =="Portugal"| Country ==  "UK" | Country =="Finland" | Country =="Hungary"| Country =="Switzerland"| Country== "Netherlands", "Europe",
                                              if_else(Country == "China" | Country =="Turkey"| Country =="India"| Country =="Indonesia"| Country =="Viet Nam"| Country =="Japan"| Country =="Malaysia"| Country == "Israel"|Country =="Philippines", "Asia",
                                                      if_else(Country =="Cameroon"| Country =="Egypt"| Country == "Ghana"| Country == "Nigeria"| Country =="South Africa"| Country =="Kenya"| Country =="Malawi" | Country =="Uganda"|Country == "Ethiopia"| Country =="Sao Tome and Principe"| Country =="Benin"| Country == "Zambia", "Africa",
                                                              if_else(Country =="Costa Rica"| Country == "Panama"| Country =="Guatemala"| Country =="Nicaragua", "North America",
                                                                      if_else(Country =="Dominican Republic" | Country == "Jamaica", "North America",
                                                                              if_else(Country == "New Zealand"|Country =="Australia", "Oceania",
                                                                                      Country)))))))))
  
  
  
  # Moderator preparation ####
  data <- data %>%
    mutate(Crop_FAO_extra = ifelse(Crop_ann_pen_T == "Perennial" & Crop_woodiness_T %in% c("Tree","Liana","Shrub"), paste0(Crop_FAO_T," (","Perennial Woody)"),
                                   ifelse(Crop_ann_pen_T == "Perennial", paste0(Crop_FAO_T," (Perennial ",Crop_woodiness_T,")"),paste0(Crop_FAO_T," (",Crop_ann_pen_T," ",Crop_woodiness_T,")")))) %>%
    mutate(Crop_FAO_extra = ifelse(Crop_FAO_extra %in% c("Fodder (Annual Herb)","Mixed/other (Annual Herb)"),"Other (Annual Herb)",
                                   ifelse(Crop_FAO_extra %in% c("Nuts (Perennial Woody)","Stimulants (Perennial Woody)"),"Nuts/Stimulants (Perennial Woody)", Crop_FAO_extra)))
  
  #table(data$Crop_FAO_extra) 
  
  data %>% filter(Crop_FAO_extra == "Cereals (Annual Herb)" & System_T == "Agroforestry") %>% select(ID,Crop_C,System_C,System_raw_C,Crop_T, System_T,System_raw_T)
  
  #table(data$Crop_ann_pen_C,data$Crop_woodiness_C)
  #table(data$Crop_FAO_C)
  data <- data %>%
    mutate(Crop_FAO_C = ifelse(Crop_FAO_C %in% c("Nuts","Stimulants"),"Nuts/Stimulants",Crop_FAO_C))
  data <- data %>%
    mutate(Crop_FAO_T = ifelse(Crop_FAO_T %in% c("Nuts","Stimulants"),"Nuts/Stimulants",Crop_FAO_T))
  
  data <- data %>%
    mutate(Crop_FAO_C_extra = ifelse(Crop_ann_pen_C == "Perennial" & Crop_woodiness_C %in% c("Liana","Shrub"), paste0(Crop_FAO_C," (","Perennial Shrub/Liana)"),
                                     ifelse(Crop_ann_pen_C == "Perennial", paste0(Crop_FAO_C," (Perennial ",Crop_woodiness_C,")"),paste0(Crop_FAO_C," (",Crop_ann_pen_C," ",Crop_woodiness_C,")")))) %>%
    mutate(Crop_FAO_C_extra = ifelse(Crop_FAO_C_extra %in% c("Fodder (Annual Herb)","Fodder (Mixed or nd Herb)"),"Fodder (Annual Herb)",
                                     ifelse(Crop_FAO_C_extra %in% c("Nuts (Perennial Woody)","Stimulants (Perennial Woody)"),"Nuts/Stimulants (Perennial Woody)", Crop_FAO_C_extra)))
  
  data <- data %>%
    mutate(Crop_FAO_CT = ifelse(Crop_FAO_C==Crop_FAO_T,Crop_FAO_C,
                                ifelse(Crop_FAO_C == "Cereals","Cereals - Other" ,"Mixed")))
  
  data <- data %>%
    mutate(Crop_FAO_CT_extra = ifelse(Crop_ann_pen_C ==  "Perennial" & Crop_ann_pen_T == "Perennial" & Crop_woodiness_C %in% c("Liana","Shrub") & Crop_woodiness_T %in% c("Liana","Shrub"), paste0(Crop_FAO_CT," (","Perennial Shrub/Liana)"),
                                      ifelse(Crop_ann_pen_C == "Perennial"& Crop_ann_pen_T == "Perennial", paste0(Crop_FAO_CT," (Perennial ",Crop_woodiness_C,":",Crop_woodiness_T, ")"),
                                             ifelse(Crop_ann_pen_C =="Annual" & Crop_ann_pen_T == "Annual" & Crop_woodiness_C == Crop_woodiness_T, paste0(Crop_FAO_CT," (Annual ",Crop_woodiness_C, ")"),"Other"))))
  data <- data %>% 
    mutate(Crop_FAO_CT_extra = ifelse(Crop_FAO_CT_extra == "Fodder (Perennial Herb:Herb)","Fodder (Perennial Herb)",
                                      ifelse(Crop_FAO_CT_extra =="Fruits (Perennial Tree:Tree)","Fruits (Perennial Tree)",
                                             ifelse(Crop_FAO_CT_extra =="Nuts/Stimulants (Perennial Tree:Tree)","Nuts/Stimulants (Perennial Tree)",
                                                    ifelse(Crop_FAO_CT_extra %in% c("Fodder (Annual Herb)","Mixed/other (Annual Herb)","Mixed (Annual Herb)"),"Other",Crop_FAO_CT_extra)))))

  #check <- data %>% group_by(Crop_T,Yield_measure_group) %>% mutate(Yield_measure_n = n()) %>% arrange(Yield_measure_n)
  #ggplot(check, aes(x = Crop_T, y = Yield_measure_n,colour=Yield_measure_group)) +geom_segment(aes(x = Crop_T, y = 0, xend = Crop_T, yend = Yield_measure_n)) +geom_point(shape=16,size=3)+coord_flip()
  
  data <- data %>%
    mutate(System_T_action = ifelse(System_T %in% c("Agroforestry"),"add trees",
                                    ifelse(System_T %in% c("Crop rotation","Cultivar mixture","Intercropping"),"add annual crop species/varieties",
                                           ifelse(System_T == "Embedded natural","add natural habitat",
                                                  ifelse(System_T =="Associated plant species","add cover, trap or green manure crops",
                                                         ifelse(System_T == "Combined practices","integrate crop-animal production",System_T)))))) %>%
    mutate(Crop_ann_pen_C_action = paste0(Crop_ann_pen_C," crop: ",System_T_action)) %>%
    mutate(Crop_woodiness_C_action = ifelse(Crop_woodiness_C %in% c("Tree","Liana","Shrub"),paste0("Woody crop: ",System_T_action),paste0(Crop_woodiness_C," crop: ",System_T_action)))
  
  data <- data %>%
    mutate(Crop_type_CT = ifelse(Crop_ann_pen_C ==  "Perennial" & Crop_ann_pen_T == "Perennial"& Crop_woodiness_C == Crop_woodiness_T, paste0("Perennial ",Crop_woodiness_C),
                                 ifelse(Crop_ann_pen_C =="Annual" & Crop_ann_pen_T == "Annual" & Crop_woodiness_C == Crop_woodiness_T, paste0("Annual ",Crop_woodiness_C),"Other"))) %>%
    mutate(Crop_type_CT = ifelse(Crop_type_CT %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type_CT))
  
  #table(data$Crop_type_CT)
  #data[data$Crop_type_CT %in% c("Other"),c("ID","Crop_C","Crop_T","Crop_ann_pen_C","Crop_ann_pen_T", "Crop_woodiness_C","Crop_woodiness_T", "Crop_type_CT","Crop_FAO_C","Crop_FAO_T","System_details_C", "System_details_T")]
  
  # ID 643 uses perennial alfalfa as control but harvests it annually in the experiment, so classified here as an annual
  # ID 740 uses cassava as control which is a perennial shrub producing edible roots after 8-9 months and grown as an annual generally.
  # Harvested from 36 days to 576 days in the experiment, so classified here as an annual
  data <- data %>%
    mutate(Crop_type_CT = ifelse(Crop_type_CT == "Other" & Crop_T == "Wheat-Alfalfa","Annual Herb",
                                 ifelse(Crop_type_CT == "Other" & Crop_T =="Maize-Cassava","Annual Shrub",Crop_type_CT)))
  #table(data$Crop_type_CT)
  
  data <- data %>%
    mutate(Crop_type_C = paste0(Crop_ann_pen_C," ",Crop_woodiness_C)) %>%
    mutate(Crop_type_C = ifelse(Crop_type_C %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type_C))
  
  #table(data$Crop_type_C)         
  data <- data %>%
    mutate(Crop_type_C = ifelse(Crop_type_C == "Mixed or nd Herb" & Crop_T == "Wheat-Alfalfa","Annual Herb",Crop_type_C))
  
  data <- data %>% mutate(Region.Name = as.factor(Region.Name),Biome = BIOME_NAME)
  names(data)
  data  <- data %>% mutate(Biome = ifelse(as.character(BIOME_NAME) %in% c("Tropical & Subtropical Coniferous Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Moist Broadleaf Forests"),"Tropical & Subtropical Forests",as.character(BIOME_NAME))) %>%
    mutate(Biome = ifelse(as.character(BIOME_NAME) %in% c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests"),"Temperate Forests",Biome))
  
  # Fix agrochemical and biome labels and make taxa_simp variable ####
  data <- data  %>% mutate(Agrochem_CT = as.character(Agrochem_CT)) %>%
    mutate(Agrochem_CT = ifelse(Agrochem_CT =="Pesticides and/or fertilizers","Agrochemicals",
                                ifelse(Agrochem_CT =="No pesticides and/or fertilizers","No agrochemicals",as.character(Agrochem_CT)))) 
  data <- data %>%  mutate(Agrochem_CT = factor(Agrochem_CT,levels=c("Agrochemicals","No agrochemicals","Mixed","No data")))
  
  data <- data %>% mutate(Biome = as.character(Biome)) %>% 
    mutate(Biome = ifelse(Biome=="Tropical & Subtropical Grasslands, Savannas & Shrublands","Tropical & Subtropical Non-forest",
                          ifelse(Biome == "Temperate Grasslands, Savannas & Shrublands","Temperate Non-forest",
                                 ifelse(Biome =="Mediterranean Forests, Woodlands & Scrub","Mediterranean",
                                        ifelse(Biome=="Montane Grasslands & Shrublands","Montane",
                                               ifelse(Biome=="Deserts & Xeric Shrublands" ,"Deserts",
                                                      ifelse(Biome =="Boreal Forests/Taiga", "Boreal",Biome)))))))
  data <- data %>% mutate(Biome = factor(Biome))
  
  
  data <- data %>% mutate(Biome_simp = ifelse(Biome %in% c("Boreal","Deserts","Montane"),"Other",as.character(Biome)))

  data <- data %>% mutate(Taxa_group_simp = ifelse(Taxa_group %in% c("Annelids & millipedes","Arachnids","Arthropods","Coleoptera","Hymenoptera","Insect (other)","Lepidoptera"),"Invertebrates",ifelse(Taxa_group %in% c("Birds","Mammals"),"Vertebrates",as.character(Taxa_group))))
  data$Taxa_group_simp <- fct_relevel(data$Taxa_group_simp,"Invertebrates")

  # Remove unnecessary columns####
  
  if(outcome %in% c("bio_yield")){
    data <- data %>% 
      mutate(se_Y = sqrt(vi_Y),
             se_B = sqrt(vi_B)) %>%
      mutate(ci.lb_Y=yi_Y-1.96*se_Y,ci.ub_Y=yi_Y+1.96*se_Y,ci.lb_B=yi_B-1.96*se_B,ci.ub_B=yi_B+1.96*se_B) %>%
      mutate(Synergies_sig = ifelse(yi_B >0 & ci.lb_B>0 & yi_Y >0 & ci.lb_Y>0,"Win B-Win Y",
                                    ifelse(yi_B >0 & ci.lb_B >0 & (yi_Y <0 & ci.ub_Y<0),"Win B-Lose Y",
                                           ifelse(yi_B >0 & ci.lb_B >0,"Win B-Equiv Y",
                                                  ifelse(yi_B <0 & ci.ub_B <0 & yi_Y>0 & ci.lb_Y>0,"Lose B-Win Y",
                                                         ifelse(yi_Y >0 & ci.lb_Y>0,"Equiv B-Win Y",
                                                                ifelse((yi_B<0 & ci.ub_B<0) & (yi_Y<0 & ci.ub_Y<0),"Lose B-Lose Y",
                                                                       ifelse((yi_B<0 & ci.ub_B<0) & (yi_Y<=0 & ci.ub_Y>0),"Lose B-Equiv Y",
                                                                              ifelse((yi_B<=0 & ci.ub_B>0) & (yi_Y<0 & ci.ub_Y<0),"Equiv B-Lose Y","Equiv B-Equiv Y"))))))))) %>%
      mutate(Synergies_withEquiv = ifelse(yi_B >0 & yi_Y >0,"Win B-Win Y",
                                          ifelse(yi_B>0 & yi_Y<0,"Win B-Lose Y",
                                                 ifelse(yi_B>0 & yi_Y==0,"Win B-Equiv Y",
                                                        ifelse(yi_B<0 & yi_Y >0,"Lose B-Win Y",
                                                               ifelse(yi_B<0 & yi_Y == 0,"Lose B-Equiv Y",
                                                                      ifelse(yi_B==0 & yi_Y>0,"Equiv B-Win Y",
                                                                             ifelse(yi_B==0 & yi_Y<0,"Equiv B-Lose Y",
                                                                                    ifelse(yi_B==0 & yi_Y==0,"Equiv B-Equiv Y",
                                                                                           ifelse(yi_B<0 & yi_Y<0,"Lose B-Lose Y","check")))))))))) %>%
      mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ,  "Win B-Win Y"),"Win B-Win Y",
                                         ifelse(Synergies_sig %in% c("Lose B-Lose Y" ,"Lose B-Equiv Y" , "Equiv B-Lose Y" ), "Lose B-Lose Y",Synergies_sig))) %>%
      mutate(Synergies_sig_binary = ifelse(((yi_B >0 & ci.lb_B>0) | (yi_B <0 & ci.ub_B<0)) & ((yi_Y >0 & ci.lb_Y>0)|(yi_Y<0 & ci.ub_Y<0)),1,0 ))
    
    
    data <- data %>% select(ID,Effect_ID,Comparison_ID_C,Comparison_ID_T,Experiment_stage,
                            B_value_C,B_SD_C,B_N_C,B_value_T,B_SD_T,B_N_T,
                            Yield_value_C,Yield_SD_C,Yield_N_C,Yield_value_T,Yield_SD_T,Yield_N_T,
                            yi_Y,vi_Y,se_Y,ci.lb_Y,ci.ub_Y, yi_Y_pc, 
                            yi_B,vi_B,se_B,ci.lb_B,ci.ub_B,yi_B_pc,
                            Synergies,Synergies_sig,
                            System_C,System_details_C,System_T,System_details_T,System_T_action, 
                            Crop_C,Crop_T,Crop_FAO_C,Crop_type_C,Crop_ann_pen_C,Crop_woodiness_C, Crop_FAO_T,Crop_ann_pen_T,Crop_woodiness_T, 
                            Crop_type_CT,Crop_FAO_CT,#Crop_FAO_CT_extra,Crop_ann_pen_C_action,Crop_woodiness_C_action,
                            Agrochem_CT,Pesticide_CT,Fertiliser_CT,Tillage_CT,Time_state_C,Time_state_T,Farm_size,
                            Taxa_group,Taxa_class,Taxa_order,Taxa_details, #Taxa_group_map,
                            B_measure,B_measure_group,B_ground,Pest_group,
                            Yield_measure,Yield_measure_group,
                            DevelopmentStatus,Region.Name,Sub.region.Name, Continent,BIOME_NAME,Biome,Biome_simp,
                            Lat_original_C,Long_original_C,Lat_original_T,Long_original_T,Lat_C,Long_C,Lat_T,Long_T,
                            Country,
                            # ID_C_bio,ID_T_bio,ID_C_yield,ID_T_yield,ID_CT_bio
                            Validity_yield_overall,Validity_biodiversity_overall,
                            Validity_location,
                            Validity_N_yield_C,Validity_N_yield_T,
                            Validity_N_biodiversity_C,Validity_N_biodiversity_T,
                            Validity_time_C,Validity_time_T,
                            Authors,Title,Year)
    
  }
  
  return(data)
}

# run function and remove simplified other controls (unless system is embedded natural): 

d_bioyield <- data_formating(data=d_full,"d_bioyield",run="LER",zeros=0,outcome="bio_yield") %>% filter(!(System_C == "Simplified other" & System_T != "Embedded natural")) 

# Prepare data for sensitivity analyses ####

data_sensitivity = function(data,name,run="bioyield"){
  
  if(run=="bioyield"){
    d_validity.time.high <- data %>%
      filter(Validity_time_C==1 & Validity_time_T == 1)
    
    d_validity.location.high <- data %>%
      filter(Validity_location ==1)
    
    d_validity.N.high <- data %>%
      filter(Validity_N_yield_C==1 & Validity_N_yield_T == 1 & Validity_N_biodiversity_C ==1 & Validity_N_biodiversity_T==1)
    
    data <- data %>%
      mutate(Validity_time_C = ifelse(Validity_time_C=="nd",NA,Validity_time_C),
             Validity_time_T = ifelse(Validity_time_T=="nd",NA,Validity_time_T),
             Validity_location = ifelse(Validity_location =="nd",NA,Validity_location)) %>%
      mutate(Validity_location = as.numeric(Validity_location)) %>%
      mutate(Validity_time_CT = as.numeric(Validity_time_C) + as.numeric(Validity_time_T),
             Validity_N_yield_CT = Validity_N_yield_C+ Validity_N_yield_T,
             Validity_N_bio_CT = Validity_N_biodiversity_C + Validity_N_biodiversity_T) %>%
      mutate(Validity_time_CT = ifelse(is.na(Validity_time_CT),Validity_time_CT,
                                       ifelse(Validity_time_CT == 2,1,0)),
             Validity_N_yield_CT = ifelse(is.na(Validity_N_yield_CT),Validity_N_yield_CT,
                                          ifelse(Validity_N_yield_CT==2,1,0)),
             Validity_N_bio_CT = ifelse(is.na(Validity_N_bio_CT),Validity_N_bio_CT,
                                        ifelse(Validity_N_bio_CT == 2,1,0))) %>%
      mutate(Validity_yield_overall = Validity_time_CT + Validity_N_yield_CT + Validity_location,
             Validity_biodiversity_overall = Validity_time_CT + Validity_N_bio_CT + Validity_location) %>%
      mutate(Validity_yield_overall = ifelse(is.na(Validity_yield_overall),Validity_yield_overall,
                                             ifelse(Validity_yield_overall>=2,1,0)),
             Validity_biodiversity_overall = ifelse(is.na(Validity_biodiversity_overall),Validity_biodiversity_overall,
                                                    ifelse(Validity_biodiversity_overall>=2,1,0)))
    
    d_validity.high <- data %>% filter(Validity_yield_overall ==1 & Validity_biodiversity_overall ==1)
  }
  else{
    d_validity.time.high <- data %>%
      filter(Validity_time_C==1 & Validity_time_T == 1)
    
    d_validity.location.high <- data %>%
      filter(Validity_location ==1)
    
    d_validity.N.high <- data %>%
      filter(Validity_N_yield_C==1 & Validity_N_yield_T == 1)
    
    data <- data %>%
      mutate(Validity_time_C = ifelse(Validity_time_C=="nd",NA,Validity_time_C),
             Validity_time_T = ifelse(Validity_time_T=="nd",NA,Validity_time_T),
             Validity_location = ifelse(Validity_location =="nd",NA,Validity_location)) %>%
      mutate(Validity_location = as.numeric(Validity_location)) %>%
      mutate(Validity_time_CT = as.numeric(Validity_time_C) + as.numeric(Validity_time_T),
             Validity_N_yield_CT = Validity_N_yield_C+ Validity_N_yield_T) %>%
      mutate(Validity_time_CT = ifelse(is.na(Validity_time_CT),Validity_time_CT,
                                       ifelse(Validity_time_CT == 2,1,0)),
             Validity_N_yield_CT = ifelse(is.na(Validity_N_yield_CT),Validity_N_yield_CT,
                                          ifelse(Validity_N_yield_CT==2,1,0))) %>%
      mutate(Validity_yield_overall = Validity_time_CT + Validity_N_yield_CT + Validity_location) %>%
      mutate(Validity_yield_overall = ifelse(is.na(Validity_yield_overall),Validity_yield_overall,
                                             ifelse(Validity_yield_overall>=2,1,0)))
    
    d_validity.high <- data %>% filter(Validity_yield_overall ==1)
  }
  
  assign(paste0(name,"_validity.high"),d_validity.high,envir = .GlobalEnv )
}

data_sensitivity(data=d_bioyield,name="d_bioyield")

# Export data ####
write.xlsx(list("d_bioyield" = d_bioyield,
                "d_bioyield_validity" = d_bioyield_validity.high),
           file=paste0(outpath,"Bio_yields_effects.xlsx"), overwrite=TRUE)

# Export list of articles included in trade-off analysis 
articles_included_tradeoff <- d_bioyield %>% group_by_at(c("ID","Country")) %>% summarise(n_cases = n(),n_biodiversity = n()) %>%
  left_join(
    d_bioyield %>% select(ID,Comparison_ID_C,Comparison_ID_T,Yield_value_T,Yield_SD_T,Yield_N_T,Crop_C) %>% unique() %>% 
      group_by_at("ID") %>% summarise(n_yield = n()),by="ID") %>% left_join(
        articles %>% select(ID,Authors,Year,Title,Source.title, DOI)) %>% ungroup() %>%
  select(Authors,Year,Title,Source.title, DOI,Country,n_cases)#n_biodiversity,n_yield)

write.xlsx(list("articles_included" = articles_included_tradeoff), file=paste0(outpath,"Bio_yields_paper_articles_included.xlsx"))

