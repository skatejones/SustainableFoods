### Explore trade-offs between biodiversity and yield in diversified versus simplified farming systems ####
### Data prep ####
# Authors: SJ

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
library(dmetar) # https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/fitting-a-three-level-model.html
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

wd <- readline() #at the prompt, copy and paste your filepath and press enter
D:\02_Bioversity\30_SustainableFoods\R
setwd(wd)

run <- "DIVSIM_yields"
outpath <- "./Results_tradeoffs_wLER/"

col.intervention = c("navy","forestgreen", "seagreen","gold","orange", "lightblue","purple")
col.synergies = c("steelblue","gold","orange","forestgreen")
col.synergies = c("#225ea8","#41b6c4","#a1dab4","#74c476") # "#ffffcc, "#238b45"

# Import data ####
# Info on ID columns:
# C_ID represents rows using the same control across treatments (crop diversity practices), e.g. monoculture versus agroforestry versus crop rotation, where insect richness is the outcome measure
# C_T_ID represents rows with same control and treatment, so rows sharing this ID represent repeated measures (e.g. multiple timesteps in sampling period) 
# Effect_ID = unique effects 
d_bio_divsim <- read.xlsx("Bio_yields_meta_datasets_split.xlsx",sheet="Data_divsim") %>% 
  mutate(Merge_ID = substr(Merge_ID,1,50))

d_bio_divnat <- read.xlsx("Bio_yields_meta_datasets_split.xlsx",sheet="Data_divnat") %>% 
  mutate(Merge_ID = substr(Merge_ID,1,50))


d_full <- read.xlsx("Bio_yields_meta_datasets_split.xlsx",sheet="Data_divsim_yields") %>% 
  mutate(Merge_ID = substr(Merge_ID,1,50))

bioclim <- read.dbf("./GIS/Biodiversity_yield_data_all_T_ecoreg_v2.dbf") %>% 
  mutate(ID = as.character(ID),Location=as.character(Location)) %>%
  mutate(Merge_ID = substr(Merge_ID,1,50))

unsd <- read.csv("D:/02_Bioversity/24_ABD_index/Science/ABDI_Tool/data/UNSD/UNSD Methodology.csv") %>%
  mutate(Country.or.Area = ifelse(M49.Code %in% c(" Hong Kong Special Administrative Region"," Macao Special Administrative Region"),
                                  M49.Code,Country.or.Area)) # Do this so china is not repeated on three rows

# Make functions used in code ####
transf = function(x){
  return((exp(x)-1)*100)
}

# Format data for analysis ####
data = d_full
data = d_bio_divsim
data_formating = function(data,name,run=0,zeros=0,outcome="bio_yield"){
  
  if(outcome %in% c("yield","bio_yield")){
    # remove rows that have no yield data
    data <- data %>% filter(Yield_data==1)
  }
  
  # Without adjusting zero mean or SD to allow their inclusion (so these are excluded).
  # This is to check how diff the results are if keeping zeros
  # instead of adding 0.0001 to all zero values
  if(zeros=="keep_zeros"){
    data = data
  }
  if(zeros!="keep_zeros"){
    if(outcome %in% c("bio","bio_yield")){
      # Add 0.0001 to zeros to allow their inclusion
      data$B_value_C[data$B_value_C == 0] <- 0.0001
      data$B_SD_C[data$B_SD_C == 0] <- 0.0001
      data$B_value_T[data$B_value_T == 0] <- 0.0001
      data$B_SD_T[data$B_SD_T == 0] <- 0.0001
    }
    if(outcome %in% c("yield","bio_yield")){
      data$Yield_value_C[data$Yield_value_C == 0] <- 0.0001
      data$Yield_SD_C[data$Yield_SD_C == 0] <- 0.0001
      data$Yield_value_T[data$Yield_value_T == 0] <- 0.0001
      data$Yield_SD_T[data$Yield_SD_T == 0] <- 0.0001
    }
  }
  
  data <- data %>% left_join(bioclim[,c("Merge_ID", "BIOME_NUM", "BIOME_NAME","REALM")],by=c("Merge_ID")) %>% 
    unique()
  #table(data$BIOME_NAME) # 10 biomes with data points (includes non-LER data)
  #unique(data$BIOME_NAME)
  #check <- data %>% filter(is.na(BIOME_NAME)) %>% select(Merge_ID,Lat_C,Long_C,Lat_T,Long_T,Country,BIOME_NAME)
  
  # Not all the biodiversity metrics are present in the yield-bio dataset - check classes
  #addmargins(table(data$B_measure_group))
  #freq(data$B_measure)
  if(outcome %in% c("bio_yield","bio")){
    data <- data %>% mutate(B_measure_group = ifelse(B_measure %in% c("Activity-density","Abundance","Vistitation frequency"),"Abundance",
                                                     ifelse(B_measure %in% c("Species Richness"),"Richness",
                                                            ifelse(B_measure %in%  c("Species Eveness ","Shannon Eveness Index", "Shannon Index" , "Shannon Index " ,   "Shannon-Wiener Index"),"Shannon diversity","Other")))) 
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
    
    tiff(paste0(outpath,"Funnel plot biodiversity ",name," n=",nrow(data),".tif"),width=6,height=4,units="in",res=120)
    metafor::funnel(x=data$yi_B,vi=data$vi_B,yaxis="sei",atransf=exp,
                    steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=1,pch=20,col="black")
    dev.off()
    
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
    
    #addmargins(table(data$Yield_measure_group))
    #freq(data$Yield_measure_group)
    #addmargins(table(data$System_T,data$Yield_measure_group))
    
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
    
    tiff(paste0(outpath,"Funnel plot yields ",name," n=",nrow(data),".tif"),width=6,height=4,units="in",res=120)
    metafor::funnel(x=data$yi_Y,vi=data$vi_Y,yaxis="sei",atransf=exp,steps=7,ylim=c(0,1.2),at=log(c(0.08,1,10,100)),digits=1,pch=20,col="black")
    dev.off()
    
  }
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
  
  #sort(unique(data$Fertiliser_CT))
  #freq(data$Fertiliser_CT)
  #data$Fertiliser_CT <- factor(data$Fertiliser_CT,levels=c("No","Mixed","Yes", "nd"))
  
  #table(data$Pesticide_C,data$Pesticide_T)
  #table(data$Pesticide_CT)
  #table(data$Fertiliser_C,data$Fertiliser_T)
  #table(data$Fertiliser_chem_C,data$Fertiliser_chem_T)
  #table(data$Fertiliser_CT)
  #addmargins(table(data$Fertiliser_CT,data$Pesticide_CT))
  
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
  #freq(data$Agrochem_CT)
  #addmargins(table(data$Agrochem_CT,data$Pesticide_CT))
  #addmargins(table(data$Agrochem_CT,data$Fertiliser_CT))
  
  data <- data %>% mutate(Crop_type = paste0(Crop_ann_pen_T," ",Crop_woodiness_T))
  data <- data %>% mutate(Crop_type = ifelse(Crop_type %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type))
  
  if(outcome %in% c("bio_yield","bio")){
    data <- data %>% mutate(Taxa_group = ifelse(Taxa_class%in% c("Diplopoda","Annelida"),"Annelids & millipedes",
                                                ifelse(Taxa_class == "Arachnida","Arachnids",
                                                       ifelse(Taxa_class %in% c("Tree","Herb"),"Plants",
                                                              ifelse(Taxa_order %in% c("Coleoptera","Hymenoptera","Lepidoptera"),Taxa_order,Taxa_class)))))%>%
      mutate(Taxa_group = ifelse(Taxa_group =="Insect","Insect (other)",Taxa_group))
    
    #freq(data$Taxa_group)
    #freq(data$Taxa_phylum)
    #freq(data$Taxa_class)
    #table(data$Taxa_order,data$Taxa_group)
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
    mutate(Crop_FAO_C_extra = ifelse(Crop_ann_pen_C == "Perennial" & Crop_woodiness_C %in% c("Liana","Shrub"), paste0(Crop_FAO_C," (","Perennial Shrub/Liana)"),
                                     ifelse(Crop_ann_pen_C == "Perennial", paste0(Crop_FAO_C," (Perennial ",Crop_woodiness_C,")"),paste0(Crop_FAO_C," (",Crop_ann_pen_C," ",Crop_woodiness_C,")")))) %>%
    mutate(Crop_FAO_C_extra = ifelse(Crop_FAO_C_extra %in% c("Fodder (Annual Herb)","Fodder (Mixed or nd Herb)"),"Fodder (Annual Herb)",
                                     ifelse(Crop_FAO_C_extra %in% c("Nuts (Perennial Woody)","Stimulants (Perennial Woody)"),"Nuts/Stimulants (Perennial Woody)", Crop_FAO_C_extra)))
  
  data <- data %>%
    mutate(Crop_FAO_CT = ifelse(Crop_FAO_C==Crop_FAO_T,Crop_FAO_C,
                                ifelse(Crop_FAO_C == "Cereals","Cereals - Other" ,"Mixed"))) %>%
    mutate(Crop_FAO_CT = ifelse(Crop_FAO_CT %in% c("Stimulants","Nuts"),"Nuts/Stimulants",Crop_FAO_CT))
  #table(data$Crop_FAO_C)
  #table(data$Crop_FAO_C,data$Crop_FAO_T)
  #table(data$Crop_FAO_C,data$System_T)
  #table(data$Crop_ann_pen_C,data$System_T)
  #table(data$Crop_C,data$System_T)
  #table(data$Crop_T,data$System_T)
  #table(data$Crop_FAO_C,data$Crop_T)
  #table(data$Crop_FAO_CT) 
  #table(data$Crop_FAO_CT,data$Crop_T) 
  
  data <- data %>%
    mutate(Crop_FAO_CT_extra = ifelse(Crop_ann_pen_C ==  "Perennial" & Crop_ann_pen_T == "Perennial" & Crop_woodiness_C %in% c("Liana","Shrub") & Crop_woodiness_T %in% c("Liana","Shrub"), paste0(Crop_FAO_CT," (","Perennial Shrub/Liana)"),
                                      ifelse(Crop_ann_pen_C == "Perennial"& Crop_ann_pen_T == "Perennial", paste0(Crop_FAO_CT," (Perennial ",Crop_woodiness_C,":",Crop_woodiness_T, ")"),
                                             ifelse(Crop_ann_pen_C =="Annual" & Crop_ann_pen_T == "Annual" & Crop_woodiness_C == Crop_woodiness_T, paste0(Crop_FAO_CT," (Annual ",Crop_woodiness_C, ")"),"Other"))))
  data <- data %>% 
    mutate(Crop_FAO_CT_extra = ifelse(Crop_FAO_CT_extra == "Fodder (Perennial Herb:Herb)","Fodder (Perennial Herb)",
                                      ifelse(Crop_FAO_CT_extra =="Fruits (Perennial Tree:Tree)","Fruits (Perennial Tree)",
                                             ifelse(Crop_FAO_CT_extra =="Nuts/Stimulants (Perennial Tree:Tree)","Nuts/Stimulants (Perennial Tree)",
                                                    ifelse(Crop_FAO_CT_extra %in% c("Fodder (Annual Herb)","Mixed/other (Annual Herb)","Mixed (Annual Herb)"),"Other",Crop_FAO_CT_extra)))))
  #table(data$Crop_FAO_CT_extra) 
  #table(data$Crop_FAO_CT_extra) 
  #table(data$Crop_FAO_CT_extra,data$Yield_measure_group) 
  #table(data$Crop_T,data$System_T)
  #ggplot(data,aes(x=Crop_T,fill=Yield_measure_group))+geom_bar()+coord_flip()
  
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
  
  #table(data$System_T_action)
  #table(data$Crop_ann_pen_C_action)
  #table(data$Crop_woodiness_C_action)
  
  #check <- data %>% filter(System_T %in% c("Embedded natural")) %>% select(ID,Crop_C,Crop_T,crops_all_scientific_T, System_raw_C,System_raw_T,System_details_T)
  #check <- data %>% filter(System_T %in% c("Associated plant species")) %>% select(ID,Crop_C,Crop_T,crops_all_scientific_T, System_raw_C,System_raw_T,System_details_T)
  
  data <- data %>%
    mutate(Crop_type_CT = ifelse(Crop_ann_pen_C ==  "Perennial" & Crop_ann_pen_T == "Perennial"& Crop_woodiness_C == Crop_woodiness_T, paste0("Perennial ",Crop_woodiness_C),
                                 ifelse(Crop_ann_pen_C =="Annual" & Crop_ann_pen_T == "Annual" & Crop_woodiness_C == Crop_woodiness_T, paste0("Annual ",Crop_woodiness_C),"Other"))) %>%
    mutate(Crop_type_CT = ifelse(Crop_type_CT %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type_CT))
  
  #table(data$Crop_type_CT)
  #data[data$Crop_type_CT %in% c("Other"),c("ID","Crop_C","Crop_T","Crop_ann_pen_C","Crop_ann_pen_T", "Crop_woodiness_C","Crop_woodiness_T", "Crop_type_CT","Crop_FAO_C","Crop_FAO_T","System_details_C", "System_details_T")]
  
  # ID 643 uses perennial alfalfa as control but harvests it annually in the experiment, so classified here as an annual
  # ID 740 uses cassava as control whichis a perennial shrub producing edible roots after 8-9 months and grown as an annual generally.
  # Harvested from 36 days to 576 days in the experiment, so classified here as an annual
  data <- data %>%
    mutate(Crop_type_CT = ifelse(Crop_type_CT == "Other" & Crop_T == "Wheat-Alfalfa","Annual Herb",
                                 ifelse(Crop_type_CT == "Other" & Crop_T =="Maize-Cassava","Annual Shrub",Crop_type_CT)))
  #table(data$Crop_type_CT)
  
  data <- data %>%
    mutate(Crop_type_C = paste0(Crop_ann_pen_C," ",Crop_woodiness_C)) %>%
    mutate(Crop_type_C = ifelse(Crop_type_C %in% c("Perennial Liana","Perennial Shrub"),"Perennial Shrub/Liana",Crop_type_C))
  
  #table(data$Crop_type_C)         
  #check <- data %>% filter(Crop_type_C %in% c("Mixed or nd Herb","Perennial Herb")) %>% select(ID,Crop_C,Crop_T,crops_all_scientific_T, System_raw_C,System_details_C,System_raw_T,System_details_T)
  #check <- data %>% filter(Crop_C %in% c("Cotton")) %>% select(ID,Crop_C,Crop_T,crops_all_scientific_T, System_raw_C,System_details_C,System_raw_T,System_details_T,Crop_ann_pen_C,Crop_woodiness_C)
  data <- data %>%
    mutate(Crop_type_C = ifelse(Crop_type_C == "Mixed or nd Herb" & Crop_T == "Wheat-Alfalfa","Annual Herb",Crop_type_C))
  
  #table(data$Crop_type_CT,data$Crop_type_C)
  #table(data$Agrochem_CT) # no data = 99 (wLER) or 194 (full dataset)
  
  # Add potential random effects columns ####
  # make unique identifier for shared controls or treatments
  #data <- data %>% 
  #  mutate(ID_C_bio = paste(ID,Comparison_ID_C,Taxa,Functional_group,B_measure,yi_B,sep="_"),
  #         ID_T_bio = paste(ID,Comparison_ID_T,Taxa,Functional_group, B_measure,yi_B,sep="_"),
  #         ID_C_yield = paste(ID,Comparison_ID_C,Yield_measure,yi_Y,sep="_"),
  #         ID_T_yield = paste(ID,Comparison_ID_T,Yield_measure,yi_Y,sep="_"))
  
  # make unique identifier for C-T pairs repeatedly sampled for biodiversity
  #data <- data %>% mutate(ID_CT_bio = paste(ID,Comparison_ID_C,Comparison_ID_T,System_raw_C,System_raw_T,Taxa,Functional_group, B_measure,sep="_"))
  
  data <- data %>% mutate(Region.Name = as.factor(Region.Name),Biome = BIOME_NAME)
  
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
  if(outcome %in% c("yield")){
    
    data <- data %>% 
      mutate(se_Y = sqrt(vi_Y)) %>% mutate(ci.lb_Y=yi_Y-1.96*se_Y,ci.ub_Y=yi_Y+1.96*se_Y)
    
    data <- data %>% select(ID,Comparison_ID_C,Comparison_ID_T,
                            yi_Y,vi_Y,se_Y,ci.lb_Y,ci.ub_Y, 
                            #Crop_type_CT,Crop_FAO_CT,Crop_FAO_CT_extra,Crop_ann_pen_C_action,Crop_woodiness_C_action,
                            Crop_C,Crop_T,Crop_FAO_C,Crop_type_C,Crop_ann_pen_C,Crop_woodiness_C, 
                            System_C,System_details_C,System_T,System_T_action,System_details_T, 
                            Agrochem_CT,Pesticide_CT,Fertiliser_CT,
                            Yield_measure_group,
                            DevelopmentStatus,Region.Name,Sub.region.Name, Continent,BIOME_NAME,Biome,
                            Lat_original_C,Long_original_C,Lat_original_T,Long_original_T,Lat_C,Long_C,Lat_T,Long_T,
                            Country,
                            #ID_C_yield,ID_T_yield,
                            Validity_yield_overall,Validity_biodiversity_overall,
                            Validity_location,
                            Validity_N_yield_C,Validity_N_yield_T,
                            Validity_time_C,Validity_time_T) %>% unique()
    
  }
  if(outcome %in% c("bio")){
    data <- data %>% 
      mutate(se_B = sqrt(vi_B)) %>% mutate(ci.lb_B=yi_B-1.96*se_B,ci.ub_B=yi_B+1.96*se_B)
    
    data <- data %>% select(ID,Effect_ID, Comparison_ID_C,Comparison_ID_T,
                            yi_B,vi_B,se_B,ci.lb_B,ci.ub_B,
                            #Crop_type_CT,Crop_FAO_CT,Crop_FAO_CT_extra,Crop_ann_pen_C_action,Crop_woodiness_C_action,
                            Crop_C,Crop_T,Crop_FAO_C,Crop_type_C,Crop_ann_pen_C,Crop_woodiness_C, 
                            System_C,System_details_C,System_T,System_T_action,System_details_T, 
                            Agrochem_CT,Pesticide_CT,Fertiliser_CT,
                            Taxa_group,Taxa_class,Taxa_order,Taxa_details,#Taxa_group_map
                            B_measure_group,B_ground,Pest_group,
                            DevelopmentStatus,Region.Name,Sub.region.Name, Continent,BIOME_NAME,Biome,
                            Lat_original_C,Long_original_C,Lat_original_T,Long_original_T,Lat_C,Long_C,Lat_T,Long_T,
                            Country,
                            #ID_C_bio,ID_T_bio,ID_CT_bio
                            Validity_biodiversity_overall,
                            Validity_location,
                            Validity_N_biodiversity_C,Validity_N_biodiversity_T,
                            Validity_time_C,Validity_time_T)
    
  }
  
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
      mutate(Synergies_sig_simp = ifelse(Synergies_sig %in% c("Win B-Equiv Y"  , "Equiv B-Win Y" ,  "Win B-Win Y"),"Win B-Win Y",
                                         ifelse(Synergies_sig %in% c("Lose B-Lose Y" ,"Lose B-Equiv Y" , "Equiv B-Lose Y" ), "Lose B-Lose Y",Synergies_sig))) %>%
      mutate(Synergies_sig_binary = ifelse(((yi_B >0 & ci.lb_B>0) | (yi_B <0 & ci.ub_B<0)) & ((yi_Y >0 & ci.lb_Y>0)|(yi_Y<0 & ci.ub_Y<0)),1,0 ))
    
    
    data <- data %>% select(ID,Effect_ID,Comparison_ID_C,Comparison_ID_T,
                            yi_Y,vi_Y,se_Y,ci.lb_Y,ci.ub_Y, yi_Y_pc, 
                            yi_B,vi_B,se_B,ci.lb_B,ci.ub_B,yi_B_pc,
                            Synergies,Synergies_sig,
                            #Crop_type_CT,Crop_FAO_CT,Crop_FAO_CT_extra, Crop_ann_pen_C_action,Crop_woodiness_C_action,
                            Crop_C,Crop_T,Crop_FAO_C,Crop_type_C,Crop_ann_pen_C,Crop_woodiness_C, 
                            System_C,System_details_C,System_T,System_T_action,System_details_T, 
                            Agrochem_CT,Pesticide_CT,Fertiliser_CT,
                            Taxa_group,Taxa_class,Taxa_order,Taxa_details, #Taxa_group_map,
                            B_measure_group,B_ground,Pest_group,
                            Yield_measure_group,
                            DevelopmentStatus,Region.Name,Sub.region.Name, Continent,BIOME_NAME,Biome,
                            Lat_original_C,Long_original_C,Lat_original_T,Long_original_T,Lat_C,Long_C,Lat_T,Long_T,
                            Country,
                            # ID_C_bio,ID_T_bio,ID_C_yield,ID_T_yield,ID_CT_bio
                            Validity_yield_overall,Validity_biodiversity_overall,
                            Validity_location,
                            Validity_N_yield_C,Validity_N_yield_T,
                            Validity_N_biodiversity_C,Validity_N_biodiversity_T,
                            Validity_time_C,Validity_time_T)
    
  }
  
  
  return(data)
}

d_bioyield <- data_formating(data=d_full,"d_bioyield",run="LER",zeros=0,outcome="bio_yield")
d_bioyield_LERplus <- data_formating(data=d_full,"d_bioyield_LERplus",run=0,zeros=0,outcome="bio_yield")
d_bioyield_zeros <- data_formating(data=d_full,"d_bioyield_zeros",run="LER",zeros="keep_zeros",outcome="bio_yield")

d_bio_divsim <- data_formating(data=d_bio_divsim,name="d_bio_divsim",run=0,zeros=0,outcome="bio")
d_bio_divnat <- data_formating(data=d_bio_divnat,name="d_bio_divnat",run=0,zeros=0,outcome="bio")

d_yield <- data_formating(data=d_full,"d_yield",run="LER",zeros=0,outcome="yield")
d_yield <- d_yield %>% mutate(Effect_ID = 1:nrow(d_yield)) %>% mutate(Effect_ID = as.character(Effect_ID))
#id <- d_bioyield[,c("Effect_ID","ID","ID_C_yield","ID_T_yield","Yield_measure_group")] %>% unique()
#d_yield <- d_yield %>% left_join(id,by=c("ID","ID_C_yield","ID_T_yield","Yield_measure_group"))
d_yield_LERplus <- data_formating(data=d_full,name="d_yield_LERplus",run=0,zeros=0,outcome="yield")
d_yield_LERplus <- d_yield_LERplus %>% mutate(Effect_ID = 1:nrow(d_yield_LERplus)) %>% mutate(Effect_ID = as.character(Effect_ID))

d_yield_zeros <- data_formating(data=d_full,"d_yield_zeros",run="LER",zeros="keep_zeros",outcome="yield")
d_yield_zeros <- d_yield_zeros %>% mutate(Effect_ID = 1:nrow(d_yield_zeros)) %>% mutate(Effect_ID = as.character(Effect_ID))


write.xlsx(list("d_bio_divsim" = d_bio_divsim,
                "d_bio_divnat" = d_bio_divnat,
                "d_yield" = d_yield),
           #d_validity.location.high,d_validity.N.high,d_validity.time.high),
           file=paste0(outpath,"Bio_yields_effects_website.xlsx"), overwrite=TRUE)

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
  assign(paste0(name,"_validity.N.high"),d_validity.N.high,envir = .GlobalEnv )
  assign(paste0(name,"_validity.location.high"),d_validity.location.high,envir = .GlobalEnv )
  assign(paste0(name,"_validity.time.high"),d_validity.time.high,envir = .GlobalEnv )
}

data_sensitivity(data=d_bioyield,name="d_bioyield")
data_sensitivity(data=d_yield,name="d_yield",run=0)

# Export data ####
write.xlsx(list("d_full" = d_full,
                "d_bioyield" = d_bioyield,
                "d_bioyield_LERplus" = d_bioyield_LERplus,
                "d_bioyield_zeros" = d_bioyield_zeros,
                "d_bioyield_validity" = d_bioyield_validity.high,
                "d_yield" = d_yield, 
                "d_yield_zeros" = d_yield_zeros,
                "d_yield_validity" = d_yield_validity.high),
           #d_validity.location.high,d_validity.N.high,d_validity.time.high),
           file=paste0(outpath,"Bio_yields_effects.xlsx"), overwrite=TRUE)


### Bubble chart of N effect sizes per outcome, diversification system and crop ####
d <- d_bioyield

ggplot(d,aes(x=Taxa_group,fill=Synergies))+
  geom_bar(position="dodge",width=0.6)+
  scale_fill_manual(values=col.synergies)+
  scale_colour_manual(values=c("black","grey70","green"))+
  ylab("Number of effect sizes")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#ggsave(paste0(outpath,"Fig yield biodiversity synergies by taxa group"),device="tiff",dpi=150)

ggplot(d,aes(x=Continent,fill=Synergies))+
  geom_bar(position="dodge",width=0.6)+
  scale_fill_manual(values=col.synergies)+
  scale_colour_manual(values=c("black","grey70","green"))+
  ylab("Number of effect sizes")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

d_freq <- data.frame(table(System_T=d$System_T,Crop_FAO_T=d$Crop_FAO_T,Crop_T=d$Crop_T,Synergies=d$Synergies)) %>% filter(Freq>0)

ggplot(d_freq,aes(x=Synergies,y=reorder(Crop_T,desc(Crop_T)),size=Freq,colour=System_T))+
  geom_point(alpha=0.7)+
  scale_size(range = c(1.6, 7), name="Number of effect sizes")+
  ylab("")+xlab("")+
  scale_colour_manual(values=col.intervention,name="Intervention")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))+
  facet_grid(rows=vars(Crop_FAO_T), 
             #switch="x",
             #strip.position = "right", 
             space="free",
             scales = "free")+
  theme(panel.spacing = unit(0, "lines"),
        legend.title=element_text(face="bold",size=10),
        legend.text=element_text(size=10),
        strip.text.y=element_text(angle=0,hjust=0,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave(paste0(outpath,"Fig yield biodiversity synergies by crop and intervention.tif"),device="tiff",dpi=150)

### Scatterplot with four quadrants showing outcomes per diversification system ####
labels <- labels_d_bioyield
g <- ggplot(d_bioyield,aes(x=yi_B,y=yi_Y,colour=System_T))+
  geom_point(size=2,alpha=0.5)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  ylab("Log RR for Yield")+
  xlab("Log RR for Biodiversity")+
  #ylim(c(-5,5))+
  #xlim(c(-5,5))+
  geom_text(aes(x=10,y=10,label=labels[which(labels$Synergies=="Win B-Win Y"),c("Freq_pc")]),inherit.aes=FALSE,colour="forestgreen")+
  geom_text(aes(x=10,y=-10,label=labels[which(labels$Synergies=="Win B-Lose Y"),c("Freq_pc")]),inherit.aes=FALSE,colour="black")+
  geom_text(aes(x=-10,y=10,label=labels[which(labels$Synergies=="Lose B-Win Y"),c("Freq_pc")]),inherit.aes=FALSE,colour="black")+
  geom_text(aes(x=-10,y=-10,label=labels[which(labels$Synergies=="Lose B-Lose Y"),c("Freq_pc")]),inherit.aes=FALSE,colour="red")+
  scale_colour_brewer(type="qual",palette="Set1",name="")+
  #scale_fill_manual(values=col.intervention,name="Intervention")+
  #scale_colour_manual(values=c("grey","black","orange"),name="")+
  #scale_shape_manual(name="")+
  theme_bw()+
  guides(shape=guide_legend(title="",title.position = "top"))
g

ggsave(paste0(outpath,"Fig yield biodiversity synergies all yield effects n=",nrow(d),".tif"),device="tiff",dpi=150)
ggsave(paste0(outpath,"Fig yield biodiversity synergies all yield effects n=",nrow(d),".pdf"),device="pdf",dpi=150)
#png(paste0("Fig yield biodiversity synergies all yield effects n=",nrow(d), ".png"), width=900,height=500,res=110)

table(d_bioyield$Synergies, d_bioyield$Synergies_sig)
table(d_bioyield$Synergies_sig)
g <- ggplot(d_bioyield,aes(x=yi_B,y=yi_Y,colour=Synergies_sig))+
  geom_point(size=2,alpha=0.5)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  ylab("Log RR for Yield")+
  xlab("Log RR for Biodiversity")+
  #ylim(c(-5,5))+
  #xlim(c(-5,5))+
  scale_colour_brewer(type="qual",palette="Set1",name="")+
  #scale_fill_manual(values=col.intervention,name="Intervention")+
  #scale_colour_manual(values=c("grey","black","orange"),name="")+
  #scale_shape_manual(name="")+
  theme_bw()+
  guides(shape=guide_legend(title="",title.position = "top"))
g

# Check for collinearity between variables ####
# Using VIF ###
# https://www.statology.org/variance-inflation-factor-r/ 
# if the largest value is >10, indicates collinearity (Bowerman and O'Connell 1990 in https://www.bookdown.org/chua/ber642_advanced_regression/multinomial-logistic-regression.html#checking-assumptionl-multicollinearity)
model_Y_full = rma.mv(yi_Y, vi_Y, mods=~System_T+Crop_FAO_CT+Crop_type_CT+Agrochem_CT+BIOME_NAME, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, three VIFs are >10 but all <17 (cereals - other, BIOME deserts and zeric shrublands, Biome Tropical grasslands). 
# Intercropping 9.8, associated plant species 9.6.
# Probable interaction between System_T and Crop_FAO_CT.

model_Y_full = rma.mv(yi_Y, vi_Y, mods=~System_T+Crop_FAO_extra+Agrochem_CT+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, all VIFs are <7

model_Y_full = rma.mv(yi_Y, vi_Y, mods=~System_T+Crop_FAO_CT_extra+Agrochem_CT+BIOME_NAME, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, all VIFs are

vif.rma(model_Y_full,digits=2)
vif <- capture.output(vif.rma(model_Y_full,digits=2))
#vif <- data.frame(vif.rma(model_Y_full,digits=2)) # not working!
vif <- data.frame(vif)
vif$covariate <- row.names(vif)
colnames(vif)[1] <-"value"
barplot(vif[,c("value")], main = "VIF Values", horiz = TRUE, col = "lightblue")
abline(v = 10, lwd = 3, lty = 2)
g <- ggplot(vif,aes(y=reorder(covariate,value),x=value))+
  geom_col(fill="grey60",width=0.8)+scale_x_continuous(expand=c(0,0),limits=c(0,7),breaks=seq(0,7,1))+labs(y="Variables",x="VIF")+ 
  theme_bw()+theme(text=element_text(size=10),panel.grid.major.x=element_line(size=0.5,colour="grey80"),panel.border=element_blank())
g
png(paste0(outpath,"vif yield.png"),width=6,height=4,units="in", res=300)
g
dev.off()

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_FAO_extra+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, two VIFs are >10, though only just (intercropping, use of pesticides and fertilisers) 

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, all VIFs are <9

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_FAO_extra+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, all VIFs are <10

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_type_CT+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, all VIFs are <10
# this one makes theoretical sense to me, except maybe change Continent to Biome

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_type_CT+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])
# with this model (excluding agrochem = no data), one VIF is >10 (pesticide and fertiliser use, followed by intercropping)

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_type_CT+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])
# with this model, one VIF is >10 (pesticide and fertiliser use, followed by intercropping, no pesticides, associated plant species)

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_FAO_CT_extra+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d)
# with this model, two VIFs are >10 (intercropping, pesticide and fertiliser use)

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_type_CT:Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])
# with this model, one VIF is >10 but only just (intercropping)
# this one makes the most theoretical sense to me, except maybe change Continent to Biome

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_type_CT:Agrochem_CT+Taxa_group+B_measure_group+BIOME_NAME, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])

# with this model, several vifs very high

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_FAO_CT_extra:Agrochem_CT+Taxa_group+B_measure_group+BIOME_NAME, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T:Crop_type_CT+Crop_FAO_CT:Agrochem_CT+Taxa_group+B_measure_group+BIOME_NAME, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])
# several vifs very high in taxa group, biome, system T....

table(d$Crop_type_CT,d$Agrochem_CT) # nearly half of cells are empty...
table(d$Crop_FAO_CT_extra,d$Agrochem_CT)
table(d$Crop_FAO_CT,d$Agrochem_CT)

model_B_full = rma.mv(yi_B, vi_B, mods=~System_T+Crop_FAO_CT_extra+Agrochem_CT+Taxa_group+B_measure_group+Continent, 
                      random=list(~1|Effect_ID,~1|ID), intercept=TRUE,method = "ML",test="t", tdist=TRUE, data=d[which(d$Agrochem_CT != "No data"),])
# model cannot run

vif.rma(model_B_full,digits=2)
vif <- data.frame(vif.rma(model_B_full,digits=2))
vif$covariate <- row.names(vif)
colnames(vif)[1] <-"value"
barplot(vif[,c("value")], main = "VIF Values", horiz = TRUE, col = "lightblue")
abline(v = 10, lwd = 3, lty = 2)
g <- ggplot(vif,aes(y=reorder(covariate,value),x=value))+
  geom_col(fill="grey60",width=0.8)+scale_x_continuous(expand=c(0,0),limits=c(0,12),breaks=seq(0,12,1))+labs(y="Variables",x="VIF")+ 
  theme_bw()+theme(text=element_text(size=10),panel.grid.major.x=element_line(size=0.5,colour="grey80"),panel.border=element_blank())
g


# Using Cramer's V ###
# 0 means no association, 1 means complete association
d <- d_bioyield
cramersv = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}

chisq.test(x=d$BIOME_NAME,y=d$Region.Name) 
cramersv(x=d$BIOME_NAME,y=d$Region.Name) 

chisq.test(x=d$System_T,y=d$Crop_FAO_C) 
cramersv(x=d$System_T,y=d$Crop_FAO_C) 
table(x=d$System_T,y=d$Crop_FAO_C)

chisq.test(x=d$System_T,y=d$Crop_type_C) 
cramersv(x=d$System_T,y=d$Crop_type_C) 

chisq.test(x=d$System_T,y=d$Agrochem_CT) 
cramersv(x=d$System_T,y=d$Agrochem_CT) 

chisq.test(x=d$Agrochem_CT,y=d$Crop_type_C) 
cramersv(x=d$Agrochem_CT,y=d$Crop_type_C)

chisq.test(x=d$Crop_type_C,y=d$Crop_FAO_C) 
cramersv(x=d$Crop_type_C,y=d$Crop_FAO_C) 
table(x=d$Crop_type_C,y=d$Crop_FAO_C) # strongly associated

chisq.test(x=d_yield$Crop_type_C,y=d_yield$Crop_FAO_C) 
cramersv(x=d_yield$Crop_type_C,y=d_yield$Crop_FAO_C) 
table(x=d_yield$Crop_type_C,y=d_yield$Crop_FAO_C) # strongly associated

chisq.test(x=d$Crop_FAO_C,y=d$Crop_woodiness_C) 
cramersv(x=d$Crop_FAO_C,y=d$Crop_woodiness_C)  # strongly associated
table(x=d$Crop_FAO_C,y=d$Crop_woodiness_C)

chisq.test(x=d$Crop_ann_pen_C,y=d$Crop_woodiness_C) 
cramersv(x=d$Crop_ann_pen_C,y=d$Crop_woodiness_C)  # moderately associated

chisq.test(x=d$System_T,y=d$Agrochem_CT) 
cramersv(x=d$System_T,y=d$Agrochem_CT)
table(x=d$System_T,y=d$Agrochem_CT)

chisq.test(x=d$Crop_FAO_C,y=d$Pesticide_CT) # 
cramersv(x=d$Crop_FAO_C,y=d$Pesticide_CT) #  

chisq.test(x=d$Crop_FAO_C,y=d$Fertiliser_CT) # 
cramersv(x=d$Crop_FAO_C,y=d$Fertiliser_CT) #

chisq.test(x=d$Crop_FAO_C,y=d$Agrochem_CT) # 
cramersv(x=d$Crop_FAO_C,y=d$Agrochem_CT) #

chisq.test(x=d$B_measure_group,y=d$Taxa_group) #
cramersv(x=d$B_measure_group,y=d$Taxa_group) # 
table(x=d$B_measure_group,y=d$Taxa_group)

chisq.test(x=d$Pesticide_CT,y=d$Fertiliser_CT) # 
cramersv(x=d$Pesticide_CT,y=d$Fertiliser_CT) #  
addmargins(table(x=d$Pesticide_CT,y=d$Fertiliser_CT))

chisq.test(x=d$DevelopmentStatus,y=d$Region.Name) # 
cramersv(x=d$DevelopmentStatus,y=d$Region.Name) # strongly associated
addmargins(table(x=d$DevelopmentStatus,y=d$Region.Name))