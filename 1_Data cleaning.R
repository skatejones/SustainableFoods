# Code to clean and standardise data extracted from primary articles
# for use in a meta-analysis of diversified farming effects on biodiversity and yields

library(cairoDevice)
library(ggplot2)
library(funModeling) 
library(tidyverse) 
library(Hmisc)
library(reshape2)
library(data.table)
library(gridExtra)
#install.packages("BiocManager")
#BiocManager::install("EBImage") 
library(EBImage)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(magrittr) # for piping (%>%)
library(stringr) # for pattern searching
library(viridis)
library(rgeos)
library(rworldxtra)
library(rworldmap)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gtable)
library(maps)

# Data preparation --------------------------------------------------------

### Import data ####

data <- read.xlsx("../Biodiversity/Meta-analysis biodiversity data 2021-06-11.xlsx",
                  sheet="DATABASE LOCAL",skipEmptyRows = TRUE)

articles <- read.xlsx("../Biodiversity/Meta-analysis ABD-biodiversity article list 2021-02-24.xlsx",
                      sheet="Articles_local",skipEmptyRows = TRUE)

# Clean article data for export and get stats for prisma flowchart ####
names(articles)
articles <- articles %>% 
  dplyr::rename("Exclusion_reason" = "Database.local.-.exclusion.reason", 
         "Inclusion_yes_no" = "Database.local.-.included?.(Yes/No/To.do)")

sort(unique(articles$Exclusion_reason))
articles <- articles %>% 
  mutate(Inclusion_yes_no = ifelse(Inclusion_yes_no == "yes","Yes",
                                   ifelse(Inclusion_yes_no %in% c("no","No"),"No",Inclusion_yes_no)),
         Exclusion_reason = ifelse(Inclusion_yes_no %in% c("No"),Exclusion_reason,NA))%>%
  mutate(Exclusion_reasion_pico = ifelse(Exclusion_reason %in% c("Compare forest plantations" ,"Compare natural vs. monoculture","Compare one plant species in different settings","Management system comparison only" ,"Natural land comparison only"),"Unsuitable intervention",
                ifelse(Exclusion_reason %in% c("Cropping system description only","Irrelevant ","No comparison"  ),"Unsuitable population",
                       ifelse(Exclusion_reason %in% c("Unsuitable reference system used"),"Unsuitable comparator",
                              ifelse(Exclusion_reason %in% c("Different Output ","Effect on yield only"   ),"Unsuitable outcomes",
                                     ifelse(Exclusion_reason %in% c( "It's a meta-analysis or review" ,"No conducted in the field" ,"Qualitative ","Secondary data","Unpublished studies" ),"Unsuitable context",
                                            ifelse(Exclusion_reason %in% c("Author contacted by AS 20/01/2021 to confirm meaning of error bars. The author did not reply. Decided to include assuming error bars are SE, based on information reported in Table 4"),NA,Exclusion_reason)))))))
sort(unique(articles$Exclusion_reasion_pico))
addmargins(table(articles$Exclusion_reasion_pico,articles$Inclusion_yes_no))
addmargins(table(articles$Inclusion_yes_no))

#sort(unique(articles$Article_source))
df_status(articles)
articles <- articles %>%
  mutate(ID = as.character(ID),
         #Title = str_to_lower(Title, locale = "en")
         Article_source =as.character(Article_source),
         Article_source = if_else(str_detect(Article_source, paste(c("From references","Beillouin"), collapse = "|"),negate=FALSE), "Previous meta-analysis", 
                                  if_else(str_detect(Article_source, paste(c("Scopus", "Web of Science", "Web of science","WoS"),collapse= "|"),negate=FALSE), "Scopus or WoS",
                                                  Article_source)),
         Year = round(as.numeric(Year),0),
         Volume = round(as.numeric(Volume),0),
         Issue = round(as.numeric(Issue),0),
         Page_start = round(as.numeric(Page_start),0),
         Page_end = round(as.numeric(Page_end),0),
         Page_count = round(as.numeric(Page_count),0)) %>%
  select(ID,Article_source,Inclusion_yes_no,Exclusion_reasion_pico,Exclusion_reason,
         Authors,Year, Title,Source.title,Volume,Issue,Art_No,Page_start,Page_end,Page_count,
         DOI,Link,ISSN,ISBN)

addmargins(table(articles$Article_source))
ggplot(articles,aes(x=Year))+
  geom_bar(colour="white")

# Export articles for paper 
write.csv(articles,paste0("A List of articles reviewed.csv"),row.names=F)

# clean main data ####
names(data)
data <- data[c(4:nrow(data)),] # drop extra header rows

#names(data)
data <- data[!is.na(data$ID),]
data[data=="ND"] <- "nd"
data[data=="check"] <- NA
data[data=="Check"] <- NA
data[data=="CHECK"] <- NA
data[data=="na"] <- NA
data[data=="NA"] <- NA
data[data==""] <- NA
data[data=="."] <- NA
data[data=="N/A"] <- NA

# exclude rows identified as having no comparison during data entry
data %>%
  filter((System %in% "NA - exclude")) %>%
  select(ID,System,System_raw)
data %>%
  filter((Comparison_class %in% "NA - exclude")) %>%
  select(ID,System,System_raw,Comparison_class)

# check where there are missing lat-longs
data %>%
  filter(is.na(Lat)|is.na(Long)|Lat %in% ("nd") | Long %in% c("nd")) %>%
  select(ID,Country)
data <- data %>%
  filter(!(System %in% "NA - exclude")) %>%
  filter(!(Comparison_class %in% "NA - exclude")) %>%
  mutate(Lat = gsub("^(.*?),.*", "\\1", Lat), # when several entries, take the first
         Lat = gsub("^(.*?);.*", "\\1", Lat),
         Long = gsub("^(.*?),.*", "\\1", Long),
         Long = gsub("^(.*?);.*", "\\1", Long)) %>%
  mutate(Lat = ifelse(Lat =="nd",NA,Lat),
         Long = ifelse(Long =="nd",NA,Long))
nrow(data[which(is.na(data$Lat)),]) # 198
nrow(data[which(is.na(data$Long)),]) # 198

sort(unique(data$Long))
sort(unique(data$Lat))

# format columns
sort(unique(data$Experiment_stage))

data<- data %>%
  ##To transform factors to numeric, and lower case##  
  mutate(ID=as.character(round(as.numeric(ID),0)),
         Comparison_ID = as.character(Comparison_ID),
         Experiment_stage = as.character(round(as.numeric(Experiment_stage),0)),
         B_value = as.numeric(as.character(B_value)),
         B_error_measure = as.character(B_error_measure),
         B_error_value = as.numeric(as.character(B_error_value)),
         B_error_range_l = as.numeric(as.character(B_error_range_l)),
         B_error_range_u = as.numeric(as.character(B_error_range_u)),
         B_SD = as.numeric(as.character(B_SD)), 
         B_N = as.numeric(as.character(B_N)),
         Yield_measure = as.character(Yield_measure),
         Yield_value = as.numeric(as.character(Yield_value)),
         Yield_error_measure = as.character(Yield_error_measure),
         Yield_error_value = as.numeric(as.character(Yield_error_value)),
         Yield_error_range_l = as.numeric(as.character(Yield_error_range_l)),
         Yield_error_range_u = as.numeric(as.character(Yield_error_range_u)),
         Yield_SD = as.numeric(as.character(Yield_SD)),
         Yield_N = as.numeric(as.character(Yield_N)),
         Crop1_C_yield = as.numeric(Crop1_C_yield),
         Crop1_C_SD = as.numeric(Crop1_C_SD),
         Crop1_T_yield = as.numeric(Crop1_T_yield),
         Crop1_T_SD = as.numeric(Crop1_T_SD),
         Crop2_C_yield = as.numeric(Crop2_C_yield),
         Crop2_C_SD = as.numeric(Crop2_C_SD),
         Crop2_T_yield = as.numeric(Crop2_T_yield),
         Crop2_T_SD = as.numeric(Crop2_T_SD),
         LER_value = as.numeric(LER_value),
         LER_SD = as.numeric(LER_SD),
         LER_N = as.numeric(LER_N),
         Lat= as.numeric(as.character(Lat)),
         Long = as.numeric(as.character(Long)),
         Country = as.character(Country),
         Functional_group = as.character(Functional_group),
         Taxa_class = as.character(Taxa_class),
         Crop = as.character(Crop),
         System_raw = as.character(System_raw),
         System = as.factor(as.character(System))) %>%
  mutate(Yield_measure = ifelse(Yield_measure =="nd",NA,Yield_measure))

# Make comparison class column 
sort(unique(data$System))
data <- data %>%
  mutate(Comparison_class = ifelse(System %in% c("Natural","Abandoned farmland"),"Natural",
                                   ifelse(System %in% c("Monoculture","Simplified other"),"Simplified","Diversified")))
### Join the article source information and year of publication   
setDT(data)
setDT(articles)
data <- articles[data,on="ID"]

### Remove unnecessary columns
data <- data %>%
  # select all columns we will include in database for data paper
  select(ID, Comparison_ID, Comparison_class,Experiment_stage,
         Crop,	Crop_FAO,Crop_species_richness, Crop_variety_richness, Crop_ann_pen,	Crop_woodiness,
         crops_all_common,	crops_all_scientific,	crops_all_scientific_level, # exclude these?
         System_raw, System_details, System, 
         Taxa, Taxa_details,Taxa_phylum,Taxa_class,Taxa_order, Functional_group, B_ground, 
         #taxa_common,	taxa_scientific,	taxa_scientific_level, # exclude these?
         Farm_size,	Farm_context,	Fertiliser,	Fertiliser_chem,	Pesticide,	Pesticide_quantity,	Soil_management,	
         Time_state,	Study_length,	Landscape_context,	Sampling_unit,
         B_measure,B_value, B_error_measure, B_error_value,B_error_range_l,B_error_range_u, B_SD, B_N, 
         Yield_measure,Yield_value,Yield_error_measure,Yield_error_value, Yield_error_range_l,Yield_error_range_u,Yield_SD,Yield_N,
         Crop1_C_yield,Crop1_C_SD,Crop1_T_yield,Crop1_T_SD,Crop2_C_yield,Crop2_C_SD,Crop2_T_yield,Crop2_T_SD,LER_value,LER_SD,LER_N,
         Location, Country,  Lat, Long, Data_entry, Data_validation, Notes,Notes_yield,
         Authors,Title,Year) %>%
  mutate(Notes = ifelse(!is.na(Notes_yield),paste0(Notes, " ",Notes_yield),Notes)) %>%
  select(-c(Notes_yield))

# see where there is missing data
colSums(sapply(data, is.na)) 

# Remove rows with missing biodiversity means, SDs, or sample size
setDT(data)
data <- data[complete.cases(data[,c("B_value","B_error_value","B_N")])] # none removed

## Get % missing data
nrow(data[Soil_management == "nd" | is.na(Soil_management),])/nrow(data)*100 #74.5%
nrow(data[Time_state == "nd" | is.na(Time_state),])/nrow(data)*100 # 17.5%
nrow(data[Farm_size == "nd" | is.na(Farm_size),])/nrow(data)*100 # 61.6%
nrow(data[Data_validation  == "nd" | is.na(Data_validation ),])/nrow(data)*100 #0%

## summarise values
#describe(data)
#df_status(data)
#freq(data)
# Important!: Should not be any "#N/A" =so check in database why automated classification didn't work

### Check through frequency charts and clean data again to fix spelling errors or inconsistencies
sort(unique(data$Country))
data$Country <- as.character(data$Country)
data$Country[data$Country %in% c("Republic of Benin")] <- "Benin" 
data$Country[data$Country %in% c("Brazil ")] <- "Brazil" 
data$Country[data$Country %in% c("Canada ")] <- "Canada" 
data$Country[data$Country %in% c("Costa Rica ")] <- "Costa Rica" 
data$Country[data$Country %in% c("China ")] <- "China" 
data$Country[data$Country %in% c("Germany ")] <- "Germany" 
data$Country[data$Country %in% c("Indonesia ")] <- "Indonesia" 
data$Country[data$Country %in% c("Kenia")] <- "Kenya" 
data$Country[data$Country %in% c("Malaysia ")] <- "Malaysia"
data$Country[data$Country %in% c("Mexico ")] <- "Mexico"
data$Country[data$Country %in% c("Nigeria ")] <- "Nigeria" 
data$Country[data$Country %in% c("Sao Tome ","Sao Tome")] <- "Sao Tome and Principe" 
data$Country[data$Country %in% c("South Africa ")] <- "South Africa"
data$Country[data$Country %in% c("Swiss")] <- "Switzerland"
data$Country[data$Country %in% c("Sweeden")] <- "Sweden" 
data$Country[data$Country %in% c("Vietnam")] <- "Viet Nam" 
data$Country[data$Country %in% c("England","Scotland","Wales","UK")] <- "United Kingdom" 
data$Country[data$Country %in% c("California, USA","Washington, USA","Washington, USA","USA")] <- "United States of America" 
sort(unique(data$Country))

# add country centroid as lat-long for rows with missing lat-long
# get a data.frame with centroids
centroids <- getMap(resolution="high")
centroids <- gCentroid(centroids, byid=TRUE)
centroids <- as.data.frame(centroids)
centroids$Country <- row.names(centroids)
data <- data %>% 
  left_join(centroids,by="Country")# %>%

unique(data[which(is.na(data$Lat)),c("ID", "Country","Lat","x","y")])

data <- data %>%  
  mutate(Lat_original = Lat,
         Long_original = Long,
         Lat = ifelse(is.na(Lat),y,Lat),
         Long = ifelse(is.na(Long),x,Long))

sort(unique(data$Crop))
data$Crop <- as.character(data$Crop)
data$Crop[data$Crop %in% c("Green bean","Beans, green","French bean","Haricot bean")] <- "Green bean"
#data$Crop[data$Crop %in% c("Lima bean")] <- "Beans, dry"
data$Crop[data$Crop %in% c("Black walnut", "Walnut ")] <- "Walnut"
data$Crop[data$Crop %in% c("Broccolli")] <- "Broccoli"
data$Crop[data$Crop %in% c("Brussels sprout")] <- "Brussel sprout"
data$Crop[data$Crop %in% c("Cabbage, mustard")] <- "Cabbage-mustard" 
data$Crop[data$Crop %in% c("Cereal","Cereals","Wheat, maize"  , "Maize and Wheat"  )] <- "Cereals" 
data$Crop[data$Crop %in% c("Cocoa ","Cacao")] <- "Cocoa" 
data$Crop[data$Crop %in% c("Coffee/Cocoa" )] <- "Coffee" 
data$Crop[data$Crop %in% c("forest", "Forest ","Forest  ","Forest verges","Primary forest","Woodland", "Trees, closed canopy","Trees, open canopy")] <- "Forest"
data$Crop[data$Crop %in% c("Field boundary")] <- "Field margin"
data$Crop[data$Crop %in% c("Fruit" )] <- "Fruit trees" 
data$Crop[data$Crop %in% c( "Grassland: combined mowing and grazing, Winter wheat, Oats+peas","Prairie herbaceous as biofuel")] <- "Grassland"
data$Crop[data$Crop %in% c( "Grassland","Grass","Smooth bromegrass")] <- "Grassland"
data$Crop <- ifelse(data$Comparison_class %in% c("Simplified","Diversified") & data$Crop %in% c("Grassland", "Forage"),"Pasture",data$Crop)
data$Crop[data$Crop %in% c("Grapevines","Grapevine","Vineyard","vineyard")] <- "Grape" 
data$Crop[data$Crop %in% c("Curly-leaf kale")] <- "Kale" 
data$Crop[data$Crop %in% c("Maize ","Maize  ", "Maize, Hedgerow" )] <- "Maize"
data$Crop[data$Crop %in% c( "Maize, Mustard" , "Maize, Potato", "Maize, Sesame"  , "Maize, Sweet potato","Corn","Maize, cassava" )] <- "Maize-other"
data$Crop[data$Crop %in% c( "Maize, Cowpea" ,   "Maize, Desmodium" ,  "Maize, Desmodium, Napier grass" , 
                            "Maize, Faba bean" , "Maize, haricot bean", "Maize, Haricot bean","Maize, lablab",  "Maize, Leucaena")] <- "Maize-legume"
data$Crop[data$Crop %in% c("Oenothera biennis")] <- "Common evening primrose"
data$Crop[data$Crop %in% c("Palm Oil ","Palm Oil")] <- "Palm oil"
data$Crop[data$Crop %in% c("Root crops","Mixed tropical")] <- "Mixed"
data$Crop[data$Crop %in% c("Natural land")] <- "Unspecified"
data$Crop[data$Crop %in% c("Spring Pea")] <- "Pea"
data$Crop[data$Crop %in% c("Sorghum, Sesame","Sorghum, Sweet potato")] <- "Sorghum-other"
data$Crop[data$Crop %in% c("Sorghum, Cowpea","Sorghum, Desmodium","Sorghum, Haricot bean")] <- "Sorghum-legume"
data$Crop[data$Crop %in% c("Sunflowers")] <- "Sunflower"
data$Crop[data$Crop %in% c("Rubber plantation" )] <- "Rubber"
data$Crop[data$Crop %in% c("Soy Beans" , "Soybean, grasses")] <- "Soybean"
data$Crop[data$Crop %in% c("Wheat\n","Winter Wheat ","Winter Wheat", "Winter wheat")] <- "Wheat" 
data$Crop[data$Crop %in% c("Wheat, garlic","Wheat, oilseed rape","Canola, wheat")] <- "Wheat-other" 
data$Crop[data$Crop %in% c("Wheat, Alfalfa","Wheat, pea")] <- "Wheat-legume" 
data$Crop[data$Crop %in% c("No specified" ,"Perennial energy crop" ,"Perrenials and annuals ","Annual crops","Agroforestry")] <- "nd" 
table(data[which(is.na(data$Crop)),]$Comparison_class)
check <- data[which(is.na(data$Crop)),]
data$Crop[is.na(data$Crop)] <- "nd" 
sort(unique(data$Crop))

data$Crop <- gsub(".*, .*", "Mixed", data$Crop) # change everything that has more than one crop separated by ',' to 'Mixed'
data$Crop <- gsub(".* or .*", "Mixed", data$Crop) # change everything that has more than one crop separated by 'or' to 'Mixed'
data$Crop <- gsub(".* and .*", "Mixed", data$Crop) # change everything that has more than one crop separated by 'and' to 'Mixed'
data$Crop[data$Crop %in% c("nd" ,"Mixed","Mixed tropical","Prairie herbaceous as biofuel")] <- "Mixed or nd" 
sort(unique(data$Crop))

ggplot(data,aes(y=Crop))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3)+
  theme_bw()
sort(names(data))

sort(unique(data$Crop_FAO))
data$Crop_FAO[data$Crop_FAO %in% c("NATURAL")] <- NA
data$Crop_FAO[data$Crop_FAO %in% c("MIXED", "OTHER","nd","Other or  nd")] <- "Other or nd"
data$Crop_FAO[data$Crop %in% c("Sunflower","Common evening primrose","Palm oil")] <- "6 - OIL-BEARING CROPS AND DERIVED PRODUCTS" 
data$Crop_FAO[data$Crop %in% c("Maize")] <- "1 - CEREALS AND CEREAL PRODUCTS"   
data$Crop_FAO[data$Crop %in% c("Mixed or nd","Other or nd","Other or  nd")] <- "Other or nd"   
data$Crop_FAO[data$Crop %in% c("Grape")] <- "8 - FRUITS AND DERIVED PRODUCTS"   
data$Crop_FAO[data$Crop %in% c("Leucaena")] <- "11 - FODDER CROPS AND PRODUCTS"   

ggplot(data,aes(y=Crop_FAO))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 

sort(unique(data$System)) 
sort(unique(data[data$System=="Natural",]$System_raw))
data$System_raw[data$System_raw %in% c("Forest", "Forest ", "Forest reserve","Natural forest","Primary Forest","Primary Forest ","Semi-natural forest","Secondary Forest", "Remnant","Secondary Succesional Habitat ")] <- "Forest" 
data$System_raw[data$System_raw %in% c("Fynbos","Shrubland"  )] <- "Shrubland" 
data$System_raw[data$System_raw %in% c("Grassland - natural","Grassland - semi-natural"  )] <- "Grassland" 
data$System_raw[data$System_raw %in% c("Wetlands"  )] <- "Wetland" 
ggplot(data[data$System=="Natural",],aes(y=System_raw))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 
ggplot(data[data$System=="Abandoned farmland",],aes(y=System_raw))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 
check <- unique(data[data$System =="Abandoned farmland",c("System","System_details","System_raw")])
unique(check$System_details)
data$System_raw[data$System == "Abandoned farmland" & !(data$System_raw %in% c("Forest"  ))] <-  "Mixed"

ggplot(data,aes(y=System))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 
data$System <- as.character(data$System)
#data$System[data$System %in% c("Embedded natural" )] <- "Diversified other" 
data$System[data$System %in% c("Abandoned farmland")] <- "Abandoned"

sort(unique(data$Pesticide))
data$Pesticide <- as.character(data$Pesticide)
data$Pesticide[data$Pesticide %in% c( "Fungicide" ,
                                     "Fungicide, insecticide",
                                     "Sometimes",
                                     "YES","Yes ","yes")] <- "Yes" 
data$Pesticide[data$Pesticide %in% c("NO","No ","no")] <- "No"
sort(unique(data$Pesticide))
data$Pesticide <- ifelse(data$Pesticide %in% c("Yes","No"),data$Pesticide,"nd")

sort(unique(data$Fertiliser))
data$Fertiliser <- as.character(data$Fertiliser)
data$Fertiliser[data$Fertiliser %in% c("Probably yes", "Part","Sometimes", "Mixed",
                                       "YES","Yes ","varies","YES\n","yes")] <- "Yes" 
data$Fertiliser[data$Fertiliser %in% c("NO","No ","no")] <- "No"
sort(unique(data$Fertiliser))
data$Fertiliser <- ifelse(data$Fertiliser %in% c("Yes","No"),data$Fertiliser,"nd")

sort(unique(data$Soil_management))
data$Soil_management <- as.character(data$Soil_management)
data$Soil_management[data$Soil_management %in% c("No Tillage","No Tillage ","No tillage","no tillage","No plough","no","no\n")] <- "No Tillage"
data$Soil_management[data$Soil_management %in% c("Yes","yes", "Tillage ", "tillage", "Reduced Tillage", "Reduced Tillage (10cm)","Chisel Plough","Tillage (Varrying levels) ", "Plough","Plough tillage","Plowing, Harrowing, Liming ")] <- "Tillage"
data$Soil_management[data$Soil_management %in%  c("Annual vegetation clearance and burning","slash and burn","Burning")] <- "Burning"
data$Soil_management[data$Soil_management %in%  c("nd\n","nd")] <- "nd"
sort(unique(data$Soil_management))
ggplot(data,aes(y=Soil_management))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 

# only 8 rows = burning, so remove
#data$Soil_management <- ifelse(data$Soil_management %in% c("No Tillage", "Tillage"),data$Soil_management,"nd")

ggplot(data,aes(y=Taxa_order,fill=Taxa_phylum))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+30,label=..count..),size=3) 
data$Taxa_order[data$Taxa_order %in% c("Amphibia","Amphibian")] <- "Amphibians"
data$Taxa_order[data$Taxa_order %in% c("Reptile")] <- "Reptiles"

ggplot(data,aes(y=Taxa_class,fill=Taxa_phylum))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+10,label=..count..),size=3) 
data$Taxa_group <- data$Taxa_order
sort(unique(data$Taxa_group))
data$Taxa_group[data$Taxa_group == "Aves"] <- "Birds"
data$Taxa_group[data$Taxa_group %in% c("Arachnida","Araneae")] <- "Arachnids"
data$Taxa_group[data$Taxa_group == "Chiroptera"] <- "Bats"
data$Taxa_group[data$Taxa_group %in% c("Nematoda","Annelida")] <- "Nematodes and annelids"
data$Taxa_group[data$Taxa_group %in% c("Lagomorpha")] <- "Mammals"
data$Taxa_group <- ifelse(data$Taxa_group %in% c("Mammals"),data$Taxa_class,data$Taxa_group)
data$Taxa_group[data$Taxa_group %in% c("Mammals")] <- "Mammals other"
#sort(unique(data[data$Taxa_group == "Mammals other",]$Taxa_details))
data$Taxa_group[data$Taxa_group %in% c("Arthropods",
                                       "Insect n.s","Insects ns/oth",
                                       "Collembola",
                                       "Diptera",
                                       "Dermaptera",
                                       "Homoptera",
                                       "Isoptera",
                                       "Isopoda",
                                       "Orthoptera",
                                       "Neuroptera",
                                       "Neoptera",
                                       "Thysanoptera")] <-"Arthropods other"
data$Taxa_group[data$Taxa_group %in% c("Plant vascular non-woody" ,
                                       "Plant nonvascular" , 
                                       "Plants n.s.","Plants n.s",
                                       "Plant vascular n.s.")] <-"Plants other"
data$Taxa_group[data$Taxa_group %in% c("Plant vascular woody")] <-"Plants woody"
data$Taxa_group[data$Taxa_group %in% c("Lepidoptera")] <-"Butterflies, moths"
data$Taxa_group[data$Taxa_group %in% c("Hymenoptera")] <-"Sawflies, wasps, bees, ants"
data$Taxa_group[data$Taxa_group %in% c("Hemiptera")] <-"True bugs"
#sort(unique(data[data$Taxa_group == "Coleoptera",]$Taxa_details))
data$Taxa_group[data$Taxa_group %in% c("Coleoptera")] <-"Beetles"
#sort(unique(data[data$Taxa_group == "Myriapoda",]$Taxa_details))
data$Taxa_group[data$Taxa_group %in% c("Myriapoda")] <-"Millipedes, centipedes, symphylans"
sort(unique(data$Taxa_group))
ggplot(data,aes(y=Taxa_group,fill=Taxa_phylum))+geom_bar()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_text(stat='count', aes(x=..count..+10,label=..count..),size=3) 

# Arthropods:
# mostly aphids, true bugs, beetles, butterflies or moths, ladybirds, thrips
# Aphids, subfamily Aphidoidea, order Hemiptera
# Butterflies, suborder Rhopalocera, order Lepidoptera
# Coleoptera = beetles, order Coleoptera
# Chrysopa = lacewigs, order Neuroptera
# Dermaptera = earwigs, order Dermaptera
# Diptera = flies, order Diptera
# Hemiptera = true bug, order Hemiptera
# Homoptera = aphids, suborder of order Hemiptera
# Hoverflies, family Syrphidae, order Diptera
# Hymenoptera = sawflies, wasps, bees, ants, order Hymenoptera
# Isoptera = termites, infraorder of order Blattodea
# Lepidoptera = butterflies and moths, order Lepidoptera
# Moths, unranked Heterocera, order Lepidoptera
# Neoptera, includes several orders of winged insects
# Neuroptera = lacewig, order Neuroptera
# Psocoptera = barkflies, order Psocoptera
# Thysanoptera = thrips (<1mm winged insects), order Thysanoptera

sort(unique(data$Functional_group))
data$Functional_group <- as.character(data$Functional_group)
data$Functional_group[data$Functional_group %in% c("decomposer", "decomposers","Decomposers","Descomposers" )] <-"Decomposers"
data$Functional_group[data$Functional_group %in% c("Frugivore"  , "frugivores"  , "Frugivores","Frugivorous"  )] <-"Frugivores"
data$Functional_group[data$Functional_group %in% c("granivore","granivores"  )] <-"Granivores"
data$Functional_group[data$Functional_group %in% c("Herbivore" , "herbivores" ,"Herbivores")] <-"Herbivores"
data$Functional_group[data$Functional_group %in% c("insectivores")] <-"Insectivores"
data$Functional_group[data$Functional_group %in% c("omnivores","omnivore","Omnivore")] <-"Omnivores"
data$Functional_group[data$Functional_group %in% c( "Other" , "Others" ,"others","Fugivores/Insectivores","mixed","Mixed","nd")] <-"Other"
data$Functional_group[data$Functional_group %in% c( "pest" , "Pest")] <-"Pests"
data$Functional_group[data$Functional_group %in% c("Pest control","Pest parasitoid", "Predator", "predators" ,"Predators","Pest enemy","natural enemies")] <-"Natural enemies"
data$Functional_group[data$Functional_group %in% c("Pollinator","pollinator")] <-"Pollinators"
data$Functional_group[data$Functional_group %in% c("weed","Weed","weeds")] <-"Weeds"
data$Functional_group[data$Functional_group %in% c("Plant","Plants")] <-"Autotrophs" # since plants is not a functional group
data$Functional_group[data$Functional_group %in% c("Overall carabids")] <-"Carnivores"
funModeling::freq(data$Functional_group) 

data$Pest_group <- as.character(data$Functional_group)
data$Pest_group[data$Pest_group %in% c("Pests","Weeds")] <-"Pests"
data$Pest_group <- ifelse(data$Pest_group %in% c("Pests"),data$Pest_group, "Other")
funModeling::freq(data$Pest_group)

funModeling::freq(data$Data_entry)
funModeling::freq(data$Data_validation)

# calculate SD again here for full tracability
# Equations to calculate the SD from SE, SEM, inter quartil-range and CI##

### Equations to calculate the SD from SE, IC and IQR ####
# Equation to calculate the SD from SE (Higgins & Green 2011)(a= B_error_value; b= N_samples)
# http://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
# a = standard error, b = number of samples
SE_SD <- function (a, b) {  
  result<- a * sqrt(b)
  return(result)
}

##Equation to calculate the SD from M_IQR (Hozo et al., 2005) 
# a= Number of samples; b= lower IQR; c= median; d= upper IQR
M_IQR_SD<- function (a, b,c,d) {  
  result<- sqrt(((a + 1)/(48 * a*((a-1)^2))) * (((a^2) + 3) * ((b - (2*c) + d)^2) + (4* (a^2)) * ((d - b)^2)))
  return(result)
}

##Equation to calculate the SD from CI (Higgins & Green 2011) (a= N_samples; b= B_error_value)
##http://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
# a = number of samples, b = CI range
CI_SD<- function (a, b) {  
  # qt(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
  result<- sqrt(a) * (b/((qt((1-(0.05/2)), (a - 1)))*2)) 
  return(result)
}


sort(unique(data$B_error_measure))
sort(unique(data$Yield_error_measure))

data <- data %>% 
  mutate(B_SD = if_else(B_error_measure %in% c("Standard error" ,"standard error", "Square error of the Mean"), SE_SD(B_error_value, B_N),
                                      if_else(B_error_measure == "median,IQR", M_IQR_SD( B_N, B_error_range_l, B_error_value, B_error_range_u),
                                              if_else(B_error_measure == "Confidence intervals", CI_SD(B_N,B_error_value),
                                                      (B_error_value))))) %>%
  mutate(Yield_SD = if_else(Yield_error_measure %in% c("Standard error", "standard error" ,"Square error of the Mean"), SE_SD(Yield_error_value, Yield_N),
                            if_else(Yield_error_measure == "median,IQR", M_IQR_SD(Yield_N, Yield_error_range_l, Yield_error_value, Yield_error_range_u),
                                    if_else(Yield_error_measure %in% c("Confidence intervals","Confidence interval"), CI_SD(Yield_N,Yield_error_value),
                                            (Yield_error_value))))) %>%
  # Prioritise LER for intercropping and other metrics (kg/ha or g/plant etc) for everything else
  mutate(Yield_measure = ifelse(System == "Intercropping" & !(is.na(LER_value)) & LER_value !="nd","LER",ifelse(System=="Intercropping" & (is.na(Yield_value) | Yield_measure =="nd"),NA,Yield_measure)))  %>%
  mutate(Yield_value = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",LER_value,Yield_value),
         Yield_SD = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",LER_SD,Yield_SD),
         Yield_N = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",LER_N,Yield_N),
         Yield_error_measure = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",NA,Yield_error_measure), 
         Yield_error_value = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",NA,Yield_error_value), 
         Yield_error_range_l = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",NA,Yield_error_range_l), 
         Yield_error_range_u = ifelse(System == "Intercropping" & !(is.na(LER_value))& LER_value !="nd",NA,Yield_error_range_u)) %>%
  select(-c(x,y, 
            #Yield_error_range_l,Yield_error_range_u,B_error_range_l,B_error_range_u,
            #B_error_measure,Yield_error_measure,B_error_value,Yield_error_value,
            Data_entry,Data_validation))

#data <- data %>%
#  select(-c(Yield_value_2,Yield_SD_2,Yield_N_2, Yield_error_range_l_2,Yield_error_range_u_2,Yield_error_measure_2,Yield_error_value_2,Yield_measure_2))

check <- data[data$System == "Intercropping",]
addmargins(table(check$Yield_measure_1))
table(check$System,check$Yield_measure)

names(data)
#data <- data[,c(1:18,50,19:43,48,49,44:47,51)]
#data <- data[,c(1:51,56,57,52:55,58,59)]
# export cleaned database 

write.csv(select(data,-c("Pest_group","Taxa_group")), paste0("B_biodiversity_local_data_cleaned_",format(Sys.time(),"%Y%m%d"),".csv"),row.names=FALSE)

### Match controls and treatments within each study ####

# select columns we will include in merge or effect size database
data_meta <- data %>%
  select(ID, Experiment_stage, Comparison_ID, Comparison_class, 
         Crop,	Crop_FAO, Crop_ann_pen,	Crop_woodiness,crops_all_common, crops_all_scientific,crops_all_scientific_level,
         System_raw, System_details, System,
         Taxa, Taxa_details, Taxa_group, Taxa_class, Taxa_order, Taxa_phylum, Functional_group,B_ground,Pest_group,
         Fertiliser,	Fertiliser_chem,	Pesticide,Pesticide_quantity,		Soil_management,
         Farm_size, Time_state, Study_length, Farm_context, Landscape_context, Sampling_unit,
         B_error_measure, B_error_value, B_error_range_l, B_error_range_u,
         B_measure,B_value, B_SD,B_N, 
         Yield_measure,Yield_value, Yield_SD,Yield_N,
         Yield_error_measure, Yield_error_value, Yield_error_range_l, Yield_error_range_u,
         Location, Country,  Lat, Long, Lat_original,Long_original,
         #Data_entry, Data_validation, 
         Authors, Title, Year,Notes)

data_meta <- data_meta[which(!is.na(data_meta$Country)),] # none removed
data_meta <- unique(data_meta) # 4765

# Match comparators and interventions to calculate effect sizes  ####

# Merge_ID specifies what parameters have to match between the control and treatment for them to be compared
data_meta <- data_meta %>%
  mutate(Merge_ID = paste(ID,
                       Experiment_stage,
                       Taxa_group,
                       Taxa_details,
                       Functional_group, 
                       B_measure,
                       #Yield_measure, # don't use this because the LER entries mess things up
                       Country, # doesn't change anything so could exclude
                       sep="_")) 
length(unique(data_meta$Merge_ID)) # 1445

setDT(data_meta)
data_C <- data_meta[which(data_meta$Comparison_class %in% c("Simplified","Natural")),]
data_T <- data_meta[which(data_meta$Comparison_class=="Diversified"),]

# 
names(data_meta)
data_meta <- merge.data.frame(data_C[,-c( "Taxa" ,"Taxa_details", "Taxa_group","Taxa_order","Taxa_class","Taxa_phylum", "Functional_group","B_ground","Pest_group" ,
                                         "B_measure",
                                         "Farm_size","Farm_context", "Landscape_context","Location", "Country", 
                                         "Authors","Title","Year","Notes")],
                              # remove columns that are definitely the same in C and T data 
                              # otherwise they will be duplicated 
                              data_T[,-c("ID","Experiment_stage")],  
                         by.x="Merge_ID",by.y="Merge_ID",all.x=T,all.y=T,suffixes=c("_C","_T"))
names(data_meta)

data_meta <- unique(data_meta) # none removed
length(unique(data_meta$ID)) #237 articles, 4076 effects after fixing IDs 651 and 754
nrow(data_meta)

#write.xlsx(list(data_meta), sheetName = c("data_meta"), file="data_check.xlsx")

# Check that all rows have a comparison (if not, check why in raw data) 
nrow(data_meta[(is.na(data_meta$ID)),])
sort(unique(data_meta[which(is.na(data_meta$B_value_T)),]$ID))
sort(unique(data_meta[which(is.na(data_meta$B_value_C)),]$ID))
sort(unique(data_meta[which(is.na(data_meta$Yield_value_T)),]$ID)) # no yields for many where there were only biodiversity outcomes
sort(unique(data_meta[which(is.na(data_meta$Yield_value_C)),]$ID)) # no yields when only natural habitat as control (142, 568)

# Set yield_measure_C to LER if yield_value_T is LER (but yield_value_C and other value columns to NA)
# Set yield_measure_C to NA if yield_measure_T is NA
data_meta <- data_meta %>%
  # set to NA yield measure in controls that use LER to measure yields (which is already calculated and recorded in the treatment)
  mutate(Yield_measure_C = ifelse(Yield_measure_T == "LER" |  is.na(Yield_measure_T) | is.na(Yield_value_C) | is.na(Yield_value_T),NA,Yield_measure_C)) %>%
  # set to NA yield values in treatments with no matching control data
  mutate(Yield_value_T = ifelse(is.na(Yield_value_C) & Yield_measure_T != "LER",NA,Yield_value_T)) %>%
  # set to NA yield values in controls with no matching treatment data, or that use LER
  mutate(Yield_value_C = ifelse(is.na(Yield_value_T) | Yield_measure_T == "LER" | is.na(Yield_measure_C),NA,Yield_value_C)) %>%
  # set to NA yield variance and sample size if there is no yield value data
  mutate(Yield_N_C = ifelse(is.na(Yield_value_C),NA,Yield_N_C),
         Yield_SD_C = ifelse(is.na(Yield_value_C),NA,Yield_SD_C),
         Yield_error_measure_C = ifelse(is.na(Yield_value_C)|Yield_error_measure_C=="nd",NA,Yield_error_measure_C),
         Yield_error_range_l_C = ifelse(is.na(Yield_value_C),NA,Yield_error_range_l_C),
         Yield_error_range_u_C = ifelse(is.na(Yield_value_C),NA,Yield_error_range_u_C),
         Yield_error_value_C = ifelse(is.na(Yield_value_C),NA,Yield_error_value_C),
         Yield_N_T = ifelse(is.na(Yield_value_T),NA,Yield_N_T),
         Yield_SD_T = ifelse(is.na(Yield_value_T),NA,Yield_SD_T),
         Yield_error_measure_T = ifelse(is.na(Yield_value_T)|Yield_error_measure_T=="nd",NA,Yield_error_measure_T),
         Yield_error_range_l_T = ifelse(is.na(Yield_value_T),NA,Yield_error_range_l_T),
         Yield_error_range_u_T = ifelse(is.na(Yield_value_T),NA,Yield_error_range_u_T),
         Yield_error_value_T = ifelse(is.na(Yield_value_T),NA,Yield_error_value_T))

# check for rows with yield means but no SD or sample size
#check <- data_meta[which(!is.na(data_meta$Yield_value_T) & (data_meta$Yield_N_T =="nd" | is.na(data_meta$Yield_N_T))),]
#check <- data_meta[which(!is.na(data_meta$Yield_value_T) & (data_meta$Yield_SD_T =="nd" | is.na(data_meta$Yield_SD_T))),]

### Compute validity values (0 low quality, 1 high quality, nd or NA based on Sanchez et al. 2021) ####
data_meta <- data_meta %>%
  mutate(Yield_data = ifelse(!is.na(Yield_value_T)&!is.na(Yield_SD_T) & !is.na(Yield_N_T) & 
                               Yield_value_T !="nd"& Yield_SD_T !="nd" & Yield_N_T !="nd",1,0),
         Yield_problems = ifelse((Yield_measure_C == Yield_measure_T) | Yield_measure_T =="LER",0,1))

nrow(data_meta[which(data_meta$Yield_problems ==1),]) #none

data_meta <- data_meta %>%
  mutate(Yield_measure = Yield_measure_T) %>%
  select(-c(Yield_measure_C,Yield_measure_T,Yield_problems))

table(data_meta$Yield_data) #1222 effect sizes

data_meta <- data_meta %>%
  mutate(Year = round(as.numeric(Year),0),
         Validity_N_biodiversity_C = ifelse(B_N_C <5,0,1),
         Validity_N_biodiversity_T = ifelse(B_N_T <5,0,1),
         Validity_N_yield_C = ifelse(Yield_data==0,NA,
                                     ifelse(Yield_measure=="LER"&Yield_N_T<5,0,
                                            ifelse(Yield_measure=="LER"&Yield_N_T>=5,1,
                                                   ifelse(Yield_N_C <5,0,1)))),
         Validity_N_yield_T = ifelse(is.na(Yield_N_T|Yield_data==0),NA,
                                     ifelse(Yield_N_T <5,0,1)),
         Validity_time_C = ifelse(Time_state_C == "nd",Time_state_C,
                                ifelse(str_detect(Time_state_C, paste(">"),negate=FALSE),1,
                                       ifelse(str_detect(Time_state_C, paste("<"),negate=FALSE),0,1))),
         Validity_time_T = ifelse(Time_state_T == "nd",Time_state_T,
                                  ifelse(str_detect(Time_state_T, paste(">"),negate=FALSE),1,
                                         ifelse(str_detect(Time_state_T, paste("<"),negate=FALSE),0,1))),
         # for location use lat-long to 2.d.p. which means sites are within 1.1 km of each other
         Validity_location = ifelse(is.na(data_meta$Lat_original_C) | is.na(data_meta$Lat_original_T) | is.na(data_meta$Long_original_C) | is.na(data_meta$Long_original_T),"nd",
                                      ifelse(paste(round(data_meta$Lat_C,2),round(data_meta$Long_C,2),sep=",")==paste(round(data_meta$Lat_T,2),round(data_meta$Long_T,2),sep=","),1,0)))

df_status(data_meta)

data_meta <- data_meta %>%
  mutate(Validity_biodiversity_overall = ifelse(Validity_N_biodiversity_C == 1 & Validity_N_biodiversity_T == 1 &
                                                  Validity_time_C == "1" & Validity_time_T == "1" & Validity_location =="1","1",
                                                ifelse(Validity_time_C == "nd" | Validity_time_T =="nd" | Validity_location =="nd","nd","0")))
addmargins(table(data_meta$Validity_biodiversity_overall))

data_meta <- data_meta %>%
  mutate(Validity_yield_overall = ifelse(Yield_data==0, NA,
                                         ifelse(Validity_N_yield_C == 1 & Validity_N_yield_T == 1 &
                                                    Validity_time_C == "1" & Validity_time_T == "1" & Validity_location =="1","1",
                                                  ifelse((Validity_time_C == "nd") | (Validity_time_T =="nd") | (Validity_location =="nd"),"nd","0"))))
addmargins(table(data_meta$Validity_yield_overall))

addmargins(table(data_meta$System_T, data_meta$Validity_yield_overall))

names(data_meta)
data_meta <- data_meta[c(1:53,63:88,54:62,91:94, 89:90, 95:109)]

# Export
write.csv(select(data_meta,-c("Pest_group","Taxa_group","Yield_data","Merge_ID")),paste0("C_biodiversity_local_data_meta_",format(Sys.time(),"%Y%m%d"), ".csv"),row.names=FALSE)

## Make figure showing results of validity assessment ####
validity <- data_meta %>%
  # set to NA the validity check for yield N in treatments against natural controls
  #mutate(Validity_N_yield_T = ifelse(is.na(Validity_N_yield_C),NA,Validity_N_yield_T),
  #       Validity_N_yield_C = ifelse(is.na(Validity_N_yield_T),NA,Validity_N_yield_C),
  #       Validity_yield_overall = ifelse(is.na(Validity_N_yield_C),NA,Validity_yield_overall)) %>%
  reshape2::melt(id.vars = c("System_T"),
                                      measure.vars=c("Validity_N_biodiversity_C",
                                                     "Validity_N_yield_C",
                                                     "Validity_time_C",          
                                                     "Validity_N_biodiversity_T",    
                                                     "Validity_N_yield_T", 
                                                     "Validity_time_T" ,             
                                                     "Validity_location",
                                                     "Validity_biodiversity_overall",
                                                     "Validity_yield_overall"),
                            variable.name="Criteria",
                            value.name="Bias")
validity$Count <- 1
addmargins(table(validity$Criteria,validity$Bias))
addmargins(table(validity$System_T))

validity <- validity %>%
  group_by(System_T,Criteria,Bias) %>%
  summarise(Count=sum(Count))

validity <- reshape2::dcast(validity,System_T+Criteria~Bias,value.var="Count",sum)
validity <- validity %>%
  dplyr::rename(High_quality = `1`,
                Low_quality = `0`,
                Unknown_quality = nd,
                Not_applicable = `NA`) %>%
  mutate(ES_count = High_quality+Low_quality+Unknown_quality) %>%
  mutate(High_quality_pc = round(High_quality/ES_count*100,1),
         Low_quality_pc = round(Low_quality/ES_count*100,1),
         Unknown_quality_pc = round(Unknown_quality/ES_count*100,1))

validity_total <- validity %>%
  group_by(Criteria) %>%
  summarise(ES_count=sum(ES_count),
            High_quality = sum(High_quality),
            Low_quality = sum(Low_quality),
            Unknown_quality = sum(Unknown_quality),
            Not_applicable = sum(Not_applicable)) %>%
  mutate(System_T = "All diversity practices",
         High_quality_pc = round(High_quality/ES_count*100,1),
         Low_quality_pc = round(Low_quality/ES_count*100,1),
         Unknown_quality_pc = round(Unknown_quality/ES_count*100,1))
names(validity)
names(validity_total)

validity_total <- validity_total[,c("System_T",  "Criteria", "Low_quality", "High_quality",      
                                    "Unknown_quality" ,   "Not_applicable" ,    "ES_count"  , 
                                    "High_quality_pc"   ,"Low_quality_pc"  ,   "Unknown_quality_pc")] 

names(validity)
names(validity_total)
validity <- rbindlist(list(validity,validity_total),use.names=FALSE)

sort(unique(validity$Criteria))
validity$Criteria <- factor(validity$Criteria,
                               levels=c( "Validity_N_biodiversity_C" , "Validity_N_biodiversity_T",
                                         "Validity_N_yield_C" , "Validity_N_yield_T" ,
                                          "Validity_time_C" , "Validity_time_T" ,
                                         "Validity_location",
                                         "Validity_yield_overall",
                                         "Validity_biodiversity_overall" ),
                                                                  
                               labels=c( "Validity_N_biodiversity_C" = "\u2265 5 biodiversity samples for comparator",
                                         "Validity_N_biodiversity_T" = "\u2265 5 biodiversity samples for intervention",
                                          "Validity_N_yield_C" = "\u2265 5 yield samples for comparator" ,
                                        "Validity_N_yield_T" = "\u2265 5 yield samples for intervention",
                                         "Validity_time_C" = "\u2265 1yr in current state for comparator",
                                         "Validity_time_T" = "\u2265 1yr in current state for intervention",
                                         "Validity_location" = "Similar site characteristics for comparator and intervention",
                                        "Validity_yield_overall" = "All criteria met for yield effects",
                                        "Validity_biodiversity_overall" = "All criteria met for biodiversity effects"))

sort(unique(validity$System_T))
validity$System_T <- factor(validity$System_T, levels=c("Agroforestry","Cover crops","Crop rotation","Cultivar mixture",
                                                        "Diversified other","Embedded natural", "Intercropping","All diversity practices")) 
validity$Check <- validity$High_quality_pc + validity$Low_quality_pc + validity$Unknown_quality_pc

addmargins(table(validity$System_T,validity$ES_count))

g <- ggplot(validity,aes(System_T,Criteria,fill=High_quality_pc))+
  geom_tile()+
  geom_text(aes(label=High_quality_pc))+
  ylab("")+xlab("")+
  scale_fill_gradient(low="white", high="forestgreen",name="Number of\neffect sizes\nmeeting criteria (%)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
g

tiff("Fig 4_Validity assessment.tif",width=1000,height=500,units="px",compression="lzw",res=100)
g
dev.off()

png("Fig 4_Validity assessment.png",width=900,height=500,res=100)
g
dev.off()

validity.high <- validity %>%
  dcast(Criteria~System_T,value.var=c("High_quality_pc"))

validity.low <- validity %>%
  dcast(Criteria~System_T,value.var=c("Low_quality_pc"))

validity.unknown <- validity %>%
  dcast(Criteria~System_T,value.var=c("Unknown_quality_pc"))

validity.table <- merge.data.frame(validity.high,validity.low,
                                   by.x="Criteria",by.y="Criteria",suffixes=c(" High"," Low"))
validity.table <- merge.data.frame(validity.table,validity.unknown,
                                   by.x="Criteria",by.y="Criteria",suffixes=c(""," Unknown"))
names(validity.table)
#validity.table <- validity.table[,c(1,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14,21,8,15,22)]

write.xlsx(validity.table, paste0("Table_validity_assessment_",format(Sys.time(),"%Y%m%d"),".xlsx"),row.names=FALSE)

## Make table showing number of effects per category ####
data_meta %>%
  #group_by(System_C) %>% 
  summarise(ES_biodiversity = sum(!is.na(B_value_T)),
            Articles_biodiversity = n_distinct(ID))
data_meta[(data_meta$Yield_data == 1),] %>%
  summarise(ES_yield = sum(Yield_data == 1),
            Articles_yield = n_distinct(ID))
table(data_meta$System_C)
nrow(articles)
freq(data_meta[data_meta$System_C %in% c("Simplified other","Monoculture"),]$Crop_C)

data_meta <- data_meta %>%
  mutate(Taxa_group_map = ifelse(Taxa_group %in% c("Birds", "Amphibians","Mammals","Reptiles"),Taxa_group,
                                 ifelse(Taxa_group %in% c("Mammals other","Bats"),"Mammals", Taxa_phylum)),
         Functional_group = ifelse(Functional_group == "Weeds","Pests",Functional_group))


freq(data_meta$System_C)
freq(data_meta$Country)
freq(data_meta$System_T)
freq(data_meta$Taxa_phylum)
freq(data_meta$Taxa_group_map)
freq(data_meta$Functional_group)
round(addmargins(table(data_meta$Taxa_phylum,data_meta$Taxa_group)/4108),2)
freq(data_meta$Taxa_class)
freq(data_meta$Taxa_order)
freq(data_meta$Crop_T)
freq(data_meta$Crop_FAO_T)
freq(data_meta$Yield_data)

#check <- data_meta[!(is.na(data_meta$Yield_value_T)),]
# Make figure of crops per country
#check <- data_meta[which(data_meta$Crop_FAO_T %in% c("Other or  nd" ,"Other/nd")),]
#check <- unique(check[,c("Crop_T","Crop_FAO_T")])

#check <- data_meta[which(data_meta$Comparison_class_C == "Natural" & is.na(data_meta$Crop_FAO_C)),]

data_meta <- data_meta %>%
  mutate(Crop_FAO_C = ifelse(Comparison_class_C == "Natural","Non-crop (Natural)",Crop_FAO_C),
         Crop_ann_pen_C = ifelse(Comparison_class_C == "Natural","Perennial",Crop_ann_pen_C),
         Crop_woodiness_C = ifelse(Comparison_class_C == "Natural" & Crop_C %in% c("Forest"),"Tree",Crop_woodiness_C),
         Crop_woodiness_C = ifelse(Comparison_class_C == "Natural" & Crop_C %in% c("Shrubland","Mangrove"),"Shrub",Crop_woodiness_C),
         Crop_woodiness_C = ifelse(Comparison_class_C == "Natural" & Crop_C %in% c("Wetland"),"Mixed",Crop_woodiness_C))

#d <- data.frame(table(System=data_meta$System_T,Crop=data_meta$Crop_T,Country=data_meta$Country,FAO_crop=data_meta$Crop_FAO_T))
d <- data.frame(table(System=data_meta$System_T,Crop=data_meta$Crop_T,Outcome=data_meta$Yield_data,FAO_crop=data_meta$Crop_FAO_T))
d <- d[(d$Freq!=0),]

d <- d %>%
  mutate(Country = as.character(Country)) %>%
  mutate(Continent= if_else(Country == "Argentina"| Country =="Brazil" | Country =="Ecuador" | Country =="Uruguay" |  Country == "Peru"| Country =="Colombia", "South America",
                            if_else(Country %in% c("Canada","USA", "United States of America","Mexico"), "North America", 
                                    if_else(Country %in% c("Germany","Belgium","France", "Swiss", "Sweden", "Poland","Spain","Italy","Portugal", "United Kingdom", "Finland","Hungary","Switzerland","Netherlands"), "Europe",
                                            if_else(Country == "China" | Country =="Turkey"| Country =="India"| Country =="Indonesia"| Country =="Viet Nam"| Country =="Japan"| Country =="Malaysia"| Country == "Israel"|Country =="Philippines", "Asia",
                                                    if_else(Country =="Cameroon"| Country =="Egypt"| Country == "Ghana"| Country == "Nigeria"| Country =="South Africa"| Country =="Kenya"| Country =="Malawi" | Country =="Uganda"|Country == "Ethiopia"| Country =="Sao Tome and Principe"| Country =="Benin"| Country == "Zambia", "Africa",
                                                            if_else(Country =="Costa Rica"| Country == "Panama"| Country =="Guatemala"| Country =="Nicaragua", "Central America",
                                                                    if_else(Country =="Dominican Republic" | Country == "Jamaica", "Central America",
                                                                            if_else(Country == "New Zealand"|Country =="Australia", "Oceania",
                                                                                    Country)))))))))


addmargins(table(d$Continent))
addmargins(table(d$Outcome))

d <- d %>% 
  mutate(FAO_crop_N_total = sum(Freq)) %>%
  group_by(FAO_crop) %>%
  mutate(FAO_crop_N = sum(Freq),
         FAO_crop_pc = round(sum(Freq)/FAO_crop_N_total*100,1)) %>%
  ungroup()

sort(unique(d$FAO_crop))
d <- d %>%
  mutate(FAO_crop = factor(FAO_crop,levels=c("1 - CEREALS AND CEREAL PRODUCTS",                    
                                             "2 - ROOTS AND TUBERS AND DERIVED PRODUCTS"  ,        
                                             "3 - SUGAR CROPS AND SWEETENERS AND DERIVED PRODUCTS",
                                             "4 - PULSES AND DERIVED PRODUCTS"     ,               
                                             "5 - NUTS AND DERIVED PRODUCTS"   ,                   
                                             "6 - OIL-BEARING CROPS AND DERIVED PRODUCTS"  ,       
                                             "7 - VEGETABLES AND DERIVED PRODUCTS" ,               
                                             "8 - FRUITS AND DERIVED PRODUCTS"   ,                 
                                             "9 - FIBRES OF VEGETAL AND ANIMAL ORIGIN"   ,
                                              #"10 - SPICES"  ,                                      
                                             "11 - FODDER CROPS AND PRODUCTS"  ,                  
                                             "12 - STIMULANT CROPS AND DERIVED PRODUCTS"  ,        
                                             "13 - TOBACCO AND RUBBER AND OTHER CROPS", 
                                            "Other or nd" 
                                            #"Non-crop (Natural)"
                                            ),
                           labels=c("Cereals",
                                    "Roots & tubers",
                                    "Sugar",
                                    "Pulses",
                                    "Nuts",
                                    "Oil-bearing crops",
                                    "Vegetables",
                                    "Fruits",
                                    "Fibres",
                                    #"Spices",
                                    "Fodder",
                                    "Stimulants",
                                    "Rubber",
                                    "Other/nd"
                                    #"N/A (natural habitat)"
                                    )))

d <- d %>%
  mutate(FAO_crop_label = paste0(FAO_crop," (",FAO_crop_N,", ", FAO_crop_pc,"%)"),
         Comparison_type = ifelse(FAO_crop == "N/A (natural habitat)","Natural habitat","Simplified farming system"))%>%
  mutate(FAO_crop_label = ifelse(FAO_crop_pc<0.1,paste0(FAO_crop," (",FAO_crop_N,", <0.1%)"),FAO_crop_label),
         FAO_crop_label = factor(FAO_crop_label,levels=unique(FAO_crop_label[order(Comparison_type,FAO_crop_label)])))
  #     Continent = factor(Continent,levels=c("Africa","South America","Central America","North America","Asia","Europe","Oceania")))

#d <- d %>%
#  mutate(Continent = factor(Continent),
#         Country = factor(Country,levels=unique(Country[order(Country,Continent)])),
#         Crop = factor(Crop,levels=unique(Crop[order(FAO_crop)])))

write.csv(d,"Biodiversity local data crop-system-outcome freq.csv",row.names = FALSE)

length(unique(d$Crop))

g <- ggplot(d,aes(x=FAO_crop,y=reorder(Country,desc(Country)),size=Freq,colour=System))+
  geom_point(alpha=0.7)+
  scale_size(range = c(1.6, 7), name="Number of effect sizes")+
  ylab("")+xlab("")+
  #scale_colour_viridis(discrete=TRUE, option="D") +
  scale_fill_manual(values="BrBG",breaks=c(0,25,50,100,150))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
g

png("Crop group-country-treatment bubble plot.png", width=800,height=1000,unit="px",res=110)
g
dev.off()

g <- ggplot(d,aes(x=System,y=reorder(Crop,desc(Crop)),size=Freq,colour=Outcome))+
  geom_point(alpha=0.7)+
  scale_size(range = c(1.7, 7), name="Number of comparisons")+
  ylab("")+xlab("")+
 # scale_colour_viridis(discrete=TRUE, option="D",direction=-1) +
  #scale_colour_manual(values=c("purple","seagreen","orange", "gold","lightblue","steelblue","navy"),name="Region")+
  #scale_colour_manual(values=c("black","grey60"),labels=c("Biodiversity","Biodiversity & yield"), name="Outcomes measured")+
  scale_colour_manual(values=c("purple","gold"),labels=c("Biodiversity","Biodiversity & yield"), name="Outcomes measured")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))+
  facet_grid(rows=vars(FAO_crop_label), 
             #switch="x",
             #strip.position = "right", 
             space="free",
             scales = "free")+
  theme(panel.spacing = unit(0, "lines"),
        strip.text.y=element_text(angle=0,hjust=0,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside")
g

#png("Crop-continent-treatment bubble plot_gold.png", width=1000,height=1300,unit="px",res=120)
png("Crop-outcome-treatment bubble plot_gold.png", width=1000,height=1300,unit="px",res=120)
g
dev.off()

#tiff("Crop-continent-treatment bubble plot_gold.tif", width=1000,height=1300,unit="px",compression="lzw", res=120)
tiff("Crop-outcome-treatment bubble plot_gold.tif", width=1000,height=1300,unit="px",compression="lzw", res=120)
g
dev.off()

### Taxa group x Functional group frequency table
d_taxa <- data.table(addmargins(table(Taxa_group_map=data_meta$Taxa_group_map,Functional_group=data_meta$Functional_group)))

d_taxa <- d_taxa %>% setDT %>% dcast.data.table(Taxa_group_map~Functional_group)

d_taxa.studies <- unique(data_meta[,c("ID","Taxa_group_map","Functional_group")])  %>%
  group_by(Taxa_group_map,Functional_group) %>% tally() %>% rename(n_studies=n) %>%
  setDT %>% dcast.data.table(Taxa_group_map~Functional_group)

d_taxa.studies <- unique(data_meta[,c("ID","Taxa_group_map","Functional_group")]) 
d_taxa.studies <- data.table(addmargins(table(Taxa_group_map=d_taxa.studies$Taxa_group_map,Functional_group=d_taxa.studies$Functional_group)))
d_taxa.studies <- d_taxa.studies %>% setDT %>% dcast.data.table(Taxa_group_map~Functional_group)

d.taxa.yields <- data_meta[which(data_meta$Yield_data ==1),]
d.taxa.yields <- data.table(addmargins(table(Taxa_group_map=d.taxa.yields$Taxa_group_map,Functional_group=d.taxa.yields$Functional_group)))
d.taxa.yields <- d.taxa.yields %>% setDT %>% dcast.data.table(Taxa_group_map~Functional_group)

d_taxa_cbind <- inner_join(d_taxa,d_taxa.studies,by=c("Taxa_group_map"), suffix=c(".x",".y")) %>%
  mutate(Autotrophs = ifelse(Autotrophs.y!=0, paste0(Autotrophs.x," (",Autotrophs.y,")"),Autotrophs.x),
         Decomposers = ifelse(Decomposers.y!=0, paste0(Decomposers.x," (",Decomposers.y,")"),Decomposers.x),
         Frugivores = ifelse(Frugivores.y!=0, paste0(Frugivores.x," (",Frugivores.y,")"),Frugivores.x),
         Granivores = ifelse(Granivores.y!=0, paste0(Granivores.x," (",Granivores.y,")"),Granivores.x),
         Herbivores = ifelse(Herbivores.y!=0, paste0(Herbivores.x," (",Herbivores.y,")"),Herbivores.x),
         Insectivores = ifelse(Insectivores.y!=0, paste0(Insectivores.x," (",Insectivores.y,")"),Insectivores.x),
         `Natural enemies` = ifelse(`Natural enemies.y`!=0, paste0(`Natural enemies.x`," (",`Natural enemies.y`,")"),`Natural enemies.x`),
         Omnivores = ifelse(Omnivores.y!=0, paste0(Omnivores.x," (",Omnivores.y,")"),Omnivores.x),
         Other = ifelse(Other.y!=0, paste0(Other.x," (",Other.y,")"),Other.x),
         Pests = ifelse(Pests.y!=0, paste0(Pests.x," (",Pests.y,")"),Pests.x),
         Pollinators = ifelse(Pollinators.y!=0, paste0(Pollinators.x," (",Pollinators.y,")"),Pollinators.x),
         Total = ifelse(Sum.y!=0,paste0(Sum.x," (",Sum.y,")"),Sum.x)) %>%
  select(Taxa_group_map,Autotrophs, Decomposers ,Frugivores,  Granivores,  Herbivores , Insectivores, `Natural enemies` , 
         Omnivores,Other,Pests,Pollinators,Total) %>%
  mutate(Taxa_group_map = ifelse(Taxa_group_map=="Sum","Total",Taxa_group_map))

write.csv(d_taxa_cbind,"Taxa group by function.csv",row.names=F)

d_taxa_cbind <- full_join(d_taxa, d.taxa.yields,by=c("Taxa_group_map"), suffix=c(".x",".y"))%>%
  mutate(Autotrophs = ifelse(Autotrophs.y!=0 & !is.na(Autotrophs.y), paste0(Autotrophs.x," (",Autotrophs.y,")"),Autotrophs.x),
         Decomposers = ifelse(Decomposers.y!=0 & !is.na(Autotrophs.y), paste0(Decomposers.x," (",Decomposers.y,")"),Decomposers.x),
         Frugivores = ifelse(Frugivores.y!=0 & !is.na(Frugivores.y), paste0(Frugivores.x," (",Frugivores.y,")"),Frugivores.x),
         Granivores = ifelse(Granivores.y!=0 & !is.na(Granivores.y), paste0(Granivores.x," (",Granivores.y,")"),Granivores.x),
         Herbivores = ifelse(Herbivores.y!=0 & !is.na(Herbivores.y), paste0(Herbivores.x," (",Herbivores.y,")"),Herbivores.x),
         Insectivores = Insectivores,
         `Natural enemies` = ifelse(`Natural enemies.y`!=0 & !is.na(`Natural enemies.y`), paste0(`Natural enemies.x`," (",`Natural enemies.y`,")"),`Natural enemies.x`),
         Omnivores = ifelse(Omnivores.y!=0 & !is.na(Omnivores.y), paste0(Omnivores.x," (",Omnivores.y,")"),Omnivores.x),
         Other = ifelse(Other.y!=0 & !is.na(Other.y), paste0(Other.x," (",Other.y,")"),Other.x),
         Pests = ifelse(Pests.y!=0 & !is.na(Pests.y), paste0(Pests.x," (",Pests.y,")"),Pests.x),
         Pollinators = ifelse(Pollinators.y!=0 & !is.na(Pollinators.y), paste0(Pollinators.x," (",Pollinators.y,")"),Pollinators.x),
         Total = ifelse(Sum.y!=0 & !is.na(Sum.y),paste0(Sum.x," (",Sum.y,")"),Sum.x)) %>%
  select(Taxa_group_map,Autotrophs, Decomposers ,Frugivores,  Granivores,  Herbivores , Insectivores, `Natural enemies` , 
         Omnivores,Other,Pests,Pollinators,Total) %>%
  mutate(Taxa_group_map = ifelse(Taxa_group_map=="Sum","Total",Taxa_group_map))

write.csv(d_taxa_cbind,"Taxa group by function_yields.csv",row.names=F)

# JOINT CODE ENDS HERE !!! ####
### prepare subgroups - SPECIFIC TO SARAHs analysis only####

# biodiversity measure ####
sort(unique(data_meta$B_measure))
data_meta$B_measure_group <- as.character(data_meta$B_measure)
data_meta$B_measure_group <- ifelse(data_meta$B_measure_group %in% c("Activity-density","Abundance","Abundanceweighted mean CWM","Vistitation frequency"),"Abundance",data_meta$B_measure_group) 
data_meta$B_measure_group <- ifelse(data_meta$B_measure_group %in% c("Species Richness", "Rareified Species Richness","Jack-knife species richness","Number of orders"),"Richness",data_meta$B_measure_group) 
data_meta$B_measure_group <- ifelse(data_meta$B_measure_group %in%  c("Species Eveness ","Shannon Eveness Index", "Shannon Index" , "Shannon Index " ,   "Shannon-Wiener Index"),"Other",data_meta$B_measure_group) 
#data_meta$B_measure_group <- ifelse(data_meta$B_measure_group %in%  c("Simpson's reciprocal index","Simpson Index "),"Simpsons",data_meta$B_measure_group) 
data_meta$B_measure_group <- ifelse(data_meta$B_measure_group %in%  c("Chao1 Index", "Pielou","Pielou Index" , "Fisher alpha", "Margalef Index " , "Jaccard similarity index",
                                                                      "Simpson","Simpson Index", "Simpson's reciprocal index","Simpson Index ","Trophic Diversity Index","Colonization Percent",
                                                                      "Colonization Percent ","Area under disease progress curve (AUDPC)"),"Other",data_meta$B_measure_group) 
sort(unique(data_meta$B_measure_group))
data_meta$B_measure_group <- factor(data_meta$B_measure_group,levels=c("Abundance","Richness", "Other"),labels=c("Abundance","Richness", "Other"))

# taxa groups ####
unique(data_meta$Taxa_group)
table(data_meta$Taxa_group)
data_meta$Taxa_group <- factor(data_meta$Taxa_group)

# FAO crop groups ####
sort(unique(data_meta$Crop_FAO_C))
funModeling::freq(data_meta$Crop_FAO_C)
# merge classes with very small sizes (n<10)
# pulses = 4
# spices = 8
# sugar = 2
data_meta$Crop_FAO_C[data_meta$Crop_FAO_C %in% c("10 - SPICES" , 
                                                 "4 - PULSES AND DERIVED PRODUCTS",
                                       "Other or nd",
                                       "3 - SUGAR CROPS AND SWEETENERS AND DERIVED PRODUCTS")] <- "OTHER"

data_meta$Crop_FAO_C[is.na(data_meta$Crop_FAO_C)] <- "OTHER"
data_meta$Crop_FAO_C[data_meta$Crop_FAO_C %in% c("nd","Other or  nd")] <- "OTHER"
funModeling::freq(data_meta$Crop_FAO_C)
table(data_meta$Crop_FAO_C)
data_meta$Crop_FAO_C <- factor(data_meta$Crop_FAO_C , 
                          levels = c( "1 - CEREALS AND CEREAL PRODUCTS"          , 
                                      "2 - ROOTS AND TUBERS AND DERIVED PRODUCTS", 
                                      "5 - NUTS AND DERIVED PRODUCTS"    ,         
                                      "6 - OIL-BEARING CROPS AND DERIVED PRODUCTS", 
                                      "7 - VEGETABLES AND DERIVED PRODUCTS" ,   
                                      "8 - FRUITS AND DERIVED PRODUCTS"       , 
                                      "9 - FIBRES OF VEGETAL AND ANIMAL ORIGIN"  , 
                                      "11 - FODDER CROPS AND PRODUCTS"    ,      
                                      "12 - STIMULANT CROPS AND DERIVED PRODUCTS" ,
                                       "OTHER",
                                      "Non-crop (Natural)"),
                          labels=c( "Cereals"          , 
                                    "Roots & tubers",
                                    "Nuts"    ,         
                                    "Oil-bearing crops", 
                                    "Vegetables" ,   
                                    "Fruits"       , 
                                    "Fibres"  , 
                                    "Fodder"    ,      
                                    "Stimulants" ,
                                    "Mixed/other",
                                    "Non-crop (Natural)"))

sort(unique(data_meta$Crop_FAO_T))
funModeling::freq(data_meta$Crop_FAO_T)
# merge classes with very small sizes (n<10)
# pulses = 4
# rubber = 4
# sugar = 2
data_meta$Crop_FAO_T[data_meta$Crop_FAO_T %in% c("4 - PULSES AND DERIVED PRODUCTS",
                                                 "13 - TOBACCO AND RUBBER AND OTHER CROPS",
                                                 "Other or nd",
                                                 "3 - SUGAR CROPS AND SWEETENERS AND DERIVED PRODUCTS")] <- "OTHER"

data_meta$Crop_FAO_T[is.na(data_meta$Crop_FAO_T)] <- "OTHER"
data_meta$Crop_FAO_T[data_meta$Crop_FAO_T %in% c("nd","Other or  nd")] <- "OTHER"
funModeling::freq(data_meta$Crop_FAO_T)
table(data_meta$Crop_FAO_T)
data_meta$Crop_FAO_T <- factor(data_meta$Crop_FAO_T , 
                               levels = c( "1 - CEREALS AND CEREAL PRODUCTS"          , 
                                           "2 - ROOTS AND TUBERS AND DERIVED PRODUCTS", 
                                           "5 - NUTS AND DERIVED PRODUCTS"    ,         
                                           "6 - OIL-BEARING CROPS AND DERIVED PRODUCTS", 
                                           "7 - VEGETABLES AND DERIVED PRODUCTS" ,   
                                           "8 - FRUITS AND DERIVED PRODUCTS"       , 
                                           "9 - FIBRES OF VEGETAL AND ANIMAL ORIGIN"  , 
                                           "11 - FODDER CROPS AND PRODUCTS"    ,      
                                           "12 - STIMULANT CROPS AND DERIVED PRODUCTS" ,
                                           "OTHER"),
                               labels=c( "Cereals"          , 
                                         "Roots & tubers",
                                         "Nuts"    ,         
                                         "Oil-bearing crops", 
                                         "Vegetables" ,   
                                         "Fruits"       , 
                                         "Fibres"  , 
                                         "Fodder"    ,      
                                         "Stimulants" ,
                                         "Mixed/other"))
# Pesticide C-T ####
data_meta$Pesticide_CT <- paste0(as.character(data_meta$Pesticide_C),":",as.character(data_meta$Pesticide_T))
sort(unique(data_meta$Pesticide_CT))
data_meta$Pesticide_CT[data_meta$Pesticide_CT %in% c("nd:nd","nd:NA")] <- "nd"
data_meta$Pesticide_CT[data_meta$Pesticide_CT %in% c("nd:No","No:nd","No:No","No:NA")] <- "No"
data_meta$Pesticide_CT[data_meta$Pesticide_CT %in% c("nd:Yes","Yes:nd","Yes:Yes")] <- "Yes"
data_meta$Pesticide_CT[data_meta$Pesticide_CT %in% c("No:Yes","Yes:No")] <- "Mixed"
sort(unique(data_meta$Pesticide_CT))
funModeling::freq(data_meta$Pesticide_CT)
data_meta$Pesticide_CT <- factor(data_meta$Pesticide_CT,levels=c("No","Mixed","Yes", "nd"))

# Tillage C-T ####
data_meta$Tillage_CT <- paste0(as.character(data_meta$Soil_management_C),":",as.character(data_meta$Soil_management_T))
sort(unique(data_meta$Tillage_CT))
data_meta$Tillage_CT[data_meta$Tillage_CT %in% c("nd:nd","Burning:Burning","nd:NA","NA:NA","NA:nd")] <- "nd"
data_meta$Tillage_CT[data_meta$Tillage_CT %in% c("nd:No Tillage","No Tillage:nd","No Tillage:No Tillage","NA:No Tillage")] <- "No"
data_meta$Tillage_CT[data_meta$Tillage_CT %in% c("nd:Tillage","Tillage:nd","Tillage:Tillage","NA:Tillage")] <- "Yes"
data_meta$Tillage_CT[data_meta$Tillage_CT %in% c("No Tillage:Tillage","Tillage:No Tillage")] <- "Mixed"
sort(unique(data_meta$Tillage_CT))
funModeling::freq(data_meta$Tillage_CT)
data_meta$Tillage_CT <- factor(data_meta$Tillage_CT,levels=c("No","Mixed","Yes", "nd"))

# Fertiliser C-T ####
data_meta$Fertiliser_CT <- paste0(as.character(data_meta$Fertiliser_C),":",as.character(data_meta$Fertiliser_T))
sort(unique(data_meta$Fertiliser_CT))
data_meta$Fertiliser_CT[data_meta$Fertiliser_CT %in% c("nd:nd")] <- "nd"
data_meta$Fertiliser_CT[data_meta$Fertiliser_CT %in% c("nd:No","No:nd","No:No","No:NA")] <- "No"
data_meta$Fertiliser_CT[data_meta$Fertiliser_CT %in% c("nd:Yes","Yes:nd","Yes:Yes","Yes:NA")] <- "Yes"
data_meta$Fertiliser_CT[data_meta$Fertiliser_CT %in% c("No:Yes","Yes:No")] <- "Mixed"
sort(unique(data_meta$Fertiliser_CT))
funModeling::freq(data_meta$Fertiliser_CT)
data_meta$Fertiliser_CT <- factor(data_meta$Fertiliser_CT,levels=c("No","Mixed","Yes", "nd"))

# Diversity practice ####
unique(data_meta$System_T)
funModeling::freq(data_meta$System_T)
data_meta$System_T <- factor(data_meta$System_T,levels=c("Agroforestry",
                                                 "Crop rotation",
                                                 "Cover crops",
                                                 "Intercropping",
                                                 "Embedded natural",
                                                 "Diversified other"))

# Crop habitat C-T ####
#data_meta$Crop_ann_pen_CT<- paste0(as.character(data_meta$Crop_ann_pen_C),":",as.character(data_meta$Crop_ann_pen_T))
#sort(unique(data_meta$Crop_ann_pen_CT))
#data_meta$Crop_ann_pen_CT[data_meta$Crop_ann_pen_CT %in% c("nd:nd")] <- "nd"
#data_meta$Crop_woodiness_CT <- paste0(as.character(data_meta$Crop_woodiness_C),":",as.character(data_meta$Crop_woodiness_T))
#sort(unique(data_meta$Crop_woodiness_CT))
sort(unique(data_meta$Crop_ann_pen_T))
funModeling::freq(data_meta$Crop_ann_pen_T)
data_meta$Crop_ann_pen_T <- as.character(data_meta$Crop_ann_pen_T)
data_meta[data_meta$Crop_ann_pen_T == "Biennual",]$Crop_T
data_meta$Crop_ann_pen_T[data_meta$Crop_ann_pen_T %in% c("Biennual")] <- "Annual"
data_meta$Crop_ann_pen_T[data_meta$Crop_ann_pen_T %in% c("Unknown","Mixed")] <-  "Mixed or ns"
data_meta$Crop_ann_pen_T[is.na(data_meta$Crop_ann_pen_T)] <- "Mixed or ns"
funModeling::freq(data_meta$Crop_ann_pen_T)

sort(unique(data_meta$Crop_ann_pen_C))
check <- data_meta[which(data_meta$Crop_ann_pen_C =="#N/A"),]
data_meta$Crop_ann_pen_C <- as.character(data_meta$Crop_ann_pen_C)
funModeling::freq(data_meta$Crop_ann_pen_C)
data_meta[data_meta$Crop_ann_pen_C == "Biennual",]$Crop_T
data_meta$Crop_ann_pen_C[data_meta$Crop_ann_pen_C %in% c("Biennual")] <- "Annual"
data_meta$Crop_ann_pen_C[data_meta$Crop_ann_pen_C %in% c("Unknown","Mixed","#N/A")] <-  "Mixed or ns"
data_meta$Crop_ann_pen_C[is.na(data_meta$Crop_ann_pen_C)] <- "Mixed or ns"
funModeling::freq(data_meta$Crop_ann_pen_C)

sort(unique(data_meta$Crop_woodiness_T))
data_meta$Crop_woodiness_T <- as.character(data_meta$Crop_woodiness_T)
funModeling::freq(data_meta$Crop_woodiness_T)
data_meta$Crop_woodiness_T[data_meta$Crop_woodiness_T %in% c("Unknown","Mixed")] <-  "Mixed or ns"
data_meta$Crop_woodiness_T[is.na(data_meta$Crop_woodiness_T)] <- "Mixed or ns"
funModeling::freq(data_meta$Crop_woodiness_T)

sort(unique(data_meta$Crop_woodiness_C))
data_meta$Crop_woodiness_C <- as.character(data_meta$Crop_woodiness_C)
funModeling::freq(data_meta$Crop_woodiness_C)
data_meta$Crop_woodiness_C[data_meta$Crop_woodiness_C %in% c("Unknown","Mixed")] <-  "Mixed or ns"
data_meta$Crop_woodiness_C[is.na(data_meta$Crop_woodiness_C)] <- "Mixed or ns"
funModeling::freq(data_meta$Crop_woodiness_C)

sort(unique(data_meta$Crop_woodiness_T))
data_meta$Crop_woodiness_T <- as.character(data_meta$Crop_woodiness_T)
funModeling::freq(data_meta$Crop_woodiness_T)
data_meta$Crop_woodiness_T[data_meta$Crop_woodiness_T %in% c("Unknown","Mixed")] <-  "Mixed or ns"
data_meta$Crop_woodiness_T[is.na(data_meta$Crop_woodiness_T)] <- "Mixed or ns"
funModeling::freq(data_meta$Crop_woodiness_T)

### Prepare random effects identifiers ####
# same study = data_meta$ID 
length(unique(data_meta$ID)) 

# same control
data_meta$C_ID <- factor(paste0(data_meta$ID,data_meta$Comparison_ID_C,data_meta$Taxa_group,data_meta$Taxa_details,data_meta$B_measure_group,data_meta$Fertiliser_CT, data_meta$Pesticide_CT,data_meta$Crop_FAO_C))
length(unique(data_meta$C_ID)) #1286

# same control, same treatment (so multiple rows should represent repeated measures)
data_meta$C_T_ID <- factor(paste0(data_meta$C_ID,data_meta$Comparison_ID_T))
length(unique(data_meta$C_T_ID)) # 2240

# export distributions for all data_meta 
funModeling::freq(data_meta, path_out = "./AllComparisons") # saves plots to working directory

# NOTE:
# Control, # 2 classes, can have multiple rows representing diff controls in same class, e.g. maize monoculture, wheat monoculture both in monoculture class
# Treatment, # 6 classes, can have multiple rows representing diff treatments within same class, e.g. hedgerow, flower strip both in embedded natural class
# Taxa_group, # 14 classes, can have multipel rows represting diff taxa within same class, e.g. wasp, bee both in wasp and bee class
# B_measure, # 4 classes, can have multiple rows representing diff measures within same class, e.g. two different evenness metrics
# Crop_FAO_C) # 8 classes, can have multiple rows represtenting diff crops within same class, e.g. maize, wheat both in cereal class

# Split - Simplified systems as control ####
sort(unique(data_meta$System_C))
data.divsim <- data_meta[which(data_meta$System_C %in% c("Monoculture","Simplified other")),]
funModeling::freq(data.divsim, path_out = "./Diversified_Simplified") # saves plots to working directory
freq(data.divsim$Functional_group)
#data.divsim$Functional_group <- as.character(data.divsim$Functional_group)
#data.divsim$Taxa_pest_group <- paste0(data.divsim$Taxa_group," - ",data.divsim$Functional_group)
#sort(unique(data.divsim$Pest_group))
funModeling::freq(data.divsim$Pest_group)
#data.divsim$Pest_group <- ifelse(data.divsim$Pest_group %in% c("Arachnids - Pests","Arachnids - Other"),"Arachnids - mainly Pests",data.divsim$Taxa_pest_group)
funModeling::freq(data.divsim$Taxa_group)
# group taxa with n<10
# reptiles = 1
# plants woody = 6
data.divsim <- data.divsim %>%
  mutate(Taxa_group = as.character(Taxa_group),
         Taxa_group = ifelse(Taxa_group %in% c("Reptiles","Amphibians"), "Amphibians, reptiles",
                             ifelse(Taxa_group %in% c("Plants woody", "Plants other"),"Plants",Taxa_group)))
funModeling::freq(data.divsim$Taxa_group)

length(unique(data.divsim$C_ID)) # 932
length(unique(data.divsim$C_T_ID)) # 1632 meaning that some controls are compared with multiple treatments

# format experiment stage column which represents timesteps in repeated measures entries
unique(data.divsim$Experiment_stage)
data.divsim$Experiment_stage[(data.divsim$Experiment_stage) %in% c(NA)] <- 1
data.divsim$Experiment_stage <- as.numeric(data.divsim$Experiment_stage)
unique(data.divsim$Experiment_stage)

# unique effects 
# NOTE: Merge_ID column is not unique because same ID is given for 
# treatments compared against different controls, or 
# controls compared against multiple treatments 
data.divsim <- unique(data.divsim) 
data.divsim$Effect_ID <- factor(1:nrow(data.divsim))
length(unique(data.divsim$Effect_ID)) # 3381

# export to be able to share full original data with others
write.csv(data.divsim,paste0("MA biodiversity data_divsim.csv"),row.names=F)
#write.csv(data.divsim[data.divsim$Pest_group=="Other",],"MA biodiversity data_divsim_nopests.csv",row.names=F)

# Split - Natural systems as control ####
unique(data_meta$System_C)
data.divnat <-  data_meta[which(data_meta$System_C %in% c("Abandoned","Natural")),]
funModeling::freq(data.divnat, path_out = "./Diversified_Natural") # saves plots to working directory
funModeling::freq(data.divnat$Pest_group) # only 11 pests and weeds...
funModeling::freq(data.divnat$Taxa_group) 
# merge or remove taxa classes with very small sizes (n<10)
# nematodes and annelids = 2
# fungi = 8
data.divnat <- data.divnat %>%
  mutate(Taxa_group = as.character(Taxa_group),
         Taxa_group = ifelse(Taxa_group %in% c("Fungi","Bacteria"),"Bacteria, fungi",Taxa_group))%>%
  filter(Taxa_group != "Nematodes and annelids")

unique(data.divnat$Experiment_stage)
data.divnat$Experiment_stage[(data.divnat$Experiment_stage) %in% c(NA)] <- 1
data.divnat$Experiment_stage <- as.numeric(data.divnat$Experiment_stage)
unique(data.divnat$Experiment_stage)

data.divnat <- unique(data.divnat) 
data.divnat$Effect_ID <- factor(1:nrow(data.divnat))
length(unique(data.divnat$Effect_ID)) #730

names(data.divnat)

# export to be able to share full original data with Marc and others
write.csv(data.divnat,paste0("MA biodiversity data_divnat.csv"),row.names=F)
#write.csv(data.divnat[data$Pest_group == "Other",],paste0("MA biodiversity data_divnat_nopests.csv"),row.names=F)


### Make maps ####
#https://www.datanovia.com/en/blog/how-to-create-a-map-using-ggplot2/
##https://datacarpentry.org/r-raster-vector-geospatial/10-vector-csv-to-shapefile-in-r/
#https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
#Label: https://www.datanovia.com/en/blog/how-to-change-ggplot-labels/
#Colors: https://color.adobe.com/search?q=tree

# Map showing outcome measured
#data_meta <- read.csv(paste0("C_biodiversity_local_data_meta_",format(Sys.time(),"%Y%m%d"), ".csv"),row.names=FALSE)

map.data <- data_meta %>%
  mutate(Country = as.character(Country)) %>%
  mutate(Country = ifelse(Country %in% c("United States of America"),"USA",
                                         ifelse(Country %in% c("United Kingdom"),"UK",Country))) %>%
  mutate(Continent= if_else(Country == "Argentina"| Country =="Brazil" | Country =="Ecuador" | Country =="Uruguay" |  Country == "Peru"| Country =="Colombia", "South America",
                            if_else(Country == "Canada"|Country =="USA" |Country =="Mexico", "North America", 
                                    if_else(Country =="Germany"| Country =="Belgium"|Country == "France"|Country == "Swiss"| Country == "Sweden"| Country == "Poland"| Country =="Spain"| Country =="Italy"| Country =="Portugal"| Country ==  "UK" | Country =="Finland" | Country =="Hungary"| Country =="Switzerland"| Country== "Netherlands", "Europe",
                                            if_else(Country == "China" | Country =="Turkey"| Country =="India"| Country =="Indonesia"| Country =="Viet Nam"| Country =="Japan"| Country =="Malaysia"| Country == "Israel"|Country =="Philippines", "Asia",
                                                    if_else(Country =="Cameroon"| Country =="Egypt"| Country == "Ghana"| Country == "Nigeria"| Country =="South Africa"| Country =="Kenya"| Country =="Malawi" | Country =="Uganda"|Country == "Ethiopia"| Country =="Sao Tome and Principe"| Country =="Benin"| Country == "Zambia", "Africa",
                                                            if_else(Country =="Costa Rica"| Country == "Panama"| Country =="Guatemala"| Country =="Nicaragua", "North America",
                                                                    if_else(Country =="Dominican Republic" | Country == "Jamaica", "North America",
                                                                            if_else(Country == "New Zealand"|Country =="Australia", "Oceania",
                                                                                    Country)))))))),
         Taxa_mapgroup = Taxa_group) %>%
         #Taxa_mapgroup = ifelse(Taxa_phylum=="Chordata",Taxa_class,Taxa_phylum)) %>%
  #mutate(Taxa_mapgroup = ifelse(Taxa_mapgroup == "Amphibian","Amphibians",
  #                              ifelse(Taxa_mapgroup == "Annelida","Annelids",
  #                                     ifelse(Taxa_mapgroup == "Arthropoda","Arthropods",
  #                                            ifelse(Taxa_mapgroup == "Nematoda","Nematodes",
  #                                                  ifelse(Taxa_mapgroup == "Plantae","Plants",Taxa_mapgroup)))))) %>%
  #group_by(ID, Taxa_phylum, Taxa_mapgroup, Country, Continent,Lat_T, Long_T) %>% tally()
  group_by(ID, Yield_data, Country, Continent,Lat_T, Long_T) %>% tally() %>% rename(n_effects = n)

map.data.studies <- unique(data_meta[,c("ID","Country")]) %>% 
  mutate(Country = as.character(Country),
         Country = ifelse(Country %in% c("United States of America"),"USA",
                          ifelse(Country %in% c("United Kingdom"),"UK",Country))) %>%
  group_by(Country) %>% tally() %>% rename(n_studies=n)

unique(map.data$Yield_data)

map.data.country <- map.data %>%
  group_by(Country,Continent) %>%
  summarise(n_effects=sum(n_effects))
map.data.country <- map.data.country %>%
  left_join(map.data.studies)

ggplot(map.data.country,aes(x=n_effects,fill=n_studies))+
  geom_bar()

sort(unique(map.data.country$n_effects))
sort(unique(map.data.country$n_studies))

map.data.country <- map.data.country %>%
  mutate(n_effects_classes = ifelse(n_effects<25,"1-25",
                                    ifelse(n_effects<51,"26-50",
                                           ifelse(n_effects<101,"51-100",
                                                  ifelse(n_effects<201,"101-200",
                                                 ifelse(n_effects<501,"201-500","898")))))) %>%
  mutate(n_studies_classes = ifelse(n_studies<2,"1",
                                           ifelse(n_studies<11,"2-10",
                                                  ifelse(n_studies<31,"11-30","52")))) %>%
  mutate(n_effects_classes = factor(n_effects_classes,levels=c("1-25","26-50","51-100","101-200","201-500","898")),
         n_studies_classes = factor(n_studies_classes,levels=c("1","2-10","11-30","52")))


world_map <- ggplot2::map_data("world")%>%filter(region != "Antarctica")

map.data.joined <- left_join(world_map,map.data.country,by=c("region"="Country")) 
map.data.joined <- map.data.joined %>%
  mutate(n_effects_classes = as.character(n_effects_classes),
         n_studies_classes = as.character(n_studies_classes)) %>%
  mutate(n_effects_classes = ifelse(is.na(n_effects_classes),"0",n_effects_classes),
         n_studies_classes = ifelse(is.na(n_studies_classes),"0",n_studies_classes))# %>%

sort(unique(map.data.joined$n_effects_classes))
sort(unique(map.data.joined$n_studies_classes))

map.data.joined <- map.data.joined %>%
  mutate(n_effects_classes = factor(n_effects_classes,levels=c("0", "1-25","26-50","51-100","101-200","201-500","918")),
         n_studies_classes = factor(n_studies_classes,levels=c("0", "1","2-10","11-30","52")))


map.data <- map.data %>%
  mutate(Yield_data = factor(as.character(Yield_data),levels=c("1","0"),labels=c("Biodiversity & yield", "Biodiversity")))

legend.practices <- ggplot() +
  geom_point(data = map.data, mapping = aes(x=Long_T, y=Lat_T, color = Yield_data,size=n_effects), 
             shape=16,alpha=0.7,show.legend = TRUE)+
  scale_color_manual(values=c("purple","gold"))+
  #scale_color_manual(values=c("purple","blue","lightgreen","forestgreen","orange","gold"))+
  #scale_color_brewer(palette="RdYlBu",name="Diversified farming\nintervention")+
  #scale_color_viridis_d(name="Diversified farming\nintervention")+
  scale_size(range = c(1.6, 7), name="Number of\ncomparisons")+
  #scale_color_brewer(palette="RdYlBu",name="Diversified farming\nsystem")+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "left",
    legend.direction = "vertical",
    legend.justification = c("center", "bottom"),
    legend.text = element_text(size =8),
    legend.title=element_text(size=10),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank())+
  guides(colour = guide_legend(title = "Outcomes measured", title.position = "top",
                             title.hjust = 0,  
                             label.hjust = 0))
legend.practices 

legend.practices<- gtable_filter(ggplot_gtable(ggplot_build(legend.practices)), "guide-box") 

g <- ggplot() +
  geom_polygon(data = map.data.joined, aes(x = long, y = lat, group = group, fill=n_studies_classes), color = "grey")+
  geom_polygon(fill=NA, color = "grey")+
  geom_point(data = map.data, mapping = aes(x=Long_T, y=Lat_T, color = Yield_data,size=n_effects), 
             shape=16,show.legend = FALSE)+
  coord_fixed(1.2)+
  scale_color_manual(values=c("purple","gold"))+
  #scale_color_manual(values=c("purple","blue","lightgreen","forestgreen","orange","gold"))+
  #scale_color_brewer(palette="RdYlBu",name="Diversified farming\nsystem")+
  #scale_color_viridis_d(name="Diversified farming\nsystem")+
  #scale_color_brewer(palette="RdYlBu",name="Diversified farming\nsystem")+
  scale_fill_manual(values=c("white","grey30","grey40","grey60","grey80"),na.value="white")+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c("centre", "bottom"),
    legend.text = element_text(size =8),
    legend.title=element_text(size=10),
    legend.key=element_blank(),
    legend.background = element_blank(),
    panel.background = element_blank())+
  guides(fill = guide_legend(title = "Number of articles", title.position = "top",
                              title.hjust = 0.4,  
                              #ncol = 8,
                              label.hjust = 0))
g

plot <- g +  annotation_custom(grob = legend.practices, xmin = -150, xmax = -150, ymin = -90, ymax = -90)

grid.newpage()
grid.draw(plot)
dev.off()

tiff("Fig map outcome and effect N.tif",width=1100,height=700,res=120)
grid.newpage()
grid.draw(plot)
dev.off()

# Map showing taxa group

map.data <- data_meta %>%
  mutate(Country = as.character(Country)) %>%
  mutate(Country = ifelse(Country %in% c("United States of America"),"USA",
                          ifelse(Country %in% c("United Kingdom"),"UK",Country))) %>%
  mutate(Continent= if_else(Country == "Argentina"| Country =="Brazil" | Country =="Ecuador" | Country =="Uruguay" |  Country == "Peru"| Country =="Colombia", "South America",
                            if_else(Country == "Canada"|Country =="USA" |Country =="Mexico", "North America", 
                                    if_else(Country =="Germany"| Country =="Belgium"|Country == "France"|Country == "Swiss"| Country == "Sweden"| Country == "Poland"| Country =="Spain"| Country =="Italy"| Country =="Portugal"| Country ==  "UK" | Country =="Finland" | Country =="Hungary"| Country =="Switzerland"| Country== "Netherlands", "Europe",
                                            if_else(Country == "China" | Country =="Turkey"| Country =="India"| Country =="Indonesia"| Country =="Viet Nam"| Country =="Japan"| Country =="Malaysia"| Country == "Israel"|Country =="Philippines", "Asia",
                                                    if_else(Country =="Cameroon"| Country =="Egypt"| Country == "Ghana"| Country == "Nigeria"| Country =="South Africa"| Country =="Kenya"| Country =="Malawi" | Country =="Uganda"|Country == "Ethiopia"| Country =="Sao Tome and Principe"| Country =="Benin"| Country == "Zambia", "Africa",
                                                            if_else(Country =="Costa Rica"| Country == "Panama"| Country =="Guatemala"| Country =="Nicaragua", "North America",
                                                                    if_else(Country =="Dominican Republic" | Country == "Jamaica", "North America",
                                                                            if_else(Country == "New Zealand"|Country =="Australia", "Oceania",
                                                                                    Country)))))))),
         Taxa_mapgroup = ifelse(Taxa_phylum=="Chordata",Taxa_class,Taxa_phylum)) %>%
  mutate(Taxa_mapgroup = ifelse(Taxa_mapgroup == "Amphibian","Amphibians",
                                ifelse(Taxa_mapgroup == "Annelida","Annelids",
                                       ifelse(Taxa_mapgroup == "Arthropoda","Arthropods",
                                              ifelse(Taxa_mapgroup == "Nematoda","Nematodes",
                                                    ifelse(Taxa_mapgroup == "Plantae","Plants",Taxa_mapgroup)))))) %>%
  group_by(ID, Taxa_phylum, Taxa_mapgroup, Country, Continent,Lat_T, Long_T) %>% tally() %>% rename(n_effects = n)

map.data.studies <- unique(data_meta[,c("ID","Country")]) %>% 
  mutate(Country = as.character(Country),
         Country = ifelse(Country %in% c("United States of America"),"USA",
                          ifelse(Country %in% c("United Kingdom"),"UK",Country))) %>%
  group_by(Country) %>% tally() %>% rename(n_studies=n)

unique(map.data$Taxa_mapgroup)

map.data.country <- map.data %>%
  group_by(Country,Continent) %>%
  summarise(n_effects=sum(n_effects))
map.data.country <- map.data.country %>%
  left_join(map.data.studies)

ggplot(map.data.country,aes(x=n_effects))+
  geom_histogram()

map.data.country <- map.data.country %>%
  mutate(n_effects_classes = ifelse(n_effects<25,"1-25",
                                    ifelse(n_effects<51,"26-50",
                                           ifelse(n_effects<101,"51-100",
                                                  ifelse(n_effects<201,"101-200",
                                                         ifelse(n_effects<501,"201-500","918")))))) %>%
  mutate(n_studies_classes = ifelse(n_studies<2,"1",
                                    ifelse(n_studies<11,"2-10",
                                           ifelse(n_studies<31,"11-30","52")))) %>%
  mutate(n_effects_classes = factor(n_effects_classes,levels=c("1-25","26-50","51-100","101-200","201-500","918")),
         n_studies_classes = factor(n_studies_classes,levels=c("1","2-10","11-30","52")))


world_map <- ggplot2::map_data("world")%>%filter(region != "Antarctica")

map.data.joined <- left_join(world_map,map.data.country,by=c("region"="Country"))

legend.taxa <- ggplot() +
  geom_point(data = map.data, mapping = aes(x=Long_T, y=Lat_T, color = Taxa_mapgroup,size=n_effects), 
             shape=16,alpha=0.7,show.legend = TRUE)+
  scale_color_manual(values=c("#40004b", "#762a83", "#bf812d","#8c510a", "#dfc27d", "gold", 
                              "#d9f0d3","#a6dba0","#5aae61","#1b7837", "#00441b"))+
  scale_size(range = c(1.6, 7), name="Number of\neffect sizes")+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "left",
    legend.direction = "vertical",
    legend.justification = c("center", "bottom"),
    legend.text = element_text(size =8),
    legend.title=element_text(size=10),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank())+
  guides(colour = guide_legend(title = "Taxonomic group", title.position = "top",
                               title.hjust = 0,  
                               label.hjust = 0))
legend.taxa 
legend.taxa<- gtable_filter(ggplot_gtable(ggplot_build(legend.taxa)), "guide-box") 

g <- ggplot() +
  geom_polygon(data = map.data.joined, aes(x = long, y = lat, group = group, fill=n_studies_classes), color = "grey")+
  geom_polygon(fill=NA, color = "grey")+
  geom_point(data = map.data, mapping = aes(x=Long_T, y=Lat_T, color = Taxa_mapgroup,size=n_effects), 
             shape=16,show.legend = FALSE)+
  coord_fixed(1.2)+
  scale_color_manual(values=c("#40004b", "#762a83", "#bf812d","#8c510a", "#dfc27d", "gold", 
                              "#d9f0d3","#a6dba0","#5aae61","#1b7837", "#00441b"))+
  #scale_color_brewer(palette="PRGn",name="Taxonomic group")+
  scale_size(range = c(1.6, 7), name="Number of\neffect sizes")+
  scale_fill_brewer(palette="Blues",na.value="white")+
 # scale_fill_manual(values=c("grey30","grey40","grey60","grey80"),na.value="white")+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c("centre", "bottom"),
    legend.text = element_text(size =8),
    legend.title=element_text(size=10),
    legend.key=element_blank(),
    legend.background = element_blank(),
    panel.background = element_blank())+
  guides(fill = guide_legend(title = "Number of articles", title.position = "top",
                             title.hjust = 0.4,  
                             #ncol = 8,
                             label.hjust = 0))
g

plot <- g +  annotation_custom(grob = legend.taxa, xmin = -150, xmax = -150, ymin = -90, ymax = -90)

grid.newpage()
grid.draw(plot)

png("Fig map taxa and effect N.png",width=2200,height=1500,res=200)
grid.newpage()
grid.draw(plot)
dev.off()

# Various data checks

freq(data_meta$Country)
freq(data_meta$Taxa_group)
freq(data_meta$Taxa_phylum)
length(unique(data_meta$Country))
length(unique(data_meta$Taxa_order))
