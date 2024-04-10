# notes ----
#Summarize station-level benthic invert mean CPUE to append to snow crab condition datasets
  #to use as a covariate for modeling 

# Author: EJF
# last updated: 2024/3/13 with 2023 groundfish data 

#Note: This script uses a groundfish dataset that is generated via the Snow Crab 
  #ESP indicator development script (see gf_data_pull.R, which queries directly from Racebase)

# load ----
library(tidyverse)
library(mgcv)

# data ----

#Load groundfish data queried directly from Racebase (see gf_data_pull.R script)
benthic <- read_csv("./data/gf_cpue_timeseries.csv")

#Calculate mean CPUE (in kg/km^2) for each guild across years 
benthic %>%
  filter(YEAR > 2018,
    !(SPECIES_CODE %in% c(68560, 68580, 69322, 69323))) %>% #remove commercial crab species 
  group_by(CRUISE, YEAR, STATION, LATITUDE_DD_START, LONGITUDE_DD_START) %>%
  summarise(Gersemia_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(41201:41221)], na.rm = T),
            Pennatulacea_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(42000:42999)], na.rm = T),
            Actinaria_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(43000:43999)], na.rm = T),
            Polychaeta_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(50000:59099)], na.rm = T),
            Barnacles_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(65100:65211)], na.rm = T),
            Shrimps_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(66000:66912)], na.rm = T),
            Crabs_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(68000:69599)], na.rm = T),
            Gastropods_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(71000:73999)], na.rm = T),
            Bivalves_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(74000:75799)], na.rm = T),
            Asteroidea_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(80000:82499)], na.rm = T),
            Echinoidea_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(82500:82729)], na.rm = T),
            Ophiuroidea_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(83000:84999)], na.rm = T),
            Holothuroidea_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(85000:85999)], na.rm = T),
            Porifera_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(91000:91999)], na.rm = T),
            Bryozoans_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(95000:95499)], na.rm = T),
            Ascidians_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(98000:99909)], na.rm = T),
            Total_Benthic_cpue = sum(CPUE_KGKM2[SPECIES_CODE %in% c(41201:99909)], na.rm = T)) %>%
            rename_with(tolower) %>%
            rename(gis_station = station) -> benthic_cpue
  
write.csv(benthic_cpue, file = "./output/benthic_cpue.csv")

#load output for figures
benthic_cpue <- read_csv("./output/benthic_cpue.csv")

#Spatial maps 

#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#total benthic cpue
benthic_cpue %>% 
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = longitude_dd_start, y = latitude_dd_start, size=total_benthic_cpue), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  facet_wrap(~year)
#missing stations in 2023 due to unfinished NBS grid 



