#GOALS: 
# 1) Append haul data to snow crab biometrics datasheet
#2) Add a column for maturity using clutch codes/chela heights
    #Distribution based cutlines for males published in 2019 tech memo
#3) Calculate CPUE at each station
#4) Join haul data with fatty acid data 
 
# Author: Erin Fedewa
# last updated: 2/20/23

# load ----
library(tidyverse)

#############################
#Append maturity to biometrics data
bio_dat <- read.csv("./data/2019_2022 biometrics data.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
bio_dat %>%
  mutate(CW = as.numeric(CW),
      maturity = case_when((Sex == 2 & CH_CC > 0) ~ 1,
                              (Sex == 2 & CH_CC == 0) ~ 0,
                              (Sex == 1 & log(CH_CC) < -2.20640 + 1.13523 * log(CW))| (Sex == 1 & CW < 50) ~ 0,
                              (Sex == 1 & log(CH_CC) > -2.20640 + 1.13523 * log(CW)) ~ 1)) -> cond_mat

#############################
#Append EBS & NBS haul data 
ebs_haul <- read.csv("./data/haul_opilio.csv")
nbs_haul <- read.csv("./data/haul_opilio_nbs.csv")

#Combine ebs and nbs haul files 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  filter(cruise %in% c(201901, 201902, 202101, 202102, 202201, 202202),
         haul_type==3) %>%
  select(vessel, cruise, haul, mid_latitude, mid_longitude, gis_station, 
         bottom_depth, gear_temperature) %>%
  distinct() -> snow_haul

#Join haul and biometric datasets 
cond_mat %>% 
  rename_with(tolower) %>%
  left_join(snow_haul, by = c("cruise", "vessel", "haul")) -> mat_haul

#Add in CPUE data for each station 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  filter(cruise %in% c(201901, 201902, 202101, 202102, 202201, 202202),
         haul_type==3) %>%
  group_by(cruise, gis_station, area_swept) %>% 
  summarise(cpue = sum(sampling_factor, na.rm = T) / mean(area_swept)) %>%
  right_join(mat_haul, by = c("cruise", "gis_station")) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}")))  ->  snow_cpue

#Add in BSIERP sampling regions associated with each station
  #read in lookup table
regions <- read.csv("./data/regions_lookup.csv")
  #join
snow_cpue %>%
  left_join(regions, by="gis_station") -> snow_cpue_area

##################################################
#Joining Lipid Lab 2/9/23 FA data (2019 and 2021 only)

lipid <- read.csv("./data/2019_2022 FA data.csv")
colnames(lipid)<-gsub("X","",colnames(lipid))

#Data wrangling
lipid %>%
  pivot_longer(!vial_id, names_to= "id", values_to = "data") %>% 
  pivot_wider(names_from ="vial_id", values_from="data") %>%
  mutate(vial_id = gsub(".","-",id, fixed = TRUE)) %>%
  select(-id, -Order_processed, -Instd_Vial) %>%
  mutate(year = as.numeric(year)) -> lipid.dat

#Create % Weight FA Master by joining to haul data 
lipid.dat %>%
  select(contains(c("_percWT","vial_id","year"))) %>%
  full_join(snow_cpue_area, by=c("vial_id", "year")) %>%
  write_csv(file="./data/percWT_FA_master.csv")

#Create FA per WWT Master by joining to haul data
lipid.dat %>%
  select(contains(c("_perWWT","vial_id","year"))) %>%
  full_join(snow_cpue_area, by=c("vial_id","year")) %>%
  write_csv(file="./data/perWWT_FA_master.csv")

#Create Total FA Master by joining to haul data
lipid.dat %>%
  select(-contains(c("_perWWT","_percWT"))) %>%
  full_join(snow_cpue_area, by=c("vial_id","year")) %>%
  #calculate additional WWT:DWT/FA metrics
  mutate(DWT_WWT = hepato_dwt/hepato_wwt,
         Perc_DWT = DWT_WWT*100,
         Total_FA = as.numeric(Total_FA_Conc)/DWT_WWT,
         WWT_DWT = hepato_wwt/hepato_dwt) %>%
write_csv(file="./data/total_FA_master.csv")



