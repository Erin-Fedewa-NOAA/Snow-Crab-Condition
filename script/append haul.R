#GOALS: 
# 1) Append haul data to snow crab biometrics datasheet
#2) Add a column for maturity using clutch codes/chela heights
    #Distribution based cutlines for males published in 2019 tech memo
#3) Calculate snow crab CPUE at each station
#4) Append station-level benthic invert densities
#5) Join haul data with fatty acid data 

#Follow ups: Change chela NBS vrs EBS script section to recognize 01 vrs 02 string 
 
# Author: Erin Fedewa
# last updated: 10/10/23

# load ----
library(tidyverse)

#############################
#Append maturity to biometrics data
bio_dat <- read.csv("./data/2019_2023 biometrics data.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
bio_dat %>%
  mutate(CW = as.numeric(CW),
      maturity = case_when((Sex == 2 & CH_CC > 0) ~ 1, #mature female (EBS & NBS)
                              (Sex == 2 & CH_CC == 0) ~ 0, #immature female (EBS & NBS)
                              #EBS male cutlines 
                              (Sex == 1 & Cruise %in% c(201901,202101,202201,202301) & 
                                 log(CH_CC) < -2.20640 + 1.13523 * log(CW))| (Sex == 1 & CW < 50 &
                                  Cruise %in% c(201901,202101,202201,202301)) ~ 0, #immature male EBS
                              (Sex == 1 & Cruise %in% c(201901,202101,202201,202301) &
                                 log(CH_CC) >= -2.20640 + 1.13523 * log(CW)) ~ 1, #mature male EBS
                              #NBS male cutlines
                              (Sex == 1 & Cruise %in% c(201902,202102,202202,202302) & 
                                 log(CH_CC) < -1.916947 + 1.070620 * log(CW))| (Sex == 1 & CW < 40 &
                                  Cruise %in% c(201902,202102,202202,202302)) ~ 0, #immature male NBS
                           (Sex == 1 & Cruise %in% c(201902,202102,202202,202302) &
                              log(CH_CC) >= -1.916947 + 1.070620 * log(CW)) ~ 1)) -> cond_mat #mature male NBS
                                                            

#############
#Append EBS & NBS haul data 
ebs_haul <- read.csv("./data/haul_opilio.csv")
  unique(ebs_haul$CRUISE) 
nbs_haul <- read.csv("./data/haul_opilio_nbs.csv")
  unique(nbs_haul$CRUISE)

#Combine ebs and nbs haul files 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  filter(cruise %in% c(201901, 201902, 202101, 202102, 202201, 202202, 202301, 202302),
         haul_type==3) %>%
  select(vessel, cruise, haul, mid_latitude, mid_longitude, gis_station, 
         bottom_depth, gear_temperature) %>%
  distinct() -> snow_haul

#Join haul and biometric datasets 
cond_mat %>% 
  rename_with(tolower) %>%
  left_join(snow_haul, by = c("cruise", "vessel", "haul")) -> mat_haul

#################################
#Add in station-level covariates of interest for bayesian modeling

#Station-level snow crab CPUE (model covariate- competition/density dependence)
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  filter(cruise %in% c(201901, 201902, 202101, 202102, 202201, 202202, 202301, 202302),
         haul_type==3) %>%
  group_by(cruise, gis_station, area_swept) %>% 
  summarise(cpue = sum(sampling_factor, na.rm = T) / mean(area_swept)) %>%
  right_join(mat_haul, by = c("cruise", "gis_station")) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}")))  ->  snow_cpue

#Add in benthic invert CPUE data for each station (model covariate- prey quantity)
  #Pulling from a different script, which generates the CPUE estimates

source("./script/benthic_invert.R")

benthic_cpue %>%
  select(cruise,year,gis_station,total_benthic_cpue) %>% 
  right_join(snow_cpue) -> snow_invert_cpue

#Add in BSIERP and sampling regions associated with each station
  #read in lookup table
regions <- read.csv("./data/regions_lookup.csv")
  #join
snow_invert_cpue %>%
  left_join(regions, by="gis_station") -> snow_cpue_area

##################################################
#Joining Lipid Lab 2/9/23 FA data (no 2023 data yet)

lipid <- read.csv("./data/2019_2022 FA data.csv", na.strings="")
colnames(lipid)<-gsub("X","",colnames(lipid))

#Data wrangling
lipid %>%
  pivot_longer(!vial_id, names_to= "id", values_to = "data", values_transform = as.numeric) %>% 
  pivot_wider(names_from ="vial_id", values_from="data") %>%
  mutate(vial_id = gsub(".","-",id, fixed = TRUE)) %>%
  select(-id, -Lost_sample) %>%
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

#Create FA per DWT Master by joining to haul data
lipid.dat %>%
  select(contains(c("_perDWT","vial_id","year"))) %>%
  full_join(snow_cpue_area, by=c("vial_id","year")) %>%
  write_csv(file="./data/perDWT_FA_master.csv")

#Create Total FA Master by joining to haul data
lipid.dat %>%
  select(-contains(c("_perWWT","_percWT", "_perDWT"))) %>%
  full_join(snow_cpue_area, by=c("vial_id","year")) %>%
  #calculate additional WWT:DWT/FA metrics
  mutate(DWT_WWT = hepato_dwt/hepato_wwt,
         Perc_DWT = DWT_WWT*100,
         Total_FA = as.numeric(Total_FA_Conc_WWT)/DWT_WWT,
         WWT_DWT = hepato_wwt/hepato_dwt) %>%
write_csv(file="./data/total_FA_master.csv")



