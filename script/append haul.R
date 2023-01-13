#GOALS: 
# 1) Append haul data to master condition datasheet
#2) Add a column for maturity using clutch codes/chela heights
    #Distribution based cutlines for males published in 2019 tech memo
#3) Calculate CPUE at each station
 
# Author: Erin Fedewa
# last updated: 1/13/23

# load ----
library(tidyverse)

#############################
#Append maturity 
condition_master <- read.csv("./data/Snow Condition Master.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
condition_master %>%
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

#Join haul and condition datasets 
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

##################################################

#Write new master csv                               
write_csv(snow_cpue, file="./data/condition_haul_master.csv")
