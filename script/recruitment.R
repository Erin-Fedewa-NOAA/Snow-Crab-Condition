#Goals
  # Develop incides of abundance for recruits (newshell mature males and females)
  #Calculate annual means for snow crab condition
  #Assess relationship between recruitment and condition the year prior

# Author: Erin Fedewa
# last updated: 2/27/23

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(gganimate)

#condition data
condition_master <- read.csv("./data/total_FA_master.csv")

#EBS & NBS haul & strata data 
ebs_haul <- read.csv("./data/haul_opilio.csv")
ebs_strata <- read.csv("./data/strata_opilio.csv")

nbs_haul <- read.csv("./data/haul_opilio_nbs.csv")
nbs_strata <- read.csv("./data/strata_opilio_nbs.csv")

#Functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

##############################################
#Calculate Mean condition by LME/sex - add BSIERP region too? 

#data wrangling 
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),
         maturity != 1) %>%
  group_by(lme, sex, year) %>%
  summarise(avg_Total_FA = mean(Total_FA, na.rm=T),
            ssd = sd(Total_FA, na.rm = TRUE),
            count = n()) %>%
  mutate(se = ssd / sqrt(count),
         lower_ci = lower_ci(avg_Total_FA, se, count),
         upper_ci = upper_ci(avg_Total_FA, se, count)) -> condition

#calculate EBS mature female and mature male abundance the following year 
ebs_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(HAUL_TYPE == 3, 
         YEAR > 2018, 
         SEX %in% c(1,2)) %>%
  mutate(MAT_SEX = case_when((SEX == 2 & CLUTCH_SIZE > 0 & SHELL_CONDITION %in% c(0:2)) ~ "Mature Female",
                             (SEX == 1 & WIDTH_1MM >= 102 & SHELL_CONDITION %in% c(0:2)) ~ "Mature Male")) %>%
  filter(!is.na(MAT_SEX)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, MAT_SEX) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(ebs_haul %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3,
                      YEAR > 2018) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue_cnt = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  
  #join to stratum
  left_join(ebs_strata %>%
              select(STATION_ID, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR > 2018) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA, MAT_SEX) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  group_by(YEAR, MAT_SEX) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6) %>%
  filter(!is.na(MAT_SEX)) %>%
  mutate(lag_year = YEAR -1) -> ebs_abundance

#so you'd have one plot, with both sexes, both regions 


#NBS
  #what size class is next up from size calculating- and will these crab move south?







#############################################