#Data exploration plots & data summaries for total FA dataset

# Author: Erin Fedewa
# last updated: 2/27/23

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(gganimate)

condition_master <- read.csv("./data/total_FA_master.csv")

####################################
#SAMPLE SIZES AND SPATIAL EXTENT

#plot sampling locations by year 
#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#Sample size by year plot
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, size=n_crab), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  facet_wrap(~year)

# Sample sizes by year 
condition_master %>%
  group_by(year,lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(lme), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "Sample size")

# Sample sizes by sampling region 
condition_master %>%
  group_by(year,region, lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(region), y= n), stat='identity') +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Bering Sea Region", y = "Sample size")

#Sample sizes by maturity
condition_master %>%
  filter(maturity != "NA") %>%
  group_by(year,maturity) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") #22 mature males in dataset

#Sample sizes by sex
condition_master %>%
  filter(sex != "NA") %>%
  group_by(year,sex) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(sex), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Sex", y = "Sample size")

#Size composition sampled by region/yr
condition_master %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  group_by(year) %>%
  ggplot() +
  geom_histogram(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")

#Pull very large females that appear to be Tanner crab and re-plot
condition_master %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>% #large females likely tanners
  group_by(year) %>%
  ggplot() +
  geom_histogram(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")
#will need to factor size into models! 

#Size range sampled across years by region
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>%
  group_by(lme, sex) %>%
  summarize(avg_cw = mean(cw, na.rm=T), 
            max_cw = max(cw, na.rm=T), 
            min_cw = min(cw, na.rm=T))

#############################################
#RXN BETWEEN CONDITION METRICS 

#data wrangling 
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) -> new.dat

#% DWT vrs total FA
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA, color=factor(year), label=vial_id)) +
  geom_point() +
  geom_text(hjust=.7, vjust=-.5) +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA (mg/g DWT)")
#two major outliers from 2019 may need a closer look 

#% DWT vrs total FA concentration
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA_Conc, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA per WWT (mg FA/g WWT)")
#two major outliers from 2019 may need a closer look 

#Crab weight vrs % DWT by sex
new.dat %>%
  ggplot(aes(crab_wgt, Perc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Crab weight (g)", y = "% DWT in hepatopancreas") +
  facet_wrap(~sex, scales = "free_x")
#Lots of indvidual variability! Might not expect a strong size/weight
  # rxn with FA's because growth is incremental and we only sampled one cohort? 

#Crab size vrs % DWT by sex
new.dat %>%
  ggplot(aes(cw, Perc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Carapace width (mm)", y = "% DWT in hepatopancreas") +
  facet_wrap(~sex, scales = "free_x")

# Condition factor K vrs % DWT
new.dat %>%
mutate(K=crab_wgt/(cw^3)*100000) %>%
  ggplot(aes(K, Perc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw()  +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Fultons K Condition Factor", y = "% DWT in hepatopancreas") +
  facet_wrap(~sex, scales = "free_x")
#I believe there's a better way to do this though...follow up

#############################################
#SPATIAL/INTERANNUAL VARIATION IN CONDITION METRICS 
  #Using %DWT for now as proxy until 2022 data comes in

#%DWT FA by lme and year
new.dat %>%
  ggplot(aes(factor(year), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")

#%DWT by region and year
new.dat %>%
  ggplot(aes(factor(region), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")

#%DWT by lme, year and sex
new.dat %>%
  ggplot(aes(factor(sex), Perc_DWT)) +
  geom_boxplot() +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
#EBS large males in 2019 are the anomaly 

#%DWT x size by year and region
new.dat %>%
  ggplot(aes(cw, Perc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  labs(x= "Carapace width (mm)", y = "% DWT in hepatopancreas") +
  facet_wrap(~lme, scales = "free_x")

#Mean total FA by size bin
  new.dat %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>%
  mutate(size_bin = cut(cw, breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110))) %>%
  group_by(size_bin, year, sex) %>%
  summarise(Avg_DWT = mean(Perc_DWT, na.rm=T)) %>%
  filter(size_bin != "NA") %>%
  ggplot(aes(as.factor(size_bin), Avg_DWT)) +
  geom_col() +
  theme_bw() +
  facet_grid(sex~year) +
  labs(x= "", y = "% DWT Hepatopancreas") 

################################################
#RELATIONSHIPS WITH COVARIATES
  
#range of observed temperature data by year 
condition_master %>%
  group_by(year, region, gis_station) %>%
  summarise(temperature = mean(gear_temperature)) %>%
  ggplot(aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year)

#Depth
condition_master %>%
    group_by(year, region, gis_station) %>%
    summarise(depth = mean(bottom_depth)) %>%
    ggplot(aes(depth)) +
    geom_histogram(bins = 12, fill = "dark grey", color = "black") +
    facet_wrap(~year)

#Cpue
condition_master %>%
  group_by(year, region, gis_station) %>%
  summarise(cpue = mean(cpue)) %>%
  ggplot(aes(cpue)) +
  geom_histogram(fill = "dark grey", color = "black") +
  facet_wrap(~year)

#Plot explanatory variables as predictors of %DWT by year/station 
new.dat %>%
  group_by(year, lme, gis_station) %>%
  summarise(size = mean(cw), 
            temperature = mean(gear_temperature),
            CPUE = mean(cpue^0.25), #fourth root transform
            avg_percDWT = mean(Perc_DWT)) -> plot

#Mean size-at-station vrs %DWT 
ggplot(plot, aes(size, avg_percDWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam") 

#Temp-at-station vrs %DWT
ggplot(plot, aes(temperature, avg_percDWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam")

#CPUE-at-station vrs %DWT
ggplot(plot, aes(CPUE, avg_percDWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam")

#Maybe we should split out by sex?

################################################
#SPATIAL PLOTS

#Avg total FA by station/year
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(avg_FA=mean(Total_FA)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, size=avg_FA), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  facet_wrap(~year)

##########################################
##TO NOTE:
#In 2019, there are 23 crab samples that have no associated FA data. Labels were lost 
  #during shipment of these samples so FA were not run
#Note that there are a few BCS visually positive crab noted. These crab were not extreme outliers 
  #so were not removed from further analyses 
#Despite protocols specifying immature males, there are 22 mature males in the dataset. These should be
  #removed from further analyses to control for ontogeny.
#The following vial ID's should be removed from further analyses. Crab appear to either be mature females
  #or female tanner crab based on size: 2019-65, 2019-67, 2019-68, 2019-71, 2019-66


