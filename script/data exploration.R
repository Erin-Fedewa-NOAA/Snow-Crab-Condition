#Data exploration plots & data summaries for total FA dataset

# Author: Erin Fedewa
# last updated: 1/5/24

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(gganimate)
library(viridis)
library(ggridges)
library(RColorBrewer)

condition_master <- read.csv("./data/total_FA_master.csv")

####################################
#SAMPLE SIZES AND SPATIAL EXTENT

#plot sampling locations by year 
#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#Sample size by year map
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, size=n_crab), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  facet_wrap(~year)
  ggsave("./figures/data exploration/n_year.png", dpi=300)

# Sample sizes by year plot
condition_master %>%
  group_by(year,lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(lme), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "Sample size")
ggsave("./figures/data exploration/n_year2.png", dpi=300)

# Sample sizes by BSIERP region 
condition_master %>%
  group_by(year,bsierp_region, lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(bsierp_region), y= n), stat='identity') +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "BSIERP Region", y = "Sample size")
ggsave("./figures/data exploration/n_bsierpregion.png", dpi=300)

#Sample sizes by maturity
condition_master %>%
  filter(maturity != "NA") %>%
  group_by(lme,year,maturity) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") #45 mature males in dataset

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
ggsave("./figures/data exploration/n_sex.png", dpi=300)

#Size composition sampled by region/yr
condition_master %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  group_by(year) %>%
  ggplot() +
  geom_density(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
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
  geom_density(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")
ggsave("./figures/data exploration/size_comp.png", dpi=300)
#will need to factor size into models! 

#Size range sampled across years by region
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>%
  group_by(lme, year, sex) %>%
  summarize(avg_cw = mean(cw, na.rm=T), 
            max_cw = max(cw, na.rm=T), 
            min_cw = min(cw, na.rm=T))

#############################################
#RXN BETWEEN CONDITION METRICS 

#data wrangling 
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),
         maturity != 1) -> new.dat

#% DWT vrs total FA
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA, color=factor(year), label=vial_id)) +
  geom_point() +
  geom_text(hjust=.7, vjust=-.5) +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA (mg/g DWT)")

#% DWT vrs total FA concentration
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA_Conc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA per DWT (mg FA/g WWT)") +
  theme(legend.title=element_blank())
ggsave("./figures/data exploration/DWTvFA.png", dpi=300)

#Crab weight at size vrs % DWT by sex
new.dat %>%
  mutate(lw = crab_wgt/cw) %>%
  ggplot(aes(crab_wgt, lw, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Crab weight at size", y = "% DWT in hepatopancreas") +
  facet_wrap(~sex, scales = "free_x")

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
  
#%DWT by lme and year
new.dat %>%
  ggplot(aes(factor(year), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
ggsave("./figures/data exploration/DWT_year.png", dpi=300)

#total FA concentration by lme and year
new.dat %>%
  ggplot(aes(factor(year), Total_FA)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "Total FA")
ggsave("./figures/data exploration/TotalFA_year.png", dpi=300)

#Bar plot
lme_names <- as_labeller(c("EBS" = "Eastern Bering Sea",
    "NBS" = "Northern Bering Sea"))

new.dat %>%
  ggplot(aes(factor(year), Perc_DWT, fill=lme)) +
  geom_point(aes(color=lme),stat="summary", size=4) +
  geom_errorbar(stat="summary", 
                 colour="darkgray", size=.6) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 11)) +
  facet_wrap(~lme, labeller = lme_names)

#density plot, EBS only 
new.dat %>%
  filter(lme == "EBS") %>%
  ggplot(aes(Perc_DWT, factor(year))) +
  geom_density_ridges(aes(fill=factor(year)), scale=2,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x)) +
  theme_ridges(center=TRUE) +
  scale_fill_brewer() +
  labs(y= "", x = "Snow Crab Energetic Condition (%DWT)") +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 11))
  

#%DWT by region and year
new.dat %>%
  ggplot(aes(factor(bsierp_region), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
 # labs(x= "", y = "Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/DWT_bsierpregion.png", dpi=300)

#%DWT by lme, year and sex
new.dat %>%
  ggplot(aes(factor(sex), Perc_DWT)) +
  geom_boxplot() +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
ggsave("./figures/data exploration/DWT_year_lme.png", dpi=300)
#EBS large males in 2019 are the anomaly - note that two outliers for %DWT figure
  #above are also large outliers here too (2019-207 and 2019-212)

#%DWT x size by year and region
new.dat %>%
  ggplot(aes(cw, Perc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  labs(x= "Carapace width (mm)", y = "% DWT in hepatopancreas") +
  facet_wrap(~lme, scales = "free_x")

#Mean %DWT by size bin
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
  labs(x= "Carapace width size bin (mm)", y = "% DWT in hepatopancreas")
  ggsave("./figures/data exploration/DWT_year_size.png", dpi=300)

################################################
#RELATIONSHIPS WITH COVARIATES
  
#range of observed temperature data by year 
condition_master %>%
  group_by(year, bsierp_region, gis_station) %>%
  summarise(temperature = mean(gear_temperature)) %>%
  ggplot(aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Bottom Temperature (C)", y = "Count")
  ggsave("./figures/data exploration/temp.png", dpi=300)

#Depth
condition_master %>%
    group_by(year, bsierp_region, gis_station) %>%
    summarise(depth = mean(bottom_depth)) %>%
    ggplot(aes(depth)) +
    geom_histogram(bins = 12, fill = "dark grey", color = "black") +
    facet_wrap(~year) +
    theme_bw() +
    labs(x= "Depth (m)", y = "Count")

#Crab Cpue
condition_master %>%
  group_by(year, bsierp_region, gis_station) %>%
  summarise(cpue = mean(cpue)) %>%
  ggplot(aes(cpue)) +
  geom_histogram(fill = "dark grey", color = "black") +
  facet_wrap(~year)+
  theme_bw() +
  labs(x= "Snow Crab Density", y = "Count")

#Benthic invert Cpue
condition_master %>%
  group_by(year, bsierp_region, gis_station) %>%
  summarise(invert_cpue = mean(total_benthic_cpue)) %>%
  ggplot(aes(invert_cpue)) +
  geom_histogram(fill = "dark grey", color = "black") +
  facet_wrap(~year)+
  theme_bw() +
  labs(x= "Benthic Invert Density", y = "Count")

#Plot explanatory variables as predictors of % DWT by year/station
new.dat %>%
  group_by(year, lme, gis_station) %>%
  summarise(size = mean(cw), 
            temperature = mean(gear_temperature),
            CPUE = mean(cpue^0.25), #fourth root transform
            invert = mean(total_benthic_cpue^0.25),
            avg_Perc_DWT = mean(Perc_DWT)) -> plot

#Mean size-at-station vrs %DWT
ggplot(plot, aes(size, avg_Perc_DWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Mean carapace width at station (mm)", y="% DWT in hepatopancreas")
ggsave("./figures/data exploration/stationxsizexDWT.png", dpi=300)

#Temp-at-station vrs %DWT
ggplot(plot, aes(temperature, avg_Perc_DWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Mean temperature at station (C)", y="% DWT in hepatopancreas")
ggsave("./figures/data exploration/stationxtempxDWT.png", dpi=300)

#Temp by year vrs %DWT
new.dat %>%
  group_by(year, lme) %>%
  summarise(temperature_annual = mean(gear_temperature),
            avg_Perc_DWT_annual = mean(Perc_DWT)) %>%
ggplot(aes(temperature_annual, avg_Perc_DWT_annual)) +
  geom_point() + 
  theme_bw() +
  labs(x="Mean temperature(C)", y="% DWT in hepatopancreas")
ggsave("./figures/data exploration/yearxtempxDWT.png", dpi=300)

#CPUE-at-station vrs %DWT
ggplot(plot, aes(CPUE, avg_Perc_DWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Snow crab density at station", y="% DWT in hepatopancreas")
ggsave("./figures/data exploration/stationxcpuexDWT.png", dpi=300)

#Benthic invert CPUE-at-station vrs %DWT
ggplot(plot, aes(invert, avg_Perc_DWT)) +
  geom_point() + 
  facet_grid(lme~year) +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Benthic invert density at station", y="% DWT in hepatopancreas")
ggsave("./figures/data exploration/stationxinvertcpuexDWT.png", dpi=300)

#Maybe we should split out by sex?

################################################
#SPATIAL PLOTS

#Avg %DWT by station/year
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(avg_dwt=mean(Perc_DWT)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=avg_dwt))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  scale_color_viridis() +
  theme_bw() +
  facet_wrap(~year) +
  labs(y= "Latitude", color = "Energetic Condition\n(% DWT)")
ggsave("./figures/data exploration/avgDWT_map.png", dpi=300)

#Follow up on this- move into survey grid shapefiles- potentially use an 
  #IDW approach? 

##########################################
##TO NOTE:
#In 2019, there are 23 crab samples that have no associated FA data. Labels were lost 
  #during shipment of these samples so FA were not run
#Note that there are a few BCS visually positive crab noted. These crab were not extreme outliers 
  #so were not removed from further analyses 
#Despite protocols specifying immature males, there are mature males in the dataset. These should be
  #removed from further analyses to control for ontogeny.
#The following vial ID's should be removed from further analyses. Crab appear to either be mature females
  #or female tanner crab based on size: 2019-65, 2019-67, 2019-68, 2019-71, 2019-66
#200 samples from 2023 were prioritized for fatty acids, and the rest were only measured for DWT/WWT


