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
library(broom)


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

# Sample sizes by year - total FA WWT
#Note that in 2023, not all samples were prioritized for FA analysis
condition_master %>%
  filter(Total_FA_Conc_WWT > 0) %>%
  group_by(year,lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(lme), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "Sample size")

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
#SPATIAL/INTERANNUAL VARIATION IN CONDITION METRICS 

#data wrangling - use this dataset for all further analyses!!
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),#likely tanners
         maturity != 1) -> new.dat
  
#%DWT by lme and year
new.dat %>%
  ggplot(aes(factor(year), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
ggsave("./figures/data exploration/DWT_year.png", dpi=300)

#total FA concentration (DWT) by lme and year
new.dat %>%
  ggplot(aes(factor(year), Total_FA_Conc_DWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "Total FA per DWT")
ggsave("./figures/data exploration/TotalFA_year.png", dpi=300)

#total FA concentration (WWT) by lme and year
new.dat %>%
  ggplot(aes(factor(year), Total_FA_Conc_WWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "Total FA per WWT")

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
  #2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
  #used as response variable in all further analyses. 

#Bar plot
lme_names <- as_labeller(c("EBS" = "Eastern Bering Sea",
    "NBS" = "Northern Bering Sea"))

new.dat %>%
  ggplot(aes(factor(year), Total_FA_Conc_WWT, fill=lme)) +
  geom_bar(aes(color=lme),stat="summary") +
  geom_errorbar(stat="summary", 
                 colour="darkgray", size=.6) +
  theme_bw() +
  labs(x= "", y = "Total FA Concentration (mg FA/g WWT)") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 11)) +
  facet_wrap(~lme, labeller = lme_names)

#density plot 
new.dat %>%
  #filter(lme == "EBS") %>%
  ggplot(aes(Total_FA_Conc_WWT, factor(year))) +
  geom_density_ridges(aes(fill=factor(year)), scale=2,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x),
                      rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2,
                      alpha = 0.7) +
  theme_ridges(center=TRUE) +
  scale_fill_brewer() +
  labs(y= "", x = "Snow Crab Energetic Condition (Total FA per WWT)") +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 11)) +
  facet_wrap(~lme, labeller = lme_names)
  
#FA by region and year
new.dat %>%
  ggplot(aes(factor(bsierp_region), Total_FA_Conc_WWT)) +
  geom_boxplot() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/DWT_bsierpregion.png", dpi=300)

#FA by lme, year and sex
new.dat %>%
  ggplot(aes(factor(sex), Total_FA_Conc_WWT)) +
  geom_boxplot() +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "", y = "Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/FA_year_lme.png", dpi=300)
#EBS large males in 2019 are the anomaly 

#FA x size by year and region
new.dat %>%
  ggplot(aes(cw, Total_FA_Conc_WWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  labs(x= "Carapace width (mm)", y = "Total FA per WWT (mg FA/g WWT)") +
  facet_wrap(~lme, scales = "free_x")

#FA by size bin
  new.dat %>%
    mutate(size_bin = cut(cw, breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110))) %>%
    group_by(size_bin, year, sex) %>%
    summarise(Avg_FA = mean(Total_FA_Conc_WWT, na.rm=T)) %>%
    filter(size_bin != "NA") %>%
  ggplot(aes(as.factor(size_bin), Avg_FA)) +
  geom_col() +
  theme_bw() +
  facet_grid(sex~year) +
  labs(x= "Carapace width size bin (mm)", y = "Total FA per WWT (mg FA/g WWT)")
  ggsave("./figures/data exploration/FA_year_size.png", dpi=300)

################################################
#RELATIONSHIPS WITH COVARIATES
  
#range of observed temperature data by year 
  new.dat  %>%
  group_by(year, lme, gis_station) %>%
  summarise(temperature = mean(gear_temperature)) %>%
  ggplot(aes(temperature,factor(year))) +
    geom_density_ridges(aes(fill=factor(year)), scale=2,
                        quantile_lines=TRUE,
                        quantile_fun=function(x,...)mean(x),
                        rel_min_height = 0.01, jittered_points = TRUE,
                        position = position_points_jitter(width = 0.5, height = 0),
                        point_shape = "|", point_size = 2,
                        alpha = 0.7) +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "Bottom Temperature (C)", y = "Count")
  ggsave("./figures/data exploration/temp.png", dpi=300)

#Depth
new.dat %>%
    group_by(year, lme, gis_station) %>%
    summarise(depth = mean(bottom_depth)) %>%
    ggplot(aes(depth)) +
    geom_histogram(bins = 12, fill = "dark grey", color = "black") +
    facet_grid(vars(year), vars(lme)) +
    theme_bw() +
    labs(x= "Depth (m)", y = "Count")

#Crab Cpue
new.dat %>%
  group_by(year, lme, gis_station) %>%
  summarise(cpue = mean(cpue)) %>%
  ggplot(aes(cpue)) +
  geom_density(fill = "dark grey", color = "black") +
  facet_grid(vars(year), vars(lme)) +
  theme_bw() +
  labs(x= "Snow Crab Density", y = "Count")

#Benthic invert Cpue
new.dat %>%
  group_by(year, lme, gis_station) %>%
  summarise(invert_cpue = mean(total_benthic_cpue)) %>%
  ggplot(aes(invert_cpue)) +
  geom_histogram(fill = "dark grey", color = "black") +
  facet_grid(vars(year), vars(lme)) +
  theme_bw() +
  labs(x= "Benthic Invert Density", y = "Count")

#Plot explanatory variables as predictors of Total FA by year/station
new.dat %>%
  group_by(year, lme, gis_station) %>%
  summarise(size = mean(cw, na.rm=T), 
            temperature = mean(gear_temperature, na.rm=T),
            CPUE = mean(cpue^0.25, na.rm=T), #fourth root transform
            invert = mean(total_benthic_cpue^0.25, na.rm=T),
            avg_FA = mean(Total_FA_Conc_WWT, na.rm=T)) -> plot

#Mean size-at-station vrs FA
ggplot(plot, aes(size, avg_FA, color=as.factor(year))) +
  geom_point() + 
  facet_wrap(~lme) +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Mean carapace width at station (mm)", y="Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/stationxsizexFA.png", dpi=300)

#Temp-at-station vrs FA
ggplot(plot, aes(temperature, avg_FA, color=as.factor(year))) +
  facet_wrap(~lme) +
  geom_point() + 
  geom_smooth(method = "gam", knots = 3) +
  theme_bw() +
  labs(x="Mean temperature at station (C)", y="Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/stationxtempxFA.png", dpi=300)

#CPUE-at-station vrs %DWT
ggplot(plot, aes(CPUE, avg_FA, color=as.factor(year))) +
  facet_wrap(~lme) +
  geom_point() + 
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x="Snow crab density at station", y="Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/stationxcpuexFA.png", dpi=300)

#Benthic invert CPUE-at-station vrs %DWT
ggplot(plot, aes(invert, avg_FA, color=as.factor(year))) +
  facet_wrap(~lme) +
  geom_point() + 
  geom_smooth(method = "gam", knots = 3) +
  theme_bw() +
  labs(x="Benthic invert density at station", y="Total FA per WWT (mg FA/g WWT)")
ggsave("./figures/data exploration/stationxinvertcpuexFA.png", dpi=300)

################################################
#SPATIAL PLOTS - response and covariates 

#Avg Total FA by station/year
new.dat %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(avg_FA=mean(Total_FA_Conc_WWT, na.rm=T)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=avg_FA))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  scale_color_viridis() +
  theme_bw() +
  facet_wrap(~year) +
  labs(y= "Latitude", color = "Energetic Condition\n(Total FA per WWT)")
ggsave("./figures/data exploration/avgWWT_map.png", dpi=300)

#temperature by station/year
new.dat %>% 
  group_by(year, lme, gis_station, mid_latitude, mid_longitude) %>%
  summarise(temperature = mean(gear_temperature)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=temperature))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  scale_color_viridis() +
  theme_bw() +
  facet_wrap(~year) +
  labs(y= "Latitude", color = "Bottom Temperature (C)")

#bentic invert cpue by station/yr
new.dat %>% 
  group_by(year, lme, gis_station, mid_latitude, mid_longitude) %>%
  summarise(invert = mean(total_benthic_cpue^0.25, na.rm=T)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=invert))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  scale_color_viridis() +
  theme_bw() +
  facet_wrap(~year) +
  labs(y= "Latitude", color = "Invert CPUE (kg/km2)")

#snow crab cpue by station/yr
new.dat %>% 
  group_by(year, lme, gis_station, mid_latitude, mid_longitude) %>%
  summarise(CPUE = mean(cpue^0.25, na.rm=T)) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=CPUE))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  scale_color_viridis() +
  theme_bw() +
  facet_wrap(~year) +
  labs(y= "Latitude", color = "Snow Crab Density")

##########################################
##TAKE AWAYS TO NOTE:
#In 2019, there are 23 crab samples that have no associated FA data. Labels were lost 
  #during shipment of these samples so FA were not run
#Note that there are a few BCS visually positive crab noted. These crab were not extreme outliers 
  #so were not removed from further analyses 
#Despite protocols specifying immature males, there are mature males in the dataset. These should be
  #removed from further analyses to control for ontogeny.
#The following vial ID's should be removed from further analyses. Crab appear to either be mature females
  #or female tanner crab based on size: 2019-65, 2019-67, 2019-68, 2019-71, 2019-66
#200 samples from 2023 were prioritized for fatty acids, and the rest were only measured for DWT/WWT
  #125 for EBS, 75 for NBS


