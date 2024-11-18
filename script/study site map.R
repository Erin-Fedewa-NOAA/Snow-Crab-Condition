#Script to create Fig 1 in manuscript
  #1) ERA 5 ice cover processing to add ice extent boundary to map
  #2) Map of sampling locations
  #3) EBS/NBS population trajectories
  #4) Temperature and density bar plots

#Authors: E. Fedewa, E. Ryznar 

#To downlaod and process ice cover data from ERA 5 monthly averaged data
  # 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview
  # 2) Click on "Download data" tab
  # 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”
  # 4) For ice cover, click on the “Other” drop down, and then check the “Sea-ice cover” box. 
  # 5) Select years you’d like the data to cover. Sometimes including the most recent year with earlier years 
  #    in your selection results in an error, which can be resolved by downloading recent year data separately.
  # 6) Select months you’d like the data to cover (Jan-Apr in this case) 
  # 7) Select geographical area for the data. Subsetting to a specific region speeds up the data querying/download.
  #    For the data below, the region has been North 66°, West -175°, South 52°, and East -155°. 
  # 8) Select NetCDF(experimental) as the data format, and submit form to query/download data. 

# **NOTE** 
  # For requests which return a mixture of ERA5 and ERA5T data  
  # (such as for data from the 1st of the month), instantaneous variables (e.g temperature) come from ERA5T (which has 'experiment version'  of 5) 
  # while accumulated variables (fluxes, precipitation) come from both datasets with the following structure:
  # 00-06 UTC on 1 day of the month from ERA5 (expver 1)
  # 07-23 UTC on 1 day of the month (and the following dates up to 5 day from present) from ERA5T (expver 5)
  
  # When these data are converted to netCDF a new dimension is created called expver containing 1 and 5. 
    #Moreover, a single time coordinate is used which covers the entire requested period.

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(tidync)
library(sf)
library(terra)
library(akgfmaps)
library(ggridges)
library(patchwork)
library(hrbrthemes)
library(ggtext)
library(ggpubr)


#install.packages("remotes")
#remotes::install_github("afsc-gap-products/akgfmaps")

### LOAD PLOTTING COLORS -------------------------------------------------------
my_colors <- c("#D55E00","#9ECAE1", "#4292C6", "#084594")
options(scipen = 999)

### PROCESS ICE DATA ----------------------------------------------------------
  
data <- tidync("./Data/ERA5_ice_1975_2024.nc") %>% 
  hyper_tibble() %>% 
  separate(date, into = c("Year", "Month", "Time"), sep = c(4,6)) %>%
  select(-Time) %>%
  mutate(Year = as.numeric(Year),
         Month = as.numeric(Month))

## SET COORDINATE REFERENCE SYSTEMS (CRS) --------------------------------------

in.crs <- "+proj=longlat +datum=NAD83" #CRS is in lat/lon
map.crs <- "EPSG:3338" # final crs for mapping/plotting: Alaska Albers

## LOAD SHELLFISH ASSESSMENT PROGRAM GEODATABASE -------------------------------

survey_gdb <- "./Data/SAP_layers" 
survey_strata <- terra::vect(survey_gdb, layer = "EBS.NBS_surveyarea")

## LOAD CONDITION AND SURVEY DATA ---------------------------------------------

condition_master <- read.csv("./data/total_FA_master.csv")
ebs_haul <- read.csv("./data/crabhaul_opilio.csv")
nbs_haul <- read.csv("./data/crabhaul_opilio_nbs.csv")
ebs_strata <- read.csv("./data/crabstrata_opilio.csv")
nbs_strata <- read.csv("./data/crabstrata_opilio_nbs.csv")

## LOAD ALASKA REGION LAYERS (FROM AKGFMAPS) -----------------------------------

ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area
ebs_survey_areas$survey_name <- c("Eastern Bering Sea", "Northern Bering Sea")

# Survey areas plot
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = ebs_survey_areas, mapping = aes(fill = survey_name)) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_fill_viridis_d(name = "Survey") +
  theme_bw()

# EBS/NBS shelf Survey grid plot
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = ebs_layers$survey.grid, fill = NA) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  theme_bw()

### PANEL A MAP PLOT -----------------------------------------------------------

#Filter for sea ice extent threshold 
  #We'll use 15% as our threshold for sea ice extent from estimates of sea ice 
  #concentration (e.g. concentration >= 0.15 is ice-covered)- per NSIDC
data %>%
  filter(Year > 2018, Year < 2024, Month == 03, Year != 2020) %>%
  filter(siconc >= 0.15) %>%
  #there are some ice-covered areas in the AI that really drive the southern 
    #bound of the hull so we'll eliminate these few points manually 
  filter(latitude > 55.5) %>%
  rename(year = Year) %>%
#transform sea ice data into spatial data frame
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
  sf::st_transform(sf::st_crs(map.crs)) %>%
  vect(.) %>%
  mask(., survey_strata) %>% #this bounds ice extent data to survey region
  sf::st_as_sf() -> ice_extent  

#Transform crab data into spatial data frame
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  # Convert lat/long to an sf object
  st_as_sf(coords = c("mid_longitude", "mid_latitude"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat

#Use the spatial data frame to generate a convex hull around the data extent
ice_hull  <-  st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(ice_extent))), 
                          dist = 15000), dTolerance = 5000)

#Add hull, ice extent layer and crab data to map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
#add hull for sea ice extent
  #geom_sf(data = ice_hull,
          #fill = NA,
          #color = alpha("red", 0.85),
          #linewidth = 1) +
#add ice extent
  geom_sf(data=ice_extent , aes(), color = "#9ECAE1", alpha = 0.25 ) +
#add crab sampling
  geom_sf(data=crab_dat, aes(size = n_crab), color = "grey30") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  labs(x="", y="", size = expression(paste("Snow crab \n samples"))) +
  theme(legend.position = "bottom") +
  facet_wrap(~year) -> map

### PANEL B ------------------------------------------------------------------
#EBS and NBS abundance timeseries 

#calculate EBS abundance timeseries 
ebs_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(HAUL_TYPE == 3,
         YEAR >= 1988) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(ebs_haul %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3,
                      YEAR >= 1988) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue_cnt = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  
  #join to stratum
  left_join(ebs_strata %>%
              select(STATION_ID, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR >= 1988) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            N_CPUE = n(),
            VAR_CPUE = (var(cpue_cnt)*(TOTAL_AREA^2))/N_CPUE,
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  distinct() %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6,
            SD_CPUE = sqrt(sum(VAR_CPUE)),
            ABUNDANCE_CI = (1.96*(SD_CPUE))/1e6) %>%
  mutate(lme = rep("Eastern Bering Sea")) -> ebs_abundance

#----------------------
#calculate NBS abundance timeseries 
nbs_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR >= 1988) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(nbs_haul %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(YEAR >= 1988) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue_cnt = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  
  #join to stratum
  left_join(nbs_strata %>%
              select(GIS_STATION, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR >= 1988) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            N_CPUE = n(),
            VAR_CPUE = (var(cpue_cnt)*(TOTAL_AREA^2))/N_CPUE,
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  distinct() %>%
  filter(STRATUM != "NA") %>% #NAs are issues documented in 11/17 SAP email! 
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6,
            SD_CPUE = sqrt(sum(VAR_CPUE)),
            ABUNDANCE_CI = (1.96*(SD_CPUE))/1e6) %>%
  mutate(lme = rep("Northern Bering Sea")) -> nbs_abundance

#Combine EBS and NBS datasets
ebs_abundance %>%
  bind_rows(nbs_abundance) -> plot_abun

#Plot with lines
plot_abun %>%
  filter (YEAR >= 1995) %>%
ggplot() +
  #geom_point(data = subset(plot_abun, lme == "Eastern Bering Sea"), 
             #aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 2) +
  geom_line(data = subset(plot_abun, lme == "Eastern Bering Sea"), 
            aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 1, alpha = 7) +
  geom_ribbon(data = subset(plot_abun, lme == "Eastern Bering Sea"),
                aes(ymin = ABUNDANCE_MIL - ABUNDANCE_CI, 
                    ymax = ABUNDANCE_MIL + ABUNDANCE_CI, x=YEAR),
              fill = "#9ECAE1", alpha = 0.2) +
  geom_point(data = subset(plot_abun, lme == "Northern Bering Sea"), 
             aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 1.5) +
  geom_errorbar(data = subset(plot_abun, lme == "Northern Bering Sea"),
              aes(ymin = ABUNDANCE_MIL - ABUNDANCE_CI, 
                  ymax = ABUNDANCE_MIL + ABUNDANCE_CI, x=YEAR),
              color = "#74C476", alpha = .7) +
  scale_color_manual(values = c("#4292C6", "#41AB5D")) +
  theme_bw() +
  labs(y="Snow Crab \nAbundance (millions)", x="") +
  theme(legend.position="top") +
  theme(axis.title.y = element_text(size=10, color = "grey30")) +
  theme(legend.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=2.5, shape=19))) -> abun_plot
  

#Stacked bar plot
plot_abun %>%
  filter (YEAR >= 1992) %>%
  ggplot(aes(fill=lme, y=ABUNDANCE_MIL, x=YEAR)) + 
  geom_col(position = position_dodge(preserve = 'single')) +
  theme_bw() +
  scale_fill_manual(values = c("#4292C6", "#084594")) +
  labs(y="Snow Crab Abundance (millions)", x="")

### PANELS C AND D -----------------------------------------------------------
#Barplots of annual mean cpue and temp for EBS and NBS

#cpue data wrangling
condition_master %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  mutate(lme = recode(lme, EBS = "Eastern Bering Sea", NBS = "Northern Bering Sea")) %>%
  distinct(year, lme, gis_station, cpue) %>%
  group_by(year, lme) %>%
  summarize(mean_cpue = mean(cpue, na.rm = T)/1000, #converting to thous crab/nmi2
            sd_cpue = sd(cpue, na.rm = T)/1000,
            n_cpue = n()) %>%
  mutate(se_cpue = sd_cpue / sqrt(n_cpue),
         lower.ci = mean_cpue - qt(1 - (0.05 / 2), n_cpue - 1) * se_cpue,
         upper.ci = mean_cpue + qt(1 - (0.05 / 2), n_cpue - 1) * se_cpue,
         year = as.factor(year)) -> cpue.dat
#and plot
ggplot(cpue.dat, aes(year, mean_cpue)) +
  geom_col(aes(fill = ordered(year)), size=3) +
  geom_errorbar(aes(year, ymin=mean_cpue - se_cpue, ymax=mean_cpue + se_cpue), 
                width=0.3, size=0.5, color = "grey40") +
  labs(y = "Mean Snow Crab<br>Density (num/nmi<sup> 2</sup>)", x= "") +
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme, scales = "free_y") +
  geom_vline(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 1, y=155, label = "Mid-collapse"),
            size = 2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 3, y=155, label = "Post-collapse"),
            size = 2, color = "#0072B2") +
  theme_ipsum(axis_title_size = 12, axis_text_size =10) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(colour="grey30", size = 11)) +
  theme(strip.text = element_text(hjust = .5)) +
  theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5)) +
  theme(panel.grid.major.x = element_blank())  -> mean_cpue_plot

#temperature data wrangling
condition_master %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  mutate(lme = recode(lme, EBS = "Eastern Bering Sea", NBS = "Northern Bering Sea")) %>%
  distinct(year, lme, gis_station, gear_temperature) %>%
  group_by(year, lme) %>%
  summarize(mean_temp = mean(gear_temperature, na.rm = T),
            sd_temp = sd(gear_temperature, na.rm = T),
            n_temp = n()) %>%
  mutate(se_temp = sd_temp / sqrt(n_temp),
         lower.ci = mean_temp - qt(1 - (0.05 / 2), n_temp - 1) * se_temp,
         upper.ci = mean_temp + qt(1 - (0.05 / 2), n_temp - 1) * se_temp,
         year = as.factor(year)) -> temp.dat
#plot  
ggplot(temp.dat, aes(year, mean_temp)) +
  geom_col(aes(fill = ordered(year)), size=3) +
  geom_errorbar(aes(year, ymin = ifelse(mean_temp - se_temp < 0, 0, mean_temp - se_temp), 
                    ymax=mean_temp + se_temp),
                width=0.3, size=0.5, color = "grey40") +
  labs(y = expression("Mean Bottom \nTemperature " ( degree~C)), x = "") +
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme) +
  geom_vline(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 1, y=3.3, label = "Mid-collapse"),
            size = 2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 3, y=3.3, label = "Post-collapse"),
            size = 2, color = "#0072B2") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 12, axis_text_size =10) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(colour="grey30", hjust = 0.5, size = 11)) +
  theme(strip.text = element_text(hjust = .5)) +
  theme(panel.grid.major.x = element_blank())  -> mean_temp_plot

### COMBINE PANELS AND SAVE FIGURE --------------------------------------------

#Figure 1 for ms: combined abun, map, temp and cpue plot

abun_plot + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position  = c(.05, .95))
ggsave("./figures/Fig1a.png", height=2.5 , width=4, units="in")

map + plot_annotation(tag_levels = list('b')) & 
  theme(plot.tag.position  = c(.1, 1)) 
ggsave("./figures/Fig1b.png", height=6 , width=6, units="in")

mean_cpue_plot / plot_spacer() / mean_temp_plot  + plot_layout(heights = c(6, -4 , 6)) +
  plot_annotation(tag_levels = list(c('c', 'd'))) 
ggsave("./figures/Fig1cd.png", height=6 , width=6, units="in")


-----------------------------------------------------------------------------
#Bonus Figs: And now pdfs by year

#temperature
condition_master  %>%
  filter(lme %in% c("EBS", "NBS")) %>%
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
  scale_fill_manual(values = my_colors) +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "Bottom Temperature (C)", y = "Count") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) 

#cpue 
condition_master  %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  group_by(year, lme, gis_station) %>%
  summarise(cpue = mean(cpue)) %>%
  ggplot(aes(cpue,factor(year))) +
  geom_density_ridges(aes(fill=factor(year)), scale=2,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x),
                      rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2,
                      alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "Snow Crab Density", y = "Count") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) 


condition_master %>%
     group_by(year, lme) %>%
  summarise(mean_cpue = mean(cpue), 
            mean_temp = mean(gear_temperature)) %>%
  ggplot( ) +
    geom_point(aes(mean_cpue, mean_temp, color = lme)) +
  geom_line(aes(mean_cpue, mean_temp, color = lme)) +
  geom_text(label = year)

condition_master %>%
  mutate(fourth.root.cpue = as.numeric(cpue^0.25)) %>%
  ggplot( ) +
  geom_point(aes(fourth.root.cpue, gear_temperature, color = lme)) +
  geom_line(aes(fourth.root.cpue, gear_temperature, color = lme)) +
  geom_text(label = year)
  



















   
  