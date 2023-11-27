#Change in CPUE plot (2018 vrs 2021)

# load ----
library(tidyverse)
library(sf)
library(ggmap)

#data
ebs_haul <- read.csv("./data/haul_opilio.csv")

ebs_haul %>%
  rename_with(tolower) %>%
  filter(cruise %in% c(201801, 202101),
         haul_type==3) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  group_by(year, gis_station, area_swept) %>% 
  summarise(cpue = sum(sampling_factor, na.rm = T) / mean(area_swept)) %>%
  ungroup() -> dat

dat %>%
  dplyr::select(-area_swept) %>%
  pivot_wider(names_from = year, values_from = cpue) %>%
  mutate(perc_change_cpue = (((`2021` - `2018`)/`2018`)*100)) %>%
  mutate_if(is.numeric, ~ replace_na(., 0) %>% 
              replace(., is.infinite(.), 100)) %>%
  dplyr::select(-c(`2018`,`2021`)) -> cpue

ebs_haul %>%
  rename_with(tolower) %>%
  filter(cruise == 202101,
         haul_type==3) %>%
  dplyr::select(gis_station, mid_latitude, mid_longitude) %>%
  distinct(gis_station, mid_latitude, mid_longitude) -> lat

cpue %>%
  full_join(lat) -> plot

#plot % change in cpue
#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#plot
plot %>% 
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group=group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, fill=as.numeric(perc_change_cpue)))+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 64)) +
  theme_bw() 

#The scale isn't showing negative values b/c data are skewed 
#Calucalte breakpoints or scale (discrete colors vrs continuous)

plot %>% print(n=500)
