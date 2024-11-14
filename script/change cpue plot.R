#Change in station-level mean snow crab CPUE plot (2018 vrs 2021)
  #Plot created for AMSS talk 

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(scales) 

#data
ebs_haul <- read.csv("./data/haul_opilio.csv")

ebs_haul %>%
  rename_with(tolower) %>%
  filter(cruise %in% c(201801, 202101),
         haul_type==3) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  group_by(year, gis_station, area_swept) %>% 
  summarise(cpue = sum(sampling_factor, na.rm = T) / mean(area_swept)) %>%
  ungroup() %>%
  dplyr::select(-area_swept) %>%
  mutate(cpue = cpue + 1) %>%
  mutate(log_cpue = log(cpue))-> dat

#calculate % change in CPUE
dat %>%
  pivot_wider(names_from = year, values_from = cpue) %>%
  mutate(perc_change_cpue = (((`2021` - `2018`)/`2018`)*100)) %>%
  mutate_if(is.numeric, ~ replace_na(., 0) %>% 
              replace(., is.infinite(.), 100)) %>%
  dplyr::select(-c(`2018`,`2021`)) -> cpue
  
#The scale isn't showing negative values b/c data are skewed 
#Calculate breakpoints or scale (discrete colors vrs continuous), or log transform
  
#Extract lat and longs for stations
ebs_haul %>%
  rename_with(tolower) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  filter(cruise %in% c(201801, 202101),
         haul_type==3) %>%
  dplyr::select(year, gis_station, mid_latitude, mid_longitude) %>%
  distinct(year, gis_station, mid_latitude, mid_longitude) -> lat

dat %>%
  full_join(lat) -> plot

#plot change in cpue as 2018 vrs 2021 side by side 
#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#plot % change in cpue
plot %>% 
  filter(cpue > 1) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group=group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, color=as.numeric(cpue), size=as.numeric(cpue)))+
  coord_quickmap(xlim = c(-178.5, -158), ylim = c(55, 62.5)) +
  theme_bw() +
  facet_wrap(~year) +
  xlab("") + ylab("") + 
  scale_color_continuous(name = "Snow Crab Density", labels = comma,
    limits=c(1, 2850000), breaks=seq(1, 2850000, by=570000)) +
  guides(color= guide_legend(), size=guide_legend()) +
  scale_size_continuous(name = "Snow Crab Density", labels = comma,
      limits=c(1, 2850000), breaks=seq(1, 2850000, by=570000))
ggsave("./figures/data exploration/cpue.png", dpi=300)


