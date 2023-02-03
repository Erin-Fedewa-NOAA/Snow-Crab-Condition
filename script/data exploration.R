#Data exploration plots & data summaries 

# Author: Erin Fedewa
# last updated: 1/13/23

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(gganimate)

condition_master <- read.csv("./data/condition_haul_master.csv")

####################################
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

#2022 only
condition_master %>% 
  filter(year == 2022) %>%
  group_by(mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, size=n_crab), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  labs(x="", y="", size="no. of crab")

# Table of sample sizes by sampling region 
condition_master %>%
  group_by(year,region, lme) %>%
  count() %>%
  print(n=50) -> plot

#Plot of sample sizes
plot %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(region), y= n), stat='identity') +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Bering Sea Region", y = "Sample size")

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

