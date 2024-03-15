# Accessing EBS survey shapefiles using the akgfmaps package
# Created by Sean Rohan <sean.rohan@noaa.gov>
# March 13, 2024

install.packages("remotes")
remotes::install_github("afsc-gap-products/akgfmaps")

library(akgfmaps)

condition_master <- read.csv("./data/total_FA_master.csv")

# Plotting layers
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

slope_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "EPSG:3338")

survey_areas <- dplyr::bind_rows(ebs_layers$survey.area,
                                 slope_layers$survey.area)

survey_areas$survey_name <- c("EBS Shelf", "NBS", "EBS Slope")

# Survey areas
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = survey_areas, mapping = aes(fill = survey_name)) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_fill_viridis_d(name = "Survey") +
  theme_bw()


# EBS/NBS shelf Survey grid
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = ebs_layers$survey.grid, fill = NA) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  theme_bw()

#Now how to overlay this on grid above?
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








# EBS/NBS shelf survey strata
ggplot() +
  geom_sf(data = ebs_layers$survey.strata,
          mapping = aes(fill = factor(Stratum))) +
  geom_sf_text(data = sf::st_centroid(ebs_layers$survey.strata),
               mapping = aes(label = Stratum)) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_fill_viridis_d(name = "Stratum") +
  theme_bw()
