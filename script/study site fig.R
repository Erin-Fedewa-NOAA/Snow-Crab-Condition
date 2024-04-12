# Goal: Create Fig 1 study area/cpue/temp plots

#NOTE: Accessing EBS survey shapefiles using the akgfmaps package

library(akgfmaps)
library(ggridges)
library(patchwork)
library(sf)
library(hrbrthemes)

#load data 
condition_master <- read.csv("./data/total_FA_master.csv")

#install.packages("remotes")
#remotes::install_github("afsc-gap-products/akgfmaps")

#colors
my_colors <- c("#D55E00","#9ECAE1", "#4292C6", "#084594")
options(scipen = 999)

#################################

# Plotting layers
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

#Transform crab data 
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
# Convert lat/long to an sf object
  st_as_sf(coords = c("mid_longitude", "mid_latitude"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat
  
#And now plot survey grid with crab data  
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add crab layers
  geom_sf(data=crab_dat, aes(size = n_crab), color = "#9ECAE1") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  labs(x="", y="", size = expression(paste("Snow crab \n samples"))) +
  facet_wrap(~year) -> map

###################################
#(B) and (C) Panel Figures:
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
  ylab("Mean Snow Crab Density") + xlab("") +
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme, scales = "free_y") +
  geom_vline(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 1, y=155, label = "Mid-collapse"),
    size = 2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Eastern Bering Sea"), aes(x = 3, y=155, label = "Post-collapse"),
    size = 2, color = "#0072B2") +
  geom_vline(data = subset(cpue.dat, lme == "Northern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Northern Bering Sea"), aes(x = 1, y=515, label = "Mid-collapse"),
            size = 2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Northern Bering Sea"), aes(x = 3, y=515, label = "Post-collapse"),
            size = 2, color = "#0072B2") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 12, axis_text_size =10) +
  theme(legend.position="none") +
    theme(strip.text = element_text(colour = "grey40", hjust = .5)) +
  theme(panel.grid.major.x = element_blank()) -> mean_cpue_plot

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
    ylab("Mean Bottom Temperature") + xlab("") +
    scale_fill_manual(values=my_colors) +
    facet_wrap(~lme) +
    #geom_vline(aes(xintercept = 1.5), linetype="dashed") +
    #geom_text(aes(x = 1, y=3.3, label = "Marine heat wave"),
              #size = 2.4, color = "#D55E00") +
    #geom_text(aes(x = 3, y=3.3, label = "Post-heat wave"),
              #size = 2.4, color = "#0072B2") +
    theme_ipsum(axis_title_just = "cc", axis_title_size = 12, axis_text_size =10) +
    theme(legend.position="none") +
    theme(strip.text = element_text(colour = "grey40", hjust = .5)) +
    theme(panel.grid.major.x = element_blank())  -> mean_temp_plot
  
#Figure 1 for ms: combined map, temp and cpue plot
mean_cpue_plot / plot_spacer() / mean_temp_plot  + plot_layout(heights = c(6, -3 , 6)) +
  plot_annotation(tag_levels = list(c('b', 'c'))) -> b_c_plot

map + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position  = c(.1, .95)) -> a_plot

ggarrange(a_plot, b_c_plot, ncol=2, nrow=1)
ggsave("./figures/Fig1.png", height=8 , width=8, units="in")

###############
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





