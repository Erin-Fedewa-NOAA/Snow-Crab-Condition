#Create metadata files for snow crab condition dataset

#Author: Erin Fedewa
#Metadata created on 2/2/23

library(dataspice)
library(tidyverse)

# Data 
scdat <- read.csv("./data/condition_haul_master.csv")

######################################

#Create metadata folder and files 
create_spice(dir="./")

#Data creator details 
edit_creators()

#Access details: where data can be accessed 
prep_access()
edit_access()

#Biblio details: metadata (title, spatiotemporal coverage)
range(scdat$year) 
range(scdat$mid_longitude, na.rm=T)
range(scdat$mid_latitude, na.rm=T)

edit_biblio()

#Attribute details: details about the variables in the dataset
prep_attributes()
edit_attributes()

#Create metadata file 
write_spice()

jsonlite::read_json(here::here("data", "metadata", "dataspice.json")) %>% 
  listviewer::jsonedit()

#BuildReadMe Site 
build_site()
