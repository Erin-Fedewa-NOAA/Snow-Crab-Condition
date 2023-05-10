library(dplyr)


#################
### data prep ###
#################
dir("internal_data_jens")
g<- readRDS('internal_data_jens/glob_bloom_type_DECISION_tree_data_feb2023.RDS')

g$peak_timing_all_log[g$gl_gap_sizeTS>32]<-NA
g<- g[!is.na(g$gl_peak_timing_all_log),] # crop data to those

table(g$gl_type) # removing ice_cycle. Erin don't worry about this / just a subset of estimates that were tricky - ice coming and going. I am removing those for your index
g<-g[!g$gl_type=='ice_cycle',]



########################################
### assigning spatial BSIERP regions ###
########################################
# Erin I started with the BSIERP regions. I can share the shape file with you - but for now I just stored it in gitignore.
# You can also just swap the shape file yourself

bsregions <- rgdal::readOGR("internal_data_jens/BSIERP-Regions", layer="BSIERP_regions_2012")
x <- g$gl_lon
y  <- g$gl_lat

#Now you have to define the coordinate system. First, put the x and y data into a
#data frame. By doing this you are creating an object of the class SpatialPoints
d <- data.frame (lon=x, lat=y)
sp::coordinates(d) <- c("lon", "lat")
sp::proj4string(d) <- sp::CRS("+init=epsg:4269")

#this is the NAD83 projection used by federal agencies for lat-long and what was
#used to make the shapefile. Now you have to define a new coordinate system for your
#data because the shapefile may or may not be in lat-long. So, you need to see what
#the CRS definitions are for the shapefile. Since you read it in already, type:
CRSargs <-bsregions@proj4string

#this will return the coordinate system string. Now you can convert your lat-lon
#file to this projection
CRS.new <- sp::CRS(CRSargs@projargs)

#Use the spTransform function from sp package to convert your data
d_new <- sp::spTransform(d, CRS.new)

#Now your data should be in the proper coordinates to compare to the shapefile.
#Finally, we have to compare the data to the shapefile and place the points within
#the polygons
q <- sp::over(d_new, as(bsregions, "SpatialPolygons"))

#This will create a new variable, q, that assigns each row of data set d.new to a
#particular polygon found in the polygons file.
#We can append this variable to the  data set
g$bsregion <- NULL
g$bsregion <- as.character(q)
head(g)
tail(g)

### subsetting data / saving for Erin project ### 
s <- subset(g,select = c(year, gl_lon,   gl_lat , gl_depth ,gl_peak_timing_all_log, gl_cum_anomSST ,gl_mean_ice_frac, gl_ice_retr_roll15,gl_type ,bsregion))

head(s)
# renaming columns to make it easier to follow #
s<-s %>% rename(longitude = gl_lon)
s<-s %>% rename(latitude = gl_lat)
s<-s %>% rename(bottom_depth = gl_depth)

s<-s %>% rename(peak_timing  = gl_peak_timing_all_log)
s<-s %>% rename(SSTcumanom45 = gl_cum_anomSST) # april may sea surf temperature cummulative anomaly 
s<-s %>% rename(spring_ice_fraction = gl_mean_ice_frac ) # consistency of ice cover
s<-s %>% rename(ice_retreat_timing = gl_ice_retr_roll15  )
s<-s %>% rename(bloom_type  = gl_type )

head(s)
write.csv(s,"data/bloom_type_data_JMN.csv")


###
### Erin you can start here for the bloom type analyses / data setup
###

s<-read.csv("data/bloom_type_data_JMN.csv")
head(s)

s<-s[!is.na(s$bsregion),]




