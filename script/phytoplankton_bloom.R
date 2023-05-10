library(dplyr)
library(ggplot2)
library(cmocean) # just for colors - you can pick a different palette if you want 
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
s<-s[!is.na(s$bsregion),]
head(s)
# getting the crab years #

sub_crab<-c(2019,2021,2022) # picking areas from all mooring areas
s2<- s[s$year %in% sub_crab, ]



aggS<-aggregate(ice_retreat_timing ~longitude*latitude*bsregion,data=s,mean) 
# note here - because there are NA estimates for bloom type (which I filtered out) - the lat and long wobble a bit. That is why this map looks a bit strange. You can safely ignore that
head(aggS)
# loading map 
data_map<-map_data("world2") # map_data from ggmap mapping package. Slow so I leave it out for now.
# note map is in 0-360 longitude +
ggplot()+
  coord_equal(xlim=c(180,205),ylim=c(54,66),ratio = 1.8)+
  geom_polygon(data = data_map, aes(x=long, y = lat, group = group),colour="black", fill="darkgrey")+ # map blanked out
  #geom_tile(data=aggS,aes(x=longitude+360 ,y=latitude,fill=bsregion))#+
  geom_point(data=aggS,aes(x=longitude+360 ,y=latitude,color=as.factor(bsregion)))#+
#scale_fill_gradientn(colours = (cmocean('thermal')(200)),name = "") +
  #geom_point(data=d3,aes(x=longitude,y=latitude),color='white')+
  #geom_point(data=d3,aes(x=effective_contour_longitude ,y=effective_contour_latitude ),color='orange') #+



windows(20,10)
ggplot()+
  coord_equal(xlim=c(180,205),ylim=c(54,66),ratio = 1.8)+
  geom_polygon(data = data_map, aes(x=long, y = lat, group = group),colour="black", fill="darkgrey")+ # map blanked out
  #geom_tile(data=aggS,aes(x=longitude+360 ,y=latitude,fill=bsregion))#+
  geom_point(data=s2,aes(x=longitude+360 ,y=latitude,color=as.factor(bloom_type)),size=5)+
  facet_wrap(.~year,ncol=3)

head(s2)
# ice retreat timing (which link to bloom type)

windows(20,10)
ggplot()+
  coord_equal(xlim=c(180,205),ylim=c(54,66),ratio = 1.8)+
  geom_polygon(data = data_map, aes(x=long, y = lat, group = group),colour="black", fill="darkgrey")+ # map blanked out
  geom_point(data=s2,aes(x=longitude+360 ,y=latitude,color=(ice_retreat_timing )),size=5)+
  scale_color_gradientn(colours = (cmocean('haline')(200)),name = "") +
  facet_wrap(.~year,ncol=3)


### summary of bloom type (percent ice associated vs open water). 
### Hunt et al. 2011 provides a good overview of why a change in bloom type might influence fish. Same could apply for crab
### Happy to explain more
### my only concern here is that we split (by bs region) into relatively few number of obs per region - so the percent calc is a based on 5-10 values. We can talk
head(s)
bl_reg <- s %>% group_by(bsregion,year) %>% summarize (nb_ice = sum(bloom_type=="ice_full"),
                                                       nb_open = sum(bloom_type=="ice_free") ) 
bl_reg$total<- bl_reg$nb_ice+bl_reg$nb_open
bl_reg$perc_open<- (bl_reg$nb_open/bl_reg$total)*100 # this would be your indicator 

head(bl_reg)
### looking at a single region ###
### region 13 and 15 - were difficult ot estimate - so not a lot of esimates (5, 6 years). I would avoid using those. 
bl_reg13<-bl_reg[bl_reg$bsregion==13,]


windows(22,20)
ggplot()+
  geom_point(data=bl_reg,aes(x=year,y=perc_open),size=4,col='dodgerblue')+
  geom_line(data=bl_reg,aes(x=year,y=perc_open),size=2,col='dodgerblue')+
  facet_wrap(.~bsregion,ncol=5)



bl_reg_erin<- bl_reg[bl_reg$year %in% sub_crab, ]

write.csv(bl_reg_erin,"data/final_perc_open_water_2019_2022.csv")


