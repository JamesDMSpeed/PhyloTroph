###Testing

library(raster)
library(rgdal)
library(sf)
library(fasterize)
library(ggplot2)
library(leaflet)
library(tmap)

#Reading shapefile and layer into spatial vector object
mamN <- readOGR("C:/Phylotroph/Data/Ranges/Mammals")
mamN1 <- mamN

#Transformation between datum(s) and conversion between projections, crs = coordinate reference system
mamTrans1<-spTransform(mamN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)

#Creating a RasterLayer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s<-raster(e, resolution=20000, crs=(mamTrans1))
mamSF <- st_as_sf(mamTrans1)

#one layer per species (sf, raster)
stackedDistributionMaps<-fasterize(mamSF,s,by="BINOMIAL")
plot(stackedDistributionMaps)


#merge polygons per species for species richness map, unioning geometries 
speciesPoly<-aggregate(mamTrans1,by="BINOMIAL")
speciesPoly2<-st_as_sf(speciesPoly)

#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
plot(richnessMap)

#Seems like the richnessmap is including areas outside the coastline of Norway? Would want to
#crop it. Or Are the pale areas part of the land? 

#assess the species list: mamN@data$BINOMIAL
#summary(richnessMap) - 
#bbox(richnessMap)- explore the bounding area of any spatial object
#proj4string(richnessMap)- explore the projection system of any spatial object

#class(mamN)
#shape <- read_sf(dsn = "C:/Phylotroph/Data/Ranges/Mammals", layer = "NorwayNew_")
#Amph <-read_sf(dsn = "C:/Phylotroph/Data/Ranges/Amphibians", layer = "amphibian_Norway_IUCN")
#Just get a lot of error messages 
#Amphibians<-readOGR("C:/Phylotroph/Data/Ranges/Amphibians")

##leaflet/tmap
#qtm(mamN)
mammap <- qtm(richnessMap)



#Get map of Norway
norway<-getData('GADM',level=0,country='NOR')
#Create clipping polygon
norwayTrans <- spTransform(norway,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
