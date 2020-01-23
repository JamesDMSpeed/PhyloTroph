###Testing

library(raster)
library(rgdal)
library(sf)
library(fasterize)
library(ggplot2)
library(leaflet)
library(tmap)
library(dplyr)

#Problems with reading in files
repN <- readOGR("Data/Ranges/Reptiles", "NorwayReptiles")
amphN <- readOGR ('Data/Ranges/Amphibians','amphibian_Norway_IUCN')

#Reading shapefile and layer into spatial vector object
mamN <- readOGR("Data/Ranges/Mammals")#Should be read from repository project directory, not local drive
mamN1 <- mamN

List <- c('Pusa hispida','Cystophora cristata','Halichoerus grypus','Phoca vitulina','Erignathus barbatus')
#mamN1 <- filter(mamN1@data$BINOMIAL, !(BINOMIAL %in% List))
mamN1 <-mamN1[mamN1@data$BINOMIAL != c('Halichoerus grypus','Phoca vitulina'),]
mamN1 <-mamN1[mamN1@data$BINOMIAL != c('Cystophora cristata', 'Erignathus barbatus'),]
mamN1 <-mamN1[mamN1@data$BINOMIAL != 'Pusa hispida',]

#Transformation between datum(s) and conversion between projections, crs = coordinate reference system
mamTrans1<-spTransform(mamN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)

#Creating a RasterLayer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s<-raster(e, resolution=20000, crs=(mamTrans1))
mamSF <- st_as_sf(mamTrans1)
#Filter out species from list in the sf-object
#mamSF <- filter(mamSF, !(BINOMIAL %in% List))

#one layer per species (sf, raster)
stackedDistributionMaps<-fasterize(mamSF,s,by="BINOMIAL")
plot(stackedDistributionMaps)

#merge polygons per species for species richness map, unioning geometries 
speciesPoly<-aggregate(mamTrans1,by="BINOMIAL")
speciesPoly2<-st_as_sf(speciesPoly)

#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
plot(richnessMap)

#Get map of Norway and check the projection
norway<-getData('GADM',level=0,country='NOR')
plot(norway)
#change projetction to the same as mammal-map
norwayTrans <- spTransform(norway,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
extent(norwayTrans)
extent(richnessMap) 

###TESTING PACKAGE TMAP###
qtm(richnessMap)
mammap <- qtm(richnessMap, fill = "BINOMIAL")

#Stacking tmap elements

#assess the species list: mamN@data$BINOMIAL
#summary(richnessMap) - 
#bbox(richnessMap)- explore the bounding area of any spatial object
#proj4string(richnessMap)- explore the projection system of any spatial object
#class(mamN)
#Terrmam <- filter(mamN1@data, !(BINOMIAL %in% List))
