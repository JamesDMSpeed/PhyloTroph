###Testing


library(raster)
library(rgdal)
library(sf)
library(fasterize)

#Reading shapefile and layer into spatial vector object
mamN <- readOGR("C:/Phylotroph/Data/Ranges/Mammals")
mamN1 <- mamN

#Object with coordinates transformed to the new coordinate reference system
mamTrans1<-spTransform(mamN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)

#description of a Coordinate Reference System (map projection). 
#Setting the resolution of a Raster object
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
richnessMap
plot(richnessMap)


