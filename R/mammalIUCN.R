
library(raster)
library(rgdal)
library(sf)
library(fasterize)


mamN <- readOGR("NorwayNew_.shp")
mamN1 <- mamN
mamTrans1<-spTransform(mamN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)
s<-raster(e, resolution=10000, crs=(mamTrans1))
mamSF <- st_as_sf(mamTrans1)
#one layer per species
stackedDistributionMaps<-fasterize(mamSF,s,by="BINOMIAL")
plot(stackedDistributionMaps)
#merge polygons per species for species richness map
speciesPoly<-aggregate(mamTrans1,by="BINOMIAL")
speciesPoly2<-st_as_sf(speciesPoly)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
plot(richnessMap)
