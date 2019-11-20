
library(raster)
library(rgdal)
library(raster)
library(sf)
library(fasterize)

bird <- readOGR("Data/Ranges/Birds/Norway_bird.shp")
bird1 <- bird
birdTrans1<-spTransform(bird1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(birdTrans1)
s<-raster(e, resolution=10000, crs=(birdTrans1))
birdSF <- st_as_sf(birdTrans1)
#one layer per species, so here 98 distribution maps
birdstacked<-fasterize(birdSF,s,by="SCINAME")
plot(birdstacked)
#merge polygons per species for species richness map
birdPoly<-aggregate(birdTrans1,by="SCINAME")
birdPoly2<-st_as_sf(birdPoly)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
birdrichness <- fasterize(birdPoly2,s,fun="count",field="SCINAME")
plot(birdrichness)
