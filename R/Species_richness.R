### Making species richness maps

### Making species richness maps

library(pacman)
pacman::p_load(raster,
               rgdal,
               sf,
               fasterize,
               ggplot2,
               leaflet,
               tmap,
               dplyr,
               tiff,
               rgeos,update = F)

library(sp)

setwd("C:/Users/Kjers/OneDrive/Dokumenter/9.-10. semester/Master/Master_prooject")

####Reading shapefiles####
#Reading shapefile for mammal
mamN <- readOGR("Data/Ranges/Mammals")#Should be read from repository project directory, not local drive
mamN1 <- mamN

##Reptiles
#repN <- readOGR("Reptiles_","NorwayReptiles")
repN_Sillero <- readOGR("Reptiles_","reptiles_Norway_Sillero")
repN1 <- repN_Sillero

##Amphibians
#amphN <- readOGR("Amphibians_","amphibian_Norway_IUCN")
amphN_Sillero <- readOGR("Amphibians_","amphibian_Norway_Sillero")
amphN1 <- amphN_Sillero

#Plants
#plantList <- read.csv("Data/Ranges/Plants/biodiverse_results_concatenated_good100_utm32_0402_2.csv")

plantSR <- raster("Data/Ranges/Plants/PlantSpeciesRichness.tif")
plantN1 <- plantSR

#Birds
birdN <- readOGR("Birds_","Birbies_fix")
#norBirds <- levels(birdN@data$SCINAME) #No. of species
birdN1 <- birdN
birdN1 <- birdN1[-461,]
#Removing a small polygon from the species of Calidris maritima (small island)
#birdN1@data$REVIEWE <- NULL


#####SR PLANTS####
#Change crs for the plant raster (resolution is correct)
#crs(plantN1) <- "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"
#plantTrans <- projectRaster(plantN1, crs="+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
plantN2 <-plantN1
projection(plantN2) <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
richnessMap_plant <- plantN2

####SR MAMMALS####
##Subset from species from Artsdatabanken
NorTerrMam <- read.csv("Data/Rodlista_terrestrisk.csv")
NorList <- levels(NorTerrMam$Vitenskapelig_navn)
NorMam <- subset(mamN1,mamN1@data$BINOMIAL%in% NorList)

#Transformation between datum(s) and conversion between projections, crs = coordinate reference system
mamTrans1<-spTransform(NorMam,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)
#Creating a RasterLa<yer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s<-raster(e, resolution=20000, crs=(mamTrans1))
mamSF <- st_as_sf(mamTrans1)
#one layer per species (sf, raster)
stackedDistributionMaps<-fasterize(mamSF,s,by="BINOMIAL")
plot(stackedDistributionMaps)

#merge polygons per species for species richness map, unioning geometries 
speciesPoly<-aggregate(mamTrans1,by="BINOMIAL")
plot(speciesPoly)
#speciesPoly_mam <- aggregate(mamSF, by= "BINOMIAL")
#Evaluation error: TopologyException: found non-noded intersection between LINESTRING
speciesPoly2<-st_as_sf(speciesPoly)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
plot(richnessMap)

####SR AMPHIBIANS####
#Transformation between datum(s) and conversion between projections
amphTrans1<-spTransform(amphN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_amph<-extent(amphTrans1)

#make raster and convert to sf-object 
s_amph<-raster(e_amph, resolution=20000, crs=(amphTrans1))
amphSF <- st_as_sf(amphTrans1)
#one layer per species (sf, raster)
#stackedDistributionMaps_amph<-fasterize(amphSF,s_amph,by="amp")
#plot(stackedDistributionMaps_amph)
#merge polygons per species for species richness map, unioning geometries 
speciesPoly_amph<-aggregate(amphTrans1,by="amp")
speciesPoly2_amph<-st_as_sf(speciesPoly_amph)

#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_amph <- fasterize(speciesPoly2_amph,s_amph,field="amp")
plot(richnessMap_amph)
#table(getValues(richnessMap_amph))  

####SR REPTILES####
repTrans1 <-spTransform(repN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_rep<-extent(repTrans1)

#Creating a RasterLayer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s_rep<-raster(e_rep, resolution=20000, crs=(repTrans1))
repSF <- st_as_sf(repTrans1)

#merge polygons per species for species richness map, unioning geometries 
speciesPoly_rep<-aggregate(repTrans1,by="rep")
speciesPoly2_rep<-st_as_sf(speciesPoly_rep)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_rep <- fasterize(speciesPoly2_rep,s_rep,field="rep")
plot(richnessMap_rep)
y <- getValues(richnessMap_rep)

####SR BIRDS####
birdTrans1<-spTransform(birdN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_bird<-extent(birdTrans1)
s_bird<-raster(e_bird, resolution=20000, crs=(birdTrans1))
birdSF <- st_as_sf(birdTrans1)

#one layer per species (sf, raster)
stackedDistributionMaps_bird<-fasterize(birdSF,s_bird,by="SCINAME")
plot(stackedDistributionMaps_bird)

#merge polygons per species for species richness map, unioning geometries 
#ERROR!!
speciesPoly_bird<-raster::aggregate(birdTrans1,by="SCINAME") #ERROR!
speciesPoly2_bird<-st_as_sf(speciesPoly_bird)

richnessMap_bird <- fasterize(speciesPoly2_bird,s_bird,fun="count",field="SCINAME")
plot(richnessMap_bird)

####SR TOTAL####

#Make empty raster
vec <- c(-99551.21, 1340000, 6440080, 7962744)
ext <- extent(vec)
emp_ras <- raster(ext)
projection(emp_ras) <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
res(emp_ras)<- 20000

r1 <- richnessMap_plant
r2 <- richnessMap_rep
r3 <- richnessMap_amph
r4 <- richnessMap
r5 <- richnessMap_bird
r6 <- emp_ras

#Resample to empty raster
r2_rs <- resample(r2,r5)
r3_rs <- resample(r3,r5)
r4_rs <- resample(r4,r5)
r5_rs <- resample(r5,r5)
r1_rs <- resample(r1,r5)

###MAKING A RASTERSTACK
raster_result <- overlay(r1_rs,r2_rs,r3_rs,r4_rs,r5_rs, fun=function(x,y,z,v,c){return(x+y+z+v+c)})
plot(raster_result)

my.stack2 <- stack(r2,r3, RAT=T)
sp <- overlay(r2,r3, fun=function(x,y){return(x+y)})
sp_2 <- overlay(my.stack2, fun = function(x,y)(x+y))
plot(sp_2)  

#combo <- resample(richnessMap,plantSR, method="bilinear")
#total_richness <- aggregate(richness)
#raster_result <- stack(richnessMap,plantSR)
#merged_raster <- raster::merge(richnessMap,plantTrans)
#plot(merged_raster)


