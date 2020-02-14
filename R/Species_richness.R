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
               tiff,update = F)

##Reading shapefiles
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
  plantTrans <- projectRaster(plantSR, res = 20000, crs="+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
  #csr(plantSR) <- "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"

  #Birds
  birdN <- readOGR("Birds_","Birbies")
  levels(birdN@data$SCINAME) #No. of species
  birdN1 <- birdN

#SPECIES RICHNESS MAMMALS
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
  speciesPoly2<-st_as_sf(speciesPoly)
  #count number of overlapping polygons (which is 1 per species so counting gives richness)
  richnessMap <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
  plot(richnessMap)

###AMPHIBIANS
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

###REPTILES
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

##BIRDS SR
  birdTrans1<-spTransform(birdN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
  e_bird<-extent(birdTrans1)
  s_bird<-raster(e_bird, resolution=20000, crs=(birdTrans1))
  birdSF <- st_as_sf(birdTrans1)
  
  #one layer per species (sf, raster)
  stackedDistributionMaps_bird<-fasterize(birdSF,s_bird,by="SCINAME")
  plot(stackedDistributionMaps_bird)

  #merge polygons per species for species richness map, unioning geometries 
  #ERROR!!
  speciesPoly_bird<-aggregate(birdSF,by="SCINAME") #ERROR!
  speciesPoly2_bird<-st_as_sf(speciesPoly_bird)
  richnessMap_bird <- fasterize(speciesPoly2_bird,s_bird,fun="count",field="SCINAME")
  plot(richnessMap_bird)
  #birds <- row.names(as(birdTrans1, "data.frame"))
  #birds1 <- spChFIDs(birdTrans1,as.character(birdTrans1@data$SCINAME))

###MAKING A RASTERSTACK
  #combo <- resample(richnessMap,plantSR, method="bilinear")
  #total_richness <- aggregate(richness)
  #raster_result <- stack(richnessMap,plantSR)
  merged_raster <- raster::merge(richnessMap,plantTrans)
  plot(merged_raster)
