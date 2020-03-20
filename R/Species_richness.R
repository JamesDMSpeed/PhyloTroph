### Making species richness maps
###loading library####
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
               sp,
               gdalUtils,
               GISTools,
               spatial,
               spatial.tools,
               plyr,
               tidyverse,
               rgeos,update = F)

####READING SHAPEFILES####
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

#Get the right names in birdN1, so it can be alignes with traitsData
###Changing name of some birds####  
birdN1@data[["SCINAME"]] <- mapvalues((birdN1@data[["SCINAME"]]), 
                                      c("Bubo scandiacus","Calidris falcinellus","Calidris pugnax",
                                        "Hydrocoloeus minutus","Periparus ater","Poecile cinctus",
                                        "Lophophanes cristatus","Poecile montanus","Poecile palustris",
                                        "Cyanistes caeruleus","Acanthis flammea","Ardenna gravis","Ardenna grisea",
                                        "Chloris chloris","Cyanecula svecica","Dryobates minor","Hydrobates leucorhous",
                                        "Linaria cannabina","Linaria flavirostris","Lyrurus tetrix","Mareca penelope",
                                        "Mareca strepera","Spatula clypeata","Spatula querquedula","Spinus spinus"),
                                      c("Bubo scandiaca","Limicola falcinellus","Philomachus pugnax","Larus minutus",
                                        "Parus ater","Parus cinctus","Parus cristatus","Parus montanus",
                                        "Parus palustris","Parus caeruleus","Carduelis flammea","Puffinus gravis",
                                        "Puffinus griseus","Carduelis chloris","Luscinia svecica","Dendrocopos minor",
                                        "Oceanodroma leucorhoa","Carduelis cannabina","Carduelis flavirostris",
                                        "Tetrao tetrix","Anas penelope","Anas strepera","Anas clypeata",
                                        "Anas querquedula","Carduelis spinus"),warn_missing = T)
#####SR PLANTS####
GDALinfo("Data/Ranges/Plants/PlantSpeciesRichness.tif")
nlayers(plantN1)

####SR MAMMALS####
##Remove species that are not terrestrial (?)
RemoveSpecies <- c("Cystophora cristata","Erignathus barbatus",
                   "Halichoerus grypus","Phoca vitulina","Pusa hispida")
mamN2 <- subset(mamN1,!(mamN1@data$BINOMIAL%in% RemoveSpecies)) #57 species

#Transformation between datum(s) and conversion between projections
mamTrans2<-spTransform(mamN2,crs(plantSR))
#Creating a RasterLa<yer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
e_mam<-extent(mamTrans2)
s_mam<-raster(e_mam, resolution=20000, crs=(mamTrans2))
mamSF2 <- st_as_sf(mamTrans2)
#Rasterize set of polygons (sf to raster)
stackedDistributionMaps_mam<-fasterize(mamSF2,plantN1,by="BINOMIAL") #rasterbrick
#merge polygons per species for species richness map, unioning geometries 
speciesPoly_mam<-aggregate(mamTrans2,by="BINOMIAL")
speciesPoly2_mam<-st_as_sf(speciesPoly_mam)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_mam <- fasterize(speciesPoly2_mam,plantN1,fun="count",field="BINOMIAL")

####SR AMPHIBIANS####
amphTrans1<-spTransform(amphN1,crs(plantSR))
e_amph<-extent(amphTrans1)
s_amph<-raster(e_amph, resolution=20000, crs=(amphTrans1))
amphSF <- st_as_sf(amphTrans1)
stackedDistributionMaps_amph<-fasterize(amphSF,s_amph,by="amp")
speciesPoly_amph<-aggregate(amphTrans1,by="amp")
speciesPoly2_amph<-st_as_sf(speciesPoly_amph)
richnessMap_amph <- fasterize(speciesPoly2_amph,s_amph,field="amp")
plot(richnessMap_amph)

####SR REPTILES####
repTrans1 <-spTransform(repN1,crs(plantSR))
e_rep<-extent(repTrans1)
s_rep<-raster(e_rep, resolution=20000, crs=(repTrans1))
repSF <- st_as_sf(repTrans1)
stackedDistributionMaps_rep<-fasterize(repSF,s_rep,by="rep")
speciesPoly_rep<-aggregate(repTrans1,by="rep")
speciesPoly2_rep<-st_as_sf(speciesPoly_rep)
richnessMap_rep <- fasterize(speciesPoly2_rep,s_rep,field="rep")
plot(richnessMap_rep)

####SR BIRDS####
birdTrans1<-spTransform(birdN1,crs(plantSR))
e_bird<-extent(birdTrans1)
s_bird<-raster(e_bird, resolution=20000, crs=(birdTrans1))
birdSF <- st_as_sf(birdTrans1)
stackedDistributionMaps_bird<-fasterize(birdSF,s_bird,by="SCINAME")
speciesPoly_bird<-raster::aggregate(birdTrans1,by="SCINAME")
speciesPoly2_bird<-st_as_sf(speciesPoly_bird)
richnessMap_bird <- fasterize(speciesPoly2_bird,s_bird,fun="count",field="SCINAME")
plot(richnessMap_bird)

####SR TOTAL####
r1 <- plantN1
r2 <- richnessMap_rep
r3 <- richnessMap_amph
r4 <- richnessMap_mam
r5 <- richnessMap_bird
#Stack: all layers need same extent and resolution
#Make all layers have the same dimensions (number of row & columns)
#Changing the col/row might lead to loss of values associated w/ the layer
r111 <- modify_raster_margins(r1, extent_delta = c(0,0,2,0), value = NA)
r222 <- modify_raster_margins(r2, extent_delta = c(1,0,0,0), value = NA)
r333 <- modify_raster_margins(r3, extent_delta = c(1,0,0,0), value = NA)
r444 <- modify_raster_margins(r4, extent_delta = c(1,0,2,0), value = NA)
#Use plantlayer as standard extent
plant <- extent(r111)
rep <- extent(r222)
amp <- extent (r333)
mam <- extent (r444)
bir <- extent (r4)

#Change extent of all layers
r11 <- setExtent(r111,plant,keepres = T)
r22 <- setExtent(r222,plant,keepres = T)
r33 <- setExtent(r333,plant,keepres = T)
r44 <- setExtent(r444,plant,keepres = T)
r55 <- setExtent(r5,plant,keepres = T)

#Make rasterstack with all of the richness layers
my.stack = stack(r11,r22,r33,r44,r55)
names(my.stack)<-c('PlantSR','ReptileSR','AmphibianSR','MammalSR','BirdSR')
#Extract values from the rasterstack, in order to get SR within each cell
stackValdf<-raster::extract(my.stack,(1:ncell(my.stack)), df = T)
stackValdf$total <- rowSums(stackValdf[,2:6], na.rm = T, dims = 1)
View(stackValdf)
plot(stackValdf$PlantSR,stackValdf$MammalSR)

#####Load data####
#Trait data Mammals
traitDataMam <- read.table("MamFuncDat.txt", sep = '\t', header = T, fill = T)
NorMam <- traitDataMam %>%
  filter(Scientific%in% mamN2@data[["BINOMIAL"]]) #52

#Trait data Birds
traitDataBird <- read.table("BirdFuncDat.txt", sep = '\t',header = T, fill = T, quote ='')
NorwayBirds <- read_csv("NorwayBirds.csv") #250 species (251 in birdN1, but Pinguinus impennis is extinct)
NorBird <- traitDataBird %>%
  filter(Scientific %in% birdN1@data[["SCINAME"]]) #250 

#Trait data Amphibians
traitDataAmph <- read.csv("AmphiBIO_v1.csv")
NorAmphNames <- c('Bufo bufo','Lissotriton vulgaris','Pelophylax lessonae','Rana arvalis',
                  'Rana temporaria','Triturus cristatus')
NorAmph <- traitDataAmph %>%
  filter(Species%in% NorAmphNames) #6

#Reptiles
NorRepNames <- c('Angius fragilis','Coranella austriaca','Natrix natrix','Vipera berus',
                 'Lacerta vivipara')

#####MAMMALS ELTON GROUPS####
MamHerbi <- NorMam[NorMam$Diet.PlantO>=70,] #12
nameMamHerbi<- MamHerbi$Scientific

MamGrani <- NorMam[NorMam$Diet.Seed>=70,] #0 (if >= 50, Sciurus vulgaris is considered granivore)
nameMamGrani <- MamGrani$Scientific

MamCarni <- rbind((NorMam[NorMam$Diet.Vend >=50,]),
                  (NorMam[NorMam$Diet.Vect >=50,]),
                  (NorMam[NorMam$Diet.Vfish>=50,]),
                  (NorMam[NorMam$Diet.Vunk>=50]),
                  (NorMam[NorMam$Diet.Scav>=50]))#8
nameMamCarni<- MamCarni$Scientific

MamInsect <- NorMam[NorMam$Diet.Inv>=50,] #17
nameMamInsect <- MamInsect$Scientific

MamOmni <- NorMam[NorMam$Diet.Inv<50&
                    NorMam$Diet.Vend<50&
                    NorMam$Diet.Vect<50&
                    NorMam$Diet.Vfish<50&
                    NorMam$Diet.Vunk<50&
                    NorMam$Diet.Scav<50&
                    NorMam$Diet.Seed<70&
                    NorMam$Diet.PlantO<70,] #15
nameMamOmnivore <- MamOmni$Scientific

checkMam <-  NorMam%>%
  filter((!(Scientific %in% MamCarni$Scientific))&
           !(Scientific %in% MamHerbi$Scientific)&
           !(Scientific %in% MamGrani$Scientific)&
           !(Scientific %in% MamInsect$Scientific)&
           !(Scientific %in% MamOmni$Scientific)) #0

####BIRDS Elton GROUPS####
BirdHerbi <- NorBird[NorBird$Diet.PlantO>=70,] #17
nameBirdHerbi<- BirdHerbi$Scientific

BirdGrani <- NorBird[NorBird$Diet.Seed>=70,] #4
nameBirdGrani <- BirdGrani$Scientific

BirdCarni <- rbind((NorBird[NorBird$Diet.Vend >=50,]),
                   (NorBird[NorBird$Diet.Vect >=50,]),
                   (NorBird[NorBird$Diet.Vfish>=50,]),
                   (NorBird[NorBird$Diet.Vunk>=50]),
                   (NorBird[NorBird$Diet.Scav>=50])) #50
nameBirdCarni<- BirdCarni$Scientific

BirdInsect <- NorBird[NorBird$Diet.Inv>=50,] #117
nameBirdInsect <- BirdInsect$Scientific

BirdOmni <- NorBird[NorBird$Diet.Inv<50&
                      NorBird$Diet.Vend<50&
                      NorBird$Diet.Vect<50&
                      NorBird$Diet.Vfish<50&
                      NorBird$Diet.Vunk<50&
                      NorBird$Diet.Scav<50&
                      NorBird$Diet.Seed<70&
                      NorBird$Diet.PlantO<70,] #63
nameBirdOmnivore <- BirdOmni$Scientific
##Total of 251 species? 

checkBird <-  NorBird%>%
  filter((!(Scientific %in% BirdCarni$Scientific))&
           !(Scientific %in% BirdHerbi$Scientific)&
           !(Scientific %in% BirdGrani$Scientific)&
           !(Scientific %in% BirdInsect$Scientific)&
           !(Scientific %in% BirdOmni$Scientific))#0

#Check if some birds are in two categories
BirdInsect %in% BirdOmni | BirdInsect%in% BirdGrani | BirdInsect%in% BirdHerbi | BirdInsect%in% BirdCarni
BirdGrani %in% BirdInsect | BirdGrani%in% BirdOmni | BirdGrani%in% BirdHerbi | BirdGrani%in% BirdCarni
BirdHerbi %in% BirdInsect | BirdHerbi%in% BirdGrani | BirdHerbi%in% BirdOmni | BirdHerbi%in% BirdCarni
BirdCarni %in% BirdInsect | BirdCarni%in% BirdGrani | BirdCarni%in% BirdHerbi | BirdCarni%in% BirdOmni
BirdOmni %in% BirdInsect | BirdOmni%in% BirdGrani | BirdOmni%in% BirdHerbi | BirdOmni%in% BirdCarni

####AMPHIBIANS Elton GROUPS####





#####GROUPS MAMMALS/FEEDING ####
#SR herbivore mammals
herbivoreMam <- subset(mamN2, BINOMIAL %in% nameMamHerbi)
herbivoreMamTrans<-spTransform(herbivoreMam,crs(plantSR))
e_HM <-extent(herbivoreMamTrans)
s_HM<-raster(e_HM, resolution=20000, crs=(herbivoreMamTrans))
speciesPoly_HM <- aggregate(herbivoreMamTrans, by ="BINOMIAL")
speciesPoly_HM2 <- st_as_sf(speciesPoly_HM)
#adjusting
mamHerbSR <- fasterize(speciesPoly_HM2,s_HM,fun="count",field="BINOMIAL") 
mamHerbSR <- modify_raster_margins(mamHerbSR, extent_delta = c(1,0,2,0), value = NA)
mamHerbSR <- setExtent(mamHerbSR,plant,keepres = T)

#SR carnivore mammals
carnivoreMam <- subset(mamN2, BINOMIAL %in% nameMamCarni)
carnivoreMamTrans<-spTransform(carnivoreMam,crs(plantSR))
e_CM <-extent(carnivoreMamTrans)
s_CM<-raster(e_CM, resolution=20000, crs=(carnivoreMamTrans))
speciesPoly_CM <- aggregate(carnivoreMamTrans, by ="BINOMIAL")
speciesPoly_CM2 <- st_as_sf(speciesPoly_CM)
#Adjusting
mamCarnSR <- fasterize(speciesPoly_CM2,s_CM,fun="count",field="BINOMIAL")
mamCarnSR <- modify_raster_margins(mamCarnSR, extent_delta = c(2,0,3,0), value = NA)
mamCarnSR <- setExtent(mamCarnSR,plant,keepres = T)

#SR insectivores mammals
insectivoreMam <- subset(mamN2, BINOMIAL %in% nameMamInsect)
insectivoreMamTrans<-spTransform(insectivoreMam,crs(plantSR))
e_IM <-extent(insectivoreMamTrans)
s_IM<-raster(e_IM, resolution=20000, crs=(insectivoreMamTrans))
speciesPoly_IM <- aggregate(insectivoreMamTrans, by ="BINOMIAL")
speciesPoly_IM2 <- st_as_sf(speciesPoly_IM)
#Adjusting
mamInsectSR <- fasterize(speciesPoly_IM2,s_IM,fun="count",field="BINOMIAL")
mamInsectSR <- modify_raster_margins(mamInsectSR, extent_delta = c(1,0,2,0), value = NA)
mamInsectSR <- setExtent(mamInsectSR,plant,keepres = T)

#SR granivore mammals
granivoreMam <- subset(mamTrans2, BINOMIAL %in% nameMamGrani) #0

#SR omnivore mammals
omnivoreMam <- subset(mamN2, BINOMIAL %in% nameMamOmnivore)
omnivoreMamTrans<-spTransform(omnivoreMam,crs(plantSR))
e_OM <-extent(omnivoreMamTrans)
s_OM<-raster(e_OM, resolution=20000, crs=(omnivoreMamTrans))
speciesPoly_OM <- aggregate(omnivoreMamTrans, by ="BINOMIAL")
speciesPoly_OM2 <- st_as_sf(speciesPoly_OM)
#Adjusting
mamOmniSR <- fasterize(speciesPoly_OM2,s_OM,fun="count",field="BINOMIAL")
mamOmniSR <- modify_raster_margins(mamOmniSR, extent_delta = c(1,0,2,0), value = NA)
mamOmniSR <- setExtent(mamOmniSR,plant,keepres = T)

###GROUPS BIRD/FEEDING####
#SR herbivore birds
herbivoreBird <- subset(birdN1, SCINAME %in% nameBirdHerbi)
herbivoreBirdTrans<-spTransform(herbivoreBird,crs(plantSR))
e_HB <-extent(herbivoreBirdTrans)
s_HB<-raster(e_HB, resolution=20000, crs=(herbivoreBirdTrans))
speciesPoly_HB <- aggregate(herbivoreBirdTrans, by ="SCINAME")
speciesPoly_HB2 <- st_as_sf(speciesPoly_HB)
birdHerbSR <- fasterize(speciesPoly_HB2,s_HB,fun="count",field="SCINAME") 
#adjusting
birdHerbSR <- modify_raster_margins(birdHerbSR, extent_delta = c(1,0,2,0), value = NA)
birdHerbSR <- setExtent(birdHerbSR,plant,keepres = T)

#SR carnivore birds
carnivoreBird <- subset(birdN1, SCINAME %in% nameBirdCarni)
carnivoreBirdTrans<-spTransform(carnivoreBird,crs(plantSR))
e_CB <-extent(carnivoreBirdTrans)
s_CB<-raster(e_CB, resolution=20000, crs=(carnivoreBirdTrans))
speciesPoly_CB <- aggregate(carnivoreBirdTrans, by ="SCINAME")
speciesPoly_CB2 <- st_as_sf(speciesPoly_CB)
birdCarnSR <- fasterize(speciesPoly_CB2,s_CB,fun="count",field="SCINAME")
#adjusting
birdCarnSR <- setExtent(birdCarnSR,plant,keepres = T)

#SR granivore birds
granivoreBird <- subset(birdN1, SCINAME %in% nameBirdGrani) #4 species
granivoreBirdTrans <- spTransform(granivoreBird, crs(plantSR))
e_GB <- extent(granivoreBirdTrans)
s_GB <- raster(e_GB, resolution=20000, crs(granivoreBirdTrans))
speciesPoly_GB <- aggregate(granivoreBirdTrans, by ="SCINAME")
speciesPoly_GB2 <- st_as_sf(speciesPoly_GB)
birdGranSR <- fasterize(speciesPoly_GB2,s_GB,fun="count",field="SCINAME")
#Adjusting
birdGranSR <- modify_raster_margins(birdGranSR, extent_delta = c(1,0,6,0), value = NA)
birdGranSR <- setExtent(birdGranSR,plant,keepres = T)

#SR insectivore birds
insectivoreBird <- subset(birdN1, SCINAME %in% nameBirdInsect) 
insectivoreBirdTrans <- spTransform(insectivoreBird, crs(plantSR))
e_IB <- extent(insectivoreBirdTrans)
s_IB <- raster(e_IB, resolution=20000, crs(insectivoreBirdTrans))
speciesPoly_IB <- aggregate(insectivoreBirdTrans, by ="SCINAME")
speciesPoly_IB2 <- st_as_sf(speciesPoly_IB)
#Adjusting
birdInsectSR <- fasterize(speciesPoly_IB2,s_IB,fun="count",field="SCINAME")
birdInsectSR <- setExtent(birdInsectSR,plant,keepres = T)

#SR omnivore birds
omnivoreBird <- subset(birdN1, SCINAME %in% nameBirdOmnivore) 
omnivoreBirdTrans <- spTransform(omnivoreBird, crs(plantSR))
e_OB <- extent(omnivoreBirdTrans)
s_OB <- raster(e_OB, resolution=20000, crs(omnivoreBirdTrans))
speciesPoly_OB <- aggregate(omnivoreBirdTrans, by ="SCINAME")
speciesPoly_OB2 <- st_as_sf(speciesPoly_OB)
#Adjusting
birdOmnivoreSR <- fasterize(speciesPoly_OB2,s_OB,fun="count",field="SCINAME")
birdOmnivoreSR <- setExtent(birdOmnivoreSR,plant,keepres = T)

####Raster stacks####
###Stack mammals
mammalStack <- stack(mamCarnSR,mamHerbSR,mamInsectSR,mamOmniSR)
names(mammalStack)<-c('mamCarnSR','mamHerbSR','mamInsectSR','mamOmniSR')
mammalValdf<-raster::extract(mammalStack,(1:ncell(mammalStack)), df = T)
mammalValdf$total <- rowSums(mammalValdf[,2:ncol(mammalValdf)], na.rm = T, dims = 1)
View(mammalValdf)
plot(mammalStack)

#Extract species from each cell from Mammal rasterlayer
dim(stackedDistributionMaps_mam) #79,54
stackedDistributionMaps_mam2 <- modify_raster_margins(stackedDistributionMaps_mam, extent_delta = c(1,0,2,0), value = NA)
mammalPoly <- raster::extract(stackedDistributionMaps_mam2,(1:ncell(stackedDistributionMaps_mam2)), df = T)
mammalPoly$total <- rowSums(mammalPoly[,2:ncol(mammalPoly)], na.rm = T, dims = 1)
View(mammalPoly)

###Stack birds
birdStack <- stack(birdCarnSR, birdHerbSR, birdGranSR, birdInsectSR, birdOmnivoreSR)
names(birdStack)<-c('birdCarnSR','birdHerbSR','birdGranSR','birdInsectSR', 'birdOmniSR')
birdValdf<-raster::extract(birdStack,(1:ncell(birdStack)), df = T)
birdValdf$total <- rowSums(birdValdf[,2:ncol(birdValdf)], na.rm = T, dims = 1)
plot(birdStack)
View(birdValdf)

#Extract species within each cell from bird rasterlayer
dim(stackedDistributionMaps_bird) #dim = 81,55
birdPoly <- raster::extract(stackedDistributionMaps_bird,(1:ncell(stackedDistributionMaps_bird)), df = T)
birdPoly$total <- rowSums(birdPoly[,2:ncol(birdPoly)], na.rm = T, dims = 1)
View(birdPoly)


allStack <- stack(r11,birdHerbSR, birdInsectSR, birdGranSR, birdCarnSR, birdOmnivoreSR,
                  mamCarnSR,mamHerbSR,mamInsectSR,mamOmniSR)
names(allStack)<- c('plantSR','birdHerbSR','birdInsectSR','birdGranSR','birdCarnSR','birdOmnivoreSR',
                    'mamCarnSR','mamHerbSR','mamInsectSR','mamOmniSR')
allStackValdf <- raster::extract(allStack,(1:ncell(allStack)), df = T)
View(allStack)
plot(allStack)
plot(allStack$plantSR, allStack$mamHerbSR)

#####  
raster_result <- overlay(r11,r22,r33,r44,r55,fun=function(x,y,z,a,b){return(x+y+z+a+b)})
plot(raster_result)
raster_result
spplot(raster_result)
table(getValues(raster_result))
