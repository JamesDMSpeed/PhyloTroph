###TRY DIFFERENT PLANT RASTER###
####loading library####
library(raster)
library(sf)
library(rgdal)
library(fasterize)
library(ggplot2)
library(tmap)
library(dplyr)
library(spatial)
library(rgdal)
library(leaflet)
library(dplyr)
library(sp)
library(tiff)
library(gdalUtils)
library(GISTools)
library(spatial)
library(spatial.tools)
library(plyr)
library(tidyverse)
library(stats)
library(rgeos)
library(reshape)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(latticeExtra)

####READING SHAPEFILES####
#Reading shapefile for mammal
mamN <- readOGR("Data/Ranges/Mammals")#Should be read from repository project directory, not local drive
mamN1 <- mamN
##Reptiles
repN_Sillero <- readOGR("Data/Ranges/reptiles_Norway_Sillero")
repN1 <- repN_Sillero
##Amphibians
amphN_Sillero <- readOGR("Data/Ranges/amphibian_Norway_Sillero")
amphN1 <- amphN_Sillero
#Plants
#plantSR <- raster("Data/Ranges/Plants/PlantSpeciesRichness.tif")
#plantN1 <- plantSR

load("/Data/brick_native.RData") #https://ntnu.app.box.com/s/050j37osuazm4ezl9rlz9oc3srdwz9kb
#Extract species layer that is sum of all species (SR)
layer_1119 <- subset(brick_native,"index_1119")
plantN2 <- layer_1119

#Birds
birdN <- readOGR("Data/Ranges/Birds_fix")
birdN1 <- birdN
#Removing a small polygon from the species of Calidris maritima (small island)
birdN1 <- birdN1[-461,]
birdN1 <- birdN1[(birdN1@data$SCINAME != "Pinguinus impennis"),] #Removing Pinguinus impennis
#Map of norway
norway<-raster::getData('GADM',country='NOR',level=0)
norwayutm<-spTransform(norway,crs(plantN2))
norwayraster <- rasterize(norwayutm, plantN2)

##MAP OF NORWAY
#Background from GADM
map_nor<-tm_shape(norwayutm)+tm_borders(alpha = 0.5,col="black")+tm_layout(frame = F)

norwaySF <- st_as_sf(norwayutm)
coast_outline <- as(st_geometry(norwaySF), Class = "Spatial")

####BIRD PREPERATION####
#Get the right names in birdN1, so it can be alignes with traitsData
#mapvalues(x, from, to)
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

#Trait data Birds
traitDataBird <- read.table("BirdFuncDat.txt", sep = '\t',header = T, fill = T, quote ='')
#NorwayBirds <- read_csv("NorwayBirds.csv") #250 species
NorBird <- traitDataBird %>%
  filter(Scientific %in% birdN1@data[["SCINAME"]]) #250 

#Filtering out birds who are predominantly feeding pelagic  
NorBirdTerr <- NorBird[NorBird$PelagicSpecialist == 0,] #225
#Subset birdN1 to the species within NorBirdTerr
birdN2 <- subset(birdN1, birdN1$SCINAME %in% NorBirdTerr$Scientific) #225
setdiff(birdN1$SCINAME, birdN2$SCINAME) #25

####SR MAMMALS####
##Remove species that are not terrestrial from mamN1
RemoveSpecies <- c("Cystophora cristata","Erignathus barbatus",
                   "Halichoerus grypus","Phoca vitulina","Pusa hispida","Dama dama",
                   "Micromys minutus","Sorex isodon")
mamN2 <- subset(mamN1,!(mamN1@data$BINOMIAL%in% RemoveSpecies)) 
#Remove polygons with species that are extinct
mamN2 <- mamN2[mamN2@data$LEGEND != "Extinct",] #tot 49 species
mamTrans1<-spTransform(mamN2,crs(plantN2))
mamSF <- st_as_sf(mamTrans1)
#Rasterize set of polygons (sf to raster)
stackedDistributionMaps_mam<-fasterize(mamSF,plantN2,by="BINOMIAL") #rasterbrick
#Fit to the extent of plantraster (plantN2)
stackedDistributionMaps_mam[is.na(stackedDistributionMaps_mam)] <- 0
stackedDistributionMaps_mam <- mask(stackedDistributionMaps_mam, plantN2)

speciesPoly_mam<-aggregate(mamTrans1,by="BINOMIAL")
speciesPoly2_mam<-st_as_sf(speciesPoly_mam)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_mam <- fasterize(speciesPoly2_mam,plantN2,fun="count",field ="BINOMIAL")

####SR AMPHIBIANS####
amphTrans1<-spTransform(amphN1,crs(plantN2))
amphSF <- st_as_sf(amphTrans1)
stackedDistributionMaps_amph<-fasterize(amphSF,plantN2,by="amp")
speciesPoly_amph<-aggregate(amphTrans1,by="amp")
speciesPoly2_amph<-st_as_sf(speciesPoly_amph)
richnessMap_amph <- fasterize(speciesPoly2_amph,plantN2,field="amp")

####STACKED AMPHIBIANS####  
##buf_buf
stacked_bufo<-fasterize(amphSF,plantN2,by="buf_buf")#rasterbrick
#Remove buf_buf = 0
stacked_bufo <- dropLayer(stacked_bufo,1)
speciesPoly_buf <- aggregate(amphTrans1, by="buf_buf")
speciesPoly2_buf<-st_as_sf(speciesPoly_buf)
speciesPoly2_buf <- speciesPoly2_buf[2,]
#colnames(speciesPoly2_buf)[1] <- "x"

##liss_vul
stacked_liss <- fasterize(amphSF,plantN2,by="liss_vul")
stacked_liss <- dropLayer(stacked_liss,1)
speciesPoly_liss <- aggregate(amphTrans1, by="liss_vul")
speciesPoly2_liss<-st_as_sf(speciesPoly_liss)
speciesPoly2_liss <- speciesPoly2_liss[2,]
#colnames(speciesPoly2_liss)[1] <- "x"

##pelp_esc
stacked_pelp <- fasterize(amphSF,plantN2,by="pelp_esc")
stacked_pelp <- dropLayer(stacked_pelp,1)
speciesPoly_pelp <- aggregate(amphTrans1, by="pelp_esc")
speciesPoly2_pelp<-st_as_sf(speciesPoly_pelp)
speciesPoly2_pelp <- speciesPoly2_pelp[2,]
#colnames(speciesPoly2_pelp)[1] <- "x"

##ran_arv
stacked_ranA <- fasterize(amphSF,plantN2,by="ran_arv")
stacked_ranA <- dropLayer(stacked_ranA,1)
speciesPoly_ranA <- aggregate(amphTrans1, by="ran_arv")
speciesPoly2_ranA<-st_as_sf(speciesPoly_ranA)
speciesPoly2_ranA <- speciesPoly2_ranA[2,]
#colnames(speciesPoly2_ranA)[1] <- "x"

##ran_temp
stacked_ranT <- fasterize(amphSF,plantN2,by="ran_temp")
stacked_ranT <- dropLayer(stacked_ranT,1)
speciesPoly_ranT <- aggregate(amphTrans1, by="ran_temp")
speciesPoly2_ranT<-st_as_sf(speciesPoly_ranT)
speciesPoly2_ranT <- speciesPoly2_ranT[2,]
#colnames(speciesPoly2_ranT)[1] <- "x"

##tri_cris
stacked_tri <- fasterize(amphSF,plantN2,by="tri_cris")
stacked_tri <- dropLayer(stacked_tri,1)
speciesPoly_tri <- aggregate(amphTrans1, by="tri_cris")
speciesPoly2_tri<-st_as_sf(speciesPoly_tri)
speciesPoly2_tri <- speciesPoly2_tri[2,]
#colnames(speciesPoly2_tri)[1] <- "x"

stackedDistributionMaps_amph2 <- stack(stacked_bufo,stacked_liss,stacked_pelp, stacked_ranA,
                                       stacked_ranT,stacked_tri)
names(stackedDistributionMaps_amph2)<- c('buf_buf','liss_vul','pelp_esc','ran_arv','ran_temp','tri_cris')
#Fit to plant raster (plantN2)
stackedDistributionMaps_amph2[is.na(stackedDistributionMaps_amph2)] <- 0
stackedDistributionMaps_amph2 <- mask(stackedDistributionMaps_amph2, plantN2)

####SR REPTILES####
repTrans1 <-spTransform(repN1,crs(plantN2))
repSF <- st_as_sf(repTrans1)
stackedDistributionMaps_rep<-fasterize(repSF,plantN2,by="rep")
speciesPoly_rep<-aggregate(repTrans1,by="rep")
speciesPoly2_rep<-st_as_sf(speciesPoly_rep)
richnessMap_rep <- fasterize(speciesPoly2_rep,plantN2,field="rep")

####STACKED REPTILES####
##ang_sp
stacked_ang<-fasterize(repSF,plantN2,by="ang_sp")#rasterbrick
stacked_ang <- dropLayer(stacked_ang,1)
speciesPoly_ang <- aggregate(repTrans1, by="ang_sp")
speciesPoly2_ang<-st_as_sf(speciesPoly_ang)
speciesPoly2_ang <- speciesPoly2_ang[2,]
richnessMap_ang <- fasterize(speciesPoly2_ang,plantN2,field="ang_sp")

#cor_aus
stacked_cor<-fasterize(repSF,plantN2,by="cor_aus")#rasterbrick
stacked_cor <- dropLayer(stacked_cor,1)
speciesPoly_cor <- aggregate(repTrans1, by="cor_aus")
speciesPoly2_cor<-st_as_sf(speciesPoly_cor)
speciesPoly2_cor <- speciesPoly2_cor[2,]
richnessMap_cor <- fasterize(speciesPoly2_cor,plantN2,field="cor_aus")

#nat_nat
stacked_nat<-fasterize(repSF,plantN2,by="nat_nat")#rasterbrick
stacked_nat <- dropLayer(stacked_nat,1)
speciesPoly_nat <- aggregate(repTrans1, by="nat_nat")
speciesPoly2_nat<-st_as_sf(speciesPoly_nat)
speciesPoly2_nat <- speciesPoly2_nat[2,]
richnessMap_nat <- fasterize(speciesPoly2_nat,plantN2,field="nat_nat")

#vip_ber
stacked_vip<-fasterize(repSF,plantN2,by="vip_ber")#rasterbrick
stacked_vip <- dropLayer(stacked_vip,1)
speciesPoly_vip <- aggregate(repTrans1, by="vip_ber")
speciesPoly2_vip<-st_as_sf(speciesPoly_vip)
speciesPoly2_vip <- speciesPoly2_vip[2,]
richnessMap_vip <- fasterize(speciesPoly2_vip,plantN2,field="vip_ber")

#zoo_viv
stacked_zoo<-fasterize(repSF,plantN2,by="zoo_viv")#rasterbrick
stacked_zoo <- dropLayer(stacked_zoo,1)
speciesPoly_zoo <- aggregate(repTrans1, by="zoo_viv")
speciesPoly2_zoo<-st_as_sf(speciesPoly_zoo)
speciesPoly2_zoo <- speciesPoly2_zoo[2,]
richnessMap_zoo <- fasterize(speciesPoly2_zoo,plantN2,field="zoo_viv")


stackedDistributionMaps_rep2 <- stack(stacked_ang,stacked_cor,stacked_nat, stacked_vip,
                                      stacked_zoo)
names(stackedDistributionMaps_rep2)<- c('ang_sp','cor_aus','nat_nat','vip_ber','zoo_viv')
stackedDistributionMaps_rep2[is.na(stackedDistributionMaps_rep2)] <- 0
stackedDistributionMaps_rep2 <- mask(stackedDistributionMaps_rep2, plantN2)

####SR BIRDS####
birdTrans1<-spTransform(birdN2,crs(plantN2))
birdSF <- st_as_sf(birdTrans1)
stackedDistributionMaps_bird<-fasterize(birdSF,plantN2,by="SCINAME")
#Fit to plant raster
stackedDistributionMaps_bird[is.na(stackedDistributionMaps_bird)] <- 0
stackedDistributionMaps_bird <- mask(stackedDistributionMaps_bird, plantN2)

speciesPoly_bird<-raster::aggregate(birdTrans1,by="SCINAME")
speciesPoly2_bird<-st_as_sf(speciesPoly_bird)
richnessMap_bird <- fasterize(speciesPoly2_bird,plantN2,fun="count",field="SCINAME")

####SR TOTAL#### 
r1 <- plantN2

#Fit the rest of the layers to the plantraster (by masking)
r2 <- richnessMap_rep
r2[is.na(r2)] <- 0
r2 <- mask(r2, plantN2)

r3 <- richnessMap_amph
r3[is.na(r3)] <- 0
r3 <- mask(r3, plantN2)

r4 <- richnessMap_mam
r4[is.na(r4)] <- 0
r4 <- mask(r4, plantN2)

r44 <-richnessMap_mam
r44[is.na(r44)] <- 0
r44 <- mask(r44, norwayutm)

r5 <- richnessMap_bird
r5[is.na(r5)] <- 0
r5 <- mask(r5, plantN2)

#Make rasterstack with all of the richness layers#
my.stack <- stack(r1,r2,r3,r4,r5)
names(my.stack)<-c('PlantSR','ReptileSR','AmphibianSR','MammalSR','BirdSR')

#Extract values from rasterstack
stackValdf<-raster::extract(my.stack,(1:ncell(my.stack)), df = T)
stackValdf$total <- rowSums(stackValdf[,2:6], na.rm = T, dims = 1)
View(stackValdf)

####LOAD TRAIT DATA####
#Trait data Mammals
traitDataMam <- read.table("MamFuncDat.txt", sep = '\t', header = T, fill = T)
NorMam <- traitDataMam %>%
  filter(Scientific%in% mamN2@data[["BINOMIAL"]]) #50

#Trait data Amphibians
traitDataAmph <- read.csv("AmphiBIO_v1.csv")
NorAmphNames <- c('Bufo bufo','Lissotriton vulgaris','Pelophylax lessonae','Rana arvalis',
                  'Rana temporaria','Triturus cristatus')
NorAmph <- traitDataAmph %>%
  filter(Species%in% NorAmphNames) #6

####MAMMALS ELTON GROUPS####
MamHerbi <- NorMam[NorMam$Diet.PlantO>=70,] #11
nameMamHerbi<- MamHerbi$Scientific

MamGrani <- NorMam[NorMam$Diet.Seed>=70,] #0 (if >= 50, Sciurus vulgaris is considered granivore)
nameMamGrani <- MamGrani$Scientific

MamCarni <- rbind((NorMam[NorMam$Diet.Vend >=50,]),
                  (NorMam[NorMam$Diet.Vect >=50,]),
                  (NorMam[NorMam$Diet.Vfish>=50,]),
                  (NorMam[NorMam$Diet.Vunk>=50]),
                  (NorMam[NorMam$Diet.Scav>=50]))#8
nameMamCarni<- MamCarni$Scientific

MamInverti <- NorMam[NorMam$Diet.Inv>=50,] #16
nameMamInverti <- MamInverti$Scientific

MamOmni <- NorMam[NorMam$Diet.Inv<50&
                    NorMam$Diet.Vend<50&
                    NorMam$Diet.Vect<50&
                    NorMam$Diet.Vfish<50&
                    NorMam$Diet.Vunk<50&
                    NorMam$Diet.Scav<50&
                    NorMam$Diet.Seed<70&
                    NorMam$Diet.PlantO<70,] #14
nameMamOmnivore <- MamOmni$Scientific

checkMam <-  NorMam%>%
  filter((!(Scientific %in% MamCarni$Scientific))&
           !(Scientific %in% MamHerbi$Scientific)&
           !(Scientific %in% MamGrani$Scientific)&
           !(Scientific %in% MamInverti$Scientific)&
           !(Scientific %in% MamOmni$Scientific)) #0

####BIRDS ELTON GROUPS####
BirdHerbi <- NorBirdTerr[NorBirdTerr$Diet.PlantO>=70,] #17
nameBirdHerbi<- BirdHerbi$Scientific

BirdGrani <- NorBirdTerr[NorBirdTerr$Diet.Seed>=70,] #4
nameBirdGrani <- BirdGrani$Scientific

BirdCarni <- rbind((NorBirdTerr[NorBirdTerr$Diet.Vend >=50,]),
                   (NorBirdTerr[NorBirdTerr$Diet.Vect >=50,]),
                   (NorBirdTerr[NorBirdTerr$Diet.Vfish>=50,]),
                   (NorBirdTerr[NorBirdTerr$Diet.Vunk>=50]),
                   (NorBirdTerr[NorBirdTerr$Diet.Scav>=50])) #36
nameBirdCarni<- BirdCarni$Scientific

BirdInverti <- NorBirdTerr[NorBirdTerr$Diet.Inv>=50,] #113
BirdInverti <- BirdInverti[BirdInverti$Scientific != "Podiceps grisegena",] #placed in both carni and here
nameBirdInverti <- BirdInverti$Scientific

BirdOmni <- NorBirdTerr[NorBirdTerr$Diet.Inv<50&
                          NorBirdTerr$Diet.Vend<50&
                          NorBirdTerr$Diet.Vect<50&
                          NorBirdTerr$Diet.Vfish<50&
                          NorBirdTerr$Diet.Vunk<50&
                          NorBirdTerr$Diet.Scav<50&
                          NorBirdTerr$Diet.Seed<70&
                          NorBirdTerr$Diet.PlantO<70,] #55
nameBirdOmnivore <- BirdOmni$Scientific

##Total of 225 species? 
checkBird <-  NorBirdTerr%>%
  filter((!(Scientific %in% BirdCarni$Scientific))&
           !(Scientific %in% BirdHerbi$Scientific)&
           !(Scientific %in% BirdGrani$Scientific)&
           !(Scientific %in% BirdInverti$Scientific)&
           !(Scientific %in% BirdOmni$Scientific))#0

####AMPHIBIANS DIET GROUPS####
AmphInverti <- NorAmph[NorAmph$Arthro == 1,] #6
#AmphOmni <- NorAmph[NorAmph$Vert == 1 & NorAmph$Arthro == 1,] #2
#NorAmph[is.na(NorAmph$Vert)& NorAmph$Arthro == 1,]
nameAmphiInverti <- AmphInverti$Species
AmphInvertSR <- richnessMap_amph

####REPTILES DIET GROUPS####
#Ang = invertibrate, Zoo = invertibrate,
#cor = carnivore, nat = carnivore, vip = carnivore, 
RepInvertSR <- mosaic(richnessMap_ang,richnessMap_zoo,fun="sum") #2
RepCarniSR <- mosaic(richnessMap_cor, richnessMap_nat, richnessMap_vip,fun="sum")#4
nameReptInverti <- as.factor(c("Anguis fragilis","Zootoca vivipara"))
nameReptCarni <- as.factor(c("Coronella austriaca","Natrix natrix","Vipera berus"))

####GROUPS MAMMALS/FEEDING ####
#SR herbivore mammals
herbivoreMam <- subset(mamN2, BINOMIAL %in% nameMamHerbi)
herbivoreMamTrans<-spTransform(herbivoreMam,crs(plantN2))
speciesPoly_HM <- aggregate(herbivoreMamTrans, by ="BINOMIAL")
speciesPoly_HM2 <- st_as_sf(speciesPoly_HM)
mamHerbSR <- fasterize(speciesPoly_HM2,plantN2,fun="count",field="BINOMIAL") #12

#SR carnivore mammals
carnivoreMam <- subset(mamN2, BINOMIAL %in% nameMamCarni)
carnivoreMamTrans<-spTransform(carnivoreMam,crs(plantN2))
speciesPoly_CM <- aggregate(carnivoreMamTrans, by ="BINOMIAL")
speciesPoly_CM2 <- st_as_sf(speciesPoly_CM)
mamCarnSR <- fasterize(speciesPoly_CM2,plantN2,fun="count",field="BINOMIAL") #8

#SR insectivores mammals
invertivoreMam <- subset(mamN2, BINOMIAL %in% nameMamInverti)
invertivoreMamTrans<-spTransform(invertivoreMam,crs(plantN2))
speciesPoly_IM <- aggregate(invertivoreMamTrans, by ="BINOMIAL")
speciesPoly_IM2 <- st_as_sf(speciesPoly_IM)
mamInvertSR <- fasterize(speciesPoly_IM2,plantN2,fun="count",field="BINOMIAL") #18

#SR granivore mammals
granivoreMam <- subset(mamTrans1, BINOMIAL %in% nameMamGrani) #0

#SR omnivore mammals
omnivoreMam <- subset(mamN2, BINOMIAL %in% nameMamOmnivore)
omnivoreMamTrans<-spTransform(omnivoreMam,crs(plantN2))
e_OM <-extent(omnivoreMamTrans)
s_OM<-raster(e_OM, resolution=20000, crs=(omnivoreMamTrans))
speciesPoly_OM <- aggregate(omnivoreMamTrans, by ="BINOMIAL")
speciesPoly_OM2 <- st_as_sf(speciesPoly_OM)
mamOmniSR <- fasterize(speciesPoly_OM2,plantN2,fun="count",field="BINOMIAL") #15

####GROUPS BIRD/FEEDING####
#SR herbivore birds
herbivoreBird <- subset(birdN2, SCINAME %in% nameBirdHerbi)
herbivoreBirdTrans<-spTransform(herbivoreBird,crs(plantN2))
speciesPoly_HB <- aggregate(herbivoreBirdTrans, by ="SCINAME")
speciesPoly_HB2 <- st_as_sf(speciesPoly_HB)
birdHerbSR <- fasterize(speciesPoly_HB2,plantN2,fun="count",field="SCINAME")

#SR carnivore birds
carnivoreBird <- subset(birdN2, SCINAME %in% nameBirdCarni)
carnivoreBirdTrans<-spTransform(carnivoreBird,crs(plantN2))
speciesPoly_CB <- aggregate(carnivoreBirdTrans, by ="SCINAME")
speciesPoly_CB2 <- st_as_sf(speciesPoly_CB)
birdCarnSR <- fasterize(speciesPoly_CB2,plantN2,fun="count",field="SCINAME")

#SR granivore birds
granivoreBird <- subset(birdN2, SCINAME %in% nameBirdGrani)
granivoreBirdTrans <- spTransform(granivoreBird, crs(plantN2))
speciesPoly_GB <- aggregate(granivoreBirdTrans, by ="SCINAME")
speciesPoly_GB2 <- st_as_sf(speciesPoly_GB)
birdGranSR <- fasterize(speciesPoly_GB2,plantN2,fun="count",field="SCINAME")#4

#SR insectivore birds
invertivoreBird <- subset(birdN2, SCINAME %in% nameBirdInverti) 
invertivoreBirdTrans <- spTransform(invertivoreBird, crs(plantN2))
speciesPoly_IB <- aggregate(invertivoreBirdTrans, by ="SCINAME")
speciesPoly_IB2 <- st_as_sf(speciesPoly_IB)
birdInvertSR <- fasterize(speciesPoly_IB2,plantN2,fun="count",field="SCINAME")

#SR omnivore birds
omnivoreBird <- subset(birdN2, SCINAME %in% nameBirdOmnivore) 
omnivoreBirdTrans <- spTransform(omnivoreBird, crs(plantN2))
speciesPoly_OB <- aggregate(omnivoreBirdTrans, by ="SCINAME")
speciesPoly_OB2 <- st_as_sf(speciesPoly_OB)
birdOmnivoreSR <- fasterize(speciesPoly_OB2,plantN2,fun="count",field="SCINAME")


####ALLSTACKVALDF: EXTRACT AND STACK####
#Extract species from each cell from rasterlayer
#Mammals
mammalPoly <- raster::extract(stackedDistributionMaps_mam,(1:ncell(stackedDistributionMaps_mam)), df = T)
mammalPoly2 <- mammalPoly
mammalPoly$total <- rowSums(mammalPoly[,2:ncol(mammalPoly)], na.rm = T, dims = 1)

#Birds
birdPoly <- raster::extract(stackedDistributionMaps_bird,(1:ncell(stackedDistributionMaps_bird)), df = T)
birdPoly2 <- birdPoly[,-1]
birdPoly$total <- rowSums(birdPoly[,2:ncol(birdPoly)], na.rm = T, dims = 1)

###Amphibian
amphPoly <- raster::extract(stackedDistributionMaps_amph2,(1:ncell(stackedDistributionMaps_amph2)), df = T)
amphPoly2 <- amphPoly[,-1]
names(amphPoly2)<- c("Bufo.bufo","Lissotriton.vulgaris","Pelophylax.lessonae",
                     "Rana.arvalis","Rana.temporaria","Triturus.cristatus")
amphPoly$total <- rowSums(amphPoly[,2:ncol(amphPoly)], na.rm = T, dims = 1)

#Reptiles
reptPoly <- raster::extract(stackedDistributionMaps_rep2,(1:ncell(stackedDistributionMaps_rep2)), df = T)
reptPoly2 <- reptPoly[,-1]
names(reptPoly2) <- c("Anguis.fragilis","Coronella.austriaca","Natrix.natrix",
                      "Vipera.berus","Zootoca.vivipara")
reptPoly$total <- rowSums(reptPoly[,2:ncol(reptPoly)], na.rm = T, dims = 1)

##STACK ALL DIFFERENT TROPHIC GROUPS TOGETHER
allStack <- stack(RepInvertSR,RepCarniSR, AmphInvertSR,plantN2,birdHerbSR, birdInvertSR, birdGranSR, birdCarnSR, birdOmnivoreSR,
                  mamCarnSR,mamHerbSR,mamInvertSR,mamOmniSR)
#Crop to extent of Norway outline
names(allStack)<- c('repInvertSR','repCarniSR','amphInvertSR','plantSR','birdHerbSR','birdInvertSR','birdGranSR','birdCarnSR','birdOmnivoreSR',
                    'mamCarnSR','mamHerbSR','mamInvertSR','mamOmniSR')
allStack[is.na(allStack)] <- 0
allStack <- mask(allStack, plantN2)

#Extract
allStackValdf <- raster::extract(allStack,(1:ncell(allStack)), df = T)
allStackValdf$total <- rowSums(allStackValdf[,2:ncol(allStackValdf)], na.rm = T, dims = 1)

####MAKE AND EXPORT COMMUNITY MATRIX####
allPoly <- bind_cols(mammalPoly2,birdPoly2,amphPoly2,reptPoly2)
#write.table(allPoly, file = "CommunityMat_2.txt", sep = "\t",row.names = F)

####TABLE OF SPECIES: WHICH SPECIES ARE WITHIN EACH RASTERCELL(ex. plants)####
#Melt allPoly df so that you get which rastercells each species is within
allPoly_molten <- melt(allPoly,id ="ID")
allPoly_molten <- subset(allPoly_molten, value==1)
names(allPoly_molten) <-c("ID", "Scientific", "value")
levels(allPoly_molten$Scientific) <- gsub("\\."," ",levels(allPoly_molten$Scientific))

#CREATE LIST WITH SPECIES NAMES
mammal_latin <- as.character(NorMam$Scientific)
bird_latin <- as.character(NorBirdTerr$Scientific)
amph_latin <- as.character(NorAmph$Species)
rept_latin <- c("Anguis fragilis","Coronella austriaca","Natrix natrix",
                "Vipera berus","Zootoca vivipara")
Scientific <- c(mammal_latin,bird_latin,amph_latin,rept_latin) #285!

#MAKE DF W/SPECIES NAME, TAXONOMIC GROUP AND TROPHIC GROUP
table <- data.frame(Scientific)
#Make new column which categorize according to taxonomic group (mammal, bird, amph, rept)
table$taxonomic <- c(1:(length(table$Scientific)))
table$taxonomic[1:49] <- "mammal"
table$taxonomic[50:274] <- "bird"
table$taxonomic[275:280] <- "amphibian"
table$taxonomic[281:285] <- "reptile"
table$taxonomic <- as.factor(table$taxonomic)

#Make new column that categorize according to trophic group:
table$trophic <- if_else(table$Scientific%in%nameBirdInverti |
                           table$Scientific %in% nameMamInverti |
                           table$Scientific %in% nameAmphiInverti |
                           table$Scientific %in% nameReptInverti, "invertivore",
                         if_else(table$Scientific%in%nameMamCarni |
                                   table$Scientific%in%nameBirdCarni |
                                   table$Scientific%in%nameReptCarni,"carnivore",
                                 if_else(table$Scientific%in%nameBirdGrani |
                                           table$Scientific%in%nameMamGrani,"granivore",
                                         if_else(table$Scientific%in%nameMamHerbi |
                                                   table$Scientific%in%nameBirdHerbi, "herbivore",
                                                 if_else(table$Scientific%in%nameBirdOmnivore |
                                                           table$Scientific%in%nameMamOmnivore,"omnivore","NA")))))
table$trophic <- as.factor(table$trophic)

#MERGE
allPoly_merged<- merge(allPoly_molten,table, by="Scientific", all.x = T)
allPoly_merged <- allPoly_merged[,c(2,1,4,5,3)]
tableOfSpecies <- allPoly_merged[,-c(5)]

#Save and export table
#write.csv(tableOfSpecies, file="Species_table.csv")

#####SR MAPS TROPHIC GROPUS####
###SR MAPS OF HERBIVORES
#change name of birds from "SCINAME" to "BINOMIAL", in order to match mammals
speciesPoly_HB2_cop <- speciesPoly_HB2
names(speciesPoly_HB2_cop) <- c("BINOMIAL","geometry")
stackHerbi <- rbind(speciesPoly_HM2, speciesPoly_HB2_cop) #28 species
herbivoreSR <- fasterize(stackHerbi,plantN2,fun="count",field="BINOMIAL")

###SR MAPS OF CARNIVORES
#change name and class of 1 column
speciesPoly_CB2_cop <- speciesPoly_CB2
names(speciesPoly_CB2_cop) <- c("BINOMIAL","geometry")

speciesPoly2_cor_cop <- speciesPoly2_cor
speciesPoly2_cor_cop[[1]] <- "Coronella austriaca"
names(speciesPoly2_cor_cop) <- c("BINOMIAL","geometry")
speciesPoly2_cor_cop$BINOMIAL <- as.factor(speciesPoly2_cor_cop$BINOMIAL)

speciesPoly2_nat_cop <- speciesPoly2_nat
speciesPoly2_nat_cop[[1]] <- "Natrix natrix"
names(speciesPoly2_nat_cop) <- c("BINOMIAL","geometry")
speciesPoly2_nat_cop$BINOMIAL <- as.factor(speciesPoly2_nat_cop$BINOMIAL)

speciesPoly2_vip_cop <- speciesPoly2_vip
speciesPoly2_vip_cop[[1]] <- "Vipera berus"
names(speciesPoly2_vip_cop) <- c("BINOMIAL","geometry")
speciesPoly2_vip_cop$BINOMIAL <- as.factor(speciesPoly2_vip_cop$BINOMIAL)

stackCarni <- rbind(speciesPoly_CB2_cop,speciesPoly_CM2,speciesPoly2_cor_cop,#47 species
                    speciesPoly2_nat_cop,speciesPoly2_vip_cop) 
carnivoreSR <- fasterize(stackCarni,plantN2,fun="count",field="BINOMIAL")

###SR MAPS INVERTIVORES
#reptiles = zoo and ang
##ANGUIS FRAGILIS
speciesPoly2_ang_cop <- speciesPoly2_ang
speciesPoly2_ang_cop[[1]] <- "Anguis fragilis"
names(speciesPoly2_ang_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ang_cop$BINOMIAL <- as.factor(speciesPoly2_ang_cop$BINOMIAL)
##ZOOTOCA VIVIPARA
speciesPoly2_zoo_cop <- speciesPoly2_zoo
speciesPoly2_zoo_cop[[1]] <- "Zootoca vivipara"
names(speciesPoly2_zoo_cop) <- c("BINOMIAL","geometry")
speciesPoly2_zoo_cop$BINOMIAL <- as.factor(speciesPoly2_zoo_cop$BINOMIAL)
#amphibians = all species
##BUFO BUFO
speciesPoly2_buf_cop <- speciesPoly2_buf
speciesPoly2_buf_cop[[1]] <- "Bufo bufo"
names(speciesPoly2_buf_cop) <- c("BINOMIAL","geometry")
speciesPoly2_buf_cop$BINOMIAL <- as.factor(speciesPoly2_buf_cop$BINOMIAL)
##LISSOTRITON VULGARIS
speciesPoly2_liss_cop <- speciesPoly2_liss
speciesPoly2_liss_cop[[1]] <- "Lissotriton vulgaris"
names(speciesPoly2_liss_cop) <- c("BINOMIAL","geometry")
speciesPoly2_liss_cop$BINOMIAL <- as.factor(speciesPoly2_liss_cop$BINOMIAL)
##PELOPHYLAX LESSONAE 
speciesPoly2_pelp_cop <- speciesPoly2_pelp
speciesPoly2_pelp_cop[[1]] <- "Pelophylax lessonae"
names(speciesPoly2_pelp_cop) <- c("BINOMIAL","geometry")
speciesPoly2_pelp_cop$BINOMIAL <- as.factor(speciesPoly2_pelp_cop$BINOMIAL)
###RANA ARVALIS 
speciesPoly2_ranA_cop <- speciesPoly2_ranA
speciesPoly2_ranA_cop[[1]] <- "Rana arvalis"
names(speciesPoly2_ranA_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ranA_cop$BINOMIAL <- as.factor(speciesPoly2_ranA_cop$BINOMIAL)
###RANA TEMPORARIA 
speciesPoly2_ranT_cop <- speciesPoly2_ranT
speciesPoly2_ranT_cop[[1]] <- "Rana temporaria"
names(speciesPoly2_ranT_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ranT_cop$BINOMIAL <- as.factor(speciesPoly2_ranT_cop$BINOMIAL)
###TRITURUS CRISTATUS 
speciesPoly2_tri_cop <- speciesPoly2_tri
speciesPoly2_tri_cop[[1]] <- "Triturus cristatus"
names(speciesPoly2_tri_cop) <- c("BINOMIAL","geometry")
speciesPoly2_tri_cop$BINOMIAL <- as.factor(speciesPoly2_tri_cop$BINOMIAL)
#birds
speciesPoly_IB2_cop <- speciesPoly_IB2
names(speciesPoly_IB2_cop) <- c("BINOMIAL","geometry")

stackInverti<- rbind(speciesPoly2_ang_cop,speciesPoly2_zoo_cop,speciesPoly2_buf_cop, #137 species
                     speciesPoly2_liss_cop,speciesPoly2_pelp_cop,speciesPoly2_ranA_cop,
                     speciesPoly2_ranT_cop,speciesPoly2_tri_cop,speciesPoly_IB2_cop,
                     speciesPoly_IM2) 
invertivoreSR <- fasterize(stackInverti,plantN2,fun="count",field="BINOMIAL")

###SR Omnivores
speciesPoly_OB2_cop <- speciesPoly_OB2
names(speciesPoly_OB2_cop) <- c("BINOMIAL","geometry")
stackOmni <- rbind(speciesPoly_OB2_cop,speciesPoly_OM2)#69 species
omnivoreSR <- fasterize(stackOmni,plantN2,fun="count",field="BINOMIAL") 

###SR Granivore
granivoreSR <- fasterize(speciesPoly_GB2,plantN2,fun="count",field="SCINAME") #4 species 

####RASTERSTACK ALL ANIMALS####
###SR Animals; reptiles, amphibians, mammals, birds
#Copies
#Change column names to match the "BINOMIAL" of mammal layer, in order to combine together 
sfDf_rept <- speciesPoly2_rep
names(sfDf_rept) <- c("BINOMIAL","geometry")
sfDf_rept$BINOMIAL <- as.factor(sfDf_rept$BINOMIAL)
sfDf_rept <- sfDf_rept[-1,]

sfDf_amph <- speciesPoly2_amph
names(sfDf_amph) <- c("BINOMIAL","geometry")
sfDf_amph$BINOMIAL <- as.factor(sfDf_amph$BINOMIAL)
sfDf_amph <- sfDf_amph[-1,]

sfDf_mam <- speciesPoly2_mam
names(sfDf_amph) <- c("BINOMIAL","geometry")

sfDf_bird <- speciesPoly2_bird
names(sfDf_bird) <- c("BINOMIAL","geometry")
#STACK  
stackAni <- rbind(sfDf_rept,sfDf_amph,sfDf_mam,sfDf_bird)


####MAKE RICHNESS MAPS####
#COLOR PALETTE
nb.cols<- 16
mycolors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(nb.cols)
#Adjusting the legend on the sppplot
brks = seq(0,1,0.1)

###SR ANIMALS###
animalSR <- fasterize(stackAni,plantN2,fun="count",field="BINOMIAL") 
animalSR[is.na(animalSR)] <- 0
animalSR <- mask(animalSR, plantN2)

#normalize values against total SR for animals (=285)
rel_animalSR <- animalSR
values(rel_animalSR) <- (values(rel_animalSR))/285
#Plot
MAP_AnimalSR<- spplot(rel_animalSR,
                      colorkey = FALSE,
                      at = round(brks,digits = 2),
                      col.regions = mycolors,
                      main = list(label = "Mammals, reptiles, amphibians & birds", cex = 1),
                      par.settings = list(axis.line = list(col = 'transparent')),
                      sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                      contour =F)
###SR PLANTS###
rel_plantSR <- plantN2
values(rel_plantSR) <- (values(rel_plantSR)/1238) #1238, right?

MAP_PlantSR<- spplot(rel_plantSR,
                     #colorkey = FALSE,
                     at = round(brks,digits = 2),
                     col.regions = mycolors,
                     main = list(label = "Species richness", cex = 1),
                     par.settings = list(axis.line = list(col = 'transparent')),
                     sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                     contour =F)
#Change name of main to make headline for stack of all SR maps

###PLOTTING DIFFERENT TROPHIC GROUPS###
#GRANIVORES
rel_granivore <- granivoreSR
#rel_granivore[is.na(rel_granivore)] <- 0
#rel_granivore <- mask(rel_granivore, plantN2)
values(rel_granivore) <- values(rel_granivore)/4

MAP_GrSR<- spplot(rel_granivore,
                  colorkey = FALSE,
                  at = round(brks,digits = 2),
                  col.regions = mycolors,
                  main = list(label = "Granivores", cex = 1),
                  par.settings = list(axis.line = list(col = 'transparent')),
                  sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                  contour =F)

###CARNIVORES###
rel_carnivore <- carnivoreSR
rel_carnivore[is.na(rel_carnivore)] <- 0
rel_carnivore <- mask(rel_carnivore, plantN2)
values(rel_carnivore) <- values(rel_carnivore)/47

MAP_CarSR<- spplot(rel_carnivore,
                   colorkey = FALSE,
                   at = round(brks,digits = 2),
                   col.regions = mycolors,
                   main = list(label = "Carnivores", cex = 1),
                   par.settings = list(axis.line = list(col = 'transparent')),
                   sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                   contour =F)

###INVERTIVORES##
rel_invertivore <- invertivoreSR
rel_invertivore[is.na(rel_invertivore)] <- 0
rel_invertivore <- mask(rel_invertivore, plantN2)
values(rel_invertivore) <- values(rel_invertivore)/137

MAP_InvSR<- spplot(rel_invertivore,
                   colorkey = FALSE,
                   at = round(brks,digits = 2),
                   col.regions = mycolors,
                   main = list(label = "Invertivores", cex = 1),
                   par.settings = list(axis.line = list(col = 'transparent')),
                   sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                   contour =F)

###OMNIVORES###
rel_omnivore <- omnivoreSR
rel_omnivore[is.na(rel_omnivore)] <- 0
rel_omnivore <- mask(rel_omnivore, plantN2)
values(rel_omnivore) <- values(rel_omnivore)/69

MAP_OmnSR<- spplot(rel_omnivore,
                   colorkey = FALSE,
                   at = round(brks,digits = 2),
                   col.regions = mycolors,
                   main = list(label = "Omnivores", cex = 1),
                   par.settings = list(axis.line = list(col = 'transparent')),
                   sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                   contour =F)

###HERBIVORES###
rel_herbivore <- herbivoreSR
rel_herbivore[is.na(rel_herbivore)] <- 0
rel_herbivore <- mask(rel_herbivore, plantN2)
values(rel_herbivore) <- values(rel_herbivore)/28

MAP_HeSR<- spplot(rel_herbivore,
                  colorkey = FALSE,
                  at = round(brks,digits = 2),
                  col.regions = mycolors,
                  main = list(label = "Herbivores", cex = 1),
                  par.settings = list(axis.line = list(col = 'transparent')),
                  sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                  contour =F)

####PLOTTING SR THROUGH DIFFERENT TROPHIC GROUPS### 
my.SR.rasterstack <- stack(rel_plantSR,rel_granivore,rel_carnivore,rel_omnivore,rel_herbivore,
                           rel_invertivore)
names(my.SR.rasterstack) <- c('Plants','Granivores','Carnivores','Omnivores','Herbivores',
                              'Invertivores')

#SR_visu <- rasterVis::levelplot(my.SR.rasterstack,contour =F,frame = F, 
#col.regions = mycolors,scales=list(draw=FALSE ))
#tmap_arrange(plantSR_map,granivoreSR_map,carnivoreSR_map,invertivoreSR_map,omnivoreSR_map,herbivoreSR_map)

#Plants, Herbivores, Granivores, Carnivores, invertivores, omnivores
p4 <- c(MAP_PlantSR,MAP_HeSR,MAP_GrSR,MAP_CarSR, MAP_InvSR,MAP_OmnSR)


####PLOT PHYLOGENETIC DIVERSITY####
#COLOR PALETTE
nb.cols2<- 16
my.colors.green <- colorRampPalette(brewer.pal(9, "YlGn"))(nb.cols)

#ANIMALS (reptiles, amphibians, birds,mammals)
animalPD <- animalSR
values(animalPD) <- (phyloDiv_root$PD)/animal_totPd
animalPD <- mask(animalPD, plantN2)
animalPD_map <- tm_shape(animalPD)+
  tm_raster(style="cont",palette = mycolors)+
  map_nor

animalPD_map

#PLANTS
plantPD <- plantN2 
values(plantPD) <- (pd.plant.root$PD)/plant.totPd
plantPD <- mask(plantPD, plantN2)

MAP_plantPD<- spplot(plantPD,
                     colorkey = FALSE,
                     at = round(brks,digits = 2),
                     col.regions = my.colors.green,
                     main = list(label = "Phylogenetic diversity", cex = 1),
                     par.settings = list(axis.line = list(col = 'transparent')),
                     sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                     contour =F)
#Main = phylogenetic diversity, to get a main when plotting all the maps together  

#####CREATE COMMUNITY MATRICES FOR HERB,CARN,GRANI,OMNI,INVERTI####
###HERBIVORE###
#You can take the community matrix, and just crop it to the names in each group!
herbNames <- as.character(stackHerbi$BINOMIAL)
herbNames <- gsub("\\s","_",herbNames)

#setdiff(herbNames,names(commat_herbivore))
#Missing 4 species from the subsetting: 
herbNames["Anas_penelope"] <- "Mareca_penelope"
herbNames["Anas_strepera"] <- "Mareca_strepera"
herbNames["Bonasa_bonasia"] <- "Tetrastes_bonasia"
herbNames["Tetrao_tetrix"] <- "Lyrurus_tetrix"

#subset
commNames <- names(my.sample) %in% herbNames
commat_herbivore <- my.sample[commNames]
write.csv(commat_herbivore, file = "commat_herbivore.csv")

#Caluclate PD for herbivores
PD_herbivore <- pd(commat_herbivore,phylogeny,include.root = T)

herbivorePD <- herbivoreSR
values(herbivorePD) <- (PD_herbivore$PD)/(herb_totPd2)
herbivorePD <- mask(herbivorePD, plantN2)
#herbivorePD[herbivorePD == 0] <- NA

MAP_HePD<- spplot(herbivorePD,
                  colorkey = FALSE,
                  at = round(brks,digits = 2),
                  col.regions = my.colors.green,
                  main = list(label = "Herbivores", cex = 1),
                  par.settings = list(axis.line = list(col = 'transparent')),
                  sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                  contour =F)
###CARNIVORE###
carnNames <- as.character(stackCarni$BINOMIAL)
carnNames <- gsub("\\s","_",carnNames)
#Fix names
carnNames["Bubo_scandiaca"] <- "Bubo_scandiacus"
#subset
commNames.carn <- names(my.sample) %in% carnNames
commat_carnivore <- my.sample[commNames.carn]
write.csv(commat_carnivore, file = "commat_carnivore.csv")

#setdiff(names(commat_carnivore),carnNames) Missing 4 species from the subsetting: 
#Caluclate PD for herbivores
PD_carnivore <- pd(commat_carnivore,phylogeny,include.root = T)

carnivorePD <- carnivoreSR
values(carnivorePD) <- (PD_carnivore$PD)/(carn_totPd2)
carnivorePD <- mask(carnivorePD, plantN2)

MAP_CaPD<- spplot(carnivorePD,
                  colorkey = FALSE,
                  at = round(brks,digits = 2),
                  col.regions = my.colors.green,
                  main = list(label = "Carnivores", cex = 1),
                  par.settings = list(axis.line = list(col = 'transparent')),
                  sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                  contour =F)

###OMNIVORE###
omniNames <- as.character(stackOmni$BINOMIAL)
omniNames <- gsub("\\s","_",omniNames)

omniNames["Carduelis_flammea"] <- "Acanthis_flammea"
omniNames["Carduelis_spinus"] <- "Spinus_spinus"
omniNames["Parus_ater"] <- "Periparus_ater"

commNames.omn <- names(my.sample) %in% omniNames
commat_omnivore <- my.sample[commNames.omn]
write.csv(commat_omnivore, file = "commat_omnivore.csv")

PD_omnivore <- pd(commat_omnivore,phylogeny,include.root = T)

omnivorePD <- omnivoreSR
values(omnivorePD) <- (PD_omnivore$PD)/(omni_totPd)
omnivorePD <- mask(omnivorePD, plantN2)

MAP_OmnPD<- spplot(omnivorePD,
                   colorkey = FALSE,
                   at = round(brks,digits = 2),
                   col.regions = my.colors.green,
                   main = list(label = "Omnivores", cex = 1),
                   par.settings = list(axis.line = list(col = 'transparent')),
                   sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                   contour =F)

###INVERTIVORE###
invertNames <- as.character(stackInverti$BINOMIAL)
invertNames <- gsub("\\s","_",invertNames)

#Larus_minutus

invertNames["Eudromias_morinellus"] <- "Charadrius_morinellus" #Charadrius_morinellus
invertNames["Larus_minutus"] <- "Hydrocoloeus_minutus" #Hydrocoloeus_minutus
invertNames["Larus_ridibundus"] <- "Chroicocephalus_ridibundus"
invertNames["Parus_caeruleus"] <- "Cyanistes_caeruleus"
invertNames["Parus_cinctus"] <- "Poecile_cinctus"
invertNames["Parus_cristatus"] <- "Lophophanes_cristatus"
invertNames["Parus_montanus"] <- "Poecile_montanus"
invertNames["Parus_palustris"] <- "Poecile_palustris"
invertNames["Philomachus_pugnax"] <- "Calidris_pugnax"
invertNames["Phylloscopus_borealis"] <- "Seicercus_borealis"
invertNames["Saxicola_torquatus"] <- "Saxicola_torquata"

#subset
commNames.inv <- names(my.sample) %in% invertNames
commat_invertivore <- my.sample[commNames.inv]
write.csv(commat_invertivore, file = "commat_invertivore.csv")

#PLOT
PD_invertivore <- pd(commat_invertivore,phylogeny,include.root = T)
invertivorePD <- invertivoreSR
values(invertivorePD) <- (PD_invertivore$PD)/(invert_totPd)
invertivorePD <- mask(invertivorePD, plantN2)

MAP_InvPD<- spplot(invertivorePD,
                   colorkey = FALSE,
                   at = round(brks,digits = 2),
                   col.regions = my.colors.green,
                   main = list(label = "Invertivores", cex = 1),
                   par.settings = list(axis.line = list(col = 'transparent')),
                   sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                   contour =F)

###GRANIVORE###
#There are only 4 granivorous birds 
graniNames <- as.character(BirdGrani$Scientific)
graniNames <- gsub("\\s","_",graniNames)
graniNames["Emberiza_rustica"] <- "Schoeniclus_rusticus"
#Subset
commNames.grani <- names(my.sample) %in% graniNames
commat_granivore <- my.sample[commNames.grani]
write.csv(commat_granivore, file = "commat_granivore.csv")

#PLOT
PD_granivore <- pd(commat_granivore,phylogeny,include.root = T)

granivorePD <- granivoreSR
values(granivorePD) <- (PD_granivore$PD)/(grani_totPd)
granivorePD <- mask(granivorePD, plantN2)

MAP_GrPD<- spplot(granivorePD,
                  colorkey = FALSE,
                  at = round(brks,digits = 2),
                  col.regions = my.colors.green,
                  main = list(label = "Invertivores", cex = 1),
                  par.settings = list(axis.line = list(col = 'transparent')),
                  sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                  contour =F)

####RASTERSTACK ALL PD'S####  
my.rasterstack <- stack(plantPD,herbivorePD,granivorePD,omnivorePD,
                        invertivorePD, carnivorePD)
names(my.rasterstack)<- c('Plant','Herbivores','Granivores','Omnivores',
                          'Invertivores','Carnivores')  

p5 <- c(MAP_plantPD,MAP_HePD ,MAP_GrPD,MAP_CaPD, MAP_InvPD,MAP_OmnPD,
        layout=c(3,2))

####STATISTICS#####
###CORRLEATION BETWEeN PLANTS PD AND HERBIVORE PD
#All values is proportion of PD for the whole trohpic group

#Extract values from rasterstack (my.rasterstack)
rasterVal <- getValues(my.rasterstack)
rasterValdf <- as.data.frame(rasterVal)
rasterValdf <- rasterValdf[,c(1,2,3,6,5,4)]

rasterValdf <- na.omit(rasterValdf) 

###SCATTER PLOT MATRIX####
# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("r = ", r, "***")
  cex.cor <- 2
  text(0.5, 0.5, txt, cex = 2)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 1)
}
# Create the plots
pairs(rasterValdf[,c(1:2,4:6)], #exclude granivores
      lower.panel = upper.panel,
      upper.panel = panel.cor,
      cex.labels = 2)

####SCATTERPLOT CORRELATIONS####
###PLANT/HERBIVORES: Pearson's correlation: 0.35
plant.herb <- ggscatter(rasterValdf,x="Plant",y="Herbivores",xlab = "\n Plants",ylab="Herbivores \n",
                        shape = 1,
                        alpha = 0.5,
                        size = 2,
                        conf.int = T,
                        add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson",label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###PLANT/GRANIVORES: Pearson's correlation:0.44
plant.grani <- ggscatter(rasterValdf,x="Plant",y="Granivores",xlab = "\n Plants",ylab="Granivores \n",
                         shape = 1,
                         alpha = 0.5,
                         size = 2,
                         conf.int = T,
                         add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###HERBI/CARNIVORES: Pearson's correlation:0.92
herbi.carni <- ggscatter(rasterValdf,x="Herbivores",y="Carnivores",xlab = "\n Herbivores",ylab="Carnivores \n",
                         shape = 1,
                         alpha = 0.5,
                         size = 2,
                         conf.int = T,
                         add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###HERBI/INVERTIVORES: Pearson's correlation: 0.9
herbi.inverti <- ggscatter(rasterValdf,x="Herbivores",y="Invertivores",xlab = "\n Herbivores",ylab="Invertivores \n",
                           shape = 1,
                           alpha = 0.5,
                           size = 2,
                           conf.int = T,
                           add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###PLANT/INVERTIVORES: Pearson's correlation: 0.52
plant.inverti <- ggscatter(rasterValdf,x="Plant",y="Invertivores",xlab = "\n Plants",ylab="Invertivores \n",
                           shape = 1,
                           alpha = 0.5,
                           size = 2,
                           conf.int = T,
                           add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###PLANTS/OMNIVORES: Pearson's correaltion: 0.43
omni.plant <-ggscatter(rasterValdf,x="Plant",y="Omnivores",xlab = "\n Plants",ylab="Omnivores \n",
                       shape = 1,
                       alpha = 0.5,
                       size = 2,
                       conf.int = T,
                       add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###OMNIVORES/HERBIVORES: Pearson's correlation 0.95
omni.herbivores <- ggscatter(rasterValdf,x="Omnivores",y="Herbivores",xlab = "\n Omnivores",ylab="Herbivores \n",
                             shape = 1,
                             alpha = 0.5,
                             size = 2,
                             conf.int = T,
                             add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)


###GRANIVORES/OMNIVORES: Pearson's correlation 0.65
grani.omnivores <- ggscatter(rasterValdf,x="Granivores",y="Omnivores",xlab = "\n Granivores",ylab="Omnivores \n",
                             shape = 1,
                             alpha = 0.5,
                             size = 2,
                             conf.int = T,
                             add.params = list(color = "springgreen3",fill="lightgray"))+
  stat_cor(method = "pearson", size = 4.5,label.y = 1.1)+
  font("xy.text", size = 12)+
  font("xlab",size = 14)+
  font("ylab",size = 14)

###PLOT SIGNIFICANT PD#####
#color palette
simple.colors <- c('tomato2','cornsilk','steelblue3')

#Reclassify rastervalues to 'significance'-categories
reclass_df <- c(0.001,0.025,1,
                0.025,0.975,2,
                0.975,Inf,3)

reclass_m <- matrix(reclass_df,ncol = 3,byrow = TRUE)

###PLOT SIGNIFICANCE FOR ALL ANIMALS##
p_animalSR <- animalSR
values(p_animalSR) <- (randPD_ani$pd.obs.p)
p_animalSR <- mask(p_animalSR, plantN2)

animal_classified <- reclassify(p_animalSR,reclass_m,right=F)

randAnimalPD <- tm_shape(animal_classified)+
  tm_raster(style="cat",palette=simple.colors)+
  map_nor+
  tm_layout(legend.outside.position = T)

###PLANTS###
p_plantSR <- plantN2
values(p_plantSR) <- (rand_test_plant2$pd.obs.p)

plant_classified <- reclassify(p_plantSR,reclass_m,right=F)
plant_classified <- mask(plant_classified, plantN2)
names(plant_classified) <- "Plants"

randplantPD <- tm_shape(plant_classified)+
  tm_raster(style="cat",palette=simple.colors)+
  map_nor

###HERBIVORES###
p_herbiPD <- herbivoreSR
values(p_herbiPD) <- rand_test_herb$pd.obs.p

herbi_classified <- reclassify(p_herbiPD,reclass_m,right=F)
herbi_classified <- mask(herbi_classified, plantN2)
names(herbi_classified) <- "Herbivores"

randHerbiPD <- tm_shape(herbi_classified)+
  tm_raster(style="cat",palette=simple.colors)+
  map_nor

###CARNIVORES###
p_carniPD <- carnivoreSR
values(p_carniPD) <- rand_test_carni$pd.obs.p

carni_classified <- reclassify(p_carniPD,reclass_m,right=F)
carni_classified <- mask(carni_classified, plantN2)
names(carni_classified) <- "Carnivores"


randCarniPD <- tm_shape(carni_classified)+
  tm_raster(style="cat",palette = simple.colors)+
  map_nor

###INVERTIVORES###
p_invertiPD <- invertivorePD
values(p_invertiPD) <- rand_test_inverti$pd.obs.p

inverti_classified <- reclassify(p_invertiPD,reclass_m,right=F)
inverti_classified <- mask(inverti_classified, plantN2)
names(inverti_classified) <- "Invertivores"

randInvertiPD <- tm_shape(inverti_classified)+
  tm_raster(style="cat",palette = simple.colors)+
  map_nor

###OMNIVORES###
p_omniPD <- omnivorePD
values(p_omniPD) <- rand_test_omni$pd.obs.p

omni_classified <- reclassify(p_omniPD,reclass_m,right=F)
omni_classified <- mask(omni_classified, plantN2)
names(omni_classified) <- "Omnivores"


randOmniPD <- tm_shape(omni_classified)+
  tm_raster(style="cat",palette = simple.colors)+
  map_nor

###GRANIVORES###
p_graniPD <- granivorePD
values(p_graniPD) <- 0
p_graniPD <- mask(p_graniPD, plantN2)

values(p_graniPD) <- rand_test_grani$pd.obs.p

grani_classified <- reclassify(p_graniPD,reclass_m,right=F)
grani_classified <- mask(grani_classified, plantN2)
names(grani_classified) <- "Granivores"

randGraniPD <- tm_shape(grani_classified)+ 
  tm_raster(style="cat",palette = simple.colors)+ 
  map_nor

####PLOT ALL OF THEM TOGETHER
rand.p.all <-tmap_arrange(randplantPD,
                          randHerbiPD,
                          randCarniPD,
                          randInvertiPD,
                          randGraniPD,
                          randOmniPD, ncol = 3, nrow = 2)  

###ACCORDING to SES pd (pd.obs.z)####
#- Standardized effect size of PD vs. null communities
#High SES values represent cells where the mammal assemblage captures more PD
#than expected by chance, as estimated by sampling the same number
ses_animal <- animalSR
values(ses_animal) <- (randPD_ani$pd.obs.z)
ses_animal[is.na(ses_animal)] <- 0
ses_animal <- mask(ses_animal, plantN2)

SESanimalPD <- tm_shape(ses_animal)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0)+
  map_nor

###PLANTS###
ses_plants <- plantN2
values(ses_plants) <- (rand_test_plant2$pd.obs.z)
ses_plants[is.na(ses_plants)] <- 0
ses_plants <- mask(ses_plants, plantN2)

SESplantPD <- tm_shape(ses_plants)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0)+
  map_nor

###HERBIVORES###
ses_herbi <- herbivoreSR
values(ses_herbi) <- (rand_test_herb$pd.obs.z)
ses_herbi[is.na(ses_herbi)] <- 0
ses_herbi <- mask(ses_herbi, plantN2)

SESherbiPD <- tm_shape(ses_herbi)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0, title = "")+
  map_nor

###GRANIVORES###
ses_grani <- granivoreSR
values(ses_grani) <- (rand_test_grani$pd.obs.z)
ses_grani[is.na(ses_grani)] <- 0
ses_grani <- mask(ses_grani, plantN2)

SESgraniPD <- tm_shape(ses_grani)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0, title = "")+
  map_nor

###OMNIVORES##
ses_omni <- omnivoreSR
values(ses_omni) <- (rand_test_omni$pd.obs.z)
ses_omni[is.na(ses_omni)] <- 0
ses_omni <- mask(ses_omni, plantN2)

SESomniPD <- tm_shape(ses_omni)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0, title = "")+
  map_nor

###CARNIVORES###
ses_carni <- carnivoreSR
values(ses_carni) <- (rand_test_carni$pd.obs.z)
ses_carni[is.na(ses_carni)] <- 0
ses_carni <- mask(ses_carni, plantN2)

SEScarniPD <- tm_shape(ses_carni)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0, title = "")+
  map_nor

###INVERTIVORES
ses_inverti <- invertivoreSR
values(ses_inverti) <- (rand_test_inverti$pd.obs.z)
ses_inverti[is.na(ses_inverti)] <- 0
ses_inverti <- mask(ses_inverti, plantN2)

SESinvertiPD <- tm_shape(ses_inverti)+
  tm_raster(style="cont",palette = "RdYlBu",midpoint=0, title = "")+
  map_nor

###PLOT ALL OF THEM TOGETHER

ses.pd <- tmap_arrange(SESplantPD,
                       SESherbiPD,
                       SESgraniPD,
                       SESinvertiPD,
                       SESomniPD,
                       SESinvertiPD, ncol = 3, nrow = 2)

#####RESIDUALS####
#COLOR PALETTE - DIVERGING 
mycolors_res <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
brks2 = seq(-0.7,0.8,0.1) #need to change this as 

###PLANT ~ HERBIVORE REGRESSION & RESIDUALS
PDval <- getValues(my.rasterstack)
PDvaldf <- as.data.frame(PDval)
PDvaldf.1 <- PDvaldf[,1:2]
PDvaldf.1[is.na(PDvaldf.1)] <- 0

fit <- lm(Plant ~ Herbivores, data = PDvaldf.1)

PDvaldf.1$Predicted <- predict(fit)
PDvaldf.1$residuals <- residuals(fit)

#PLOT
PlHe.plot <- ggplot(PDvaldf.1, aes(x = Plant, y = Herbivores))+# Set up canvas with outcome variable on y-axis
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=Predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#MAP RESIDUALS
res_PlHr <- herbivoreSR
values(res_PlHr) <- (PDvaldf.1$residuals)
res_PlHr <- mask(res_PlHr, plantN2)

MAP_res_PlHr <- spplot(res_PlHr,
                       colorkey = list(width = 1,height = 0.5,axis.line = list(col = "black")),
                       at = round(brks2,digits = 2),
                       col.regions = mycolors_res,
                       main = list(label = "Plant/herbivores", cex = 1),
                       par.settings = list(axis.line = list(col = 'transparent')),
                       sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                       contour =F)

###HERBIVORE ~ CARNIVORES###
PDvaldf.2 <- PDvaldf[,c(2,6)]
PDvaldf.2[is.na(PDvaldf.2)] <- 0
#Fit regression
fit.2 <- lm(Herbivores ~ Carnivores, data = PDvaldf.2)
#add regression values
PDvaldf.2$predicted <- predict(fit.2)
PDvaldf.2$residuals <- residuals(fit.2)

#PLOT
HrCr.plot <- ggplot(PDvaldf.2, aes(x = Herbivores, y = Carnivores))+
  geom_smooth(method = "lm", se = FALSE, color = "lightgrey")+
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#MAP RESIDUALS
res_HrCr <- plantN2
values(res_HrCr) <- (PDvaldf.2$residuals)
res_HrCr <- mask(res_HrCr, plantN2)

MAP_res_HrCr <- spplot(res_HrCr,
                       colorkey = FALSE,
                       at = round(brks2,digits = 2),
                       col.regions = mycolors_res,
                       #main = list(label = "Herbivore/carnivores", cex = 1),
                       par.settings = list(axis.line = list(col = 'transparent')),
                       sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                       contour =F)


###PLANT~OMNIVORES###
PDvaldf.4 <- PDvaldf[,c(1,4)]
PDvaldf.4[is.na(PDvaldf.4)] <- 0
#Fit regression
fit.4 <- lm(Plant ~ Omnivores, data = PDvaldf.4)
#add regression values
PDvaldf.4$predicted <- predict(fit.4)
PDvaldf.4$residuals <- residuals(fit.4)

#PLOT
HrCr.plot <- ggplot(PDvaldf.4, aes(x = Plant, y = Omnivores))+
  geom_abline(intercept = 0.02473, slope = 0.73843)+
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#MAP RESIDUALS
res_PlOm <- plantN2
values(res_PlOm) <- (PDvaldf.4$residuals)
res_PlOm <- mask(res_PlOm, plantN2)

MAP_res_PlOm <- spplot(res_PlOm,
                       colorkey = FALSE,
                       at = round(brks2,digits = 2),
                       col.regions = mycolors_res,
                       #main = list(label = "Plant/omnivores", cex = 1),
                       par.settings = list(axis.line = list(col = 'transparent')),
                       sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                       contour =F)

###PLANT/INVERTIVORE###
PDvaldf.5 <- PDvaldf[,c(1,5)]
PDvaldf.5[is.na(PDvaldf.5)] <- 0
#Fit regression
fit.5 <- lm(Plant ~ Invertivores, data = PDvaldf.5)
#add regression values
PDvaldf.5$predicted <- predict(fit.5)
PDvaldf.5$residuals <- residuals(fit.5)

#PLOT
PlInv.plot <- ggplot(PDvaldf.5, aes(x = Plant, y = Invertivores))+
  geom_abline(intercept = 0.03007, slope = 0.83295)+
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#MAP RESIDUALS
res_PlInv <- plantN2
values(res_PlInv) <- (PDvaldf.5$residuals)
res_PlInv <- mask(res_PlInv, plantN2)

MAP_res_PlInv <- spplot(res_PlInv,
                        colorkey = FALSE,
                        at = round(brks2,digits = 2),
                        col.regions = mycolors_res,
                        #main = list(label = "Plant/invertivoers", cex = 1),
                        par.settings = list(axis.line = list(col = 'transparent')),
                        sp.layout=list("sp.polygons",norwayutm,
                                       first = F,col="gray25",alpha=0.5),
                        contour =F)

###HERBIVORES~INVERTIVORES###
PDvaldf.6 <- PDvaldf[,c(2,5)]
PDvaldf.6[is.na(PDvaldf.6)] <- 0
#Fit regression
fit.6 <- lm(Herbivores ~ Invertivores, data = PDvaldf.6)
#add regression values
PDvaldf.6$predicted <- predict(fit.6)
PDvaldf.6$residuals <- residuals(fit.6)

#PLOT
HrInv.plot <- ggplot(PDvaldf.6, aes(x = Herbivores, y = Invertivores))+
  geom_abline(intercept = -0.0007028, slope = 1.3413382)+
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#MAP RESIDUALS
res_HrInv <- plantN2
values(res_HrInv) <- (PDvaldf.6$residuals)
res_HrInv <- mask(res_HrInv, plantN2)

MAP_res_HrInv <- spplot(res_HrInv,
                        colorkey = FALSE,
                        at = round(brks2,digits = 2),
                        col.regions = mycolors_res,
                        #main = list(label = "Herbivores/invertivores", cex = 1),
                        par.settings = list(axis.line = list(col = 'transparent')),
                        sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25",alpha=0.5),
                        contour =F)

###HERBIVORES~OMNIVORES###
PDvaldf.7 <- PDvaldf[,c(2,4)]
PDvaldf.7[is.na(PDvaldf.7)] <- 0
#Fit regression
fit.7 <- lm(Herbivores ~ Omnivores, data = PDvaldf.7)
#add regression values
PDvaldf.7$predicted <- predict(fit.7)
PDvaldf.7$residuals <- residuals(fit.7)

#PLOT
HrOmn.plot <- ggplot(PDvaldf.7, aes(x = Herbivores, y = Omnivores))+
  geom_abline(intercept = 0, slope = 1.24)+
  geom_point(col= "lightsalmon")+#plot the actual points
  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#MAP RESIDUALS
res_HrOmn <- plantN2
values(res_HrOmn) <- (PDvaldf.7$residuals)
res_HrOmn <- mask(res_HrOmn, plantN2)

MAP_res_HrOmn <- spplot(res_HrOmn,
                        colorkey = FALSE,
                        at = round(brks2,digits = 2),
                        col.regions = mycolors_res,
                        #main = list(label = "Herbivores/omnivores", cex = 1),
                        par.settings = list(axis.line = list(col = 'transparent')),
                        sp.layout=list("sp.polygons",norwayutm,
                                       first = F,col="gray25",alpha=0.5),
                        contour =F)

###GRANIVORES~OMNIVORES###
#PDvaldf.8 <- PDvaldf[,c(3,4)]
#PDvaldf.8[is.na(PDvaldf.8)] <- 0
#Fit regression
#fit.8 <- lm(Granivores ~ Omnivores, data = PDvaldf.8)
#add regression values
#PDvaldf.8$predicted <- predict(fit.8)
#PDvaldf.8$residuals <- residuals(fit.8)

#PLOT
#ggplot(PDvaldf.8, aes(x = Granivores, y = Omnivores))+
# geom_abline(intercept = 0, slope = 0.705)+
#  geom_point(col= "lightsalmon")+#plot the actual points
#  geom_point(aes(y=predicted),shape = 1)+ #add the predicted values
#  theme(panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#      panel.background = element_blank(), axis.line = element_line(colour = "black"))

#MAP RESIDUALS
#res_GrOmn <- plantN2
#values(res_GrOmn) <- (PDvaldf.8$residuals)
#res_GrOmn <- mask(res_GrOmn, plantN2)

#MAP_res_GrOmn <- spplot(res_GrOmn,
#                       colorkey = FALSE,
#                        at = round(brks2,digits = 2),
#                       col.regions = mycolors_res,
#                      main = list(label = "Granivores/omnivores", cex = 1),
#                      par.settings = list(axis.line = list(col = 'transparent')),
#                      sp.layout=list("sp.polygons",norwayutm,first = F,col="gray25"),
#                      contour =F)

##MAP ALL OF THE RESIDUALS MAP TOGETHER
residualMap <- c(MAP_res_PlHr,
                 MAP_res_PlOm,
                 MAP_res_PlInv,
                 MAP_res_HrInv,
                 MAP_res_HrCr,
                 MAP_res_HrOmn,layout=c(3,3))

####CORRELATION BETWEEN SR & PD####
#Plot correlation between rel. plant PD and plant SR

AnimalPD_val <-getValues(animalPD)
AnimalSR_val <- getValues(rel_animalSR)

plantPD_val <- getValues(plantPD)
plantSR_val <- getValues(rel_plantSR)

df1 <- cbind(AnimalSR_val,AnimalPD_val,plantSR_val,plantPD_val)
df1 <- as.data.frame(df1)

plantSR_plantPD <- ggplot(df1, aes(x = plantSR_val, y = plantPD_val))+
  xlab('\n Species richness')+
  ylab('Phylogenetic diversity\n')+
  geom_point(shape=1,alpha = 0.3)+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_abline(slope = 1, intercept = 0, colour = '#FF6C90')+
  geom_text(x=0.05, y=0.8, label="a)")

aniSR_aniPD <- ggplot(df1, aes(x = AnimalSR_val, y = AnimalPD_val))+
  xlab('\n Species richness')+
  ylab('Phylogenetic diversity\n')+
  geom_point(shape=1,alpha = 0.3)+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_abline(slope = 1, intercept = 0, colour = '#FF6C90')+
  geom_text(x=0.05, y=0.83, label="b)")

cor.test(AnimalSR_val, AnimalPD_val, method=c("pearson"))
cor.test(plantSR_val, plantPD_val, method=c("pearson"))
