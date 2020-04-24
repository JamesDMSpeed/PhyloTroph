###TRY DIFFERENT PLANT RASTER###
####loading library####
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
               stats,
               rgeos,
               reshape,
               vegan,
               reshape2,update = F)

####READING SHAPEFILES####
#Reading shapefile for mammal
mamN <- readOGR("Data/Ranges/Mammals")#Should be read from repository project directory, not local drive
mamN1 <- mamN
##Reptiles
repN_Sillero <- readOGR("Data/Ranges/Reptiles_Sillero","reptiles_Norway_Sillero")
repN1 <- repN_Sillero
##Amphibians
amphN_Sillero <- readOGR("Data/Ranges/Amphibians","amphibian_Norway_Sillero")
amphN1 <- amphN_Sillero
#Plants
plantSR <- raster("Data/Ranges/Plants/PlantSpeciesRichness.tif")
plantN1 <- plantSR
#Load new raster brick with all plant species occurences
load("Data/brick_native.RData")
#Extract rasteralyer with SR
plantN2<- subset(brick_native,"index_1119")
projection(plantN2) <- crs(plantN1)
plantN2[plantN2 == 0]<- NA
#plantN2[is.na(plantN2)] <- 0
#plantN2 <- mask(plantN2, norwayutm)

#Birds
birdN <- readOGR("Birds_fix","Birbies_fix")
birdN1 <- birdN
#Removing a small polygon from the species of Calidris maritima (small island)
birdN1 <- birdN1[-461,]
birdN1 <- birdN1[(birdN1@data$SCINAME != "Pinguinus impennis"),] #Removing Pinguinus impennis
#Map of norway
norway<-getData('GADM',country='NOR',level=0)
norwayutm<-spTransform(norway,crs(plantN2))
norwayraster <- rasterize(norwayutm, plantN2)

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
traitDataBird <- read.table("Data/BirdFuncDat.txt", sep = '\t',header = T, fill = T, quote ='')
NorwayBirds <- read_csv("Data/NorwayBirds.csv") #250 species
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
traitDataMam <- read.table("Data/MamFuncDat.txt", sep = '\t', header = T, fill = T)
NorMam <- traitDataMam %>%
  filter(Scientific%in% mamN2@data[["BINOMIAL"]]) #50

#Trait data Amphibians
traitDataAmph <- read.csv("Data/Traits/AmphiBIO_v1.csv")
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
birdGranSR <- fasterize(speciesPoly_GB2,plantN2,fun="count",field="SCINAME")#5

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
names(speciesPoly2_cor_cop) <- c("BINOMIAL","geometry")
speciesPoly2_cor_cop$BINOMIAL <- as.factor(speciesPoly2_cor_cop$BINOMIAL)
speciesPoly2_nat_cop <- speciesPoly2_nat
names(speciesPoly2_nat_cop) <- c("BINOMIAL","geometry")
speciesPoly2_nat_cop$BINOMIAL <- as.factor(speciesPoly2_nat_cop$BINOMIAL)
speciesPoly2_vip_cop <- speciesPoly2_vip
names(speciesPoly2_vip_cop) <- c("BINOMIAL","geometry")
speciesPoly2_vip_cop$BINOMIAL <- as.factor(speciesPoly2_vip_cop$BINOMIAL)

stackCarni <- rbind(speciesPoly_CB2_cop,speciesPoly_CM2,speciesPoly2_cor_cop,#47 species
                    speciesPoly2_nat_cop,speciesPoly2_vip_cop) 
carnivoreSR <- fasterize(stackCarni,plantN2,fun="count",field="BINOMIAL")

###SR MAPS INVERTIVORES
#reptiles = zoo and ang
speciesPoly2_ang_cop <- speciesPoly2_ang
names(speciesPoly2_ang_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ang_cop$BINOMIAL <- as.factor(speciesPoly2_vip_cop$BINOMIAL)
speciesPoly2_zoo_cop <- speciesPoly2_zoo
names(speciesPoly2_zoo_cop) <- c("BINOMIAL","geometry")
speciesPoly2_zoo_cop$BINOMIAL <- as.factor(speciesPoly2_vip_cop$BINOMIAL)
#amphibians = all species
speciesPoly2_buf_cop <- speciesPoly2_buf
names(speciesPoly2_buf_cop) <- c("BINOMIAL","geometry")
speciesPoly2_buf_cop$BINOMIAL <- as.factor(speciesPoly2_buf_cop$BINOMIAL)
speciesPoly2_liss_cop <- speciesPoly2_liss
names(speciesPoly2_liss_cop) <- c("BINOMIAL","geometry")
speciesPoly2_liss_cop$BINOMIAL <- as.factor(speciesPoly2_liss_cop$BINOMIAL)
speciesPoly2_pelp_cop <- speciesPoly2_pelp
names(speciesPoly2_pelp_cop) <- c("BINOMIAL","geometry")
speciesPoly2_pelp_cop$BINOMIAL <- as.factor(speciesPoly2_pelp_cop$BINOMIAL)
speciesPoly2_ranA_cop <- speciesPoly2_ranA
names(speciesPoly2_ranA_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ranA_cop$BINOMIAL <- as.factor(speciesPoly2_ranA_cop$BINOMIAL)
speciesPoly2_ranT_cop <- speciesPoly2_ranT
names(speciesPoly2_ranT_cop) <- c("BINOMIAL","geometry")
speciesPoly2_ranT_cop$BINOMIAL <- as.factor(speciesPoly2_ranT_cop$BINOMIAL)
speciesPoly2_tri_cop <- speciesPoly2_tri
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

####PLOTS#####
#https://www.rdocumentation.org/packages/tmap/versions/3.0/topics/tm_raster
#https://geocompr.robinlovelace.net/adv-map.html
#help("tmap-element")
#Plot raster in general: tm_shape(raster)+tm_raster()

#MAP OF NORWAY
map_nor <- tm_shape(norwayutm)+tm_borders(alpha = 0.7)+tm_layout(frame = F)

#PLANTS
tm_shape(r1)+tm_raster(stretch.palette = T)+map_nor
#tm_shape(plantN4)+tm_raster(style="cont")+map_nor
#tm_raster(palette="OrRd") 

#REPTILES
tm_shape(r2)+tm_raster(stretch.palette = T)+map_nor

#AMPHIBIANS
tm_shape(r3)+tm_raster(stretch.palette = T)+map_nor

#MAMMALS
tm_shape(r4)+tm_raster(stretch.palette = T)+map_nor

#BIRDS
#r55.crop <- crop(r55, extent(plantN4)) #TRENG du dette? prøve med n2?
r55.crop <- mask(r55.crop, plantN2)
tm_shape(r55.crop)+tm_raster(stretch.palette = T)+map_nor

#GRANIVORES
granivoreSR_2 <- granivoreSR
granivoreSR_2[is.na(granivoreSR_2)] <- 0
granivoreSR_2 <- mask(granivoreSR_2, plantN2)
tm_shape(granivoreSR_2)+tm_raster(stretch.palette = T)+map_nor

#CARNIVORES
carnivoreSR_2 <- carnivoreSR
carnivoreSR_2[is.na(carnivoreSR_2)] <- 0
carnivoreSR_2 <- mask(carnivoreSR_2, plantN2)
tm_shape(carnivoreSR_2)+tm_raster(stretch.palette = T)+map_nor

#INVERTIVORES
invertivoreSR_2 <- invertivoreSR
invertivoreSR_2[is.na(invertivoreSR_2)] <- 0
invertivoreSR_2 <- mask(invertivoreSR_2, plantN2)
tm_shape(invertivoreSR_2)+tm_raster(stretch.palette = T)+map_nor

#OMNIVORES
omnivoreSR_2 <- omnivoreSR
omnivoreSR_2[is.na(omnivoreSR_2)] <- 0
omnivoreSR_2 <- mask(omnivoreSR_2, plantN2)
tm_shape(omnivoreSR_2)+tm_raster(stretch.palette = T)+map_nor
#tm_shape(omnivoreSR_2)+tm_raster(style = "cont",palette = "Blues")+map_nor

#HERBIVORES
herbivoreSR_2 <- herbivoreSR
herbivoreSR_2[is.na(herbivoreSR_2)] <- 0
herbivoreSR_2 <- mask(herbivoreSR_2, plantN2)
tm_shape(herbivoreSR_2)+tm_raster(stretch.palette = T)+map_nor

####HOW TO SAVE PLOTS####
#pdf("myplottt3.pdf")  ## open a pdf file for plotting
#plot(plantN3)
#plot(b, add=T)
#dev.off()  ## close pdf; you must do this!

####STATISTICS#####
plot(allStackValdf$plantSR, allStackValdf$mamHerbSR)

model_A <- lm(plantSR ~ mamHerbSR, data = allStackValdf)
summary(model_A)

model_B <-aov(plantSR ~ mamHerbSR, data = allStackValdf)
summary(model_B)