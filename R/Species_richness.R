### Making species richness maps
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
               reshape2,update = F)

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

####BIRD PREPERATION####
#Get the right names in birdN1, so it can be alignes with traitsData
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
NorwayBirds <- read_csv("NorwayBirds.csv") #250 species (251 in birdN1, but Pinguinus impennis is extinct)
NorBird <- traitDataBird %>%
  filter(Scientific %in% birdN1@data[["SCINAME"]]) #250 
#Filtering out birds who are predominantly feeding pelagic  
NorBirdTerr <- NorBird[NorBird$PelagicSpecialist == 0,]
#Subset birdN1 to the species within NorBirdTerr
birdN2 <- subset(birdN1, birdN1@data$SCINAME %in% NorBirdTerr$Scientific)

####SR PLANTS####
GDALinfo("Data/Ranges/Plants/PlantSpeciesRichness.tif")
nlayers(plantN1)

####SR MAMMALS####
##Remove species that are not terrestrial (?)
RemoveSpecies <- c("Cystophora cristata","Erignathus barbatus",
                   "Halichoerus grypus","Phoca vitulina","Pusa hispida","Dama dama",
                   "Micromys minutus")
mamN2 <- subset(mamN1,!(mamN1@data$BINOMIAL%in% RemoveSpecies)) #55 species
#Remove polygons with species that are extinct
mamN2 <- mamN2[mamN2@data$LEGEND != "Extinct",]

mamTrans2<-spTransform(mamN2,crs(plantSR))
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
e_mam2<-extent(mamTrans2)
s_mam2<-raster(e_mam2, resolution=20000, crs=(mamTrans2))
mamSF2 <- st_as_sf(mamTrans2)
#Rasterize set of polygons (sf to raster)
stackedDistributionMaps_mam2<-fasterize(mamSF2,plantN1,by="BINOMIAL") #rasterbrick
#merge polygons per species for species richness map, unioning geometries 
speciesPoly_mam2<-aggregate(mamTrans2,by="BINOMIAL")
speciesPoly2_mam2<-st_as_sf(speciesPoly_mam2)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_mam2 <- fasterize(speciesPoly2_mam2,plantN1,fun="count",field ="BINOMIAL")
#spplot(richnessMap_mam2, zlim=c(16,39))

####SR AMPHIBIANS####
amphTrans1<-spTransform(amphN1,crs(plantSR))
e_amph<-extent(amphTrans1)
s_amph<-raster(e_amph, resolution=20000, crs=(amphTrans1))
amphSF <- st_as_sf(amphTrans1)
stackedDistributionMaps_amph<-fasterize(amphSF,plantN1,by="amp")
speciesPoly_amph<-aggregate(amphTrans1,by="amp")
speciesPoly2_amph<-st_as_sf(speciesPoly_amph)
richnessMap_amph <- fasterize(speciesPoly2_amph,plantN1,field="amp")
plot(richnessMap_amph)

####STACKED AMPHIBIANS####  
##buf_buf
stacked_bufo<-fasterize(amphSF,plantN1,by="buf_buf")#rasterbrick
#Remove buf_buf = 0
stacked_bufo <- dropLayer(stacked_bufo,1)
speciesPoly_buf <- aggregate(amphTrans1, by="buf_buf")
speciesPoly2_buf<-st_as_sf(speciesPoly_buf)
speciesPoly2_buf <- speciesPoly2_buf[2,]
#colnames(speciesPoly2_buf)[1] <- "x"
richnessMap_buf <- fasterize(speciesPoly2_buf,plantN1,field="buf_buf")

##liss_vul
stacked_liss <- fasterize(amphSF,plantN1,by="liss_vul")
stacked_liss <- dropLayer(stacked_liss,1)
speciesPoly_liss <- aggregate(amphTrans1, by="liss_vul")
speciesPoly2_liss<-st_as_sf(speciesPoly_liss)
speciesPoly2_liss <- speciesPoly2_liss[2,]
#colnames(speciesPoly2_liss)[1] <- "x"
richnessMap_liss <- fasterize(speciesPoly2_liss,plantN1,field="liss_vul")

##pelp_esc
stacked_pelp <- fasterize(amphSF,plantN1,by="pelp_esc")
stacked_pelp <- dropLayer(stacked_pelp,1)
speciesPoly_pelp <- aggregate(amphTrans1, by="pelp_esc")
speciesPoly2_pelp<-st_as_sf(speciesPoly_pelp)
speciesPoly2_pelp <- speciesPoly2_pelp[2,]
#colnames(speciesPoly2_pelp)[1] <- "x"
richnessMap_pelp <- fasterize(speciesPoly2_pelp,plantN1,field="pelp_esc")

##ran_arv
stacked_ranA <- fasterize(amphSF,plantN1,by="ran_arv")
stacked_ranA <- dropLayer(stacked_ranA,1)
speciesPoly_ranA <- aggregate(amphTrans1, by="ran_arv")
speciesPoly2_ranA<-st_as_sf(speciesPoly_ranA)
speciesPoly2_ranA <- speciesPoly2_ranA[2,]
#colnames(speciesPoly2_ranA)[1] <- "x"
richnessMap_ranA <- fasterize(speciesPoly2_ranA,plantN1,field="ran_arv")

##ran_temp
stacked_ranT <- fasterize(amphSF,plantN1,by="ran_temp")
stacked_ranT <- dropLayer(stacked_ranT,1)
speciesPoly_ranT <- aggregate(amphTrans1, by="ran_temp")
speciesPoly2_ranT<-st_as_sf(speciesPoly_ranT)
speciesPoly2_ranT <- speciesPoly2_ranT[2,]
#colnames(speciesPoly2_ranT)[1] <- "x"
richnessMap_ranT <- fasterize(speciesPoly2_ranT,plantN1,field="ran_temp")

##tri_cris
stacked_tri <- fasterize(amphSF,plantN1,by="tri_cris")
stacked_tri <- dropLayer(stacked_tri,1)
speciesPoly_tri <- aggregate(amphTrans1, by="tri_cris")
speciesPoly2_tri<-st_as_sf(speciesPoly_tri)
speciesPoly2_tri <- speciesPoly2_tri[2,]
#colnames(speciesPoly2_tri)[1] <- "x"
richnessMap_tri <- fasterize(speciesPoly2_tri,plantN1,field="tri_cris")

#x <- rbind(speicesPoly2_buf2,speciesPoly2_liss,speciesPoly2_pelp,speciesPoly2_ranA,
#speciesPoly2_ranT,speciesPoly2_tri)

stackedDistributionMaps_amph2 <- stack(stacked_bufo,stacked_liss,stacked_pelp, stacked_ranA,
                                       stacked_ranT,stacked_tri)
names(stackedDistributionMaps_amph2)<- c('buf_buf','liss_vul','pelp_esc','ran_arv','ran_temp','tri_cris')
class(stackedDistributionMaps_amph2) #rasterstack instead of -brick - does that matter?

####SR REPTILES####
repTrans1 <-spTransform(repN1,crs(plantSR))
e_rep<-extent(repTrans1)
s_rep<-raster(e_rep, resolution=20000, crs=(repTrans1))
repSF <- st_as_sf(repTrans1)
stackedDistributionMaps_rep<-fasterize(repSF,plantN1,by="rep")
speciesPoly_rep<-aggregate(repTrans1,by="rep")
speciesPoly2_rep<-st_as_sf(speciesPoly_rep)
richnessMap_rep <- fasterize(speciesPoly2_rep,plantN1,field="rep")
plot(richnessMap_rep)

####STACKED REPTILES####
##ang_sp
stacked_ang<-fasterize(repSF,plantN1,by="ang_sp")#rasterbrick
stacked_ang <- dropLayer(stacked_ang,1)
speciesPoly_ang <- aggregate(repTrans1, by="ang_sp")
speciesPoly2_ang<-st_as_sf(speciesPoly_ang)
speciesPoly2_ang <- speciesPoly2_ang[2,]
richnessMap_ang <- fasterize(speciesPoly2_ang,plantN1,field="ang_sp")

#cor_aus
stacked_cor<-fasterize(repSF,plantN1,by="cor_aus")#rasterbrick
stacked_cor <- dropLayer(stacked_cor,1)
speciesPoly_cor <- aggregate(repTrans1, by="cor_aus")
speciesPoly2_cor<-st_as_sf(speciesPoly_cor)
speciesPoly2_cor <- speciesPoly2_cor[2,]
richnessMap_cor <- fasterize(speciesPoly2_cor,plantN1,field="cor_aus")

#nat_nat
stacked_nat<-fasterize(repSF,plantN1,by="nat_nat")#rasterbrick
stacked_nat <- dropLayer(stacked_nat,1)
speciesPoly_nat <- aggregate(repTrans1, by="nat_nat")
speciesPoly2_nat<-st_as_sf(speciesPoly_nat)
speciesPoly2_nat <- speciesPoly2_nat[2,]
richnessMap_nat <- fasterize(speciesPoly2_nat,plantN1,field="nat_nat")

#vip_ber
stacked_vip<-fasterize(repSF,plantN1,by="vip_ber")#rasterbrick
stacked_vip <- dropLayer(stacked_vip,1)
speciesPoly_vip <- aggregate(repTrans1, by="vip_ber")
speciesPoly2_vip<-st_as_sf(speciesPoly_vip)
speciesPoly2_vip <- speciesPoly2_vip[2,]
richnessMap_vip <- fasterize(speciesPoly2_vip,plantN1,field="vip_ber")

#zoo_viv
stacked_zoo<-fasterize(repSF,plantN1,by="zoo_viv")#rasterbrick
stacked_zoo <- dropLayer(stacked_zoo,1)
speciesPoly_zoo <- aggregate(repTrans1, by="zoo_viv")
speciesPoly2_zoo<-st_as_sf(speciesPoly_zoo)
speciesPoly2_zoo <- speciesPoly2_zoo[2,]
richnessMap_zoo <- fasterize(speciesPoly2_zoo,plantN1,field="zoo_viv")


stackedDistributionMaps_rep2 <- stack(stacked_ang,stacked_cor,stacked_nat, stacked_vip,
                                      stacked_zoo)
names(stackedDistributionMaps_rep2)<- c('ang_sp','cor_aus','nat_nat','vip_ber','zoo_viv')
class(stackedDistributionMaps_rep2) #rasterstack instead of -brick - does that matter?

####SR BIRDS####
#WITH PELAGIC SPECIALISTS
#birdTrans1<-spTransform(birdN1,crs(plantSR))
#e_bird<-extent(birdTrans1)
#s_bird<-raster(e_bird, resolution=20000, crs=(birdTrans1))
#birdSF <- st_as_sf(birdTrans1)
#stackedDistributionMaps_bird<-fasterize(birdSF,plantN1,by="SCINAME")
#speciesPoly_bird<-raster::aggregate(birdTrans1,by="SCINAME")
#speciesPoly2_bird<-st_as_sf(speciesPoly_bird)
#richnessMap_bird <- fasterize(speciesPoly2_bird,plantN1,fun="count",field="SCINAME")
#plot(richnessMap_bird)

#WITHOUT PELAGIC SPECIALISTS
birdTrans2<-spTransform(birdN2,crs(plantSR))
e_bird2<-extent(birdTrans2)
s_bird2<-raster(e_bird2, resolution=20000, crs=(birdTrans2))
birdSF2 <- st_as_sf(birdTrans2)
stackedDistributionMaps_bird2<-fasterize(birdSF2,plantN1,by="SCINAME")
speciesPoly_bird2<-raster::aggregate(birdTrans2,by="SCINAME")
speciesPoly2_bird2<-st_as_sf(speciesPoly_bird2)
richnessMap_bird2 <- fasterize(speciesPoly2_bird2,plantN1,fun="count",field="SCINAME")
plot(richnessMap_bird2) 

####SR TOTAL####
r1 <- plantN1
r2 <- richnessMap_rep
r3 <- richnessMap_amph
r4 <- richnessMap_mam2
r5 <- richnessMap_bird2

#Make rasterstack with all of the richness layers
my.stack = stack(r1,r2,r3,r4,r5)
names(my.stack)<-c('PlantSR','ReptileSR','AmphibianSR','MammalSR','BirdSR')
#Extract values from the rasterstack, in order to get SR within each cell
stackValdf<-raster::extract(my.stack,(1:ncell(my.stack)), df = T)
stackValdf$total <- rowSums(stackValdf[,2:6], na.rm = T, dims = 1)
View(stackValdf)
plot(stackValdf$PlantSR,stackValdf$MammalSR)

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

MamInverti <- NorMam[NorMam$Diet.Inv>=50,] #17
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
BirdInverti <- BirdInverti[BirdInverti$Scientific != "Podiceps grisegena",]
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

nameReptInverti <- as.factor(c("Angius fragilis","Zootoca vipera"))
nameReptCarni <- as.factor(c("Coronella austriaca","Natrix natrix","Vipera berus"))

####GROUPS MAMMALS/FEEDING ####
#SR herbivore mammals
herbivoreMam <- subset(mamN2, BINOMIAL %in% nameMamHerbi)
herbivoreMamTrans<-spTransform(herbivoreMam,crs(plantSR))
e_HM <-extent(herbivoreMamTrans)
s_HM<-raster(e_HM, resolution=20000, crs=(herbivoreMamTrans))
speciesPoly_HM <- aggregate(herbivoreMamTrans, by ="BINOMIAL")
speciesPoly_HM2 <- st_as_sf(speciesPoly_HM)
#adjusting
mamHerbSR <- fasterize(speciesPoly_HM2,plantN1,fun="count",field="BINOMIAL")

#SR carnivore mammals
carnivoreMam <- subset(mamN2, BINOMIAL %in% nameMamCarni)
carnivoreMamTrans<-spTransform(carnivoreMam,crs(plantSR))
e_CM <-extent(carnivoreMamTrans)
s_CM<-raster(e_CM, resolution=20000, crs=(carnivoreMamTrans))
speciesPoly_CM <- aggregate(carnivoreMamTrans, by ="BINOMIAL")
speciesPoly_CM2 <- st_as_sf(speciesPoly_CM)
#Adjusting
mamCarnSR <- fasterize(speciesPoly_CM2,plantN1,fun="count",field="BINOMIAL")

#SR insectivores mammals
invertivoreMam <- subset(mamN2, BINOMIAL %in% nameMamInverti)
invertivoreMamTrans<-spTransform(invertivoreMam,crs(plantSR))
e_IM <-extent(invertivoreMamTrans)
s_IM<-raster(e_IM, resolution=20000, crs=(invertivoreMamTrans))
speciesPoly_IM <- aggregate(invertivoreMamTrans, by ="BINOMIAL")
speciesPoly_IM2 <- st_as_sf(speciesPoly_IM)
#Adjusting
mamInvertSR <- fasterize(speciesPoly_IM2,plantN1,fun="count",field="BINOMIAL")

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
mamOmniSR <- fasterize(speciesPoly_OM2,plantN1,fun="count",field="BINOMIAL")

####GROUPS BIRD/FEEDING####
#SR herbivore birds
herbivoreBird <- subset(birdN1, SCINAME %in% nameBirdHerbi)
herbivoreBirdTrans<-spTransform(herbivoreBird,crs(plantSR))
e_HB <-extent(herbivoreBirdTrans)
s_HB<-raster(e_HB, resolution=20000, crs=(herbivoreBirdTrans))
speciesPoly_HB <- aggregate(herbivoreBirdTrans, by ="SCINAME")
speciesPoly_HB2 <- st_as_sf(speciesPoly_HB)
birdHerbSR <- fasterize(speciesPoly_HB2,plantN1,fun="count",field="SCINAME")

#SR carnivore birds
carnivoreBird <- subset(birdN1, SCINAME %in% nameBirdCarni)
carnivoreBirdTrans<-spTransform(carnivoreBird,crs(plantSR))
e_CB <-extent(carnivoreBirdTrans)
s_CB<-raster(e_CB, resolution=20000, crs=(carnivoreBirdTrans))
speciesPoly_CB <- aggregate(carnivoreBirdTrans, by ="SCINAME")
speciesPoly_CB2 <- st_as_sf(speciesPoly_CB)
birdCarnSR <- fasterize(speciesPoly_CB2,plantN1,fun="count",field="SCINAME")

#SR granivore birds
granivoreBird <- subset(birdN1, SCINAME %in% nameBirdGrani)
granivoreBirdTrans <- spTransform(granivoreBird, crs(plantSR))
e_GB <- extent(granivoreBirdTrans)
s_GB <- raster(e_GB, resolution=20000, crs(granivoreBirdTrans))
speciesPoly_GB <- aggregate(granivoreBirdTrans, by ="SCINAME")
speciesPoly_GB2 <- st_as_sf(speciesPoly_GB)
birdGranSR <- fasterize(speciesPoly_GB2,plantN1,fun="count",field="SCINAME")#5

#SR insectivore birds
invertivoreBird <- subset(birdN1, SCINAME %in% nameBirdInverti) 
invertivoreBirdTrans <- spTransform(invertivoreBird, crs(plantSR))
e_IB <- extent(invertivoreBirdTrans)
s_IB <- raster(e_IB, resolution=20000, crs(invertivoreBirdTrans))
speciesPoly_IB <- aggregate(invertivoreBirdTrans, by ="SCINAME")
speciesPoly_IB2 <- st_as_sf(speciesPoly_IB)
birdInvertSR <- fasterize(speciesPoly_IB2,plantN1,fun="count",field="SCINAME")

#SR omnivore birds
omnivoreBird <- subset(birdN1, SCINAME %in% nameBirdOmnivore) 
omnivoreBirdTrans <- spTransform(omnivoreBird, crs(plantSR))
e_OB <- extent(omnivoreBirdTrans)
s_OB <- raster(e_OB, resolution=20000, crs(omnivoreBirdTrans))
speciesPoly_OB <- aggregate(omnivoreBirdTrans, by ="SCINAME")
speciesPoly_OB2 <- st_as_sf(speciesPoly_OB)
birdOmnivoreSR <- fasterize(speciesPoly_OB2,plantN1,fun="count",field="SCINAME")

####RASTER STACK AND EXTRACT####
###Stack mammals
mammalStack <- stack(mamCarnSR,mamHerbSR,mamInvertSR,mamOmniSR)
names(mammalStack)<-c('mamCarnSR','mamHerbSR','mamInvertSR','mamOmniSR')
mammalValdf<-raster::extract(mammalStack,(1:ncell(mammalStack)), df = T)
mammalValdf$total <- rowSums(mammalValdf[,2:ncol(mammalValdf)], na.rm = T, dims = 1)
#Extract species from each cell from Mammal rasterlayer
mammalPoly <- raster::extract(stackedDistributionMaps_mam2,(1:ncell(stackedDistributionMaps_mam2)), df = T)
mammalPoly2 <- mammalPoly
mammalPoly$total <- rowSums(mammalPoly[,2:ncol(mammalPoly)], na.rm = T, dims = 1)

###Stack birds
birdStack <- stack(birdCarnSR, birdHerbSR, birdGranSR, birdInvertSR, birdOmnivoreSR)
names(birdStack)<-c('birdCarnSR','birdHerbSR','birdGranSR','birdInvertSR', 'birdOmniSR')
birdValdf<-raster::extract(birdStack,(1:ncell(birdStack)), df = T)
birdValdf$total <- rowSums(birdValdf[,2:ncol(birdValdf)], na.rm = T, dims = 1)
#Extract species within each cell from bird rasterlayer
birdPoly <- raster::extract(stackedDistributionMaps_bird2,(1:ncell(stackedDistributionMaps_bird2)), df = T)
birdPoly2 <- birdPoly[,-1]
birdPoly$total <- rowSums(birdPoly[,2:ncol(birdPoly)], na.rm = T, dims = 1)

###Amphibian
#Extract which amphibians are in which cell
amphPoly <- raster::extract(stackedDistributionMaps_amph2,(1:ncell(stackedDistributionMaps_amph2)), df = T)
amphPoly2 <- amphPoly[,-1]
names(amphPoly2)<- c("Bufo.bufo","Lissotriton.vulgaris","Pelophylax.lessonae",
                     "Rana.arvalis","Rana.temporaria","Triturus.cristatus")
amphPoly$total <- rowSums(amphPoly[,2:ncol(amphPoly)], na.rm = T, dims = 1)

#Reptiles
reptileStack <- stack(RepInvertSR,RepCarniSR)
names(reptileStack)<-c('repInvertSR','repCarniSR')
reptValdf<-raster::extract(reptileStack,(1:ncell(reptileStack)), df = T)
reptValdf$total <- rowSums(reptValdf[,2:ncol(reptValdf)], na.rm = T, dims = 1)
#Extract which reptiles are in which cell
reptPoly <- raster::extract(stackedDistributionMaps_rep2,(1:ncell(stackedDistributionMaps_rep2)), df = T)
reptPoly2 <- reptPoly[,-1]
names(reptPoly2) <- c("Angius.fragilis","Coronella.austriaca","Natrix.natrix",
                      "Vipera.berus","Zootoca.vipera")
reptPoly$total <- rowSums(reptPoly[,2:ncol(reptPoly)], na.rm = T, dims = 1)

###ALL
allStack <- stack(RepInvertSR,RepCarniSR, AmphInvertSR,plantN1,birdHerbSR, birdInvertSR, birdGranSR, birdCarnSR, birdOmnivoreSR,
                  mamCarnSR,mamHerbSR,mamInvertSR,mamOmniSR)
names(allStack)<- c('repInvertSR','repCarniSR','amphInvertSR','plantSR','birdHerbSR','birdInvertSR','birdGranSR','birdCarnSR','birdOmnivoreSR',
                    'mamCarnSR','mamHerbSR','mamInvertSR','mamOmniSR')
allStackValdf <- raster::extract(allStack,(1:ncell(allStack)), df = T)
allStackValdf$total <- rowSums(allStackValdf[,2:ncol(allStackValdf)], na.rm = T, dims = 1)
plot(allStack$plantSR, allStack$birdHerbSR)

####COMBINE ALL POLYGONS####
allPoly <- bind_cols(mammalPoly2,birdPoly2,amphPoly2,reptPoly2)
allPoly_molten <- melt(allPoly,id ="ID")
allPoly_molten <- subset(allPoly_molten, value==1)
names(allPoly_molten) <-c("ID", "Scientific", "value")
levels(allPoly_molten$Scientific) <- gsub("\\."," ",levels(allPoly_molten$Scientific))

####Create table with all species (ex. plants) and taxonomic/trophic grouping####
mammal_latin <- as.character(NorMam$Scientific)
bird_latin <- as.character(NorBirdTerr$Scientific)
amph_latin <- as.character(NorAmph$Species)
rept_latin <- c("Angius fragilis","Coronella austriaca","Natrix natrix",
                "Vipera berus","Zootoca vipera")
Scientific <- c(mammal_latin,bird_latin,amph_latin,rept_latin) #286
table <- data.frame(Scientific)
#write.table(table$Scientific, file = "Species_list.txt", row.names = F,col.names = F,quote = F, sep = "")

#Make new column which categorize according to taxonomic group (mammal, bird, amph, rept)
table$taxonomic <- c(1:(length(table$Scientific)))
table$taxonomic[1:50] <- "mammal"
table$taxonomic[51:275] <- "bird"
table$taxonomic[276:281] <- "amphibian"
table$taxonomic[282:286] <- "reptile"
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

####Join the table with taxonomic/trophic grouping with df of extracted species within each cell####
allPoly_merged<- merge(allPoly_molten,table, by="Scientific", all.x = T)
allPoly_merged <- allPoly_merged[,c(2,1,4,5,3)]
tableOfSpecies <- allPoly_merged[,-c(5)]

#write.csv(tableOfSpecies, file="Species_table.csv")

#####RASTER OVERLAY####  
#raster_result <- overlay(r1,r2,r3,r4,r5,fun=function(x,y,z,a,b){return(x+y+z+a+b)})
#plot(raster_result)
#raster_result
#spplot(raster_result)
#table(getValues(raster_result))
#intersect(nameBirdCarni,nameBirdInsect)
