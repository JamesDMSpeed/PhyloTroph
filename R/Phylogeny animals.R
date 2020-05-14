####LOADING LIBRARIES####
library(ape)
library(picante)
library(readr)
library(plyr)
library(raster)
library(sp)
library(ggplot2)
library(dplyr)
library(geiger)

###ANIMAL PHYLOGENY####
phylogeny <-read.tree('Data/Phylogenies/NorTerrSpecies_285.nwk') #Stored as phylo object
phylogeny #285 tips and 284 internal nodes - the tree is probably bifurcate
plot.phylo(phylogeny, show.tip.label = F, main="Norway")
#Total amount of PD
animal_totPd <- sum(phylogeny$edge.length)

#Check if the phylogeny is rooted and ultrametric
#is.rooted(phylogeny)  #T
#is.ultrametric(phylogeny) #T
#Check names in phylogeny
head(phylogeny$tip.label) #Missing M.striata. S. Isodon have been removed from SR analsysis (DD)

#Read community matrix (285 species, incl. M.striata)
my.sample <- read.table("Data/CommunityMat_2.txt",sep="\t", header = T)
names(my.sample) <- gsub("\\.","_",names(my.sample)) #M.striata included, S.isodon excluded
my.sample <- my.sample[,-1]
my.sample[is.na(my.sample)] <- 0
#SR at each site(raster cell) - should correspond to values in stackValdf
rowSums(my.sample) 

#Which species names are in allPoly, but not in the phylogeny?
(setdiff(names(my.sample),(phylogeny$tip.label))) #21
#21 species from my.sample are missing in the phylogeneis. 
names(my.sample) <- mapvalues((names(my.sample)),
                              c("Emberiza_rustica","Saxicola_torquatus",
                                "Bonasa_bonasia","Bubo_scandiaca","Philomachus_pugnax",
                                "Eudromias_morinellus","Larus_minutus","Carduelis_flammea",
                                "Carduelis_spinus","Larus_ridibundus","Parus_ater","Parus_cinctus",
                                "Parus_cristatus","Parus_montanus","Parus_palustris",
                                "Tetrao_tetrix","Anas_penelope","Anas_strepera","Parus_caeruleus",
                                "Phylloscopus_borealis"),
                              c("Schoeniclus_rusticus","Saxicola_torquata","Tetrastes_bonasia",
                                "Bubo_scandiacus","Calidris_pugnax","Charadrius_morinellus",
                                "Hydrocoloeus_minutus","Acanthis_flammea","Spinus_spinus",
                                "Chroicocephalus_ridibundus","Periparus_ater","Poecile_cinctus",
                                "Lophophanes_cristatus","Poecile_montanus","Poecile_palustris",
                                "Lyrurus_tetrix","Mareca_penelope","Mareca_strepera",
                                "Cyanistes_caeruleus","Seicercus_borealis"),warn_missing = T)
#Missing species?                                      
setdiff(names(my.sample),(phylogeny$tip.label)) #Muscicapa_striata - couldn't synonym in Timetree

#####PLANT PHYLOGENY#####
#Norway plant phylogeny (Mienna)
plantPhylo <- read.tree('Data/Phylogenies/NorwVascPlantPhylogeny.nwk') #1238 tips
plot(plantPhylo, show.tip.label = F, main = 'Norway plant phylogeny') #NOT ULTRAMETRIC
#Sum branch lengths
plant.totPd <- sum(plantPhylo$edge.length) #78.156

##Ultrametric rate smoothing - MPL method
#pp <- chronos(plantPhylo) # returned negative branch lengths
#ppp <- chronoMPL(plantPhylo) #also returned negative branch lenghts

#Load plant distributional data
load("/Data/brick_native.RData") #https://ntnu.app.box.com/s/050j37osuazm4ezl9rlz9oc3srdwz9kb

#Extract species layer that is sum of all species (SR)
layer_1119 <- subset(brick_native,"index_1119")
brick_native_new <- dropLayer(brick_native, c("index_1119")) #resulting in 1119 species layers

#Community matrix plant
comm_all <- getValues(brick_native_new)
comm_all[is.na(comm_all)]<-0

####CONVERT longlat file from Mienna TO POLYGONS (IN ORDER TO EXTRACT COMMUNITY MATRIX)
#https://datadryad.org/stash/dataset/doi:10.5061/dryad.j9kd51c7j
nor_vasc_latlong <- read_csv("Data/Ranges/nor_vasc_latlong.csv")
nor_vasc_SF <- st_as_sf(nor_vasc_latlong,coords = c("decimalLatitude","decimalLongitude"))
#assign WGS84 projection (EPSG code 4326)
st_crs(nor_vasc_SF) <- 4326 

#Convert from point geometry to polygon
#polys = st_sf(
# aggregate(
#  nor_vasc_SF$geometry,
# list(nor_vasc_SF$speciesName),
#function(g){
# st_cast(st_combine(g),"POLYGON")
#}
#)) #ERROR: polygons require at least 4 points 

#####FIX PLANTS####
#check if names match phylogeny
phydata <- match.phylo.comm(plantPhylo, comm_all) 
#17 names from cumminity matrix not in phylogeny
#136 names in phylogeny not in community matrix
#fix names that are different
colnames(comm_all) <- gsub("\\.", "-", colnames(comm_all))

colnames(comm_all)[237] <- "Celastraceae_Parnassia_palustris"
colnames(comm_all)[400] <- "LYCO_Isoetaceae_Isoetes_echinospora"
colnames(comm_all)[401] <- "LYCO_Isoetaceae_Isoetes_lacustris"
colnames(comm_all)[780] <- "Orobanchaceae_Pedicularis_sceptrum-carolinum"
colnames(comm_all)[1074] <- "Saxifragaceae_Saxifraga_osloensis"

#fix spelling in phylogeny as well
plantPhylo$tip.label[plantPhylo$tip.label=="LYCO_Isoetaceae_IsoâˆšÂ´tes_echinospora"] <- "LYCO_Isoetaceae_Isoetes_echinospora"
plantPhylo$tip.label[plantPhylo$tip.label=="LYCO_Isoetaceae_IsoâˆšÂ´tes_lacustris"] <- "LYCO_Isoetaceae_Isoetes_lacustris"
plantPhylo$tip.label[plantPhylo$tip.label=="Saxifragaceae_Saxifraga_osloâˆšÂ´nsis"] <- "Saxifragaceae_Saxifraga_osloensis"
plantPhylo$tip.label[plantPhylo$tip.label=="MONO_Poaceae_HierochloâˆšÂ´_alpina"] <- "MONO_Poaceae_Hierochloe_alpina"
plantPhylo$tip.label[plantPhylo$tip.label=="MONO_Poaceae_HierochloâˆšÂ´_hirta"] <- "MONO_Poaceae_Hierochloe_hirta"
plantPhylo$tip.label[plantPhylo$tip.label=="MONO_Poaceae_HierochloâˆšÂ´_odorata"] <- "MONO_Poaceae_Hierochloe_odorata"

#New numbers
phydata2 <- match.phylo.comm(plantPhylo, comm_all)
#1 species dropped from commat, not idn phylogeny Dryopteris_affinis)
#120 species dropped from phylogeny, not present in community data


###CALCULATIONS####
#ANIMAL
#Calculation Faith's PD by using pd function from pixante package
phyloDiv_root <- pd(my.sample,phylogeny,include.root = T)

#PLANT
pd.plant.root <- pd(phydata$comm, phydata$phy, include.root = T)

####CALCULATE PD FOR DIFFERENT TROPHIC GROUPS####
#HERBIVORES
#METHOD 1: pruning the phylogeny, and then calculate total PD of that tree
herbiTree <- phylogeny
herbMatrix <- matrix(1:28,nrow=28,ncol=1)
rownames(herbMatrix) <- c("Lemmus_lemmus","Alces_alces","Capreolus_capreolus","Rangifer_tarandus",
                          "Lepus_timidus","Arvicola_amphibius","Castor_fiber","Microtus_agrestis",
                          "Microtus_oeconomus","Myopus_schisticolor","Cervus_elaphus","Plectrophenax_nivalis",
                          "Anser_albifrons","Anser_anser","Anser_erythropus","Anser_fabalis","Tetrastes_bonasia",
                          "Branta_canadensis","Branta_leucopsis","Cygnus_columbianus","Cygnus_cygnus","Cygnus_olor",
                          "Lagopus_lagopus","Lagopus_muta","Lyrurus_tetrix","Mareca_penelope","Mareca_strepera",
                          "Tetrao_urogallus")

herbiTree <- geiger::treedata(phylogeny, herbMatrix, warnings =T)
herbPhylo <- herbiTree$phy
herb_totPd <- sum(herbPhylo$edge.length) #1259

#METHOD 2: make community matrix with only species within trophic group:
commat_herbivore <- read_csv("Data/commat_herbivore.csv",col_types = cols(X1 = col_skip()))
herbCom <- commat_herbivore[1,]
herbCom[1,] <- 1

PDherbCom <- pd(herbCom, phylogeny, include.root = T) #TOT PD = 1299
herb_totPd2 <- sum(PDherbCom$PD) #1259.41

##CARNIVORE
commat_carnivore <- read_csv("Data/commat_carnivore.csv",col_types = cols(X1 = col_skip()))
carnCom <- commat_carnivore[1,]
carnCom[1,] <- 1

PDcarncom <- pd(carnCom, phylogeny, include.root = T) #TOT PD = 2577
carn_totPd2 <- sum(PDcarncom$PD)

##OMNIVORE
commat_omnivore <- read_csv("Data/commat_omnivore.csv",col_types = cols(X1 = col_skip()))
omnCom <- commat_omnivore[1,]
omnCom[1,] <- 1

PDomncom <- pd(omnCom, phylogeny, include.root = T) #TOT PD = 2780.55
omni_totPd <- sum(PDomncom$PD)

##INVERTIVORE
commat_invertivore <- read_csv("Data/commat_invertivore.csv",col_types = cols(X1 = col_skip()))
invCom <- commat_invertivore[1,]
invCom[1,] <- 1

PDinvcom <- pd(invCom, phylogeny, include.root = T) #TOT PD = 5663.35
invert_totPd <- sum(PDinvcom$PD)

#GRANIVORE
commat_granivore <- read_csv("Data/commat_granivore.csv",col_types = cols(X1 = col_skip()))
granCom <- commat_granivore[1,]
granCom[1,] <- 1

PDgrancom <- pd(granCom, phylogeny, include.root = T) #TOT PD = 496.28
grani_totPd <- sum(PDgrancom$PD)

###RANDOMISATIONS MODELS####
#Do randomisations of SR across the raster, in order to find cells with signf. low/high PD

#RANDOMISATION MODEL #richness


#ALL ANIMALS - BIRDS,MAMMALS,REPT AND AMPH
#outputAni_tipl <- ses.pd(my.sample, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#randPD_ani <- read_csv("outputAni_tipl.csv")

#ONLY HERBIVORES  
#rand_test444 <- ses.pd(commat_herbivore, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_herb <- read_csv("rand_test444.csv")

#ONLY CARNIVORES
#rand_test_carni <- ses.pd(commat_carnivore, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_carni <- read_csv("rand_test_carni.csv")

#ONLY INVERTIVORES
#rand_test_inverti <- ses.pd(commat_invertivore, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_inverti <- read_csv("rand_test_inverti.csv")

#ONLY GRANIVORES
#rand_test_grani <- ses.pd(commat_granivore, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_grani <- read_csv("rand_test_grani.csv")

#ONLY OMNIVORES
#rand_test_omni <- ses.pd(commat_omnivore, phylogeny, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_omni <- read_csv("rand_test_omni.csv")

#PLANTS - returns to few significant areas, compared to Mieanna
#outputPlant2_tipl <- ses.pd(comm_all, phylogeny_plant,null.model = "taxa.labels",runs=999,include.root = T)
#write.csv(outputPlant2_tipl, file = "randPD_plant2.csv")
#rand_test_plant <- read_csv("randPD_plant2.csv")

#DIFFERENT COMMUNITY MATRIX AND PHYLOGENY
#rand_test222 <- ses.pd(phydata$comm, phydata$phy, null.model = "taxa.labels",runs=999,include.root = T)
#rand_test_plant2 <- read_csv("rand_test222.csv")  

#random_plants <- ses.pd(phydata$comm, phydata$phy, null.model = "richness",include.root = T)
#rand_test_plantRI <- random_plants <- read_csv("random_plants.csv") 

