####Phylogenies

library(ape)
library(picante)
library(readr)

#Norway animal phylogeny
phylogeny <-read.tree('NorTerrSpecies_285.nwk') #Stored as phylo object
phylogeny #285 tips and 284 internal nodes - the tree is probably bifurcate
plot(phylogeny,show.tip.label = F,main='Norway phylogeny')
plot.phylo(phylogeny, show.tip.label = F, type = 'fan')
#Check if the phylogeny is rooted and ultrametric
#is.rooted(phylogeny)  #T
#is.ultrametric(phylogeny) #T
head(phylogeny$edge)
head(phylogeny$tip.label)
head(phylogeny$node.label)
tail(phylogeny$node.label)

#Norway plant phylogeny

test <- read.tree('Data/Phylogenies/ArcticBorealHerbivorePhylogeny.newick')
plot(test, show.tip.label = F)
is.ultrametric(test)



plantPhylo <- read.tree('Data/Phylogenies/NorwVascPlantPhylogeny.nwk')
plot(plantPhylo, show.tip.label = F, main = 'Norway plant phylogeny')
#is.rooted(plantPhylo) TRUE
#is.ultrametric(plantPhylo)FALSE

plantPhylo_2 <- read.tree('Data/Phylogenies/Prum_merge_Genera_2016_04_22_Hackett_Stage2_10K_MCC.newick.tre')
plot(plantPhylo, show.tip.label = F)
#is.rooted(plantPhylo_2) TRUE
#is.ultrametric(plantPhylo_2) FALSE 


#Read community matrix
my.sample <- read.table("communityMat.txt",sep="\t", header = T)
names(my.sample) <- gsub("\\.","_",names(my.sample))
my.sample <- my.sample[,-1]
my.sample[is.na(my.sample)] <- 0

#Species richness at it site, should correspond to values in stackValdf
rowSums(my.sample) 

#Which species names are in allPoly, but not in the phylogeny?
x <- c(setdiff(names(my.sample),(phylogeny$tip.label))) #21
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

setdiff(names(my.sample),(phylogeny$tip.label)) #Muscicapa_striata - couldn't find alternative names
#Calculation Faith's PD by using pd function from pixante package
phyloDiv <- pd(my.sample,phylogeny,include.root = F)
phyloDiv_root <- pd(my.sample,phylogeny,include.root = T)


#plants
library(raster)
library(sp)
load("brick_native.RData")
plot(brick_native)

#remove layer "index_1119" for this, in this layer is the sum of all species (SR) 
brick_native_new <- dropLayer(brick_native, c("index_1119"))

#check if names match phylogeny
comm_all <- getValues(brick_native_new)
comm_all[is.na(comm_all)]<-0
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

phydata <- match.phylo.comm(plantPhylo, comm_all) 
#1 name from cumminity matrix not in phylogeny (not sure why, this is in the overview for the phylogeny with this name)
#120 names in phylogeny not in community matrix (these species have probably no occurrences)

pd_plant <- pd(comm_all, plantPhylo, include.root = F)
pd_plant_root <- pd(comm_all, plantPhylo, include.root = T)

