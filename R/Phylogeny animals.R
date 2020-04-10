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




