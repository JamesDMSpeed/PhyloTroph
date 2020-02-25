###Load Elton traits data
library(dplyr)

#####Load data####
traitDataMam <- read.table("MamFuncDat.txt", sep = '\t', header = T, fill = T)
View(traitDataMam)
##BIRDS
traitDataBird <- read.table("BirdFuncDat.txt", sep = '\t',header = T, fill = T, quote ='')
View(traitDataBird)

#####BIRDS####
#List of species in Norway from Handbook of the Birds of the World (?)
birdN <- readOGR("Birds_","Birbies_fix")
birdN1 <- birdN
norBirds <- levels(birdN1@data$SCINAME) #251 different bird species

#Making a subset of traitDataBird that consists only of birds found in Norway(from norBirds)
NorSpecies_birds <- traitDataBird %>%
  filter(Scientific%in% norBirds)

#Checking if all of the species from the shapefile is within Eltontraits
#####MISSING BIRDS####
#List of species found in Elton traits
eltonBirds <- as.character(NorSpecies_birds$Scientific)
#Find species BirdN1 thats not in eltonBirds, and make a list of those species
x <- setdiff(norBirds,eltonBirds)

#Read list of species that were not missing from BirdN1, but had different names
missing_birds_ <- read.csv("missing_birds_.csv")
#Make subset of those 'missing species'
NorSpecies_birds2 <- traitDataBird %>% 
  filter(Scientific %in% missing_birds_$x)

#####

NorTraitDataBird <- rbind(NorSpecies_birds,NorSpecies_birds2) #250 observations, Pinguinus impennis is missing (extint)
View(NorTraitDataBird)
