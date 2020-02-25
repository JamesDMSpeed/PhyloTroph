###Load Elton traits data
library(dplyr)
library (rgdal)
library(sp)

#####Load data####
#Mammals
traitDataMam <- read.table("Data/MamFuncDat.txt", sep = '\t', header = T, fill = T)

#Birds
traitDataBird <- read.table("Data/BirdFuncDat.txt", sep = '\t',header = T, fill = T, quote ='')

#Amphibians
traitDataAmph <- read.csv("Data/Traits/AmphiBIO_v1.csv")

#####AMPHIBIANS####
#Read list of species, extracted from shapefile (ampN) Sillero
NamesAmph <- read.csv("Data/NorwayAmphibians.csv")
NamesAmph_L <- levels(NamesAmph$x)
#Make subset from AmphiBIO made up of norwegian amphibians
NorAmph <- traitDataAmph %>%
  filter(Species %in% NamesAmph_L)


#####BIRDS####
#Read list of birds, extracted from birds shapefile
NamesBird <- read.csv("Data/NorwayBirds.csv")
NamesBird_L <- levels(NamesBird$x)
#Making a subset of traitDataBird that consists only of birds found in Norway(from norBirds)
NorBirds <- traitDataBird %>%
  filter(Scientific%in% NamesBird_L)
#There's supposed to be 251 species, based on data from birdN1, but Pinguinus impennis is extinct

##Categorize species
B.invertebrates <- NorBirds[NorBirds$Diet.Inv>=70,] #83

B.plants <- NorBirds[NorBirds$Diet.Seed>=70, ]

#& NorBirds$Diet.Nect >=70
                     #& NorBirds$Diet.Seed >=70
                     #& NorBirds$Diet.PlantO >=70,]

#Check species with not more than 70% of any diet
View(NorTerrSpecies[NorTerrSpecies$Diet.Vertebrate<=70 
                    & NorTerrSpecies$Diet.Invertebrate<=70
                    & NorTerrSpecies$Diet.Plant<=70,])#14


