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

#####MAMMALS####
NamesMam <- mamN1@data[["BINOMIAL"]]
#Make subset from traitDataMam which contains species from IUCN ranges
NorMam <- traitDataMam %>%
  filter(Scientific %in% NamesMam)

#Diets based on 3 categories from phylacine
M.DietInvert <- NorMam[NorMam$Diet.Inv>=70,]#17

M.DietVertebrate <- rbind((NorMam[NorMam$Diet.Vend >=70,]),
                          (NorMam[NorMam$Diet.Vect >=70,]),
                          (NorMam[NorMam$Diet.Vfish>=70,]),
                          (NorMam[NorMam$Diet.Vunk>=70]),
                          (NorMam[NorMam$Diet.Scav>=70])) #10

M.DietPlant <- rbind((NorMam[NorMam$Diet.PlantO >=70,]),
                     (NorMam[NorMam$Diet.Seed >=70,]),
                     (NorMam[NorMam$Diet.Nect>=70,]),
                     (NorMam[NorMam$Diet.Fruit>=70])) #12

M.DietOmni <-  NorMam%>%
  filter((!(Scientific %in% M.DietPlant$Scientific)) &
           !(Scientific %in% M.DietVertebrate$Scientific)&
           !(Scientific %in% M.DietInvert$Scientific)) #18


#####M.DIETS 10 CATEGORIES####
M.Invert <- NorMam[NorMam$Diet.Inv>=70,] #17

M.Vend <- NorMam[NorMam$Diet.Vend >=70,] #6
M.Vect <- NorMam[NorMam$Diet.Vect >=70,] #0
M.piscivore <- NorMam[NorMam$Diet.Vfish >=70,] #4
M.Vunknown <- NorMam[NorMam$Diet.Vunk >=70,]#0
M.scav <- NorMam[NorMam$Diet.Scav >=70, ] #0

M.Fruit <- NorMam[NorMam$Diet.Fruit >=70,] #0
M.Nect <- NorMam[NorMam$Diet.Nect >=70,] #0
M.Seeds <- NorMam[NorMam$Diet.Seed >=70,] #0
M.PlantO <- NorMam[NorMam$Diet.PlantO >=70,] #12

M.omnivores <- NorMam[NorMam$Diet.Inv<70&
                          NorMam$Diet.Vend<70&
                          NorMam$Diet.Vect<70&
                          NorMam$Diet.Vfish<70&
                          NorMam$Diet.Vunk<70&
                          NorMam$Diet.Scav<70&
                          NorMam$Diet.Fruit<70&
                          NorMam$Diet.Nect<70&
                          NorMam$Diet.Seed<70&
                          NorMam$Diet.PlantO<70,] #18

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

#####B.DIETS 10 CATEGORIES (ELTON)####
##Categorize species, according to 10 categories given by Elton traits (HotBotW, mostly):
B.invertebrates <- NorBirds[NorBirds$Diet.Inv>=70,] #83

B.Vend <- NorBirds[NorBirds$Diet.Vend >=70,] #19
B.Vect <- NorBirds[NorBirds$Diet.Vect >=70,] #0
B.piscivore <- NorBirds[NorBirds$Diet.Vfish >=70,] #19
B.Vunknown <- NorBirds[NorBirds$Diet.Vunk >=70,]#0
B.scav <- NorBirds[NorBirds$Diet.Scav >=70, ] #0

B.Fruit <- NorBirds[NorBirds$Diet.Fruit >=70,] #0
B.Nect <- NorBirds[NorBirds$Diet.Nect >=70,] #0
B.Seeds <- NorBirds[NorBirds$Diet.Seed >=70,] #4
B.PlantO <- NorBirds[NorBirds$Diet.PlantO >=70,] #17

Categorized <- rbind(B.invertebrates,B.Vend,B.Vect,B.piscivore,B.Vunknown,
                     B.scav,B.Fruit,B.Nect,B.Seeds,B.PlantO)

B.omnivores <- NorBirds[NorBirds$Diet.Inv<70&
                          NorBirds$Diet.Vend<70&
                          NorBirds$Diet.Vect<70&
                          NorBirds$Diet.Vfish<70&
                          NorBirds$Diet.Vunk<70&
                          NorBirds$Diet.Scav<70&
                          NorBirds$Diet.Fruit<70&
                          NorBirds$Diet.Nect<70&
                          NorBirds$Diet.Seed<70&
                          NorBirds$Diet.PlantO<70,] #108

#Checking that all species have been placed in a category 
tot_cat <- rbind(Categorized,B.omnivores)
miss <- setdiff(NamesBird$x, tot_cat$Scientific)

#####B.DIETS 4 CATEGORIES####
#Categorizing according to how the mammals are categorized in Phylacine (3 cat)
B.DietInvertebrate <- NorBirds[NorBirds$Diet.Inv>=70,] #83

B.DietVertebrate <- rbind((NorBirds[NorBirds$Diet.Vend >=70,]),
                     (NorBirds[NorBirds$Diet.Vect >=70,]),
                     (NorBirds[NorBirds$Diet.Vfish>=70,]),
                     (NorBirds[NorBirds$Diet.Vunk>=70]),
                     (NorBirds[NorBirds$Diet.Scav>=70])) #38

B.DietPlant <- rbind((NorBirds[NorBirds$Diet.PlantO >=70,]),
                     (NorBirds[NorBirds$Diet.Seed >=70,]),
                     (NorBirds[NorBirds$Diet.Nect>=70,]),
                     (NorBirds[NorBirds$Diet.Fruit>=70])) #21

B.DietOmni <-  NorBirds%>%
  filter((!(Scientific %in% B.DietPlant$Scientific)) &
           !(Scientific %in% B.DietVertebrate$Scientific)&
           !(Scientific %in% B.DietInvertebrate$Scientific)) #108

N1Species <- mamN1@data[["BINOMIAL"]]
N2Species <- NorTerrSpecies[["Binomial.1.2"]]

