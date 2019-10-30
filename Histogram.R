###MAKING SUBSET OF TraitsData w/Norwegian mammals
library(dplyr)

TraitData <- read.csv("Data/Traits/Trait_data.csv")
NorMammalsList <- read.csv("Data/NorwayMammalsfromPhylacine.csv")

#Removing the number and empty space in front of the species names in NorMammalsList
MammalNames <- gsub('[[:digit:]]+', '', NorMammalsList$x)
MammalNames <- substring(MammalNames,2)

#Making a subset of TraitData that consists only of mammals found in Norway(from MammalNames)
NorSpecies <- TraitData %>%
  filter(Binomial.1.2%in% MammalNames)

#MammalNames <- as.character(NorMammals$x)
#NorMammals <- gsub('[[:digit:]]+', '', NorMammals$x)
#NorMammals <- data.frame(NorMammals)
#colnames(NorMammals)[colnames(NorMammals)=="NorMammals"] <- "Binomial.1.2"
#NamesM <- as.character(NorMammals$Binomial.1.2)
#NamesM <-substring(NamesM,2)

#Only terrestrial species (omit also terrestrial + marine spp e.g. seals)
NorTerrSpecies<-NorSpecies[NorSpecies$Marine==0,]
dim(NorTerrSpecies)#58 species
###HISTOGRAMS###
#Frequency
hist(NorTerrSpecies$Diet.Plant,
     main = "Histogram of Plant diet",
     xlab = "Percentage of plant in diet")

#Probability density
hist(NorTerrSpecies$Diet.Plant,
     main = "Histogram of Plant diet",
     xlab = "Percentage of plant in diet",
     prob = TRUE)


hist(NorTerrSpecies$Diet.Invertebrate)

hist(NorTerrSpecies$Diet.Vertebrate)

library(MASS)
truehist(NorTerrSpecies$Diet.Plant)
lines(density(NorTerrSpecies$Diet.Plant))
