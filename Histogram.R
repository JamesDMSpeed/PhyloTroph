###MAKING SUBSET OF TraitsData w/Norwegian mammals
library(dplyr)

TraitData <- read.csv("Trait_data.csv")
NorMammalsList <- read.csv("NorwayMammalsfromPhylacine.csv")

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

###HISTOGRAMS###
#Frequency
hist(NorSpecies$Diet.Plant,
     main = "Histogram of Plant diet",
     xlab = "Percentage of plant in diet")

#Probability density
hist(NorSpecies$Diet.Plant,
     main = "Histogram of Plant diet",
     xlab = "Percentage of plant in diet",
     prob = TRUE)


hist(NorSpecies$Diet.Invertebrate)

hist(NorSpecies$Diet.Vertebrate)

library(MASS)
truehist(NorSpecies$Diet.Plant)
lines(density(NorSpecies$Diet.Plant))
