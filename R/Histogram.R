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

#Only terrestrial species (omit also seals based on taxonomy)
NorTerrSpecies<-NorSpecies[NorSpecies$Family.1.2!='Phocidae' & NorSpecies$Terrestrial==1,]
dim(NorTerrSpecies)#47 species
NorTerrSpecies[,1:3] # Note this includes humans!

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

#Check species with not more than 60% of any diet
View(NorTerrSpecies[NorTerrSpecies$Diet.Vertebrate<=60 
               & NorTerrSpecies$Diet.Invertebrate<=60
               & NorTerrSpecies$Diet.Plant<=60,]) #3
#Check species with not more than 50% of any diet
View(NorTerrSpecies[NorTerrSpecies$Diet.Vertebrate<=50 
                                & NorTerrSpecies$Diet.Invertebrate<=50
                                & NorTerrSpecies$Diet.Plant<=50,])#1
#Check species with not more than 70% of any diet
View(NorTerrSpecies[NorTerrSpecies$Diet.Vertebrate<=70 
                                & NorTerrSpecies$Diet.Invertebrate<=70
                                & NorTerrSpecies$Diet.Plant<=70,])#14


library(MASS)
truehist(NorTerrSpecies$Diet.Plant)
lines(density(NorTerrSpecies$Diet.Plant))
