#Reptiles

#data from Meiri, Shai et al. (2017), Data from: The global distribution of tetrapods reveals a need for targeted reptile conservation, Dryad, Dataset, https://doi.org/10.5061/dryad.83s7k

library(rgdal)
library(raster)

#Read reptile distribution shapefile from T drive (not pushed to GitHub)
reps<-readOGR('T:\\vm\\inh\\botanisk\\Bruker\\BOTW\\doi_10.5061_dryad.83s7k__v1\\GARD1.1_dissolved_ranges','modeled_reptiles')
reps

#Download map of Norway
norway<-getData('GADM',level=0,country='NOR')

#Crop reptile distributions to Norway outline
norreps<-reps[norway]
#Remove spp absent from Norway from level list
norreps$Binomial<-droplevels(norreps$Binomial)

#Check species list
norreps$Binomial #This list of 6 is consistant with Dolman 2018 Norske amfibier og reptiler
#Except for Lacerta agilis which Dolman does not have present in Norway

#Plot some spp
plot(norway)
plot(norreps[norreps$Binomial=='Lacerta agilis',],add=T,col=2)#Only just in Norway - can discuss this with Dag

plot(norway)
plot(norreps[norreps$Binomial=='Vipera berus',],add=T,col=2)
plot(norreps[norreps$Binomial=='Natrix natrix',],add=T,col=3)

#Write
writeOGR(norreps,'Data/Ranges/Reptiles','NorwayReptiles',driver="ESRI Shapefile")



