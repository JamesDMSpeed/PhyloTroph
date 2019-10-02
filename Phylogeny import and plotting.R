library(ape)
library(picante)


#Norway plant phylogeny

plantphylo<-read.tree('NorwVascPlantPhylogeny.nwk')
plot(plantphylo)
plot(plantphylo,show.tip.label = F,main='Norway native vascular plant phylogeny')
