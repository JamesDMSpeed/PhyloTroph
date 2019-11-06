library(ape)
library(picante)
library(rgbif)

#Norway plant phylogeny

plantphylo<-read.tree('NorwVascPlantPhylogeny.nwk')
plot(plantphylo)
plot(plantphylo,show.tip.label = F,main='Norway native vascular plant phylogeny')

#Download Norway plant data from GBIF
norway_code <- isocodes[grep("Norway", isocodes$name), "code"]
occ_count(country=norway_code,georeferenced = T)
key<-name_backbone(name='Trachaeophyta')$phylumKey
occ_count(country=norway_code,georeferenced = T,taxonKey = key)

#Read data from download of all Norwegian plants
#
occ_download_get(key="0008637-190918142434337", overwrite = TRUE) %>%
  occ_download_import



