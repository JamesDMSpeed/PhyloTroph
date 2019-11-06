#Plant data

library(raster)
library(sp)
library(rasterVis)

#Diversity maps from Ida's thesis and manuscript (in press at Journal of Biogeography)
plantdivmaps<-read.csv('Data/Ranges/Plants/biodiverse_results_concatenated_good100_utm32_0402_2.csv')

#Norway
norway<-getData('GADM',country='NOR',level=0)
norwayutm<-spTransform(norway,'+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs')
#Raster stack
plantdiv_rast<-rasterFromXYZ(plantdivmaps[,c(2:ncol(plantdivmaps))],
                             crs='+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs')
plantdiv_rast

#Plot SR and PD
levelplot(plantdiv_rast$ENDW_RICHNESS,main='Plant species richness',margin=F,scales=list(draw=F),
          par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayutm,col=grey(0.5)))
levelplot(plantdiv_rast$PD,main='Plant phylogenetic diversity',margin=F,scales=list(draw=F),
          par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayutm),col=grey(0.5))

#Write Rasters
writeRaster(plantdiv_rast$ENDW_RICHNESS,'Data/Ranges/Plants/PlantSpeciesRichness',format='GTiff')
writeRaster(plantdiv_rast$PD,'Data/Ranges/Plants/PlantPhylogeneticDiversity',format='GTiff')
