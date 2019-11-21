#Script to crop out polygons of birds outside Norway
#Save and read into R via box, avoiding large file sizes.


library(raster)
library(rgdal)
library(utils)

#Read data from BOTW From server, and crop to Norway
#This section is commented out

# #BOTW
# fc_list = ogrListLayers(path.expand('~/bio3D/BOTW/BOTW_2016.gdb'))
# print(fc_list)
# botw<-readOGR(dsn=path.expand('~/bio3D/BOTW/BOTW_2016.gdb'),layer="All_Spp")
# 
 norway<-getData('GADM',country='NOR',level=0)
# 
# bon<-botw[norway,]
# writeOGR(bon,dsn='Data/Ranges/Birds',layer='BirdsofNorway',driver='ESRI Shapefile')
# 
# #Save as a zip file
# zip('Data/Ranges/Birds/BirdsofNorway',files=list.files('Data/Ranges/Birds',pattern = 'BirdsofNorway*',full.names = T))
# #Upload this to a box, and get direct download link

#Download
birdtempfile<-tempfile()
download.file('https://ntnu.box.com/shared/static/bg0bzgzh6jokiji53gvgbxcvnzp0j8kk.zip',birdtempfile)
unzip(birdtempfile)
unlink(birdtempfile)

#Read as spatial polygons
birdsnorway<-readOGR('Data/Ranges/Birds/','BirdsofNorway')
birdsnorway
#251 species - may include species not overlapping with Norway land (only cropped to extent, not to outline)
#Also includes vagrant species - these will most likely be excluded from further analyses
#See here for further information http://datazone.birdlife.org/species/spcdistPOS

#List species
levels(birdsnorway@data$SCINAME)
tapply(birdsnorway@data$SCINAME,birdsnorway@data$ORIGIN,length)#3 introduced and 1 vagrant species

plot(birdsnorway[birdsnorway$SCINAME=='Corvus corax',])#Whole range of species overlapping with Norway. Not a problem when we come to rasterize the data.
