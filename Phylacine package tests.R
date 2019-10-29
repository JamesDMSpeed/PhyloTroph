#Testing phylacine data 

#Script modified from readme on https://github.com/MegaPast2Future/PHYLACINE_1.2 


#Loading packages
install.packages("pacman", repos="https://cloud.r-project.org")
pacman::p_load(ggplot2,
               dplyr,
               stringr,
               gridExtra,
               viridisLite,
               raster,
               rasterVis,
               rgdal,
               maptools,
               ape,
               ggtree, update = F)


#Phylogeny from phylacine from box direct link
phy<- read.nexus('https://ntnu.box.com/shared/static/rgoiukue6ywcbdfijsucg8nd34idqppi.nex')
phy
#Select a sample of these 1000 trees
set.seed(42)
phys <- phy[sample(1:1000, 50)]

# Load simple world map and subset Norway (we can use a better outline later)
norway <- getData('GADM',country='Norway',level=0)
plot(norway)

# Load trait data. Remember to always use UTF-8 encoding with PHYLACINE files to avoid weird characters 
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)
View(mam)
# Set factor levels for IUCN status. "EP" is a new status we added to designate species that went extinct in prehistory like Diprotodon  
mam$IUCN.Status.1.2 <- factor(mam$IUCN.Status.1.2, levels=c("EP", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD"))

# Subsetting current range data for Norway
# #Scripted out
# #Read in all current range data
# # Load maps of these marsupial species
# maps.current <- paste0("Data/Ranges/Current/", mam$Binomial.1.2, ".tif")
# maps.current #List of files #Check they are all there
# r.current <- stack(maps.current)#Stack up
# 
# #Transform norway to raster data crs
# norway1<-spTransform(norway,crs(r.current))
# 
# #Crop distribution maps to Norway
# norwayranges<-crop(r.current,norway1)
# 
# #Disaggregate to smaller cells
# norwayrangesD<-disaggregate(norwayranges,fact=10)
# 
# #Mask distribution maps to Norway
# norwayrangesm<-mask(norwayrangesD,norway1)
# #Subset layers with presnece within Norway
# norwayspranges<-norwayrangesm[[which(cellStats(norwayrangesm,max)>0)]]
# 
# #List mammals in Norway
# names(norwayspranges)
# #Plot richness
# plot(calc(norwayspranges,sum))
# 
##Equal area projection
#norwaymammals<-projectRaster(norwayspranges,res=10000,method='ngb',crs='+proj=laea +lat_0=62 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

#writeRaster(norwaymammals,format='GTiff',filename=paste('Data/Ranges/CurrentNorway/',names(norwaymammals)),bylayer=T)

#Read in Current distributions of Norway mammals
#Subsetted from 100km rasters - will be better to directly rasterize IUCN data to higher resolution
norwaymammalsList<-list.files('Data/Ranges/CurrentNorway',full.names=T)
norwaymammals<-stack(norwaymammalsList)
#Clean names
names(norwaymammals)<-substring(names(norwaymammals),3)

#In here need to select only terrestrial mammals

#Diversity maps
norwaymammalSR<-calc(norwaymammals,sum)
levelplot(norwaymammalSR,margin=F,scales=list(draw=FALSE))
#May be better to get the data direct from IUCN