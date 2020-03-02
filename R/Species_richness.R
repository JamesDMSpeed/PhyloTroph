### Making species richness maps

library(pacman)
pacman::p_load(raster,
               rgdal,
               sf,
               fasterize,
               ggplot2,
               leaflet,
               tmap,
               dplyr,
               tiff,
               sp,
               gdalUtils,
               rgeos,update = F)

####READING SHAPEFILES####
  #Reading shapefile for mammal
  mamN <- readOGR("Data/Ranges/Mammals")#Should be read from repository project directory, not local drive
  mamN1 <- mamN
  ##Reptiles
  #repN <- readOGR("Reptiles_","NorwayReptiles")
  repN_Sillero <- readOGR("Reptiles_","reptiles_Norway_Sillero")
  repN1 <- repN_Sillero
  ##Amphibians
  #amphN <- readOGR("Amphibians_","amphibian_Norway_IUCN")
  amphN_Sillero <- readOGR("Amphibians_","amphibian_Norway_Sillero")
  amphN1 <- amphN_Sillero
  #Plants
  #plantList <- read.csv("Data/Ranges/Plants/biodiverse_results_concatenated_good100_utm32_0402_2.csv")
  plantSR <- raster("Data/Ranges/Plants/PlantSpeciesRichness.tif")
  plantN1 <- plantSR
  #Birds
  birdN <- readOGR("Birds_","Birbies_fix")
  #norBirds <- levels(birdN@data$SCINAME) #No. of species
  birdN1 <- birdN
  birdN1 <- birdN1[-461,]
  #Removing a small polygon from the species of Calidris maritima (small island)
  #birdN1@data$REVIEWE <- NULL

#####SR PLANTS####
  #Change crs for the plant raster (resolution is correct)
  #crs(plantN1) <- "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"
  #plantTrans <- projectRaster(plantN1, crs="+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
  plantN2 <-plantN1
  projection(plantN2) <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
  richnessMap_plant <- plantN2

####SR MAMMALS####
##Subset from species from Artsdatabanken
#NorTerrMam <- read.csv("Data/Rodlista_terrestrisk.csv")
#NorList <- levels(NorTerrMam$Vitenskapelig_navn)
#NorMam <- subset(mamN1,mamN1@data$BINOMIAL%in% NorList)

#Transformation between datum(s) and conversion between projections, crs = coordinate reference system
mamTrans1<-spTransform(mamN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e<-extent(mamTrans1)
#Creating a RasterLa<yer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s<-raster(e, resolution=20000, crs=(mamTrans1))
mamSF <- st_as_sf(mamTrans1)
#one layer per species (sf, raster)
stackedDistributionMaps<-fasterize(mamSF,s,by="BINOMIAL")
plot(stackedDistributionMaps)

#merge polygons per species for species richness map, unioning geometries 
speciesPoly<-aggregate(mamTrans1,by="BINOMIAL")
#plot(speciesPoly)
#speciesPoly_mam <- aggregate(mamSF, by= "BINOMIAL")
#Evaluation error: TopologyException: found non-noded intersection between LINESTRING
speciesPoly2<-st_as_sf(speciesPoly)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_mam <- fasterize(speciesPoly2,s,fun="count",field="BINOMIAL")
plot(richnessMap_mam)
dim(richnessMap_mam)

####SR AMPHIBIANS####
#Transformation between datum(s) and conversion between projections
amphTrans1<-spTransform(amphN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_amph<-extent(amphTrans1)

#make raster and convert to sf-object 
s_amph<-raster(e_amph, resolution=20000, crs=(amphTrans1))
amphSF <- st_as_sf(amphTrans1)
#one layer per species (sf, raster)
#stackedDistributionMaps_amph<-fasterize(amphSF,s_amph,by="amp")
#plot(stackedDistributionMaps_amph)
#merge polygons per species for species richness map, unioning geometries 
speciesPoly_amph<-aggregate(amphTrans1,by="amp")
speciesPoly2_amph<-st_as_sf(speciesPoly_amph)

#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_amph <- fasterize(speciesPoly2_amph,s_amph,field="amp")
plot(richnessMap_amph)
#table(getValues(richnessMap_amph))  

####SR REPTILES####
repTrans1 <-spTransform(repN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_rep<-extent(repTrans1)

#Creating a RasterLayer object, created from an Extent object
#st_as_sf: convert object to sf-object (dataframe w/attributes and geometry)
s_rep<-raster(e_rep, resolution=20000, crs=(repTrans1))
repSF <- st_as_sf(repTrans1)

#merge polygons per species for species richness map, unioning geometries 
speciesPoly_rep<-aggregate(repTrans1,by="rep")
speciesPoly2_rep<-st_as_sf(speciesPoly_rep)
#count number of overlapping polygons (which is 1 per species so counting gives richness)
richnessMap_rep <- fasterize(speciesPoly2_rep,s_rep,field="rep")
plot(richnessMap_rep)

####SR BIRDS####
birdTrans1<-spTransform(birdN1,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
e_bird<-extent(birdTrans1)
s_bird<-raster(e_bird, resolution=20000, crs=(birdTrans1))
birdSF <- st_as_sf(birdTrans1)

#one layer per species (sf, raster)
stackedDistributionMaps_bird<-fasterize(birdSF,s_bird,by="SCINAME")
plot(stackedDistributionMaps_bird)

#merge polygons per species for species richness map, unioning geometries 
#ERROR!!
speciesPoly_bird<-raster::aggregate(birdTrans1,by="SCINAME") #ERROR!
speciesPoly2_bird<-st_as_sf(speciesPoly_bird)

richnessMap_bird <- fasterize(speciesPoly2_bird,s_bird,fun="count",field="SCINAME")
plot(richnessMap_bird)

####SR TOTAL####
r1 <- richnessMap_plant
r2 <- richnessMap_rep
r3 <- richnessMap_amph
r4 <- richnessMap_mam
r5 <- richnessMap_bird

#Change NA to NULL values, so the cells can be summarized in overlay and mosaic function
#This creates grey background to the map (since 0 = grey color), but that will be possible to fix
#By not assigning a color to value 0 when plotting
r1 <- reclassify(r1, cbind(NA, 0))
r2 <- reclassify(r2, cbind(NA, 0))
r3 <- reclassify(r3, cbind(NA, 0))
r4 <- reclassify(r4, cbind(NA, 0))
r5 <- reclassify(r5, cbind(NA, 0))

#Extending all of the rasterlayers to the dimension of the plantraster, and then aligning extent
r1 <- extend(r1, c(0,3), value = NA)
r2 <- extend(r2, c(1,0), value = NA)
r3 <- extend(r3, c(1,0), value = NA)
r4 <- extend(r4, c(1,0), value = NA)
r5 <- extend(r5, c(1,0), value = NA)

r2 <- setExtent(r2,r1, keepres = T)
r3 <- setExtent(r3,r1, keepres = T)
r4 <- setExtent(r4,r1, keepres = T)
r5 <- setExtent(r5,r1, keepres = T)

raster_result <- overlay(r1,r2,r3,r4,r5,fun=function(x,y,z,a,b){return(x+y+z+a+b)})
plot(raster_result)
raster_result
spplot(raster_result)
table(getValues(raster_result))

my.stack <- stack(r1,r2,r3,r4,r5)
plot(my.stack)
table(getValues(my.stack))

my.merge <- raster::merge(r1,r2,r3,r5, tolerance = 0)
plot (my.merge)
table(getValues(my.merge))

my.mosaic <- mosaic(r1,r2,r3,r4,r5, tolerance = 0, fun = sum)
plot(my.mosaic)
table(getValues(my.mosaic))
spplot(my.mosaic)

#####Empty raster####
#Make empty rasterlayer which cover the extents of all of the layers (r6)
# ext1 <- extent(r1)
# ext5 <- extent(r5)
# extUni <- union(ext1, ext5)
# rasUni <- raster()
# rasUni <- setExtent(rasUni, extUni)
# res(rasUni) <- 20000
# values(rasUni)<- 0
# projection(rasUni) <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
# r6 <- rasUni

#####Norway raster####
#norway <- getData('GADM', country = 'norway', level = 1)
#norwayTrans1<-spTransform(norway,crs("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"))
#e_nor<-extent(norwayTrans1)
#s_nor<-raster(e_nor, resolution=20000, crs=(norwayTrans1))
#values(s_nor)<-1

#####RESAMPLE####   
#Resample to Plant raster
# r2.plant <- resample(r2,r1,method = "ngb")
# r3.plant <- resample(r3,r1,method = "ngb")
# r4.plant <- resample(r4,r1,method = "ngb")
# r5.plant <- resample(r5,r1,method = "ngb")
# 
# resample_plant <- c(r2.plant, r3.plant, r4.plant,r5.plant)
# plant.stack <- stack(r1,r2.plant, r3.plant, r4.plant,r5.plant, fun = sum)
# plot(plant.stack)

# #Resample to empty raster
# r1.new <- resample(r1,r6,method = "ngb")
# r2.new <- resample(r2,r6,method = "ngb")
# r3.new <- resample(r3,r6,method = "ngb")
# r4.new <- resample(r4,r6,method = "ngb")
# r5.new <- resample(r5,r6,method = "ngb")
# 
# new.stack <- stack(r1.new,r2.new, r3.new, r4.new,r5.new, fun = sum)
# plot(new.stack)

#####Function 1####
reproject_align_raster<- function(rast, ref_rast=NULL, desired_origin, desired_res, desired_crs, method= "bilinear"){
  
  if (!is.null(ref_rast)) {
    desired_origin<- origin(ref_rast) #Desired origin
    desired_res<- res(ref_rast) #Desired resolution
    desired_crs<- crs(ref_rast) #Desired crs
  } #Set parameters based on ref rast if it was supplied
  if(length(desired_res)==1){
    desired_res<- rep(desired_res,2)}
  
  if(identical(crs(rast), desired_crs) & identical(origin(rast), desired_origin) & identical(desired_res, res(rast))){
    message("raster was already aligned")
    return(rast)} #Raster already aligned
  
  if(identical(crs(rast), desired_crs)){
    rast_orig_extent<- extent(rast)} else{
      rast_orig_extent<- extent(projectExtent(object = rast, crs = desired_crs))} #reproject extent if crs is not the same
  var1<- floor((rast_orig_extent@xmin - desired_origin[1])/desired_res[1])
  new_xmin<-desired_origin[1]+ desired_res[1]*var1 #Calculate new minimum x value for extent
  var2<- floor((rast_orig_extent@ymin - desired_origin[2])/desired_res[2])
  new_ymin<-desired_origin[2]+ desired_res[2]*var2 #Calculate new minimum y value for extent
  n_cols<- ceiling((rast_orig_extent@xmax-new_xmin)/desired_res[1]) #number of cols to be in output raster
  n_rows<- ceiling((rast_orig_extent@ymax-new_ymin)/desired_res[2]) #number of cols to be in output raster
  new_xmax<- new_xmin+(n_cols*desired_res[1]) #Calculate new max x value for extent
  new_ymax<- new_ymin+(n_rows*desired_res[2]) #Calculate new max y value for extent
  rast_new_template<- raster(xmn=new_xmin, xmx =new_xmax,  ymn=new_ymin, ymx= new_ymax, res=desired_res, crs= desired_crs) #Create a blank template raster to fill with desired properties
  if(!identical(desired_origin,origin(rast_new_template))){
    message("desired origin does not match output origin")
    stop()} #Throw error if origin doesn't match
  if(identical(crs(ref_rast),crs(rast))){
    rast_new<- raster::resample(x=rast, y=rast_new_template, method = method)} else{
      rast_new<- projectRaster(from=rast, to=rast_new_template, method = method)} #Use projectRaster if crs doesn't match and resample if they do
  if(!identical(desired_origin,origin(rast_new))){
    message("desired origin does not match output origin")
    stop()} #Throw error if origin doesn't match
  return(rast_new)
}

#####Function_2####
combine_rasters<- function(raster_list, ref_rast=NULL, desired_origin, desired_res, desired_crs, method= "bilinear", display_progress=TRUE){
  raster_list2<- vector("list", length = length(raster_list))
  for (i in 1:length(raster_list)) {
    if(display_progress){
      print(paste("Reprojecting", as.character(i), "of", as.character(length(raster_list))))}
    raster_list2[[i]]<- reproject_align_raster(raster_list[[i]], ref_rast=ref_rast, desired_origin=desired_origin, desired_res=desired_res, desired_crs=desired_crs, method= method)}
  for (j in 1:length(raster_list2)) {
    if(display_progress){
      print(paste("Combining raster", as.character(j), "of", as.character(length(raster_list2))))}
    if(j==1){
      combined_raster<- raster_list2[[j]]} else{
        combined_raster<- raster::merge(combined_raster, raster_list2[[j]])}}
  return(combined_raster)
}



#####Realigning and MERGE w/ function 1/2#### 
#MERGE
#By using two functions that are supposed to align the rasters, and then merge them together
des_or <- c(0,0)
des_res <- 20000
des_crs <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")

align_list <- c(r1,r2,r3,r4,r5)
#Order matters, since the method is 'merge'
align_list2 <- c(r1,r2,r3,r5)

#After aligning all of the layers have the same origin, but different extent and dimension
#Layer 2,3,4 partially overlap

##Varying between ref_rast = emp_ras/NULL and method = bilinear/ngb
#comb_1 <- combine_rasters(rast_list, ref_rast = emp_ras,des_or,des_crs,
#method = "bilinear", display_progress = F) #Decimals
#comb_3 <- combine_rasters(rast_list, ref_rast = NULL,des_or,des_crs,
#method = "bilinear", display_progress = T) #Error; attempt to replicate 'S4'
#comb_4 <- combine_rasters(rast_list, ref_rast = NULL,des_or,des_crs,
#method = "ngb", display_progress = F) #Error; attempt to replicate

comb_ras <- combine_rasters(align_list, ref_rast = r6,des_or,des_crs,
                            method = "ngb", display_progress = F) #Integers
comb_ras2 <- combine_rasters(align_list2, ref_rast = r6,des_or,des_crs,
                             method = "ngb", display_progress = F) #Integers
plot(comb_ras)
spplot(comb_ras)
plot(comb_ras2)
spplot(comb_ras2)


