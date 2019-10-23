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
data(wrld_simpl)
norway <- wrld_simpl[wrld_simpl$NAME == "Norway", ]
plot(norway)

# Load trait data. Remember to always use UTF-8 encoding with PHYLACINE files to avoid weird characters 
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)
View(mam)
# Set factor levels for IUCN status. "EP" is a new status we added to designate species that went extinct in prehistory like Diprotodon  
mam$IUCN.Status.1.2 <- factor(mam$IUCN.Status.1.2, levels=c("EP", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD"))


#Read in all current range data
# Load maps of these marsupial species
maps.current <- paste0("Data/Ranges/Current/", marsupials$Binomial.1.2, ".tif")
maps.current #List of files
r.current <- stack(maps.current)#Stack up
