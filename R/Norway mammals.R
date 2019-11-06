#Norwegian mammals

#Load packages
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

#Load phylogeny (use 50 to average)
forestphy <- read.nexus("Data/Phylogenies/Complete_phylogeny.nex")
names(forestphy) <- NULL
set.seed(42)
forestphy <- forestphy[sample(1:1000, 50)]

#Norway map
data(wrld_simpl)
norway <- wrld_simpl[wrld_simpl$NAME == "Norway", ]
plot(norway)

# Load trait data. Remember to always use UTF-8 encoding with PHYLACINE files to avoid weird characters 
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)

# Set factor levels for IUCN status. "EP" is a new status we added to designate species that went extinct in prehistory like Diprotodon  
mam$IUCN.Status.1.2 <- factor(mam$IUCN.Status.1.2, levels=c("EP", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD"))

# Load maps of all species
maps.current <- paste0("Data/Ranges/Current/", mam$Binomial.1.2, ".tif")
r.current <- stack(maps.current)
maps.pres.nat <- paste0("Data/Ranges/Present_natural/", mam$Binomial.1.2, ".tif")
r.pres.nat <- stack(maps.pres.nat)

# Project Australia map to the range map projection
norway <- spTransform(norway, crs(r.current))

# Crop range maps to just the extent of Australia for a cleaner plot
ext <- extent(norway)
r.current <- crop(r.current, ext)
r.pres.nat <- crop(r.pres.nat, ext)

# Create a blank raster of the region
blank.r <- r.current[[1]]
blank.r[] <- NA
names(blank.r) <- NA

# Load all the current raster data into a matrix for faster handling
m.current <- matrix(NA, nrow=nrow(mam), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.current) <- mam$Binomial.1.2
for(i in 1:nrow(mam)) {
  m.current[i, ] <- getValues(r.current[[i]])
}
# Current species list
sp.current <- mam$Binomial.1.2[rowSums(m.current) > 0]

# Load all the present natural raster data into a matrix for faster handling
m.pres.nat <- matrix(NA, nrow=nrow(mam), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.pres.nat) <- mam$Binomial.1.2
for(i in 1:nrow(mam)) {
  m.pres.nat[i, ] <- getValues(r.pres.nat[[i]])
}
# Present natural species list
sp.pres.nat <- mam$Binomial.1.2[rowSums(m.pres.nat) > 0]

# Create rasters of taxonomic richness
current.div <- blank.r
current.div[] <- colSums(m.current)
names(current.div) <- "Current diversity"
pres.nat.div <- blank.r
pres.nat.div[] <- colSums(m.pres.nat)
names(pres.nat.div) <- "Present natural diversity"
div <- stack(current.div, pres.nat.div)
# Change 0 to NA, for nicer plotting
div[] <- ifelse(div[] == 0, NA, div[])
```
