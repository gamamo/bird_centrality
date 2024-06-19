############################################################
#                                                          #
# Script to run ENM and calculate niche centroid distance  #
#                                                          #
############################################################


# Load packages -----------------------------------------------------------

library(devtools)
remotes::install_github("marlonecobos/ellipsenm")

install.packages("rgdal_1.6-7.tar.gz", repos = NULL, type = 'source')
 
library(rgbif)
library(ellipsenm)
library(ggplot2)
library(terra)
library(letsR)
library(MuMIn)
library(lme4)
library(tidyverse)
library(raster)
library(geodata)
library(here)
library(plyr)


# Load files from Moulatlet et al. 2023 -----------------------------------

load("spp_cen_envi2.RData")

# get spp occurrences from gbif -------------------------------------------
spocc = list()

for(i in unique(spcen_dbdf$Scientific)){
  sp_1 <- occ_search(scientificName = i, limit = 1000)
  sp_1 <- sp_1$data
  sp1_points <- dplyr::select(sp_1,decimalLongitude,decimalLatitude)
  
  #clean the data
  sp1_dups <- duplicated(sp1_points)
  sp1_nodups <- distinct(sp1_points)
  filter(sp1_nodups,decimalLongitude==0)
  filter(sp1_nodups,decimalLatitude==0)
  sp1_points_nonas  <- na.omit(sp1_nodups)
  
  #limit the distance between points to 15 km
  sp1_points_nonas = thin_data(sp1_points_nonas, longitude = "decimalLongitude",
                               latitude = "decimalLatitude", 
                               thin_distance = c(15))
  
  # prepare the data.frame
  species = rep(i, nrow(sp1_points_nonas))
  species = gsub(" ","_",species)
  occurrences = cbind(species, sp1_points_nonas)
  names(occurrences) = c("species","longitude","latitude")
  
  spocc[[i]] = occurrences
  
}

load(file="gbif_occ.RData")


# load raster layers of environmental data --------------------------------
if(F){
clim.current <- list.files(path = "C:/Users/manda/Dropbox/papers/projetos/Moulatlet - bird distribution/layers/30s"
                           , pattern='tif', full.names=TRUE) # replace this line of code by the any that you have

clim.current <- worldclim_global("bio", res=2.5, version="2.1",
                                 path="C:/Users/manda/Dropbox/papers/projetos/Moulatlet - bird distribution/layers")



# to import worldclim variables
r_terra <- rast(clim.current)

r_terra  = r_terra [[c("wc2.1_30s_bio_1","wc2.1_30s_bio_5",
            "wc2.1_30s_bio_6","wc2.1_30s_bio_12",
            "wc2.1_30s_bio_16","wc2.1_30s_bio_17")]]

r_terra  = r_terra [[c("wc2.1_30s_bio_5",
                       "wc2.1_30s_bio_6","wc2.1_30s_bio_13",
                       "wc2.1_30s_bio_14")]]
r_terra <-  scale(r_terra)
r_terra <- terra::aggregate(r_terra, fact=3, fun="mean")
writeRaster(r_terra, "climate2km_scale.tif")
}

r_terra <- rast("Clim_scaled_world.tif")
names(r_terra)
aa <- rast("climate30s_scale.tif")

r_terra  = r_terra [[c("wc2.1_10m_bio_5",
                       "wc2.1_10m_bio_6","wc2.1_10m_bio_13",
                       "wc2.1_10m_bio_14")]]

# calculate and save the convex area around PAM ---------------------------

sppoly = list()

#for (z in names(spocc)) {
  for (z in 1:length(spocc)) {
  
  print(z)
  #temp = as.data.frame(PAM_birds$Presence_and_Absence_Matrix[,c("Longitude(x)","Latitude(y)", z)])
  temp = spocc[[z]]
  #temp = temp[temp[,3]==1,]           
  #temp = temp[,c(3,1,2)]
  #colnames(temp) =c("species","longitude","latitude")
  
  sppoly[[z]] = convex_area(data = temp, longitude = "longitude",
                   latitude = "latitude", buffer_distance = 100)
  
}

save(sppoly,file="sp_poly_resu.RData")

load("sp_poly_resu.RData")


# crop rasters using the polygons and aggregate to 0.5 degree -------------

r_list=list()

for (i in 1:length(sppoly)){
#  for (i in 1:5){
  print(i)
  v <- vect(sppoly[[i]])
  rt =  terra::crop(r_terra,v,mask=TRUE,touches=T)
  r_list[[i]] = rt
  }


s <- vect(spocc[[i]],geom=c("longitude","latitude"))
plot(rt[[1]],add=T)
plot(s)
plot(sppoly[[i]],add=T)

# run the ENM models ------------------------------------------------------

cents = list()


for(i in 1:length(spocc)){
vars <- raster::stack(r_list[[i]])

#set you folder to save the products

ell_model <- ellipsoid_model(data = spocc[[i]], species = "species",
                             longitude = "longitude", latitude = "latitude",
                             raster_layers = vars,
                             method = "covmat", 
                             level = 95,
                             replicates = 1,
                             prediction = "suitability",
                             return_numeric = TRUE, format = "GTiff",
                             tolerance = 1e-60,
                             overwrite = TRUE,
                             output_directory = here("models","covmat3",unique(spocc[[i]]$species)))

ell_model@niche_volume

cents[[i]] = as.vector(ell_model@centroid)
}


save(cents,file = "centroids_sp_resu.RData")
load("centroids_sp_resu.RData")


# calculate mahalabobis distances from the centroids ----------------------
#covmat
spdist_covmat = list()

for(i in 1:length(cents)){
#for(i in 1:10){
  # check if it works for all networks, otherwise skip manually. 
  # some networks are too small in order to allow the distance calculation
  if(length(cents[[i]])!=0){
    
print(i)
sp = spocc[[i]][1,1]
sp = gsub("_"," ",sp);sp
print(sp)

todist <- spcen_dbdf[c("Scientific","PC1","long","lat")]
toext <- todist [todist$Scientific==sp,]
toextv <- vect(toext[,c(1,3,4)],geom=c("long","lat"))
toextv2 <- terra::extract(r_terra,toextv)

toextv2 = toextv2[,c(2:5)]
if(dim(toextv2)[1]<4){envidist <- 0
} else {
  
 
if(F){ 
if(length(which(duplicated(toextv2)))>0){
  for(z in which(duplicated(toextv2))){
   toextv2[z,] <- z/10000000 + toextv2[z,]
  }
}
  

if(length(which(duplicated(cov(toextv2))))>0) {
  for(z in which(duplicated(cov(toextv2)))){
    toextv2[z,] <- z/1000000 + toextv2[z,]
  }
}  
}
  
#envidist  = mahalanobis(toextv2,center=cents[[i]],cov(toextv2),inverted = F,to=1e-23)

ed <- list()  
for (u in 1:nrow(toextv2)){
  ed[[u]] <- dist(rbind(c(cents[[i]]),c(toextv2[u,])))  
}

envidist <- unlist(ed)
spdist_covmat[[i]] <- cbind(spcen_dbdf[spcen_dbdf$Scientific==sp,],envidist)

}
} 
}

# joint spp data and envidist distance

spdist_covmat_df = ldply(spdist_covmat,"data.frame")

save(spdist_covmat_df,file="sp_dists_covmat_v3_resu_euclidian.Rdata")

