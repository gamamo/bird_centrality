############################################################
#                                                          #
# Script to run ENM and calculate niche centroid distance  #
#                                                          #
############################################################


# Load packages -----------------------------------------------------------

library(devtools)
remotes::install_github("marlonecobos/ellipsenm")

library(rgbif)
library(ellipsenm)
library(ggplot2)
library(terra)
library(letsR)
library(sf)
library(plyr)
library(MuMIn)
library(lme4)
library(tidyr)
library(tidyverse)
library(reshape2)
library(raster)


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

clim.current <- list.files(path = "C:/Users/gabri/Dropbox/papers/projetos/Moulatlet - bird distribution/layers/10m"
                           , pattern='tif', full.names=TRUE)
varclimdata<-raster::stack(clim.current)

r_terra <-  rast(varclimdata)

# calculate and save the convex area around PAM ---------------------------

sppoly = list()
ptm <- proc.time()
for (i in names(spocc)) {
  
  print(i)
  temp = as.data.frame(PAM_birds$Presence_and_Absence_Matrix[,c("Longitude(x)","Latitude(y)", i)])
  temp = temp[temp[,3]==1,]           
  temp = temp[,c(3,1,2)]
  colnames(temp) =c("species","longitude","latitude")
  
  sppoly[[i]] = convex_area(data = temp, longitude = "longitude",
                   latitude = "latitude", buffer_distance = 100)
  
}
proc.time() - ptm
save(sppoly,file="sp_poly.RData")

load("sp_poly.RData")



# crop rasters using the polygons and aggregate to 0.5 degree -------------

r_list=list()
ptm <- proc.time()
for (i in 1:length(sppoly)){
  print(i)
  rt =  terra::crop(r_terra,sppoly[[i]], mask=TRUE)
  r_list[[i]] = terra::aggregate(rt,fact=3.125,fun="mean")
}

proc.time() - ptm


# run the ENM models ------------------------------------------------------

names(spocc)

cents = list()

for(i in 1:length(spocc)){

x = stack(r_list[[i]])
vars = x[[c("wc2.1_10m_bio_1","wc2.1_10m_bio_5",
            "wc2.1_10m_bio_6","wc2.1_10m_bio_12",
            "wc2.1_10m_bio_16","wc2.1_10m_bio_17")]]

setwd("C:/Users/gabri/Dropbox/papers/projetos/Moulatlet - bird distribution/models/covmat2")
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
                             output_directory = unique(spocc[[i]]$species))

cents[[i]] = as.vector(ell_model@centroid)
}


save(cents,file = "centroids_sp.RData")
load("centroids_sp.RData")




# calculate mahalabobis distances from the centroids ----------------------
#covmat
spdist_covmat = list()

for(i in 1:length(c_list_covmat)){
  # 16,41,80,129,166,238,240, no funciona
  if(length(c_list_covmat[[i]])!=0){
    
print(i)
sp = names(c_list_covmat[i])
sp = gsub("_"," ",sp);sp
print(sp)

todist = spcen_dbdf[spcen_dbdf$Scientific==sp,c("wc2.1_10m_bio_1","wc2.1_10m_bio_5",
                                               "wc2.1_10m_bio_6",
                                               "wc2.1_10m_bio_16",
                                               "wc2.1_10m_bio_17")]
envidist  = mahalanobis(todist,center=c_list_covmat[[i]],cov(todist),inverted = F,to=1e-23)

# joint spp data and mahalanobis distance

spdist_covmat[[i]] = cbind(spcen_dbdf[spcen_dbdf$Scientific==sp,],envidist)

  }
}


spdist_covmat_df = ldply(spdist_covmat,"data.frame")

save(spdist_covmat_df,file="sp_dists_covmat_v3.Rdata")



