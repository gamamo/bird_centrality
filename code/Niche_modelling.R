############################################################
#                                                          #
# Script to run ecological niche modelling and calculate   # 
# niche centroid distances                                 #
# Moulatlet et al. submitted to American Naturalist        #
#                                                          #
############################################################


# Load packages -----------------------------------------------------------
library(tidyverse) # CRAN v2.0.0
library(here)      # CRAN v1.0.1
library(devtools)  # CRAN v2.4.5
library(fs)        # CRAN v1.6.5
library(terra)     # CRAN v1.8-10
library(rgbif)     # CRAN v3.8.1
library(spThin)    # CRAN v0.2.0
library(geodata)   # CRAN v0.6-2
library(geometry)  # CRAN v0.5.1

# for further information on how to install this package, see https://github.com/marlonecobos/ellipsenm2
library(ellipsenm) # [github::marlonecobos/ellipsenm2] v0.3.5 


# PART 1: get gbif occurrences --------------------------------------------

# The code below may take a while to run. Please use the the following file
# as alternative entry point

spocc <- read_csv(here("data","gbif_occ.csv"))

# in the need to run the loop below, activate the code by replacing FALSE by TRUE when necessary

# load species list

splist <- read_csv(here("data","species_list.csv"))

if(FALSE){

spocc <-  list() # create an empity list

for(i in unique(splist$Scientific)){ # loop start
  print(i)
  sp_1 <- occ_search(scientificName = i, limit = 1000)
  sp_1 <- sp_1$data
  sp1_points <- sp_1 |> 
    dplyr::select(decimalLongitude,decimalLatitude)
  
  #clean the data by removing duplicates and NAs
  sp1_nodups <- distinct(sp1_points)
  sp1_points_nonas  <- na.omit(sp1_nodups)
  
  #rarefy the distance between points to 15 km
  
  # The function thind_data below is from the elipsenm package and this specific function
  # has not been updated with the package, so it is not working.
  # It is, however, the one originally used to generate the results of this manuscript.
  # We provide below an alternative function from the R package spThin, and our preliminary tests 
  # suggest that they may have differences (see some testing below). For the 
  # sake of reproduction, use the data file ("gbif_occ.csv") provided at
  # the beginning of this section to keep running the script.
  
  if(FALSE){
  sp1_points_nonas <-  thin_data(sp1_points_nonas, longitude = "decimalLongitude",
                               latitude = "decimalLatitude", 
                               thin_distance = c(15))
  }
  
  # rarefy with the spthin package
  sp1_points_nonas$species <- i
  sp1_points_nonas_thin <- thin(sp1_points_nonas,lat.col = "decimalLatitude", long.col ="decimalLongitude",
                           spec.col = "species",thin.par = 15,reps=1,
                           locs.thinned.list.return=T,
                           write.files=F, write.log.file=F,verbose=F)
  sp1_points_nonas_thin <- as.data.frame(sp1_points_nonas_thin)
  
  # prepare the data.frame to save it in a list
  species <-  rep(i, nrow(sp1_points_nonas_thin))
  species <-  gsub(" ","_",species)
  occurrences <-  cbind(species, sp1_points_nonas_thin)
  names(occurrences) <-  c("species","longitude","latitude")
  
  spocc[[i]] = occurrences
  
  # comparison of thinning functions with elipsenm and spthin packages
  spocc_test <- read_csv(here("data","gbif_occ.csv")) |> 
    dplyr::filter(species==gsub(" ","_",i))
  
  print(spocc_test)
  print(occurrences)
  
} # end of the loop
}


# PART 2: calculate and save the convex area around the gbif occ points ---------------------

# The code below use function that have been deprecated. Please follow the steps to create the file need
# to get to the next entry point

# Download the convex hull polygons of each species as shapefile
#list all shp files in the folder
sppoly_list <- dir_ls(here("data","convex_areas"),regexp = "*.shp")

# run the loop to import the file, transform into a SpatVector and assign the 
# species name
sppoly <-list() 
for (i in 1:length(sppoly_list)) {
  sppoly[[i]] <- vect(sppoly_list[i])
  name <- str_extract(sppoly_list[i], regex("[^/]+\\.shp$"))
  names(sppoly[[i]]) <- gsub(".shp","",name)
  print(names(sppoly[[i]]))
}

# ithe loop below calculate and save the convex area around the gbif occ points
# in the need to run the loop below, activate the code by replacing FALSE by TRUE when necessary

if(FALSE){
# load the occurrence data for each selected species
spocc <- read_csv(here("data","gbif_occ.csv"))

sppoly = list() # create an empty list

for (z in unique(spocc$species)) { #start loop
  
  print(z)
  temp <- spocc |> 
    dplyr::filter(species == z)
  
  # the function below is from the elipsenm package and the one used to generate the 
  # original results of this manuscript. But the function is not
  # available in the package update. the convex_area function calculates the convex hull and make a buffer
  # around the convex hull polygon. We provide similar functions below but they do not reproduce our
  # results precisely. We caution about use them for reproducibility.
  
  if(FALSE){
  sppoly[[z]] <-  convex_area(data = temp, longitude = "longitude",
                            latitude = "latitude", buffer_distance = 100)
  }
  # The package geometry provides an alternative function to calculate the convex hull polygons
  convex <-  convhulln(temp[,-1], options="FA")
  convex_vect <- vect(convex$hull)
  
  # the package terra provides an alternative function to calculate the buffer
  convex_vect_buffer <- terra::buffer(convex_vect,100)
  sppoly[[z]] <- convex_vect_buffer
  
} # end loop
}



# PART 3: load raster climatic layers ---------------------------------------------

# this code below may take a while to download. Load the layers using the  next line of code
climate <- geodata::worldclim_global("bio",res=0.5,path=here("data")) 

# download some example climatic layers
climate <- rast(here("data","climate_layers.tif"))

# PART 4: run the ENM models ------------------------------------------------------

# the output of the models used in the article is in the file: "distance_data.csv" inside the data folder
# Please caution that the code below uses the new version of the ellipsenm package and the results may 
# not be the same as ours. The final model output will also depend whether this code has been running with the files we provided
# or if new generated files.

# in the need to run the loop below, activate the code by replacing FALSE by TRUE when necessary
if (FALSE){
  
# to run the model we first crop the climatic layers around the convex area

# load species names

spocc <- read_csv(here("data","gbif_occ.csv")) 

# load species polygons
sppoly_list <- dir_ls(here("data","convex_areas"),regexp = "*.shp")

# run the loop to import the file, transform into a SpatVector and assign the 
# species name
sppoly <-list() 
for (i in 1:length(sppoly_list)) {
  sppoly[[i]] <- vect(sppoly_list[i])
  name <- str_extract(sppoly_list[i], regex("[^/]+\\.shp$"))
  names(sppoly[[i]]) <- gsub(".shp","",name)
  print(names(sppoly[[i]]))
}

#Firs crop the rasters using species polygons

r_list=list() # create a empty list of cropped rasters

for (i in 1:length(sppoly)){
  print(i)
  rt <-   terra::crop(climate,sppoly[[i]], mask=TRUE)
  r_list[[i]] <-terra::aggregate(rt,fact=3.125,fun="mean")
  names(r_list)[i] <- names(sppoly[[i]])
  }

cents = list() # create an empty list to save the ellipse centroids

# run the loop to run the model for each species.
# the following loop will create a folder for each species inside the folder "products".


for (k in unique(spocc$species)){

obs <- spocc |> dplyr::filter(species==k)

ell_model <- ellipsoid_model(data = obs, species = "species",
                             longitude = "longitude", latitude = "latitude",
                             raster_layers = r_list[[k]],
                             method = "covmat", 
                             level = 95,
                             replicates = 1,
                             prediction = "suitability",
                             return_numeric = TRUE, format = "GTiff",
                             tolerance = 1e-60,
                             overwrite = TRUE, 
                             output_directory = here("products",k))

cents[[i]] = as.vector(ell_model@centroid) # this list is used in PART 5
names(cents)[[i]] <- k
}


# PART 5: calculate euclidian distances from the centroids ----------------------------------
# the output of this part used in the article is in the file: "distance_data.csv" located 
# inside the data folder. Use that for running the next parts of the script.

# load ecological networks information
eni <- read_csv(here("data","network_data.csv"))

# start loop
spdist <- list() # create an empty list to save the objects
for(i in 1:length(cents)){
  
  # turn data into a matrix. This is the format necessary to calculate the euclidian distances below
  
  temp <-  unlist(cents)
  names(temp) <- NULL
  
  # some networks are too small in order to allow the distance calculation
  if(length(temp)!=0){
    
    #get the species name
    sp <- names(cents)[i]
    sp <- gsub("_"," ",sp)
    print(sp)
    
    # select the species
    toext <- eni|> 
      dplyr::filter(Scientific==sp)
  
    # convert networks lat and long into a spatvector object
    toextv <- vect(toext[,c(1,3,4)],geom=c("long","lat"))
    
    # extract network climatic values
    toextv2 <- terra::extract(climate,toextv)
    
    # calculate Euclidian distance of each network climatic position to the centroid of each 
    # species climatic model
    
    toextv2 = toextv2[,-1] # delete the column "ID"
    
    if(dim(toextv2)[1]<4){envidist <- 0 # exclude species occurring in less than four networks
    } else {
      
      ed <- list()  
      for (u in 1:nrow(toextv2)){
        ed[[u]] <- dist(rbind(c(temp),c(toextv2[u,])))  
      }
      
      envidist <- unlist(ed)
      spdist[[i]] <- cbind(eni[eni$Scientific==sp,],envidist)
      
    }
  } 
}

# convert the results to a data.frame and save it as .csv file

spdist_df <-  map_df(spdist,~.x) |> 
  write_csv(here("products","euclidian_distances.csv"))
}



