Citation: Moulatlet, G., Dattilo, W., Kissling, D., Villalobos, F. Bird species’ network centrality varies differentially across species within their climatic niches. 2025. The American Naturalist. in review.

Contact: Gabriel M. Moulatlet (mandaprogabriel@gmail.com) and Fabricio Villalobos (fabricio.villalobos@gmail.com)

Short summary of the article:Understanding how the functional role of species within seed-dispersal networks varies across geographical and climatic gradients can reveal the mechanisms driving network organization. Using data for bird species from all continents, we evaluated the variation of species’ centrality within local networks across species’ climatic niches (occupied climatic conditions) and in response to proxies of competition (number of co-existing bird species) and resource availability (number of co-existing plant species). Taken together, the variation in individual species’ centrality within climatic niches suggests the existence of areas where species achieve high centrality, which might form the substrate for evolutionary and ecological dynamics.

Responsible for collecting data and writing code: Gabriel M. Moulatlet (mandaprogabriel@gmail.com)

Folders and files: There are three folders, data (which contain the data necessary to run the R codes), code (which contain the R scripts) and outputs (where the figures and tables generated when running the scripts will be saved).

Data folder: 
	     - "species_list.csv": contains the column "Scientific" which has species names
             - "gbif_occ.csv": contains three columns. "species" is a list of species names. "latitude" and "longitude" are geographic coordinates in decimal degrees.
             - Folder "convex_areas": contains shapefiles (extension .shp) and auxiliary files (.dbf, .prj, . cpg, .shx). Each file name refers to one species.
	     - "climate_layers.tif": contains 6 climate layers downloaded from worldclim.org. For layer names, please check https://www.worldclim.org/data/worldclim21.html
	     - "network_data.csv": contains four columns. "Scientific" is the list of species names. "PC1" has the centrality values of each species in one ecological network. "latitude" and "longitude" are geographic coordinates in decimal degrees of different ecological networks where a given species occur. For further information please see Moulatlet et al. 2023. Journal of Animal Ecology.
	     - "distance_data.csv": contains seven columns. "scientific" is a list of species names. "family" is the species taxonomic family. "PC1" has the centrality values of each species in one ecological network. "nbirds" is the bird richness in each network the species participate in. "nplants" is the plant richness in each network the species participate in. "envidist" is the Euclidian distance of each network climatic location to the centroid of species climatic niche. "database" has information on the source of each network the species participate in.

Code folder:
	     - Niche_modelling.R: script with 5 parts necessary to run the niche model for each species and to calculate the variable "envidist", the main variable use in the script "Statistical Analyses".
	     - Statistical_Analyses.R: script that generates results, figures and tables.

Output folder: empty. it will be filled as the code is run.

Software version: R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Funding sources: Ciencia Básica project SEP-CONACYT CB-2017-2018 (#A1-S-34563, grant to Fabricio Villalobos). Gabriel M. Moulatet received a postdoctoral grant from this project.
