# Bird centrality

## This is the repository for the article "Bird species’ centrality in seed-dispersal networks varies within climatic niches" by Moulatlet et al.

There are two scripts in these repository that should be run in the specific order:

1) *script_enm_2024_01_10* Contains the steps for running the envelop niche models for each focal species of birds. It starts by downloading occurrences from a pre-select species list. Then, convex hull polygons are calculated and ENM are run. At the end, we calculate the mahalanobis distance from each network where the species occur to the centroid of the niche.

2) *script_mrm_2024_01_10* Contains the steps to run the multiple regression models (and model selection) for each focal species. It also has the codes to make all the figures in the manuscript.

The data necessary to run the scripts is included in the repository as well. Please note that the running the codes in the *script_enm_2024_01_10* may be computationally intensive and take several hours. As such, the intermediate files (the *.Rdata files) produced are made available in here too have loading functions wherever necessary in the scripts.
