# Bird centrality

## This is the repository for the article "Bird species’ centrality in seed-dispersal networks varies within climatic niches" by Moulatlet et al.

There are two scripts in these repository that should be run in the specific order:

1) *script_enm_2024_01_10* Contains the steps for running the envelop niche models for each focal species of birds. It starts by downloading occurrences from a pre-select species list. Then, ENM are run. At the end, it calculates the Euclidian distance from each network where the species occur to the centroid of the niche.

2) *script_mrm_mixed_2024_05_20* Contains the steps to run the mixed models for each focal species. It also has the codes to make all the figures in the manuscript.

The data necessary to run the scripts is included in the repository as well. Please note that the running the codes in the *script_enm_2024_01_10* may be computationally intensive and take several hours.

### Data description
spp_cen_envi2.RData - Data used in Moulatlet et al. 2023 (JAnimalEcology). File necessary to run the ENM.

sp_dists_covmat_v3_resu_euclidian.Rdata - Results of the ENM, used to make analysis, Figures and tables for the manuscript (Most useful for the reviewers to use with the code script_mrm_mixed_2024_05_20)

