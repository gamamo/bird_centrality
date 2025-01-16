############################################################
#                                                          #
# Script to run regression models and make the figures     #
# of the manuscript Moulatlet et al. "Bird species’ network#
# centrality varies differentially across species within   #
# their climatic niches, submitted to the                  #
# American Naturalist"                                     #
#                                                          #
############################################################

# Load packages -----------------------------------------------------------

library(tidyverse)       # CRAN v2.0.0
library(lme4)            # CRAN v1.1-35.3
library(merTools)        # CRAN v0.6.2
library(DHARMa)          # CRAN v0.4.6
library(marginaleffects) # CRAN v0.20.1
library(patchwork)       # CRAN v1.2.0 
library(sjPlot)          # CRAN v2.8.16
library(janitor)         # CRAN v2.2.0    

# Load data ---------------------------------------------------------------

data <- read_csv(here("data", "distances_data.csv"))


# create a list of species names to be added to suplementary material
speciesnames <- data
speciesnames <- speciesnames[!duplicated(speciesnames),]
speciesnames |> 
  dplyr::select(scientific,family) |> 
  distinct() |>
  rename(species=scientific) |> 
  write_csv(here("products","SpeciesNamesToSm.csv"))

# Analysis ----------------------------------------------------------------
## Make correlation between predicted variables ---------------------------

cor.test(data$nplants,data$nbirds)
cor.test(data$nplants,data$envidist)
cor.test(data$nbirds,data$envidist)


## Prepare data for the modelling -----------------------------------------

# First convert centrality values to all positive. This is done for convenience 
# when interpreting the PC1 values and ranking species from 0 (more peripheral) 
# to the highest positive value (more central), PC1 values were rescaled to non-zero, 
# positive-only values by adding up the absolute value of the minimum score plus a 
# millesimal unit to each centrality value 

datam <- data
datam$PC1 <- data$PC1 + -1*min(datam$PC1) + 0.001

# rename envidist to nichedist, thinking about the products
datam <- datam |> 
  dplyr::rename(nichedist = envidist)

# Center variables. This step is necessary prior to modelling when variables have distinct magnitude 
datam[,c("PC1","nplants","nbirds","nichedist")] <- scale(datam[,c("PC1","nplants","nbirds","nichedist")])

## Modelling --------------------------------------------------------------

# Run the global model
mI <- lmer(PC1~nichedist*nbirds*nplants+
             (nichedist*nbirds*nplants|scientific),
           datam)

# this first model had singularity issues, so we make a new model, where the cross-random terms
# were additive and not multiplicative and we used the || notation to constraint the correlation
# of the fixed terms

# Run fixed global model
mA <- lmer(PC1~nichedist*nbirds*nplants+
             (nichedist+nbirds+nplants||scientific),
           datam)

# get some statistics
summary(mA)
r.squaredGLMM(mA)
anova(mA)

#check model overdispersion. Test values larger than one indicate overdispersion
simulationOutput <- simulateResiduals(fittedModel = mA)
testDispersion(simulationOutput,type = c("PearsonChisq"),alternative = "greater")

## Evaluate the fixed effects ----------------------------------------------

# simulate fixed effects and plot. This is step is to get the effect size
feExA <- FEsim(mA, 1000)

feExA$term <- gsub("nbirds","Bird diversity",feExA$term)
feExA$term <- gsub("nplants","Plant diversity",feExA$term)
feExA$term <- gsub(":"," - ",feExA$term)
feExA$term <- gsub("nichedist","Climatic niche centroid",feExA$term)

# This code produces Figure 4 of the article
fep <- plotFEsim(feExA) +
  theme_bw()+
  labs(y="Median effect estimate",x="",title="")+
  theme(text = element_text(size=14,color = "black"),
        panel.grid = element_blank())
fep

# and save
ggsave(here("products","Figure4.jpeg"),units="cm", height = 12,width = 20)

# Sensitivity analysis - skip if necessary ---------------------------------------
if(F){ 

# it was required in the review round 2 that we
# re-run the model for subsets of species separately
# the results of these analyses were only used in the response letter, and not in the 
# manuscript
  
  envislope2 <- envislope |> dplyr::filter(direction==1)
  datam2 <- datam[datam$scientific %in% envislope2$scientific,]
  
  mI2 <- lmer(PC1~nichedist*nbirds*nplants+
                (nichedist*nbirds*nplants||scientific),
              datam2)
  summary(mI2)
  feExA <- FEsim(mI2, 1000)
  feExA$term <- gsub("nbirds","Bird diversity",feExA$term)
  feExA$term <- gsub("nplants","Plant diversity",feExA$term)
  feExA$term <- gsub(":"," - ",feExA$term)
  feExA$term <- gsub("nichedist","Climatic niche centroid",feExA$term)
  
  fep2 <- plotFEsim(feExA) +
    theme_bw()+
    labs(y="Median effect estimate",x="",title="Positive")+
    theme(text = element_text(size=14,color = "black"),
          panel.grid = element_blank())
  fep2
  
  envislope2 <- envislope |> dplyr::filter(direction==-1)
  datam2 <- datam[datam$scientific %in% envislope2$scientific,]
  mI2 <- lmer(PC1~nichedist*nbirds*nplants+
                (nichedist*nbirds*nplants||scientific),
              datam2)
  summary(mI2)
  feExA <- FEsim(mI2, 1000)
  feExA$term <- gsub("nbirds","Bird diversity",feExA$term)
  feExA$term <- gsub("nplants","Plant diversity",feExA$term)
  feExA$term <- gsub(":"," - ",feExA$term)
  feExA$term <- gsub("nichedist","Climatic niche centroid",feExA$term)
  
  fep3 <- plotFEsim(feExA) +
    theme_bw()+
    labs(y="Median effect estimate",x="",title="Negative")+
    theme(text = element_text(size=14,color = "black"),
          panel.grid = element_blank())
  fep3
  
  fep2 + fep3
  
  # End of sensitivity analysis
}

## Evaluate the random effects --------------------------------------------

# Check the original data structure. This structure
# affects how the slopes direction will be calculated, so it is necessary to look at it and,
# if necessary, we need to sort the dataframe

insight::get_data(mA)  

# Extract the estimated slope for each variable per grouping factor (i.e. species)
# Also, the slopes are classified in positive, negative or neutral based on the first and last estimated
# values. The code do that for all three variables of importance

# Variable nichedist
envislope <- mA |> slopes(type="response",variables=c("nichedist"), 
                    re.form = NA) |> 
  group_by(scientific) |>
  arrange(nichedist) |> 
  mutate(f = round(first(estimate),3),l = round(last(estimate),3)) |> 
  mutate(direction = case_when(
    l == f ~ "0",
    l < f ~ "-1",
    l > f ~ "1"
  ))

# convert the column scientific into factors
envislope$scientific <- as.factor(envislope$scientific)

# Variable bird diversity
birdslope <- mA |> slopes(type="response",variables=c("nbirds"), 
                          re.form = NA,hypothesis = "meandev") |> 
  group_by(scientific) |>
  arrange(nbirds) |> 
  mutate(f = round(first(estimate),3),l = round(last(estimate),3)) |> 
  mutate(direction = case_when(
    l == f ~ "0",
    l < f ~ "-1",
    l > f ~ "1"
  ))

# Variable plant diversity
plantslope <- mA |> slopes(type="response",variables=c("nplants"), 
                          re.form = NA) |> 
  group_by(scientific) |>
  arrange(nplants) |> 
  mutate(f = round(first(estimate),3),l = round(last(estimate),3)) |> 
  mutate(direction = case_when(
    l == f ~ "0",
    l < f ~ "-1",
    l > f ~ "1"
  ))

# join all the slopes into a single dataframe
c_slopes <- rbind(envislope,birdslope,plantslope)

# classify the slopes according to the hypothetical scenarios presented in 
# the Figure 1 of the manuscript

s3 <- c_slopes |> 
  dplyr::select(term,scientific, direction) |> 
  dplyr::distinct() |> 
  pivot_wider(names_from = term,values_from = direction)

# classification of climatic niches
s3$nicheEffect <- NA

s3[s3$nichedist < 0 , "nicheEffect"] = "slope(nichedist) <0"
s3[s3$nichedist == 0, "nicheEffect"] = "slope(nichedist) =0"
s3[s3$nichedist > 0 , "nicheEffect"] = "slope(nichedist) >0"

s3$compReso <- NA
s3$compResoDir <- NA

s3[s3$nplants == 0 & s3$nbirds == 0 , "compReso"] = "no effect"

s3[s3$nplants > 0  & s3$nbirds < 0  , "compReso"] = "opposite direction"
s3[s3$nplants > 0  & s3$nbirds < 0  , "compResoDir"] = "nplants (+)/ nbirds(-)"

s3[s3$nplants < 0  & s3$nbirds > 0   , "compReso"] = "opposite direction"
s3[s3$nplants < 0  & s3$nbirds > 0  , "compResoDir"] = "nplants (-)/ nbirds(+)"

s3[s3$nplants > 0  & s3$nbirds > 0  , "compReso"] = "same direction"
s3[s3$nplants > 0  & s3$nbirds > 0  , "compResoDir"] = "Positive"
s3[s3$nplants < 0  & s3$nbirds < 0  , "compReso"] = "same direction"
s3[s3$nplants < 0  & s3$nbirds < 0  , "compResoDir"] = "Negative"

s3[s3$nplants > 0 & s3$nbirds == 0   , "compReso"] = "single variable effect"
s3[s3$nplants > 0 & s3$nbirds == 0   , "compResoDir"] = "nplants (+)"

s3[s3$nplants == 0 & s3$nbirds < 0   , "compReso"] = "single variable effect"
s3[s3$nplants == 0 & s3$nbirds < 0   , "compResoDir"] = "nbirds (-)"

s3[s3$nplants < 0 & s3$nbirds == 0   , "compReso"] = "single variable effect"
s3[s3$nplants < 0 & s3$nbirds == 0   , "compResoDir"] = "nplants (-)"

s3[s3$nplants == 0 & s3$nbirds > 0   , "compReso"] = "single variable effect"
s3[s3$nplants == 0 & s3$nbirds > 0   , "compResoDir"] = "nbirds (+)"

s3 <- s3 |> unite("scenarios", nicheEffect:compReso,sep = " : ",remove = F)
# end of the Analyses

# Make plots and tables for the article -------------------------------------------

# plot fixed effects predictions
p1 <- plot_model(mA, type="pred", terms=c("nichedist"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Distance to the climatic centroid",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p1

p2 <- plot_model(mA, type="pred", terms=c("nbirds"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Bird diversity",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p2

p3 <- plot_model(mA, type="pred", terms=c("nplants"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Plant diversity",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p3

# combine all the plots and save. This is the figure 3
(p1 + p2 + p3)+
  plot_annotation(tag_levels = "a",tag_suffix = ")")

ggsave(here("products","Figure3.jpeg"),units = "cm",width = 40,height = 15)

# Result tables

#Table 1
s3 |> 
  tabyl(compReso,nicheEffect) |> 
  adorn_totals(c("row")) %>%
  adorn_percentages("col") %>%
  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  adorn_title(row_name="Interaction Plant-bird",
              col_name = "Niche distance effect",
              placement = "combined") %>%
  adorn_ns(position = "front") |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table1.docx"))

# Table 2
s3 |> 
  dplyr::filter(compReso == "same direction" ) |> 
    tabyl( compResoDir,nicheEffect) |> 
  adorn_totals(c("row")) %>%
  adorn_percentages("col") %>%
  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  adorn_title(row_name="Interaction Plant-bird",
              col_name = "Niche distance effect",
              placement = "combined") %>%
  adorn_ns(position = "front") |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table2.docx")) 

#Table 3
s3 |> 
  dplyr::filter(compReso == "opposite direction"  ) |> 
  tabyl( compResoDir,nicheEffect) |> 
  adorn_totals(c("row")) %>%
  adorn_percentages("col") %>%
  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  adorn_title(row_name="Interaction Plant-bird",
              col_name = "Niche distance effect",
              placement = "combined") %>%
  adorn_ns(position = "front") |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table3.docx")) 
