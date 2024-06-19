############################################################
#                                                          #
#     Script to run regression models and make figures     #
#                                                          #
############################################################

# load packages -----------------------------------------------------------

library(janitor)    
library(tidyverse)      
library(patchwork)  
library(MuMIn)      
library(lme4)
library(DHARMa)
library(merTools)
library(performance)
library(sjPlot)
library(ggeffects)
library(interactions)
library(marginaleffects)
library(psych)

# Load model outputs -----------------------------------------------------
load("sp_dists_covmat_v3_resu_euclidian.Rdata")
spcovmat <- spdist_covmat_df


# Prepare the datasets ----------------------------------------------------
# create a data object with the relevant columns
data <-  
  spcovmat |> 
  dplyr::select(Scientific,Family1, PC1, nplants, nbirds, nichedist,database) |> 
  dplyr::rename(scientific = Scientific, family = Family1)


# import traits to match with the species
traits <- read_csv("birdtraits_spmetrics_phy_v01_05_2023.csv")

# select the important traits for this study
traits <-
  traits |> 
  dplyr::select(Scientific, RangeSize_Meters, Mass,
         Trophic.Level,Wing.Length,Realm, Family1) |> 
  dplyr::rename(scientific = Scientific, range = RangeSize_Meters, 
                bodymass = Mass, trophic = Trophic.Level,
                wl= Wing.Length,realm=Realm,family=Family1)


# Get some dataset information --------------------------------------------
# the range of network numbers evaluated for each species
range(table(data$scientific))

#how many families?
length(table(data$family))

# the species names to add to suplementary material
speciesnames <- data
speciesnames <- speciesnames[!duplicated(speciesnames),]
speciesnames |> dplyr::select(scientific,family) |> distinct() |>
  rename(species=scientific) |> write_csv("SpeciesNamesToSm_resu.csv")

# Make correlation between predicted variables ----------------------------

cor.test(data$nplants,data$nbirds)
cor.test(data$nplants,data$nichedist)
cor.test(data$nbirds,data$nichedist)


# Analysis ----------------------------------------------------------------

#convert centrality values to all positive
datam <- data
datam$PC1 <- data$PC1 + -1*min(datam$PC1) + 0.001

datam <- datam |> rename(nichedist = nichedist)
#center variables
datam[,c("PC1","nplants","nbirds","nichedist")] <- scale(datam[,c("PC1","nplants","nbirds","nichedist")])


#modelling
mI <- lmer(PC1~nichedist*nbirds*nplants+
             (nichedist*nbirds*nplants|scientific),
           datam)
# this first model had singularity issues, so we make a new model, where the cross-random terms
# were additive and not multiplicative and we used the || notation to constraint the correlation
# of the fixed terms

mA <- lmer(PC1~nichedist*nbirds*nplants+
             (nichedist+nbirds+nplants||scientific),
           datam)
summary(mA)
r.squaredGLMM(mA)

#check variable effects
#fixed effects
#simulate fixed effects and plot
feExA <- FEsim(mA, 1000)
feExA$term <- gsub("nbirds","Bird diversity",feExA$term)
feExA$term <- gsub("nplants","Plant diversity",feExA$term)
feExA$term <- gsub(":"," - ",feExA$term)
feExA$term <- gsub("nichedist","Climatic niche centroid",feExA$term)

fep <- plotFEsim(feExA) +
  theme_bw()+
  labs(y="Median effect estimate",x="",title="")+
  theme(text = element_text(size=14,color = "black"),
        panel.grid = element_blank())
fep
ggsave("RESU_ESTfixed-effects.jpeg",units="cm", height = 12,width = 20)


#random effects
#simulate random effects and plot
s <- REsim(mA,oddsRatio = F,seed=123,n.sims = 1000)
plotREsim(s)

#check model overdispersion
simulationOutput <- simulateResiduals(fittedModel = mA)
testDispersion(simulationOutput,type = c("PearsonChisq"),alternative = "greater")

#Obtain the slopes for the each species

insight::get_data(mA) # this code shows the original data structure. This structure
#affects how the slopes direction will be calculated, so it is necessary to look at it and,
#if necessary, to sort the dataframe

# the code below is to extract the estimated slope for each variable per grouping factor (i.e. species)
# Also, the slopes are classifyed in positive, negative or neutral based on the first and last estimated
# values. Do that for all three variables of importance

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

# join all the new datasets into a single one
c_slopes <- rbind(envislope,birdslope,plantslope)

if(F){
#EXTRA: play with some plots for understanding the interactions
plot_slopes(mA,variables = "nichedist",by =  c("nbirds","scientific"),re.form=NA)

#select slopes
me <- predict_response(mA, terms = c("nichedist","nbirds","nplants",
                                     "scientific"), 
                       type = "random")
head(me)
plot(me, show_ci = F)
me2 <- me |> as_tibble() |> 
  dplyr::rename(nichedist=x,nbirds=group, nplants=facet,species=panel)
}

#classification ############################################

s3 <- c_slopes |> dplyr::select(term,scientific, direction) |> 
  distinct() |> 
  pivot_wider(names_from = term,values_from = direction)


#classification
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

#Make plots and tables for the article #####

#plot fixed effects predictions
p1 <- plot_model(mA, type="pred", terms=c("nichedist"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Distance to the climatic centroid",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p1

p2 <- plot_model(mA, type="pred", terms=c("nbirds"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Birds diversity",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p2

p3 <- plot_model(mA, type="pred", terms=c("nplants"),show.data = T,jitter = 0.1)+
  theme_bw()+
  labs(x="Plants diversity",title="",y="Centrality")+
  theme(text = element_text(size=16),
        panel.grid  = element_blank())
p3


(p1 + p2 + p3)+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave("RESU_PREDfixed_effects.jpeg",units = "cm",width = 40,height = 15)

# Result tables
tabyl(s3,compReso,nicheEffect) |> 
  adorn_totals(c("row")) %>%
  adorn_percentages("col") %>%
  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  adorn_title(row_name="Interaction Plant-bird",
              col_name = "Niche distance effect",
              placement = "combined") %>%
  adorn_ns(position = "front") |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = "tabyl.docx") 

s3 |> filter(compReso == "same direction"  ) |> 
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
  flextable::save_as_docx(path = "tabyl_samedir.docx") 

s3 |> filter(compReso == "opposite direction"  ) |> 
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
  flextable::save_as_docx(path = "tabyl_oppositedir.docx") 


s3 |> filter(compReso == "single variable effect"  ) |> 
  tabyl( compResoDir,nicheEffect) |> 
  adorn_totals(c("row")) %>%
  adorn_percentages("col") %>%
  adorn_pct_formatting(rounding = "half up", digits = 0) |> 
  adorn_ns(position = "front") 



#
