############################################################
#                                                          #
#     Script to run regression models and make figures     #
#                                                          #
############################################################

# load packages -----------------------------------------------------------

library(lm.beta)    # CRAN v1.7-2
library(vegan)      # CRAN v2.6-4
library(ggrepel)    # CRAN v0.9.4
library(annotater)  # CRAN v0.2.2
library(skimr)      # CRAN v2.1.5
library(readr)      # CRAN v2.1.4
library(janitor)    # CRAN v2.2.0
library(forcats)    # CRAN v1.0.0
library(tidyr)      # CRAN v1.3.0
library(patchwork)  # CRAN v1.1.3
library(ggplot2)    # CRAN v3.4.4
library(plyr)       # CRAN v1.8.9
library(dplyr)      # CRAN v1.1.3
library(tidyverse)  # CRAN v2.0.0
library(conflicted) # CRAN v1.2.0
library(paletteer)  # CRAN v1.5.0
library(MuMIn)      # CRAN v1.47.5

# Load model outputs -----------------------------------------------------
load("sp_dists_covmat_v3.Rdata")

spcovmat <- spdist_covmat_df


# Prepare the datasets ----------------------------------------------------
#delete outliers in envidist

summary(spcovmat$envidist)
spcovmat = subset(spcovmat,spcovmat$envidist>0 & spcovmat$envidist <1.690959e+10)

# create a data object and make new log-transformed columns
data <-  
  spcovmat |> 
  select(Scientific,Family1, PC1, nplants, nbirds, envidist) |> 
  dplyr::mutate(nplantslog = log(nplants),nbirdslog = log(nbirds),
         envidistlog = log(envidist),
         PC1 = 0.1 + PC1 + abs(min(PC1))) |> 
  dplyr::rename(scientific = Scientific, family = Family1)

# import traits to match with the species
traits <- read_csv("birdtraits_spmetrics_phy_v01_05_2023.csv")

# select the important traits for this study
traits <-
  traits |> 
  select(Scientific, RangeSize_Meters, Mass,
         Trophic.Level,Wing.Length,Realm, Family1) |> 
  dplyr::rename(scientific = Scientific, range = RangeSize_Meters, 
                bodymass = Mass, trophic = Trophic.Level,
                wl= Wing.Length,realm=Realm,family=Family1)


# Get some dataset information --------------------------------------------
# the range of network numbers evaluated for each species
range(table(data$scientific))

# the species names to add to suplementary material
speciesnames <- data
speciesnames <- speciesnames[!duplicated(speciesnames),]
write_csv(speciesnames,"SpeciesNamesToSm.csv")


# What are the most common families and spp?

comFam <- 
data |> 
  dplyr::mutate(occ= 1) |> 
  select(family, occ) |> 
  pivot_wider(names_from = family,
              values_from = occ,
              values_fn = sum)

comFam = as.data.frame(t(comFam))
comFam$family = rownames(comFam)

comFam |> 
  dplyr::arrange(desc(V1))

comsp <- 
  data |> 
  dplyr::mutate(occ= 1) |> 
  select(scientific, occ) |> 
  pivot_wider(names_from = scientific,
              values_from = occ,
              values_fn = sum)

comsp = as.data.frame(t(comsp))
comsp$scientific = rownames(comsp)

comsp |> 
  dplyr::arrange(desc(V1))



# Make correlation between predicted variables ----------------------------

cor.test(data$nplants,data$nbirds)
cor.test(data$nplants,data$envidist)
cor.test(data$nbirds,data$envidist)



# Analysis ----------------------------------------------------------------

#create a list to save the outputs
resu_reg_list = list()
resu_reg_list_step = list()


# STEP 1: Run Regression models for each spp. -----------------------------
# Regression are made with and without variable selection (function step)

# run the looping
for(i in unique(data$scientific)){
  
  tempsp = subset(data[data$scientific==i,])
  
  if(length(tempsp[,1])>4){
    
    print(i)
    
    # PC1 vs rich bird + envi
    a1=lm(PC1~nplants+nbirds+envidist      ,data=tempsp)
    a2=lm(PC1~nplants+nbirdslog+envidist   ,data=tempsp)
    a3=lm(PC1~nplantslog+nbirdslog+envidist,data=tempsp)
    a4=lm(PC1~nplantslog+nbirds+envidist   ,data=tempsp)
    
    templist = list(a1,a2,a3,a4)

    aic=AICc(a1,a2,a3,a4)
    aic$model = c("1","2","3","4")
    mc = aic[order(aic$AIC),]
    mc = templist[as.numeric(mc$model[1])]
    mn = unlist(mc)$call$formula
   
    # regression with all terms
    reg_rice =  summary(lm(mn,data=tempsp))
    
    rice_bind  = cbind(
      reg_rice$adj.r.squared,          # adj-rsquared
      reg_rice$r.squared,              # r-squared
      reg_rice$coefficients[c(2:4),1], # beta
      reg_rice$coefficients[c(2:4),4]  # p-value
     )
    
    std_coef =  lm.beta(lm(mn,data=tempsp))$standardized.coefficients[-1]
    
    rice_bind = as.data.frame(cbind(rice_bind,std_coef))
    
    #prepare the result table regression all variables
    resu_reg = round(as.data.frame(rbind(rice_bind)),3)
    colnames(resu_reg) = c("AdjR2","R2","b","p-value","std_b")
    resu_reg$scientific = c(rep(i,nrow((resu_reg))))
    resu_reg$family = rep(unique(data[data$scientific==i,"family"]),nrow(resu_reg))
    resu_reg$btype = c("Plant richness","Bird richness","Envidist")
    resu_reg$mean_cen = c(rep(mean(tempsp$PC1,na.rm=T),nrow(resu_reg)))
    resu_reg$max_cen = c(rep(max(tempsp$PC1,na.rm=T),nrow(resu_reg)))
    resu_reg$median_cen = c(rep(median(tempsp$PC1,na.rm=T),nrow(resu_reg)))
    resu_reg$model = as.character(mn)[3]
    
    #regression after using the "step" function
    st <- step(lm(mn,data=tempsp))
    
   if(ncol(st$model)>1){ # eliminate those species where only the intercept was signif.
    
    reg_rice2 = summary(st)
    
    rice_bind2  = cbind(
      reg_rice2$adj.r.squared,          # adj-rsquared
      reg_rice2$r.squared,              # r-squared
      reg_rice2$coefficients[c(2:length(reg_rice2$coefficients[,1])),1], # beta
      reg_rice2$coefficients[c(2:length(reg_rice2$coefficients[,1])),4]  # p-value
    )
    
    std_coef2 = lm.beta(lm(st,data=tempsp))$standardized.coefficients[-1]
    rice_bind2 = as.data.frame(cbind(rice_bind2,std_coef2))
    
    #prepare the result table regression after "step"
    resu_reg2 = round(as.data.frame(rbind(rice_bind2)),3)
    colnames(resu_reg2) = c("AdjR2","R2","b","p-value","std_b")
    resu_reg2$scientific = c(rep(i,nrow((resu_reg2))))
    resu_reg2$family = rep(unique(data[data$scientific==i,"family"]),nrow(resu_reg2))
    resu_reg2$mean_cen = c(rep(mean(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
    resu_reg2$max_cen = c(rep(max(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
    resu_reg2$median_cen = c(rep(median(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
    resu_reg2$model = as.character(st)[11]
    resu_reg2$var = NA
    
    for(x in 1:nrow(resu_reg2)){
      rn = rownames(rice_bind2)
      resu_reg2$var[x] <- rn[x]
    }
   
    } else{
     
     reg_rice2 = summary(st)
     
     rice_bind2  = cbind(
       reg_rice2$adj.r.squared,          # adj-rsquared
       reg_rice2$r.squared,              # r-squared
       reg_rice2$coefficients[1,2],      # beta
       reg_rice2$coefficients[1,4]       # p-value
     )
     std_coef2 = 0
     rice_bind2 = as.data.frame(cbind(rice_bind2,std_coef2))
     
     
     #prepare the result table regression after "step"
     resu_reg2 = round(as.data.frame(rbind(rice_bind2)),3)
     colnames(resu_reg2) = c("AdjR2","R2","b","p-value","std_b")
     resu_reg2$scientific = c(rep(i,nrow((resu_reg2))))
     resu_reg2$family = rep(unique(data[data$scientific==i,"family"]),nrow(resu_reg2))
     resu_reg2$mean_cen = c(rep(mean(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
     resu_reg2$max_cen = c(rep(max(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
     resu_reg2$median_cen = c(rep(median(tempsp$PC1,na.rm=T),nrow(resu_reg2)))
     resu_reg2$model = as.character(st)[10]
     resu_reg2$var = 1
     
     
     
    }

    
    # send it to the list
    resu_reg_list[[i]] = resu_reg
    resu_reg_list_step[[i]] = resu_reg2
  }
}

#  unlist the results of the looping

resu_reg_df = ldply(resu_reg_list)
resu_reg_step_df = ldply(resu_reg_list_step)

# edit some names in the output dataframe
resu_reg_step_df[which(resu_reg_step_df[,"var"]=="nplants"),"var"] = "Resources"
resu_reg_step_df[which(resu_reg_step_df[,"var"]=="nplantslog"),"var"] = "Resources"
resu_reg_step_df[which(resu_reg_step_df[,"var"]=="nbirds"),"var"] = "Competition"
resu_reg_step_df[which(resu_reg_step_df[,"var"]=="nbirdslog"),"var"] = "Competition"
resu_reg_step_df[which(resu_reg_step_df[,"var"]=="envidist"),"var"] = "Niche position"

resu_reg_step_df <- resu_reg_step_df |> 
  dplyr::rename(btype=var)


# filter only those that R2 are greater than 0.2  

resu_reg_plot <- 
  resu_reg_df |> 
  dplyr::filter(AdjR2>0 & R2>0.2 ) |> 
  dplyr::mutate(direction = ifelse(std_b<0,"Negative","Positive"))

resu_reg_plot_step <- 
  resu_reg_step_df |> 
  dplyr::filter(AdjR2>0 & R2>0.2) |> 
  dplyr::mutate(direction = ifelse(std_b<0,"Negative","Positive"))

# mark only those that p values lower than 0.1
resu_reg_plot$effect <- NA
resu_reg_plot_step$effect <- NA

resu_reg_plot[which(resu_reg_plot$`p-value`>0.051),"effect"] <- "No effect"
resu_reg_plot[which(is.na(resu_reg_plot$effect)),"effect"] <- 
  resu_reg_plot[which(is.na(resu_reg_plot$effect)),"direction"]

resu_reg_plot_step[which(resu_reg_plot_step$`p-value`>0.051),"effect"] <- "No effect"
resu_reg_plot_step[which(is.na(resu_reg_plot_step$effect)),"effect"] <- 
  resu_reg_plot_step[which(is.na(resu_reg_plot_step$effect)),"direction"]

# how many species significative
skim(resu_reg_plot)
skim(resu_reg_plot_step)

# change variable names
resu_reg_plot[which(resu_reg_plot[,"btype"]=="Envidist"),"btype"] = "Niche position"
resu_reg_plot[which(resu_reg_plot[,"btype"]=="Bird richness"),"btype"] = "Competition"
resu_reg_plot[which(resu_reg_plot[,"btype"]=="Plant richness"),"btype"] = "Resources"



# STEP 2: graphics --------------------------------------------------------

# Plot results for individual species (figure not includede in the final manuscript)

reg_gb = ggplot(resu_reg_plot,
                aes(x=btype,y=fct_rev(scientific),
                                  color=effect,size=std_b))+
  geom_point()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  scale_size_area()+
  ylab("")+
  xlab("")
reg_gb

ggsave("resu_reg_multi_models_v180823.jpeg",units = "cm",width=25,height = 35,dpi = 300)

reg_gb_s = ggplot(resu_reg_plot_step,
                aes(x=btype,y=fct_rev(scientific),
                    color=effect))+
  geom_point()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  scale_size_area()+
  ylab("")+
  xlab("")
reg_gb_s

ggsave("resu_reg_multi_models_STEP_v180823.jpeg",units = "cm",width=25,height = 35,dpi = 300)


# get some useful information 
# which species had negative coefficients for all variables?
a <- 
resu_reg_plot |> 
  dplyr::filter(direction=="Negative")
which(sort(table(a$scientific))==3)

a <- 
  resu_reg_plot_step |> 
  dplyr::filter(direction=="Negative")
which(sort(table(a$scientific))==3)


# which species had positive coefficients for all variables?
b <- 
  resu_reg_plot |> 
  dplyr::filter(direction=="Positive")
which(sort(table(b$scientific))==3)

b <- 
  resu_reg_plot_step |> 
  dplyr::filter(direction=="Positive")
which(sort(table(b$scientific))==3)

# which family was more representative
c <- 
  resu_reg_plot_step |> 
  distinct(scientific, .keep_all = TRUE)
sort((table(c$family)*100)/90)

resu_reg_plot_step |> 
  distinct(scientific, .keep_all = TRUE)
sort(table(c$family))


length(unique(resu_reg_plot_step$scientific))


# Scenario classification -------------------------------------------------

df1 <-  resu_reg_plot_step |> 
  select(scientific, mean_cen, max_cen, median_cen, btype, std_b) |> 
  pivot_wider(names_from = btype,
              values_from = std_b)

# add zeros to the non-significant variables
df1[which(is.na(df1$Resources)), c("Resources") ] = 0
df1[which(is.na(df1$Competition)), c("Competition") ] = 0
df1[which(is.na(df1$`Niche position`)), c("Niche position") ] = 0

# classify according to the scenarions of slide 13 in "preliminar figures. pptx"
# not used in the end
df1$scenario = NA

df1[df1$`Niche position` ==0 & df1$Resources <0  & df1$Competition >0  , "scenario"] = "a"
df1[df1$`Niche position` ==0 & df1$Resources >0  & df1$Competition <0  , "scenario"] = "a"

df1[df1$`Niche position` >0  & df1$Resources >=0  & df1$Competition <= 0 , "scenario"] = "b"
df1[df1$`Niche position` <0  & df1$Resources >=0  & df1$Competition <= 0 , "scenario"] = "b"
df1[df1$`Niche position` >0  & df1$Resources <=0  & df1$Competition >= 0 , "scenario"] = "b"
df1[df1$`Niche position` <0  & df1$Resources <=0  & df1$Competition >= 0 , "scenario"] = "b"

df1[df1$`Niche position` >=0  & df1$Resources >0  & df1$Competition > 0 , "scenario"] = "c"
df1[df1$`Niche position` <=0  & df1$Resources >0  & df1$Competition > 0 , "scenario"] = "c"
df1[df1$`Niche position` >=0  & df1$Resources <0  & df1$Competition < 0 , "scenario"] = "c"
df1[df1$`Niche position` <=0  & df1$Resources <0  & df1$Competition < 0 , "scenario"] = "c"

df1[df1$`Niche position` <0   & df1$Resources ==0 & df1$Competition == 0, "scenario"] = "d"
df1[df1$`Niche position` >0   & df1$Resources ==0 & df1$Competition == 0, "scenario"] = "d"
df1[df1$`Niche position` ==0  & df1$Resources >0  & df1$Competition == 0, "scenario"] = "d"
df1[df1$`Niche position` ==0  & df1$Resources <0  & df1$Competition == 0, "scenario"] = "d"
df1[df1$`Niche position` ==0  & df1$Resources ==0 & df1$Competition >  0, "scenario"] = "d"
df1[df1$`Niche position` ==0  & df1$Resources ==0 & df1$Competition <  0, "scenario"] = "d"


# classify according to the TABLE 1 of the manuscript
df1$scenarios_hypo = NA

df1[df1$`Niche position` <0   & df1$Resources ==0 & df1$Competition == 0 , "scenarios_hypo"] = "a"

df1[df1$`Niche position` <0   & df1$Resources >0  & df1$Competition < 0  , "scenarios_hypo"] = "b"
df1[df1$`Niche position` <0   & df1$Resources <0  & df1$Competition >0   , "scenarios_hypo"] = "b"
df1[df1$`Niche position` <0   & df1$Resources >0 & df1$Competition ==0   , "scenarios_hypo"] = "b"
df1[df1$`Niche position` <0   & df1$Resources ==0 & df1$Competition <0   , "scenarios_hypo"] = "b"
df1[df1$`Niche position` <0   & df1$Resources <0 & df1$Competition ==0   , "scenarios_hypo"] = "b"
df1[df1$`Niche position` <0   & df1$Resources ==0 & df1$Competition >0   , "scenarios_hypo"] = "b"

df1[df1$`Niche position` <0   & df1$Resources >0  & df1$Competition > 0  , "scenarios_hypo"] = "c"
df1[df1$`Niche position` <0   & df1$Resources <0  & df1$Competition < 0  , "scenarios_hypo"] = "c"

df1[df1$`Niche position` >0   & df1$Resources <=0  & df1$Competition <=0 , "scenarios_hypo"] = "d"
df1[df1$`Niche position` >0   & df1$Resources >0  & df1$Competition <= 0 , "scenarios_hypo"] = "d"
df1[df1$`Niche position` >0   & df1$Resources <=0  & df1$Competition > 0 , "scenarios_hypo"] = "d"


df1[df1$`Niche position` ==0                                             , "scenarios_hypo"] = "No niche effect"
df1[df1$`Niche position` <0   & df1$Resources <0  & df1$Competition < 0  , "scenarios_hypo"] = "Exception"


tabyl(df1, scenarios_hypo)

# merge with some traits

df1 <- left_join(df1, traits, by="scientific")
data <- left_join(data,df1, by="scientific")
data <-  left_join(data, resu_reg_plot_step[,c("scientific","direction","effect")],
                   by="scientific")


#how many significative species
length(unique(df1$scientific))

#how many species in each scenario
table(df1$scenarios_hypo)
round((table(df1$scenarios_hypo)*100)/90,1)

# Make figure 2 --------------------------------------------------------

aa <- subset(df1,df1$Resources!=0 &df1$`Niche position`!=0)
summary(lm(aa$Resources~aa$`Niche position`))

aag <- ggplot(subset(df1,df1$Resources!=0 &df1$`Niche position`!=0),
              aes(x=`Niche position`,y=Resources,
              ))+
  geom_point(size=3,color="black")+
  theme_classic()+
  #scale_color_virlabs(title= expression(paste(R^2,"=0.034",","," ","p=0.83")))+idis_c(guide = guide_legend())+
  #labs(color="Std. beta coeff. - Competition")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")+
  geom_smooth(method="lm")+
  labs(title= expression(paste(R^2,"=0.03",","," ","p=0.83")))+
  xlab("Std. beta coefficients - Niche position")+
  ylab("Std. beta coefficients - Resources")
aag

a <- subset(df1,df1$Competition!=0 &df1$`Niche position`!=0)
summary(lm(a$Competition~a$`Niche position`))

ag <- ggplot(subset(df1,df1$Competition!=0 &df1$`Niche position`!=0),
             aes(x=`Niche position`,y=Competition,
             ))+
  geom_point(size=3,color="black")+
  theme_classic()+
  #scale_color_viridis_c(guide = guide_legend())+
  #labs(color="Std. beta coeff. - Competition")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")+
  geom_smooth(method="lm")+
  labs(title= expression(paste(R^2,"=0.13",","," ","p<0.05")))+
  xlab("Std. beta coefficients - Niche position")+
  ylab("Std. beta coefficients - Competition")
ag

aag / ag+plot_annotation(tag_levels = "A",tag_suffix = ")")
ggsave("Figure2_v100124.jpeg",dpi = 300,height = 10,width = 7)


# Make Figure 3 -----------------------------------------------------------
gg <- subset(df1,df1$`Niche position`!=0)
colnames(gg)[7] = "Niche_position"
ggam <- mgcv::gam(mean_cen~s(Niche_position),data=gg)
summary(ggam)

gdc3 =ggplot(subset(df1,df1$`Niche position`!=0),
             aes(x=`Niche position`,y=`mean_cen`))+
  geom_point(size=3,color="black")+
  #geom_text_repel(aes(label = scientific),size=2,fontface="italic",
  #                max.overlaps = 20)+
  theme_classic()+
  xlab("Std. beta coefficients - Niche position")+
  ylab("Species centrality")+
  #labs(color="Std. beta coeff. - Competition")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")+
  #geom_vline(xintercept = 0,linetype="dashed")+
  geom_smooth(method = "gam",se=T,span=1)+
  labs(title= expression(paste(R^2,"=0.22",","," ","p<0.05")))
gdc3


gg2 <- subset(df1,df1$Competition!=0)
ggam2 <- mgcv::gam(mean_cen~s(Competition),data=gg2)
summary(ggam2)

gdc4 =ggplot(subset(df1,df1$Competition!=0),
             aes(x=Competition,y=`mean_cen`))+
  geom_point(size=3,color="black")+
  #geom_text_repel(aes(label = scientific),size=2,fontface="italic",
  #                max.overlaps = 20)+
  theme_classic()+
  xlab("Std. beta coefficients - Competition")+
  ylab("Species centrality")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")+
  #geom_vline(xintercept = 0,linetype="dashed")+
  geom_smooth(method = "gam",se=T,span=1)+
  labs(title= expression(paste(R^2,"=0.09",","," ","p=0.06")))
gdc4


gg3 <- subset(df1,df1$Resources!=0)
ggam3 <- mgcv::gam(mean_cen~s(Resources),data=gg3)
summary(ggam3)

gdc5 =ggplot(subset(df1,df1$Resources!=0),
             aes(x=Resources,y=`mean_cen`,
                 ))+
  geom_point(size=3,color="black")+
  #geom_text_repel(aes(label = scientific),size=2,fontface="italic",
  #                max.overlaps = 20)+
  theme_classic()+
  xlab("Std. beta coefficients - Resources")+
  ylab("Species centrality")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")+
  #geom_vline(xintercept = 0,linetype="dashed")+
  geom_smooth(method = "gam",se=T,span=1)+
  labs(title= expression(paste(R^2,"=0.01",","," ","p=0.57")))
gdc5

gdc3/gdc4/gdc5+
  plot_annotation(tag_levels = "A",tag_suffix = ")")
  
ggsave("Figure3_v100124.jpeg",units = "cm",width=18,height = 30,dpi = 300)  

# Make Figure 4 -----------------------------------------------------------
# Centrality variation between scenarios

summary(resu_reg_step_df$mean_cen)
summary(df1$mean_cen)

df3 = data[,c("scientific","PC1")]
df3=left_join(df3,df1[,c("scientific","scenarios_hypo")])
df3 = df3[!is.na(df3$scenario),]
head(df3)

#ANOVA
anov =aov(data=subset(df3,df3$scenarios_hypo!="Exception"),PC1~scenarios_hypo)
summary (anov)
TukeyHSD(anov)

bb <- resu_reg_step_df
bb <- resu_reg_step_df[,c("scientific","mean_cen")]
bb <- bb[!duplicated(bb),]


h1 <- ggplot(data=bb, aes(x=mean_cen))+
  geom_histogram(color="gray30")+
  geom_histogram(data=df1,aes(mean_cen),fill="#D82632",color="gray30")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  xlab("Mean centrality")+
  ylab("Count of species")+
  scale_x_continuous(limits = c(0,8))
h1

b1 <- ggplot(data=subset(df1,df1$scenarios_hypo!="Exception"),aes(y=mean_cen,x=scenarios_hypo))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  ylab("Mean species centrality")+
  xlab("")+
  scale_color_viridis_d()
b1

b1 + h1+
  plot_annotation(tag_levels = "A",tag_suffix = ")")
ggsave("Figure4_v100124.jpeg",units="cm",dpi=300,width = 25,height = 10)


# Make supplementary figure -----------------------------------------------

df3=left_join(resu_reg_plot_step,df1[,c("scientific","scenarios_hypo")])

cols = c("No effect" = "gray92", "Negative" ="#D82632",
         "Positive" = "#0099CC")

sa <- ggplot(subset(df3,df3$scenario=="a"),
             aes(x=btype,y=fct_rev(scientific),
                 fill=direction))+
  geom_tile()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  scale_size_area()+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = cols)
sa

sb <- ggplot(subset(df3,df3$scenario=="b"),
             aes(x=btype,y=fct_rev(scientific),
                 fill=direction))+
  geom_tile()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  #scale_size_area()+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = cols)
sb

sc <- ggplot(subset(df3,df3$scenario=="c"),
             aes(x=btype,y=fct_rev(scientific),
                 fill=direction))+
  geom_tile()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  #scale_size_area()+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = cols)
sc

sd <- ggplot(subset(df3,df3$scenario=="d"),
             aes(x=btype,y=fct_rev(scientific),
                 fill=direction))+
  geom_tile()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  #scale_size_area()+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = cols)
sd

sa + sb + sc + sd +
  plot_annotation(tag_levels = "A",tag_suffix = ")")
ggsave("FigureS1_v100124.jpeg",dpi = 300, height = 12,width = 12)

snne <- ggplot(subset(df3,df3$scenario=="No niche effect"),
             aes(x=btype,y=fct_rev(scientific),
                 fill=direction))+
  geom_tile()+
  geom_text(aes(label=std_b),color="black",size=4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_blank(),
        legend.position = "none",
        text= element_text(size=14))+
  #scale_size_area()+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = cols)
snne

snne
ggsave("FigureS2_v100124.jpeg",dpi = 300, height = 12,width = 6)



