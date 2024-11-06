## Remove leftovers
rm(list=ls())

## Needed packages
library(lme4)
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
library(lmerTest)
library(ggplot2)
library(regclass)
library(car)
library(MuMIn)
library(geodata)
library(raster)
library(FactoMineR)
library(ggcorrplot)
library(corrr)
library(factoextra)
library(paran)
library(clusterSim)
library(BBmisc)
library(arm)
library(ggeffects)
library(forestplot)
library(ggpubr)
library(sf)
library(terra)
library(dplyr)
library(spData)
library(tmap)
library(leaflet)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(merTools)
library(ggplotify)
library(patchwork)
library(ggrepel)
library(cAIC4)
library(psych)

## Save file in your directory first
## Data
data <- read.csv("./rodent.csv", header = T)
colnames(data)

## Remove rows without lat and lon
data <- data[-which(is.na(data$Lat)),]
any(is.na(data$Lon))

###############################################################################################################
###############################################################################################################
## Habitat structure 
## PCA Habitat
colnames(data)
## No NAs
length(which(is.na(data$TR1)))
length(which(is.na(data$TR2)))
length(which(is.na(data$SHR)))
length(which(is.na(data$MOS)))
length(which(is.na(data$LIC)))
length(which(is.na(data$EPL)))
length(which(is.na(data$CWD)))
length(which(is.na(data$FWD)))
length(which(is.na(data$BOU)))
length(which(is.na(data$SNA)))
length(which(is.na(data$HER)))
length(which(is.na(data$decay)))
length(which(is.na(data$stumps)))
length(which(is.na(data$heather)))
length(which(is.na(data$bilberry)))
length(which(is.na(data$lingonberry)))
length(which(is.na(data$crowberry)))
length(which(is.na(data$grass)))
length(which(is.na(data$HER)))
length(which(is.na(data$ferns)))
length(which(is.na(data$fieldlayer)))
length(which(is.na(data$manaheight)))

## data <- data[-which(is.na(data[,c(85:104,106:108)])),]
data <- data[-which(!complete.cases(data[,c(85:104,106:108)])),]

## We'll lose 59 rows
length(which(is.na(data$UMB)))

## Microbiota
which(is.na(data$gut_b_shan))
which(is.na(data$gut_f_5k_shan))
which(is.na(data$gut_f_29k_rich))
which(is.na(data$gut_f_29k_shan))

## With NAs
length(which(is.na(data$SoilpH)))
length(which(is.na(data$soil_b_rich)))
length(which(is.na(data$soil_b_shan)))
length(which(is.na(data$soil_f_rich)))
length(which(is.na(data$soil_f_shan)))

## Removing
data2 <- data[-which(is.na(data$soil_b_rich)),]
length(which(is.na(data2$soil_b_rich)))
length(which(is.na(data2$soil_b_shan)))
length(which(is.na(data2$soil_f_rich)))
length(which(is.na(data2$soil_f_shan)))

## icrobiota
which(is.na(data2$gut_f_5k_shan))
###############################################################################################################
###############################################################################################################
## Tree cover
map <- raster("./Dissimilarity_01_05_1km_uint32.tif")

## Occurence data
occ <- as.data.frame.matrix(data2[,c("Lon", "Lat")])
global.veg <- extract(map, occ)
data2$ForestCover <- global.veg
data2 <- data2[-which(is.na(data2$ForestCover)),]

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
colnames(data2)
corr_matrix <- cor(data2[,c(85:104,106:108,115)])
ggcorrplot(corr_matrix)
## Some high correlations
colnames(data2[,c(85:104,106:108,115)])
## PCA
pca.hb <- prcomp(data2[,c(85:104,106:108,115)], scale = T)
summary(pca.hb)

## Plot
fviz_pca_var(pca.hb, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
fviz_eig(pca.hb, addlabels = TRUE)

## How many axes should we use?
data_norm <- data.Normalization(data2[,c(85:104,106:108,115)], type="n1", normalization="column")
paran(data_norm, iterations=10000, quietly=FALSE, status=FALSE, all=TRUE, cfa=FALSE, graph=TRUE, color=TRUE, col=c("black","red","blue"),
      lty=c(1,2,3), lwd=1, legend=TRUE, file="",
      width=640, height=640, grdevice="png", seed=0, mat=NA, n=NA)
## We should use 6 axes: eigenvalues higher than 1, according to Horn's Parallel Analysis. However, we will use Frauke's paper directions. PC1 and PC2 explained more than 10%, thereofre, we will use only the first two axes. 

## Contribution to each axis
a <- fviz_contrib(pca.hb, "var", axes=1, xtickslab.rt=90)
plot(a,main = "Variables percentage contribution of first Principal Components")
b <- fviz_contrib(pca.hb, "var", axes=2, xtickslab.rt=90)
plot(b,main = "Variables percentage contribution of second Principal Components")
ggarrange(a,b,nrow= 2)

## Use first four axes
cor.hb.eixos <- pca.hb$rotation
eixos.pca_hb <- pca.hb$x
data2$PCHB1 <-  eixos.pca_hb[,1]
data2$PCHB2 <-  eixos.pca_hb[,2]

###############################################################################################################
###############################################################################################################
## Soil Microbiota
colnames(data2)
data2 <- data2[-which(!complete.cases(data2[,c("soil_ph", "soil_b_rich", "soil_b_shan", "soil_f_rich", "soil_f_shan")])),]
corr_matrix_soil <- cor(data2[,c("soil_ph", "soil_b_rich", "soil_b_shan", "soil_f_rich", "soil_f_shan")])
ggcorrplot(corr_matrix_soil)
## high correlations

## PCA
pca.soil <- prcomp(data2[,c("soil_ph", "soil_b_rich", "soil_b_shan", "soil_f_rich", "soil_f_shan")], scale = T)
summary(pca.soil)

## Plot
fviz_pca_var(pca.soil, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
fviz_eig(pca.soil, addlabels = TRUE)

## Contribution to each axis
a_soil <- fviz_contrib(pca.soil, "var", axes=1, xtickslab.rt=90)
plot(a_soil,main = "Variables percentage contribution of first Principal Components")
b_soil <- fviz_contrib(pca.soil, "var", axes=2, xtickslab.rt=90)
plot(b_soil,main = "Variables percentage contribution of second Principal Components")
ggarrange(a_soil,b_soil,nrow= 2)

## Use first four axes
eixos.pca_soil <- pca.soil$x
data2$PCSoil_1 <-  eixos.pca_soil[,1]
data2$PCSoil_2 <-  eixos.pca_soil[,2]


###############################################################################################################
###############################################################################################################
##########################################################
## Rodents
## Puumala orthohantavirus
colnames(data2)
data3 <- data2[-which(is.na(data2$PUUV)),]
data3$PUUV <- as.factor(data3$PUUV)
levels(data3$PUUV)
levels(data3$PUUV) <- c(0,1)
table(data3$PUUV)
## 50 infected

##########################################################
## Tick-borne
## Neoehrlichia mikurensis 
data3$NEOERL <- as.factor(data3$NEOERL)
levels(data3$NEOERL) <- c(0,1)
table(data3$NEOERL)
## 9 infected
any(is.na(data3$NEOERL))

## Babesia microriti
data3$BM <- as.factor(data3$BM)
levels(data3$BM) <- c(0,1)
table(data3$BM)
## 147 infected

## Borrelia burgdorferi
data3$BBSL <- as.factor(data3$BBSL)
levels(data3$BBSL) <- c(0,1)
table(data3$BBSL)
## 28 infected
any(is.na(data3$BBSL))

## Anaplasma
data3$ANA <- as.factor(data3$ANA)
levels(data3$ANA) <- c(0,1)
table(data3$ANA)
## 155 infected
any(is.na(data3$ANA))

##########################################################
## Fleas, lice, flies
## Bartonella
data3$BART <- as.factor(data3$BART)
levels(data3$BART) <- c(0,1)
table(data3$BART)
## 143 infected,
any(is.na(data3$BART))

###############################################################################################################
###############################################################################################################
## Treatment
#data <- data[-which(data$Treatment == "burnt_forest"),] 
table(data3$Treatment)
table(data3[which(data3$Treatment != "urban" & data3$Treatment != "suburban"),"Treatment"])
data3[which(data3$Treatment != "urban" & data3$Treatment != "suburban"),"Treatment"] <- "forest"
data3$Treatment <- as.factor(data3$Treatment)
table(data3$Treatment)

hist(log(data3$HII))
summary(lm(HII ~Treatment, data3))

## 71% of r squared
ggplot(data3, aes(x = Treatment, y = HII)) + geom_boxplot() + theme_bw()

###############################################################################################################
###############################################################################################################
table(data3$PUUV)
table(data3$BM)
table(data3$BBSL)
table(data3$ANA)
table(data3$BART)
table(data3$NEOERL)

###############################################################################################################
###############################################################################################################
##########################################   1    #############################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
colnames(data3)

## Variable distributions
hist(normalize(data3$PCHB1)) 
hist(normalize(data3$PCHB2)) 
hist(data3$PCSoil_1)
hist(data3$PCSoil_2)
hist(data3$gut_b_shan)

## Age and sex
## Host sex
any(is.na(data3$HostSex))
table(data3$HostSex)
table(data3$SpatialCode_2)
which(is.na(data3$mean_teeth_age))

## Big spleen
any(is.na(data3$big_spleen))
data3$big_spleen <- as.factor(data3$big_spleen)
levels(data3$big_spleen) <- c(0,1)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
## Maps
## Shapefiles: https://tapiquen-sig.jimdo.com/english-version/free-downloads/europe/
nordic <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")
world_points<- st_centroid(nordic)
world_points <- cbind(nordic, st_coordinates(st_centroid(nordic$geometry)))

jpeg("Map.jpeg", width = 8, height = 8, units = 'in', res = 800) 
ggplot() +
  geom_sf(data = nordic, fill= "antiquewhite") +
  coord_sf(xlim = c(10, 30), ylim = c(60, 71)) + 
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "black", fontface = "bold", check_overlap = TRUE) +
  xlab("Latitude") + ylab("Longitude") + 
  annotation_scale(location = "tl", width_hint = 0.5) + 
  annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(aes(x = data3$Lon, y = data3$Lat), color="red", size = 2)
dev.off()

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
###############################################################################################################
###############################################################################################################
## Age
hist(data3$mean_teeth_age)
dim(data3[which(data3$mean_teeth_age < 60),])

## Separated per disease
colnames(data3)

## PUUV
model <- lmer(as.integer(as.character(PUUV)) ~ BM + BBSL + ANA + NEOERL + BART + big_spleen + 
                normalize(PCHB1) + normalize(PCHB2) + 
                normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
## removed climate PC2
vif(model) 
summary(model)

## BM
model_BM <- lmer(as.integer(as.character(BM)) ~  PUUV + BBSL + ANA + NEOERL + BART + big_spleen + 
                   normalize(PCHB1) + normalize(PCHB2) + 
                   normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                   normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
vif(model_BM)
summary(model_BM)

## BART
model_BART <- lmer(as.integer(as.character(BART)) ~ PUUV + BBSL + ANA + NEOERL + BM + big_spleen + 
                     normalize(PCHB1) + normalize(PCHB2) + 
                     normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                     normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
vif(model_BART)
summary(model_BART)

## ANA
model_ANA <- lmer(as.integer(as.character(ANA)) ~ PUUV + BBSL + BART + NEOERL + BM + big_spleen + 
                    normalize(PCHB1) + normalize(PCHB2) + 
                    normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                    normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
vif(model_ANA)
summary(model_ANA)

## BBSL
model_BBSL <- lmer(as.integer(as.character(BBSL)) ~ PUUV + ANA + BART + NEOERL + BM + big_spleen + 
                     normalize(PCHB1) + normalize(PCHB2) + 
                     normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                     normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
vif(model_BBSL)
summary(model_BBSL)

## NEOERL
model_NEOERL <- lmer(as.integer(as.character(NEOERL)) ~ PUUV + ANA + BART + BBSL + BM + big_spleen + 
                       normalize(PCHB1) + normalize(PCHB2) + 
                       normalize(PCSoil_1) + normalize(PCSoil_2) + normalize(gut_b_shan) +
                       normalize(HII) + (1|HostSex) + (1|mean_teeth_age) + (1|spatial_code_2), data3, na.action = na.fail)
vif(model_NEOERL)

## Model selection
selection <- dredge(model)
selection2 <- model.avg(selection)
summary(selection2)
confint(selection2)

selection3 <- dredge(model_BM)
selection4 <- model.avg(selection3)
summary(selection4)
confint(selection4)

selection5 <- dredge(model_BART)
selection6 <- model.avg(selection5)
summary(selection6)
confint(selection6)

selection7 <- dredge(model_ANA)
selection8 <- model.avg(selection7)
summary(selection8)
confint(selection8)

selection9 <- dredge(model_BBSL)
selection10 <- model.avg(selection9)
summary(selection10)
confint(selection10)

selection11 <- dredge(model_NEOERL)
selection12 <- model.avg(selection11)
summary(selection12)
confint(selection12)

#save.image("/Users/victoraguiardesouzapenha/Desktop/AFinal.RData")
#write.csv(as.data.frame(selection),"/Users/victoraguiardesouzapenha/Desktop/PUUV.csv")
#write.csv(as.data.frame(selection3),"/Users/victoraguiardesouzapenha/Desktop/BM.csv")
#write.csv(as.data.frame(selection5),"/Users/victoraguiardesouzapenha/Desktop/BART.csv")
#write.csv(as.data.frame(selection7),"/Users/victoraguiardesouzapenha/Desktop/ANA.csv")
#write.csv(as.data.frame(selection9),"/Users/victoraguiardesouzapenha/Desktop/BBSL.csv")
#write.csv(as.data.frame(selection11), "/Users/victoraguiardesouzapenha/Desktop/NEOERL.csv")
###############################################################################################################
###############################################################################################################
###############################################################################################################
## Results
summary(selection2)
summary(selection4)
summary(selection6)
summary(selection8)
summary(selection10)
summary(selection12)

## Confidence intervals
coef_data <- as.data.frame(confint(selection2))[2,]
colnames(coef_data) <- c("low", "high")
coef_data$labeltext <- c("HII")
coef_data$Mean <- (coef_data$high + coef_data$low)/2
dim(coef_data)

first <- forestplot(
  title = "Response: PUUV",
  mean = coef_data$Mean,
  lower = coef_data$low,
  upper = coef_data$high,
  labeltext = coef_data$labeltext,
  txt_round = 2,
  zero = 0,  
  col = fpColors(box = "black", lines = "black", summary = "black")
)

###############################################################################################################
coef_data2 <- as.data.frame(confint(selection4))[c(2,3,5),]
colnames(coef_data2) <- c("low", "high")
coef_data2$labeltext <- c("Big_Spleen_Present","HII", "ANA_Present")
coef_data2$Mean <- (coef_data2$high + coef_data2$low)/2
dim(coef_data2)

second <- forestplot(
  title = "Response: BM",
  mean = coef_data2$Mean,
  lower = coef_data2$low,
  upper = coef_data2$high,
  labeltext = coef_data2$labeltext,
  txt_round = 2,
  zero = 0,  
  col = fpColors(box = "black", lines = "black", summary = "black")
)

###############################################################################################################
## Nothing for BRART

###############################################################################################################
coef_data4 <- as.data.frame(confint(selection8))[2:3,]
colnames(coef_data4) <- c("low", "high")
coef_data4$labeltext <- c("BM_Present","PCSoil_2")
coef_data4$Mean <- (coef_data4$high + coef_data4$low)/2
dim(coef_data4)

fourth <- forestplot(
  title = "Response: ANA",
  mean = coef_data4$Mean,
  lower = coef_data4$low,
  upper = coef_data4$high,
  labeltext = coef_data4$labeltext,
  txt_round = 2,
  zero = 0,  
  col = fpColors(box = "black", lines = "black", summary = "black")
)


###############################################################################################################
coef_data5 <- as.data.frame(confint(selection10))[2:3,]
colnames(coef_data5) <- c("low", "high")
coef_data5$labeltext <- c("HII","PCSoil_1")
coef_data5$Mean <- (coef_data5$high + coef_data5$low)/2
dim(coef_data5)

fifth <- forestplot(
  title = "Response: BBSL",
  mean = coef_data5$Mean,
  lower = coef_data5$low,
  upper = coef_data5$high,
  labeltext = coef_data5$labeltext,
  txt_round = 2,
  zero = 0,  
  col = fpColors(box = "black", lines = "black", summary = "black")
)

###############################################################################################################
coef_data6 <- as.data.frame(confint(selection12))[2,]
colnames(coef_data6) <- c("low", "high")
coef_data6$labeltext <- c("BBSL_present")
coef_data6$Mean <- (coef_data6$high + coef_data6$low)/2
dim(coef_data6)

sixth <- forestplot(
  title = "Response: NEOERL",
  mean = coef_data6$Mean,
  lower = coef_data6$low,
  upper = coef_data6$high,
  labeltext = coef_data6$labeltext,
  txt_round = 2,
  zero = 0,  
  col = fpColors(box = "black", lines = "black", summary = "black")
)

###############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
## Model fit
summary(selection2)
summary(selection4)
summary(selection6)
summary(selection8)
summary(selection10)
summary(selection12)

Metrics::rmse(as.integer(as.character(data3$PUUV)), predict(model))
Metrics::rmse(as.integer(as.character(data3$BM)), predict(model_BM))
Metrics::rmse(as.integer(as.character(data3$BART)), predict(model_BART))
Metrics::rmse(as.integer(as.character(data3$ANA)), predict(model_ANA))
Metrics::rmse(as.integer(as.character(data3$BBSL)), predict(model_BBSL))
Metrics::rmse(as.integer(as.character(data3$NEOERL)), predict(model_NEOERL))

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
## Plotting
coef_data$Response <- "PUUV"
coef_data2$Response <- "BM"
coef_data4$Response <- "ANA"
coef_data5$Response <- "BBSL"
coef_data6$Response <- "NEOERL"

last_one <- rbind(coef_data,coef_data2,coef_data4,coef_data5, coef_data6)
colnames(last_one)[c(1,2,4)] <- c("lower", "upper", "mean")
last_one$Response <- factor(last_one$Response)

## Reorder levels
last_one <- last_one[order(last_one$labeltext),]
unique(last_one$labeltext)
last_one$labeltext <- factor(last_one$labeltext, levels = c("ANA_Present", "BBSL_present", "BM_Present", "Big_Spleen_Present", 
                                                            "PCSoil_2", "PCSoil_1", "HII"))
# Create the forest plot
p <- ggplot(last_one, aes(x = mean, y = labeltext)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Response), height = 0.2, position = position_dodge(width = 0.5)) +
  theme_minimal() +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "95% Confidance intervals",
    y = "Variable",
    color = "Response variable"
  ) +
  scale_color_manual(values = c("blue", "red", "green", "purple","black","orange"))  +
  theme(legend.position="right") + 
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5, linewidth =1) +
  theme_classic() +theme(axis.text=element_text(size=12),
                         axis.text.y = element_text(size = 14),
                         axis.title=element_text(size=14,face="bold"))

p2 <- ggplot(last_one, aes(x = mean, y = labeltext)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "95% Confidence intervals",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_text(aes(label = Response, x = (lower + upper) / 2, y = labeltext), vjust = -1, color = "black")

p3 <- ggplot(last_one, aes(x = mean, y = labeltext)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "95% Confidence intervals",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_text_repel(aes(label = Response, x = (lower + upper) / 2, y = labeltext), 
                  direction = "y", 
                  force = 10)

p4 <- ggplot(last_one, aes(x = mean, y = labeltext)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge2(width = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "95% Confidence intervals",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_text_repel(aes(label = Response, x = (lower + upper) / 2, y = labeltext), 
                  direction = "y", 
                  force = 8)

p5 <- ggplot(last_one, aes(x = mean, y = labeltext)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "95% Confidence intervals",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_text_repel(aes(label = Response, x = (lower + upper) / 2, y = labeltext), 
                  direction = "y", 
                  force = 8)

jpeg("./Mod1_alt.jpeg", width = 10, height = 7, units = 'in', res = 800) 
p
dev.off()

## Another option would be to separate plots, instead of using the same one
p1 <- grid2grob(print(first))
p2 <- grid2grob(print(second))
p3 <- grid2grob(print(third))
p4 <- grid2grob(print(fourth))
p5 <- grid2grob(print(fifth))
p6 <- grid2grob(print(sixth))
p_both <- wrap_elements(p1) / wrap_elements(p2) / wrap_elements(p3) | wrap_elements(p4) / wrap_elements(p5) 
jpeg("./Mod1.jpeg", width = 8, height = 14, units = 'in', res = 800) 
p_both
dev.off()

## Table 2
hist(log(data3$HII))
data3$LogHII <- sqrt(data3$HII)

plt <- ggbetweenstats(
  data = data3,
  x = PUUV,
  y = LogHII, results.subtitle = F, ylab = "HII", xlab = "PUUV"
)

plt2 <- ggbetweenstats(
  data = data3,
  x = BM,
  y = LogHII, results.subtitle = F, ylab = "HII", xlab = "BM"
)

plt3 <- ggbetweenstats(
  data = data3,
  x = BBSL,
  y = LogHII, results.subtitle = F, ylab = "HII", xlab = "BBSL"
)

um <- ggplot(data3, aes(x=PUUV, y = log(HII))) + geom_boxplot() + theme_bw() + ylab("HII")
dois <- ggplot(data3, aes(x=BM, y = HII)) + geom_boxplot() + theme_bw() + ylab(NULL)
tres <- ggplot(data3, aes(x=BBSL, y = HII)) + geom_boxplot() + theme_bw() + ylab(NULL)

jpeg("./Mod2.jpeg", width = 9, height = 5, units = 'in', res = 800) 
#ggarrange(um,dois,tres,ncol=3)
ggarrange(plt,plt2,plt3, ncol = 3)
dev.off()
