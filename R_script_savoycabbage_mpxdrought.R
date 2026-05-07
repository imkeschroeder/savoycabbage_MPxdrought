#R script for graphs and models
# set up the file ####
## packages ####
library(readxl)
library(DHARMa) 
library(stats)
library(ggplot2)
library (tidyverse)
library(dbplyr)
library(broom)
library (car)
library(vegan)
library (ggthemes)

## set working repository ####
setwd("C:/Desktop/savoycabbage_MPxdrought")

## load in data ####
# biomass data
biomass <- read_excel(path = "data_savoycabbage_mpxdrought.xlsx",
                   sheet = "biomass")

# gas exchange data
GFS <- read_excel(path = "data_savoycabbage_mpxdrought.xlsx",
                  sheet = "gas_exchange")

# chlorophyll fluorescence data
elec <- read_excel(path = "data_savoycabbage_mpxdrought.xlsx",
                   sheet = "c_fluo")

# amino acid data including total concentration
aa <- read_excel(path = "data_savoycabbage_mpxdrought.xlsx",
                 sheet = "amino_acids")

# amino acid data transformed by adding 0.001 for modelling, without total concentration
aa_GLM <- read_excel(path = "data_savoycabbage_mpxdrought.xlsx",
                     sheet = "amino_acids_transformed")

# amino acid data for the NMDS
amino.acids <- read_excel(path = "data_savoycabbage_mpxdrought_NMDS.xlsx",
                          sheet = "absolute_abundance", range = cell_cols("B:T"))

# grouping information for the NMDS
amino.acids_groups <- read_excel(path = "data_savoycabbage_mpxdrought_NMDS.xlsx",
                          sheet = "groups")

# grouping information for the adonis
groups.statistic.all <- read_excel(path = "data_savoycabbage_mpxdrought_NMDS.xlsx",
                                                  sheet = "adonis")

# font and colours for the graphs
windowsFonts(`Arial` = windowsFont("Arial"))

custom_colors_8 <- c("#636363","#bdbdbd", "#fdae6b", "#a1d99b", 
                     "#9ecae1", "#e6550d", "#31a354", "#3182bd")

# plotting of the figures in the manuscript and supplement ####
# figure 2 ####
# CO2 assimilation rate A
assimilation <- ggplot(GFS, aes(x = treatment, y = A)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("A"~ "[µmol"~CO[2]~ m^-2~ s^-1*"]"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.961538, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 5), clip = "off")+
  scale_y_continuous(breaks = seq(0,5,1))
assimilation

# save the plot as an svg 
ggsave ("assimilation.svg", width = 4, height = 3)

# transpiration rate E
transpiration <- ggplot(GFS, aes(x = treatment, y = E)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("E"~ "[mmol"~H[2]*"O"~ m^-2~ s^-1*"]"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.115384, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.6), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.6,0.1))
transpiration

ggsave ("transpiration.svg", width = 4, height = 3)


# ETR
ETR <- ggplot(elec, aes(x = treatment, y = ETR)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("ETR"~ "[µmol"~ m^-2~ s^-1*"]"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = 24.230769, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(30, 60), clip = "off")+
  scale_y_continuous(breaks = seq(30,60,10))
ETR

ggsave ("ETR.svg", width = 4, height = 3)


# dry shoot biomass
shoot_dry <- ggplot(biomass, aes(x = treatment, y = shoot_dry)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Dry shoot biomass [g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.25, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 1.3), clip = "off")+
  scale_y_continuous(breaks = seq(0, 1.2, 0.3))
shoot_dry

ggsave ("shoot_dry.svg", width = 4, height = 3)


# dry root biomass
root_dry <- ggplot(biomass, aes(x = treatment, y = roots_dry)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Dry root biomass [g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.07692, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.4), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.4,0.1))
root_dry

ggsave ("root_dry.svg", width = 4, height = 3)


# root to shoot ratio
rts <- ggplot(biomass, aes(x = treatment, y = root_to_shoot_ratio)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Root-to-shoot ratio", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.153846, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.8), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.8,0.2))
rts

ggsave ("root_to_shoot.svg", width = 4, height = 3)


# figure 3 ####
## calculate NMDS ####
# prepare NMDS
# Wisconsin square root transformation
amino.acids.norm <-         
  wisconsin(sqrt(amino.acids))

# calculating distance matrix
amino.acids_distmat <- 
  vegdist(amino.acids.norm, method = "kulczynski")

# Creating easy to view matrix and writing .csv
amino.acids_distmat <- 
  as.matrix(amino.acids_distmat, labels = T)
write.csv(amino.acids_distmat, "amino.acids_distmat.csv")

# run NMDS using metaMDS
# Running NMDS in vegan (metaMDS)
amino.acids_NMS <-
  metaMDS(amino.acids_distmat,
          autotransform = FALSE,
          binary = FALSE,
          distance = "kulczynski",
          try = 1000,
          trymax = 10000,
          maxit = 10000,
          k = 2,
          PC = TRUE,
          center = TRUE)
amino.acids_NMS$stress
#stress value: 0.1380558
stressplot(amino.acids_NMS)

# calculate betadisper and adonis

# test of homogeneity of the multivariate variance spread among groups
# attach treatment to ws transformed data
data_ws <-         
  wisconsin(sqrt(amino.acids))
attach(treatments <- read_excel(path = "data_savoycabbage_mpxdrought_NMDS.xlsx",
                                sheet = 4, range = cell_cols("B:D")))

mod <- betadisper(vegdist(data_ws, method="kulczynski"), treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)
boxplot(mod)
# p < 0.001 ***

#adonis
adonis.results.all <-
  adonis2(formula = amino.acids_distmat ~ water*plastic, groups.statistic.all, permutations = 10000,
          method = "kulczynski")
print(adonis.results.all)

#adonis results
#factor water: < 0.001
#factor plastic 0.069
#factor water:plastic 0.197


## plot NMDS A ####
svg("NMDS_with_loadings.svg", width = 6, height = 5) #opens svg device

ordiplot(amino.acids_NMS, display = "sites", type= "n") 

#Create convex hulls that highlight point clusters based on grouping dataframe
ordihull(
  amino.acids_NMS,
  amino.acids_groups$group,
  display = "sites",
  draw = c("polygon"),
  col = c("#636363", "#e6550d", "#31a354", "#3182bd", 
          "#bdbdbd", "#fdae6b", "#a1d99b", "#9ecae1"),
  border = c("#636363", "#e6550d", "#31a354", "#3182bd", 
             "#bdbdbd", "#fdae6b", "#a1d99b", "#9ecae1"),
  alpha = 70,
  lty = c(1),
)

# Calculating and plotting centroids of NMDS Result

scrs <-
  scores(amino.acids_NMS, display = "sites", "species")
cent <-
  aggregate(scrs ~ group, data = amino.acids_groups, FUN = "median")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c(21, 22, 23, 24, 21, 22, 23, 24),
       col = c("black"),
       bg = c("#636363", "#e6550d", "#31a354", "#3182bd", 
              "#bdbdbd", "#fdae6b", "#a1d99b", "#9ecae1"),
       lwd = 2,
       cex = 1.4
)

# add legend with the adonis results
legend ("bottomright",
        legend =c("stress 0.138"), 
        #"adonis", 
        # bquote ("MP" ~ italic("p") ~ "= 0.069"), 
        #bquote ("W" ~ italic("p") ~ "< 0.001"), 
        # bquote ("MPxW" ~ italic("p") ~ "= 0.197")),
        bty = "n", cex = 1)

# add legend with treatment groups
legend ("topright", 
        xpd=TRUE,
        legend = c("Control w", "Control d", "PET w", "PET d", 
                   "PLA w", "PLA d", "RC w", "RC d"),
        bty = "n",
        cex = 1,
        pch = c(21, 21, 22, 22, 23, 23, 24, 24),
        col = "black",
        pt.bg = c("#636363", "#bdbdbd", "#e6550d", "#fdae6b",
                  "#31a354", "#a1d99b", "#3182bd", "#9ecae1")
)

# display loadings
en = envfit(amino.acids_NMS, amino.acids, permutations = 10000, na.rm = TRUE)

plot(en, col = "black", cex = 0.8)

dev.off()  # Closes the SVG device and saves the file

# clean up overlapping text in inkscape

## plot NMDS B ####
# Simplified NMDS with contour lines of the correlation 
# of the dry shoot biomass with the amino acid profile

contour <- read_excel(path = "data_savoycabbage_mpxdrought_NMDS.xlsx",
                      sheet = "contour_lines")


db <- contour$shoot_dry 

#farbvektoren
pchvec <- c(21, 22, 23, 24, 21, 22, 23, 24)
colvec <- c("#636363", "#e6550d", "#31a354", "#3182bd", 
            "#bdbdbd", "#fdae6b", "#a1d99b", "#9ecae1")

svg("NMDS_contour.svg", width = 6, height = 5) #opens a svg device

ordiplot(amino.acids_NMS, display = "sites", type= "n") 
with(amino.acids_groups,
     points(amino.acids_NMS,
            display = "sites",
            col = c("black"),
            pch = pchvec[group],
            bg = colvec[group],
            cex = 1))
cl_db <- ordisurf(amino.acids_NMS, db, 
                  bs = "tp",   #thin plate regression splines
                  family = gaussian(link = "identity"), #Gaussian error distribution with link identity
                  add = TRUE, #add to existing plot
                  col = "black", #line colour
                  method = "P-REML",  #restricted maximum likelihood estimation
                  lwd.cl = 1.5,
                  nlevels = 10,
                  labcex = 0.9
)


legend ("bottomleft",
        legend =c("Deviance explained: 29.2%"),
        bty = "n", cex = 1)

dev.off() #closes svg device and saves file

summary(cl_db)
# p-value < 0.001

## total amino acid concentration ####
aa_total <- ggplot(aa, aes(x = treatment, y = total)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Total amino acid concentration [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -5.384615, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(10, 90), clip = "off")+
  scale_y_continuous(breaks = seq(10,90,20))
aa_total

ggsave ("aa_total.svg", width = 4, height = 3)


# figure 4 ####
# GABA
aa_gaba <- ggplot(aa, aes(x = treatment, y = GABA)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "GABA [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -1.346153, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 7), clip = "off")+
  scale_y_continuous(breaks = seq(0, 7, 2))
aa_gaba

ggsave ("aa_gaba.svg", width = 3, height = 3)


# Hydroxyproline (hyp)
aa_hyp <- ggplot(aa, aes(x = treatment, y = Hyp)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Hydroxyproline [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.01923, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.1), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.1,0.02))
aa_hyp

ggsave ("aa_hyp.svg", width = 3, height = 3)

# Leucine (leu)
aa_leu <- ggplot(aa, aes(x = treatment, y = Leu)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Leucine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.576923, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 3), clip = "off")+
  scale_y_continuous(breaks = seq(0,3,1))
aa_leu

ggsave ("aa_leu.svg", width = 3, height = 3)

# Lysine (Lys)
aa_lys <- ggplot(aa, aes(x = treatment, y = Lys)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Lysine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.288461, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 1.5), clip = "off")+
  scale_y_continuous(breaks = seq(0,1.5,0.5))
aa_lys

ggsave ("aa_lys.svg", width = 3, height = 3)

# Methionine (Met)
aa_met <- ggplot(aa, aes(x = treatment, y = Met)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Methionine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.0576923, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.3), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.3,0.1))
aa_met

ggsave ("aa_met.svg", width = 3, height = 3)

# Phenylalanine (Phe)
aa_phe <- ggplot(aa, aes(x = treatment, y = Phe)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Phenylalanine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.480769, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 2.5), clip = "off")+
  scale_y_continuous(breaks = seq(0,2.5,0.5))
aa_phe

ggsave ("aa_phe.svg", width = 3, height = 3)

# Proline (Pro)
aa_pro <- ggplot(aa, aes(x = treatment, y = Pro)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Proline [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.76923, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 4), clip = "off")+
  scale_y_continuous(breaks = seq(0,4,1))
aa_pro

ggsave ("aa_pro.svg", width = 3, height = 3)

# Tryptophan (Trp)
aa_trp <- ggplot(aa, aes(x = treatment, y = Trp)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Tryptophan [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.288461, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 1.5), clip = "off")+
  scale_y_continuous(breaks = seq(0, 1.5, 0.5))
aa_trp

ggsave ("aa_trp.svg", width = 3, height = 3)

# Tyrosine (Tyr)
aa_tyr <- ggplot(aa, aes(x = treatment, y = Tyr)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Tyrosine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.384615, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 2), clip = "off")+
  scale_y_continuous(breaks = seq(0, 2, 0.5))
aa_tyr

ggsave ("aa_tyr.svg", width = 3, height = 3)


# figure S1 ####
# number of leaves at 49 days post sowing
nl <- ggplot(biomass, aes(x = treatment, y = number_leaves_dps49)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0.25, shape = 20, size = 2.5) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("number of leaves at 49 dps"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = 1.05, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(2, 7), clip = "off")+
  scale_y_continuous(breaks = seq(2,7,1))
nl

ggsave ("nl.svg", width = 5, height = 3)

# figure S2 ####
# stomatal conductance to water vapor gs
gh2o <- ggplot(GFS, aes(x = treatment, y = GH2O)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("g"[s]~ "[mmol"~H[2]*"O"~ m^-2~ s^-1*"]"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -11.5384, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 60), clip = "off")+
  scale_y_continuous(breaks = seq(0,60,20))
gh2o

ggsave ("gh2o.svg", width = 4, height = 3)

# intercellular CO2 concentration ci
ci <- ggplot(GFS, aes(x = treatment, y = ci)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = bquote("c"[i]~ "[ppm]"), x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -76.92, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 400), clip = "off")+
  scale_y_continuous(breaks = seq(0,400,100))
ci

ggsave ("ci.svg", width = 4, height = 3)

# figure S3 ####
# Asparagine (Asn)
aa_asn <- ggplot(aa, aes(x = treatment, y = Asn)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Asparagine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.576923, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 3), clip = "off")+
  scale_y_continuous(breaks = seq(0,3,1))
aa_asn

ggsave ("aa_asn.svg", width = 3, height = 3)

# Aspartic acid (Asp)
aa_asp <- ggplot(aa, aes(x = treatment, y = Asp)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Aspartic acid [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -6.73076, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 35), clip = "off")+
  scale_y_continuous(breaks = seq(0,35,10))
aa_asp

ggsave ("aa_asp.svg", width = 3, height = 3)

# Citrulline (Cit)
aa_cit <- ggplot(aa, aes(x = treatment, y = Cit)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Citrulline [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.048076, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 0.25), clip = "off")+
  scale_y_continuous(breaks = seq(0,0.25,0.05))
aa_cit

ggsave ("aa_cit.svg", width = 3, height = 3)

# Glutamic acid (Glu)
aa_glu <- ggplot(aa, aes(x = treatment, y = Glu)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Glutamic acid [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = 2.88461538, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(5, 16), clip = "off")+
  scale_y_continuous(breaks = seq(5,16,2.5))
aa_glu

ggsave ("aa_glu.svg", width = 3, height = 3)

# Glutamine (Gln)
aa_gln <- ggplot(aa, aes(x = treatment, y = Gln)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Glutamine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -1.6346153, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 8.5), clip = "off")+
  scale_y_continuous(breaks = seq(0,8.5,2))
aa_gln

ggsave ("aa_gln.svg", width = 3, height = 3)

# Histidine (His)
aa_his <- ggplot(aa, aes(x = treatment, y = His)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Histidine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.0153846, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0.1, 0.7), clip = "off")+
  scale_y_continuous(breaks = seq(0.1,0.7,0.1))
aa_his

ggsave ("aa_his.svg", width = 3, height = 3)

# Isoleucine (Ile)
aa_ile <- ggplot(aa, aes(x = treatment, y = Ile)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Isoleucine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.288461538, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 1.5), clip = "off")+
  scale_y_continuous(breaks = seq(0,1.5,0.5))
aa_ile

ggsave ("aa_ile.svg", width = 3, height = 3)

# Serine (Ser)
aa_ser <- ggplot(aa, aes(x = treatment, y = Ser)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Serine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -2.88461538, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 15), clip = "off")+
  scale_y_continuous(breaks = seq(0, 15, 5))
aa_ser

ggsave ("aa_ser.svg", width = 3, height = 3)

# Threonine (Thr)
aa_thr <- ggplot(aa, aes(x = treatment, y = Thr)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Threonine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.76923076, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 4), clip = "off")+
  scale_y_continuous(breaks = seq(0, 4, 1))
aa_thr

ggsave ("aa_thr.svg", width = 3, height = 3)

# Valine (Val)
aa_val <- ggplot(aa, aes(x = treatment, y = Val)) +
  geom_boxplot(aes(fill=treatment), 
               outlier.colour=NA,
               staplewidth = 0.5,
               width = 0.5)  +
  scale_x_discrete(limits=c("C","dC","PET", "dPET", "PLA", "dPLA", "RC", "dRC"),
                   labels = c("w", "d", "w", "d", "w", "d", "w", "d")) +
  geom_jitter(height = 0, width = 0, shape = 20) +
  stat_summary(fun=mean, geom="point", shape=18, size=4) +
  scale_fill_manual(values = custom_colors_8, guide = "none") +
  theme_few(base_family = "Arial",
            base_size = 13) +
  theme(axis.text = element_text(size = 12)) +
  labs(y = "Valine [µmol/g]", x = "") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -0.67307692, 
           label = c("Control", "PET", "PLA", "RC"), 
           fontface = "bold", size = 4)+
  coord_cartesian(ylim = c(0, 3.5), clip = "off")+
  scale_y_continuous(breaks = seq(0, 3.5, 1))
aa_val

ggsave ("aa_val.svg", width = 3, height = 3)

# statistics/modelling (table 1) ####
# data was first checked with Shapiro-Wilk test and Levene test
# Data was modelled with LMs or GLMs
# Model residuals were tested for normal distribution and
# variance homogeneity using diagnostic plots of the 
# DHARMa package
# p values were derived with the anova function

## CO2 assimilation rate A ####
# Shapiro-Wilk test
df.shapiro <- GFS %>%
  mutate(ratio_log = A) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(GFS$A, GFS$treatment)

# linear model
A.M1 <- lm(A ~ plastic * water, data = GFS)

Dharma_A.M1 <- simulateResiduals(A.M1, 
                                 n = 1000, 
                                 plot = "quantile")
plot(Dharma_A.M1)

testDispersion(Dharma_A.M1)
testZeroInflation(Dharma_A.M1)

r <- residuals(A.M1, type = "pearson")
hist(r)

# anova
anova(A.M1)

# MP F = 0.78; p = 0.514
# W F = 3.85; p = 0.057
# MPxW F = 1.23; p = 0.311

## transpiration rate E ####
# Shapiro-Wilk test
df.shapiro <- GFS %>%
  mutate(ratio_log = E) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(GFS$E, GFS$treatment)

# linear model
E.M1 <- lm(E ~ plastic * water, data = GFS)

Dharma_E.M1 <- simulateResiduals(E.M1, 
                                 n = 1000, 
                                 plot = "quantile")
plot(Dharma_E.M1)

testDispersion(Dharma_E.M1)
testZeroInflation(Dharma_E.M1)

r <- residuals(E.M1, type = "pearson")
hist(r)

# anova
anova(E.M1)

# MP F = 1.25; p = 0.305
# W F = 6.09; p = 0.018
# MPxW F = 1.24; p = 0.306


## ETR ####
# Shapiro-Wilk test
df.shapiro <- elec %>%
  mutate(ratio_log = ETR) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(elec$ETR, elec$treatment)

# linear model
ETR.M1 <- lm(ETR ~ plastic * water, data = elec)

Dharma_ETR.M1 <- simulateResiduals(ETR.M1, 
                                 n = 1000, 
                                 plot = "quantile")
plot(Dharma_ETR.M1)

testDispersion(Dharma_ETR.M1)
testZeroInflation(Dharma_ETR.M1)

r <- residuals(ETR.M1, type = "pearson")
hist(r)

# anova
anova(ETR.M1)

# MP F = 0.88; p = 0.458
# W F = 3.58; p = 0.063
# MPxW F = 0.30; p = 0.827

## stomatal conductance to water vapor gs ####
# Shapiro-Wilk test
df.shapiro <- GFS %>%
  mutate(ratio_log = GH2O) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(GFS$GH2O, GFS$treatment)

# linear model
G.M1 <- lm(GH2O ~ plastic * water, data = GFS)

Dharma_G.M1 <- simulateResiduals(G.M1, 
                                 n = 1000, 
                                 plot = "quantile")
plot(Dharma_G.M1)

testDispersion(Dharma_G.M1)
testZeroInflation(Dharma_G.M1)

r <- residuals(G.M1, type = "pearson")
hist(r)

# anova
anova(G.M1)

# MP F = 1.25; p = 0.305
# W F = 6.17; p = 0.017
# MPxW F = 1.25; p = 0.305

## intercellular CO2 concentration ci ####
# Shapiro-Wilk test
df.shapiro <- GFS %>%
  mutate(ratio_log = ci) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(GFS$ci, GFS$treatment)

# linear model
ci.M1 <- lm(ci ~ plastic * water, data = GFS)

Dharma_ci.M1 <- simulateResiduals(ci.M1, 
                                 n = 1000, 
                                 plot = "quantile")
plot(Dharma_ci.M1)

testDispersion(Dharma_ci.M1)
testZeroInflation(Dharma_ci.M1)

r <- residuals(ci.M1, type = "pearson")
hist(r)

# anova
anova(ci.M1)

# MP F = 1.31; p = 0.284
# W F = 0.05; p = 0.818
# MPxW F = 0.23; p = 0.873

## dry shoot biomass ####
# Shapiro-Wilk test
df.shapiro <- biomass %>%
  mutate(ratio_log = shoot_dry) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(biomass$shoot_dry, biomass$treatment)

# linear model
ds.M1 <- lm(shoot_dry ~ plastic * water, 
             data = biomass)

Dharma_ds.M1 <- simulateResiduals(ds.M1, 
                                  n = 1000, 
                                  plot = "quantile")
plot(Dharma_ds.M1)

testDispersion(Dharma_ds.M1)
testZeroInflation(Dharma_ds.M1)

r <- residuals(ds.M1, type = "pearson")
hist(r)

# anova
anova(ds.M1)

# MP F = 58.14; p < 0.001
# W F = 142.35; p < 0.001
# MPxW F = 11.99; p < 0.001

## dry root biomass ####
# Shapiro-Wilk test
df.shapiro <- biomass %>%
  mutate(ratio_log = roots_dry) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(biomass$roots_dry, biomass$treatment)

# linear model
dr.M1 <- lm(roots_dry ~ plastic * water, data = biomass)

Dharma_dr.M1 <- simulateResiduals(dr.M1, 
                                  n = 1000, 
                                  plot = "quantile")
plot(Dharma_dr.M1)

testDispersion(Dharma_dr.M1)
testZeroInflation(Dharma_dr.M1)

r <- residuals(dr.M1, type = "pearson")
hist(r)

# anova
anova(dr.M1)

# MP F = 45.91; p < 0.001
# W F = 47.61; p < 0.001
# MPxW F = 3.46; p = 0.021

## root to shoot ratio ####
# Shapiro-Wilk test
df.shapiro <- biomass %>%
  mutate(ratio_log = root_to_shoot_ratio) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(biomass$root_to_shoot_ratio, biomass$treatment)

# linear model
rts.M1 <- lm(root_to_shoot_ratio ~ plastic * water, data = biomass)

Dharma_rts.M1 <- simulateResiduals(rts.M1, n = 1000, plot = "quantile")
plot(Dharma_rts.M1)

testDispersion(Dharma_rts.M1)
testZeroInflation(Dharma_rts.M1)

r <- residuals(rts.M1, type = "pearson")
hist(r)

# anova
anova(rts.M1)

# MP F = 6.82; p < 0.001
# W F = 29.54; p < 0.001
# MPxW F = 2.92; p = 0.040

## total amino acid concentration ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = total) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$total, aa$treatment)

# linear model
taa.M1 <- lm(total ~ plastic * water, data = aa)

Dharma_taa.M1 <- simulateResiduals(taa.M1, n = 1000, plot = "quantile")
plot(Dharma_taa.M1)

testDispersion(Dharma_taa.M1)
testZeroInflation(Dharma_taa.M1)

r <- residuals(taa.M1, type = "pearson")
hist(r)

# anova
anova(taa.M1)

# MP F = 0.75; p = 0.528
# W F = 4.55; p = 0.037
# MPxW F = 0.27; p = 0.851

## asparagine (asn) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Asn) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Asn, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
asn.M3 <- glm(Asn ~ plastic * water, 
              family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_asn.M3 <- simulateResiduals(asn.M3, n = 1000, plot = "quantile")
plot(Dharma_asn.M3)

testDispersion(Dharma_asn.M3)
testZeroInflation(Dharma_asn.M3)

r <- residuals(asn.M3, type = "pearson")
hist(r)

#anova
anova(asn.M3)

# MP X2 = 0.87; p = 0.460
# W X2 = 11.62; p = 0.001
# MPxW X2 = 1.00; p = 0.398

## aspartic acid (asp) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Asp) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Asp, aa$treatment)

# linear model
asp.M1 <- lm(Asp ~ plastic * water, data = aa_GLM)

Dharma_asp.M1 <- simulateResiduals(asp.M1, n = 1000, plot = "quantile")
plot(Dharma_asp.M1)

testDispersion(Dharma_asp.M1)
testZeroInflation(Dharma_asp.M1)

r <- residuals(asp.M1, type = "pearson")
hist(r)

# anova
anova(asp.M1)

# MP F = 2.56; p = 0.062
# W F = 0; p = 0.995
# MPxW F = 0.27; p = 0.846

## gamma-aminobutyric acid (gaba) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = GABA) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$GABA, aa$treatment)

# glm with Gamma distribution and log link function
GABA.M2 <- glm(GABA ~ plastic * water, family = Gamma(link = "log"),
               data = aa_GLM)

Dharma_GABA.M2 <- simulateResiduals(GABA.M2, n = 1000, plot = "quantile")
plot(Dharma_GABA.M2)

testDispersion(Dharma_GABA.M2)
testZeroInflation(Dharma_GABA.M2)
r <- residuals(GABA.M2, type = "pearson")
hist(r)

# anova
anova(GABA.M2)

# MP X2 = 7.98; p < 0.001
# W X2 = 3.49; p = 0.066
# MPxW X2 = 4.09; p = 0.010

## glutamic acid (glu) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Glu) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Glu, aa$treatment)

# linear model
glu.M1 <- lm(Glu ~ plastic * water, data = aa_GLM)

Dharma_glu.M1 <- simulateResiduals(glu.M1, n = 1000, plot = "quantile")
plot(Dharma_glu.M1)

testDispersion(Dharma_glu.M1)
testZeroInflation(Dharma_glu.M1)

r <- residuals(glu.M1, type = "pearson")
hist(r)

# anova
anova(glu.M1)

# MP F = 1.23; p = 0.308
# W F = 4.39; p = 0.040
# MPxW F = 0.20; p = 0.896

## glutamine (gln) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Gln) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Gln, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
gln.M3 <- glm(Gln ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_gln.M3 <- simulateResiduals(gln.M3, n = 1000, plot = "quantile")
plot(Dharma_gln.M3)

testDispersion(Dharma_gln.M3)
testZeroInflation(Dharma_gln.M3)

r <- residuals(gln.M3, type = "pearson")
hist(r)

#anova
anova(gln.M3)

# MP X2 = 1.92; p = 0.136
# W X2 = 46.58; p < 0.001
# MPxW X2 = 2.50; p = 0.067

## histidine (his) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = His) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$His, aa$treatment)

# linear model
his.M1 <- lm(His ~ plastic * water, data = aa_GLM)

Dharma_his.M1 <- simulateResiduals(his.M1, n = 1000, plot = "quantile")
plot(Dharma_his.M1)

testDispersion(Dharma_his.M1)
testZeroInflation(Dharma_his.M1)

r <- residuals(his.M1, type = "pearson")
hist(r)

# anova
anova(his.M1)

# MP F = 1.16; p = 0.333
# W F = 4.40; p = 0.040
# MPxW F = 1.46; p = 0.235

## isoleucine (ile) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Ile) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Ile, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
ile.M3 <- glm(Ile ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_ile.M3 <- simulateResiduals(ile.M3, n = 1000, plot = "quantile")
plot(Dharma_ile.M3)

testDispersion(Dharma_ile.M3)
testZeroInflation(Dharma_ile.M3)

r <- residuals(ile.M3, type = "pearson")
hist(r)

#anova
anova(ile.M3)

# MP X2 = 0.73; p = 0.540
# W X2 = 16.74; p < 0.001
# MPxW X2 = 2.18; p = 0.099

## leucine (leu) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Leu) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Leu, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
leu.M3 <- glm(Leu ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_leu.M3 <- simulateResiduals(leu.M3, n = 1000, plot = "quantile")
plot(Dharma_leu.M3)

testDispersion(Dharma_leu.M3)
testZeroInflation(Dharma_leu.M3)

r <- residuals(leu.M3, type = "pearson")
hist(r)

#anova
anova(leu.M3)

# MP X2 = 2.62; p = 0.058
# W X2 = 6.18; p = 0.016
# MPxW X2 = 5.71; p = 0.002

## lysine (lys) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Lys) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Lys, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
lys.M3 <- glm(Lys ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_lys.M3 <- simulateResiduals(lys.M3, n = 1000, plot = "quantile")
plot(Dharma_lys.M3)

testDispersion(Dharma_lys.M3)
testZeroInflation(Dharma_lys.M3)

r <- residuals(lys.M3, type = "pearson")
hist(r)

#anova
anova(lys.M3)

# MP X2 = 3.51; p = 0.020
# W X2 = 8.74; p = 0.004
# MPxW X2 = 3.01; p = 0.037

## methionine (met) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Met) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Met, aa$treatment)

# linear model
met.M1 <- lm(Met ~ plastic * water, data = aa_GLM)

Dharma_met.M1 <- simulateResiduals(met.M1, n = 1000, plot = "quantile")
plot(Dharma_met.M1)

testDispersion(Dharma_met.M1)
testZeroInflation(Dharma_met.M1)

r <- residuals(met.M1, type = "pearson")
hist(r)

# anova
anova(met.M1)

# MP F = 2.87; p = 0.043
# W F = 0.26; p = 0.611
# MPxW F = 2.31; p = 0.085

## phenylalanine (phe) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Phe) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Phe, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
phe.M3 <- glm(Phe ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_phe.M3 <- simulateResiduals(phe.M3, n = 1000, plot = "quantile")
plot(Dharma_phe.M3)

testDispersion(Dharma_phe.M3)
testZeroInflation(Dharma_phe.M3)

r <- residuals(phe.M3, type = "pearson")
hist(r)

#anova
anova(phe.M3)

# MP X2 = 1.03; p = 0.384
# W X2 = 3.28; p = 0.075
# MPxW X2 = 4.67; p = 0.005

## serine (ser) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Ser) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Ser, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
ser.M3 <- glm(Ser ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_ser.M3 <- simulateResiduals(ser.M3, n = 1000, plot = "quantile")
plot(Dharma_ser.M3)

testDispersion(Dharma_ser.M3)
testZeroInflation(Dharma_ser.M3)

r <- residuals(ser.M3, type = "pearson")
hist(r)

#anova
anova(ser.M3)

# MP X2 = 2.00; p = 0.123
# W X2 = 4.81; p = 0.032
# MPxW X2 = 0.13; p = 0.939

## threonine (thr) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Thr) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Thr, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
thr.M3 <- glm(Thr ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_thr.M3 <- simulateResiduals(thr.M3, n = 1000, plot = "quantile")
plot(Dharma_thr.M3)

testDispersion(Dharma_thr.M3)
testZeroInflation(Dharma_thr.M3)

r <- residuals(thr.M3, type = "pearson")
hist(r)

#anova
anova(thr.M3)

# MP X2 = 0.61; p = 0.612
# W X2 = 9.24; p = 0.003
# MPxW X2 = 2.12; p = 0.107

## tryptophan (trp) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Trp) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Trp, aa$treatment)

# glm with Gamma distribution and log link function
trp.M2 <- glm(Trp ~ plastic * water, family = Gamma(link = "log"),
               data = aa_GLM)

Dharma_trp.M2 <- simulateResiduals(trp.M2, n = 1000, plot = "quantile")
plot(Dharma_trp.M2)

testDispersion(Dharma_trp.M2)
testZeroInflation(Dharma_trp.M2)
r <- residuals(trp.M2, type = "pearson")
hist(r)

# anova
anova(trp.M2)

# MP X2 = 1.33; p = 0.273
# W X2 = 4.51; p = 0.038
# MPxW X2 = 5.21; p = 0.003

## tyrosine (tyr) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Tyr) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Tyr, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
tyr.M3 <- glm(Tyr ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_tyr.M3 <- simulateResiduals(tyr.M3, n = 1000, plot = "quantile")
plot(Dharma_tyr.M3)

testDispersion(Dharma_tyr.M3)
testZeroInflation(Dharma_tyr.M3)

r <- residuals(tyr.M3, type = "pearson")
hist(r)

#anova
anova(tyr.M3)

# MP X2 = 1.50; p = 0.224
# W X2 = 0; p = 0.979
# MPxW X2 = 3.58; p = 0.019

## valine (val) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Val) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Val, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
val.M3 <- glm(Val ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_val.M3 <- simulateResiduals(val.M3, n = 1000, plot = "quantile")
plot(Dharma_val.M3)

testDispersion(Dharma_val.M3)
testZeroInflation(Dharma_val.M3)

r <- residuals(val.M3, type = "pearson")
hist(r)

#anova
anova(val.M3)

# MP X2 = 1.49; p = 0.226
# W X2 = 14.75; p < 0.001
# MPxW X2 = 2.60; p = 0.060

## hydroxyproline (hyp) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Hyp) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Hyp, aa$treatment)

# glm with Gamma distribution and log link function
hyp.M2 <- glm(Hyp ~ plastic * water, family = Gamma(link = "log"),
              data = aa_GLM)

Dharma_hyp.M2 <- simulateResiduals(hyp.M2, n = 1000, plot = "quantile")
plot(Dharma_hyp.M2)

testDispersion(Dharma_hyp.M2)
testZeroInflation(Dharma_hyp.M2)
r <- residuals(hyp.M2, type = "pearson")
hist(r)

# anova
anova(hyp.M2)

# MP X2 = 3.70; p = 0.016
# W X2 = 44.39; p < 0.001
# MPxW X2 = 1.15; p = 0.334

## proline (pro) ####
# Shapiro-Wilk test
df.shapiro <- aa %>%
  mutate(ratio_log = Pro) %>%
  group_by(treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ratio_log)))
df.shapiro.glance <- df.shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
df.shapiro.glance

# Levenes test
leveneTest(aa$Pro, aa$treatment)

# glm with inverse gaussian distribution and 1 mu-2 link function 
pro.M3 <- glm(Pro ~ plastic * water, family = inverse.gaussian(link = "1/mu^2"),
              data = aa_GLM)

Dharma_pro.M3 <- simulateResiduals(pro.M3, n = 1000, plot = "quantile")
plot(Dharma_pro.M3)

testDispersion(Dharma_pro.M3)
testZeroInflation(Dharma_pro.M3)

r <- residuals(pro.M3, type = "pearson")
hist(r)

#anova
anova(pro.M3)

# MP X2 = 0.12; p = 0.945
# W X2 = 43.80; p < 0.001
# MPxW X2 = 2.32; p = 0.083

# correlation tests ####
# correlate photosynthesis-related traits
## including all MP treatment groups ####
# gs vs A
cor.test(GFS$GH2O, GFS$A, method="pearson")
# p-value < 0.001

# ci vs A
cor.test(GFS$ci, GFS$A, method="pearson")
# p-value 0.832

# ETR vs A
cor.test(GFS$ETR_cor, GFS$A, method="pearson")
# 0.577

## including only single MP treatment groups ####
# subset data
GFS_C <- GFS[c(1:6, 25:30),]
GFS_PET <- GFS[c(7:12, 31:36),]
GFS_PLA <- GFS[c(13:18, 37:42),]
GFS_RC <- GFS[c(19:24, 43:48),]

## including only Control MP treatment group
# gs vs A
cor.test(GFS_C$GH2O, GFS_C$A, method="pearson")
# p-value < 0.001

# ci vs A
cor.test(GFS_C$ci, GFS_C$A, method="pearson")
# p-value 0.490

# ETR vs A
cor.test(GFS_C$ETR_cor, GFS_C$A, method="pearson")
# 0.764

## including only PET MP treatment group
# gs vs A
cor.test(GFS_PET$GH2O, GFS_PET$A, method="pearson")
# p-value < 0.001

# ci vs A
cor.test(GFS_PET$ci, GFS_PET$A, method="pearson")
# p-value 0.594

# ETR vs A
cor.test(GFS_PET$ETR_cor, GFS_PET$A, method="pearson")
# p-value 0.809

## including only PLA MP treatment group
# gs vs A
cor.test(GFS_PLA$GH2O, GFS_PLA$A, method="pearson")
# p-value < 0.001

# ci vs A
cor.test(GFS_PLA$ci, GFS_PLA$A, method="pearson")
# p-value 0.411

# ETR vs A
cor.test(GFS_PLA$ETR_cor, GFS_PLA$A, method="pearson")
# 0.023

## including only RC MP treatment group
# gs vs A
cor.test(GFS_RC$GH2O, GFS_RC$A, method="pearson")
# p-value < 0.001

# ci vs A
cor.test(GFS_RC$ci, GFS_RC$A, method="pearson")
# p-value 0.255

# ETR vs A
cor.test(GFS_RC$ETR_cor, GFS_RC$A, method="pearson")
# 0.272
