######################################
#        Clustering Analyses         #
######################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/LandGen_AridMammals/")

############
# Packages #
############

library(tess3r)
library(dartR)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(ggsn)
library(tibble)
library(scatterpie)
library(ggpubr)
library(colorblindr)

################################################
#      SET Variables AND LOAD IN DATA          #
################################################

# Paths
Out.Path <- "IBB/IBB_outputs"
gl.Path <- "SNP_Filtering/SampleClean_outputs"
str.Path <- "IBB/Data"

# Species genlights
sp.gls <- list.files(gl.Path, full.names = TRUE, pattern = ".rdata")

# Species ID
sp <- gsub(".rdata", "", gsub(paste0(gl.Path, "/RelClean.gl."), "", sp.gls))

# Species microsat str files
sp.str <- list.files(str.Path, full.names = TRUE, pattern = ".str")

# Number of Msat loci (in ordr of sp)
MSATS <- c(12, 14, 14)

# Choose maximum value of K to test, using the number of PCoA clusters plus 2
# I'm just going to test the same number to keep it easy across species
# Nothing appeared to have more than ~5 clusters max, sp I'll use 7
Kmax <- 7

################################################
# Create a loop that runs through tess         #
# for each species, then outputs K plots       #
################################################

for (i in 1:length(sp)) {

# This script takes a genlight as the input format
# I output my genlight as an R object after filtering loci and
# removing duplicated and related individuals
# The corresponding meta-data is saved in the "other" slot
# So load and rename R object

load(sp.gls[i])
gl <- sp.gl

# load in microsat data as genind
msat <- read.structure(file = sp.str[i], n.ind = nInd(gl), n.loc = MSATS[i], onerowperind = TRUE, col.lab = 1, row.marknames = 1, col.pop = 0, col.others = 0)


###################
# Package: Tess3R #
###################

# RUN FOR SNPs
# Format for Tess3R
# I will run both data sets to make sure I get the same results
Geno.Tess3R <- as.matrix(gl)
Coords.Tess3R <- as.matrix(gl@other$latlong[, c(2,1)])
Pops.Tess3R <- as.character(gl@pop)


## RUN TESS3R ##
Tess.obj <- tess3(X = Geno.Tess3R, coord = Coords.Tess3R, K = 1:Kmax, rep = 100, method = "projected.ls", ploidy = 2, openMP.core.num = 4, max.iteration = 1000, keep = "best", mask = 0.1) 

# Save Tess object
save(Tess.obj, file = paste0(Out.Path, "Tess.obj.FullDat.best100reps.K1_", Kmax, ".", sp[i], "_SNPs.R"))

# RUN FOR MICROSATS
# Format for Tess3R
# I will run both data sets to make sure I get the same results
msat <- gi2gl(msat)
Geno.Tess3R.m <- as.matrix(msat)
Coords.Tess3R.m <- as.matrix(gl@other$latlong[match(msat@ind.names, gl@ind.names), c(2, 1)])
Pops.Tess3R.m <- as.character(gl@pop[match(msat@ind.names, gl@ind.names)])


## RUN TESS3R ##
Tess.obj.m <- tess3(X = Geno.Tess3R.m, coord = Coords.Tess3R.m, K = 1:Kmax, rep = 100, method = "projected.ls", ploidy = 2, openMP.core.num = 4, max.iteration = 1000, keep = "best", mask = 0.1) 

# Save Tess object
save(Tess.obj.m, file = paste0(Out.Path, "Tess.obj.FullDat.best100reps.K1_", Kmax, ".", sp[i], "_Microsats.R"))


}



################################################
#     CREATE CROSS VALIDATION PLOTS            #
################################################

# Tess objects
sp.Tess <- list.files(Out.Path, full.names = TRUE, pattern = ".R")

# Create new species list for naming
sp.CV <- gsub(".R", "", gsub(paste0(Out.Path, "/Tess.obj.FullDat.best100reps.K1_7."), "", sp.Tess))
sp <- gsub("_Microsats", "", gsub("_SNPs", "", sp.CV))

Tess.CrossVal.combined <- data.frame(PlotOrder = integer(), K = character(),  MeanRMSE = numeric(), Min = numeric(), Max = numeric(), Species = character(), Marker = character())

##########################################
# Create a loop that uploads tess object #
# for each species, then outputs K plots #
##########################################

for (i in 1:length(sp.CV)) {
  
  print(paste0("Creating cross validation data for ", sp.CV[i]))
  
  load(sp.Tess[i])
  
  # Rename microsat tess objects to match
  if (gsub(paste0(sp[i], "_"), "", sp.CV[i]) == "Microsats") {
    Tess.obj <- Tess.obj.m
  }
  
  # Get cross validation score in format for plotting:
  
  Tess.CrossVal <- data.frame(Iteration = 1:100)
  TessCrossvalMeans <- data.frame(PlotOrder = 1:Kmax, K = paste0("K", 1:Kmax), MeanRMSE = rep(NA, Kmax), Min = rep(NA, Kmax), Max = rep(NA, Kmax))
  Kval <- 0
  
  for(Kval in 0:(length(Tess.obj)-1)){
    K <- Kval + 1
    CrossScore <- Tess.obj[[K]]$rmse
    Tess.CrossVal$temp <- CrossScore
    Mean <- 
      TessCrossvalMeans[K, "MeanRMSE"] <- mean(Tess.CrossVal$temp)
    TessCrossvalMeans[K, "Min"] <- min(Tess.CrossVal$temp)
    TessCrossvalMeans[K, "Max"] <- max(Tess.CrossVal$temp)
  }
  
  TessCrossvalMeans$Species <- sp[i]
  TessCrossvalMeans$Marker <- gsub(paste0(sp[i], "_"), "", sp.CV[i])
  
  Tess.CrossVal.combined <- rbind(Tess.CrossVal.combined, TessCrossvalMeans)
  
}

# Plot cross-validation
  
CrossValPlot.Tess <- ggplot(Tess.CrossVal.combined, aes(x=PlotOrder, y=Min)) + 
    geom_point(size = 0.8) +
    geom_line(size = 0.2) +
    geom_errorbar(aes(x = PlotOrder, ymin=Min, ymax=Max), width = 0.1, size = 0.2) +
    labs(x = "K", y = "Mean RMSE") +
    facet_wrap(Marker ~ Species, scales = "free") +
    scale_x_continuous(breaks = seq(1,Kmax,1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", 
                                      fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 5),
          axis.title.y = element_text(size = 8, margin = margin(0,10,0,0), face ="bold"),
          axis.title.x = element_text(size = 8, margin = margin(10,0,0,0), face ="bold"),
          axis.ticks = element_line(size = 0.2),
          legend.position = "none")
  
  ggsave(filename = "CrossVal.Tess.Multi.jpg",
         plot = CrossValPlot.Tess,
         device = "jpeg",
         path = Out.Path,
         scale = 1,
         width = 9.6,
         height = 5,
         units = "in",
         dpi = 300)





#################################
#      TESS BAR/PIE CHARTS      #
#################################


################
# PREPARE MAPS #
################

# Load in polygon of IBRA subregions
StudyMap.Pilb <- readOGR("Rasters_Shapefiles/PilbaraIBRA.shp")
StudyMap.WA <- readOGR("Rasters_Shapefiles/IBRA_regions_WA.shp")

# Set map limits
xlims <- c(113.7, 122.5)
ylims <- c(-24, -19.8)

# WA
StudyMap.WA.crop <- crop(StudyMap.WA, c(xlims, ylims))
StudyMap.WA.crop@data$id <- rownames(StudyMap.WA.crop@data)
StudyMap.WA.points = fortify(StudyMap.WA.crop, region = "id")
StudyMap.WA.df <- left_join(StudyMap.WA.points, StudyMap.WA.crop@data, by="id")
StudyMap.WA.df <- StudyMap.WA.df[, c(7, 9, 11, 1:2)]

# Pilbara
StudyMap.Pilb@data$id <- rownames(StudyMap.Pilb@data)
StudyMap.Pilb.points = fortify(StudyMap.Pilb, region = "id")
StudyMap.Pilb.df <- left_join(StudyMap.Pilb.points, StudyMap.Pilb@data, by="id")
StudyMap.Pilb.df <- StudyMap.Pilb.df[, c(7, 14, 12, 1:2)]

# Rename regions/sub regions so that they are ordered for plotting (so I can assign the correct colours)
StudyMap.WA.df$REG_NAME <- as.character(StudyMap.WA.df$REG_NAME)
StudyMap.WA.df <- StudyMap.WA.df[order(StudyMap.WA.df$REG_NAME), ]
Reg.name <- unique(StudyMap.WA.df$REG_NAME)

PlotOrder <- vector()
Reg <- 0

for (Reg in 0:(length(Reg.name)-1)) {
  RegNo <- Reg + 1
  PlotOrder <- c(PlotOrder, rep((RegNo), length(StudyMap.WA.df$REG_NAME[StudyMap.WA.df$REG_NAME == Reg.name[RegNo]])))
}
StudyMap.WA.df$PlotReg <- as.factor(paste0(0, PlotOrder, "_", StudyMap.WA.df$REG_NAME))
unique(StudyMap.WA.df$PlotReg)

StudyMap.Pilb.df$SUB_NAME_7 <- as.character(StudyMap.Pilb.df$SUB_NAME_7)
StudyMap.Pilb.df <- StudyMap.Pilb.df[order(StudyMap.Pilb.df$SUB_NAME_7), ]
Sub.name <- unique(StudyMap.Pilb.df$SUB_NAME)

PlotOrder_Sub <- vector()
SubReg <- 0

for (SubReg in 0:(length(Sub.name)-1)) {
  SubRegNo <- SubReg + 1
  PlotOrder_Sub <- c(PlotOrder_Sub, rep((SubRegNo), length(StudyMap.Pilb.df$SUB_NAME[StudyMap.Pilb.df$SUB_NAME == Sub.name[SubRegNo]])))
}
StudyMap.Pilb.df$PlotSubReg <- as.factor(paste((PlotOrder_Sub + (length(unique(StudyMap.Pilb.df$SUB_NAME))) + 2), StudyMap.Pilb.df$SUB_NAME, sep = "_"))
unique(StudyMap.Pilb.df$PlotSubReg)

# Create palette for background map
MapPal <- c(gray.colors(10)[seq(2,10,2)], gray.colors(10)[seq(1,10,2)])

# Plot
Map.plot <- ggplot() +
  geom_polygon(data = StudyMap.WA.df, 
               aes(x = long, y = lat, group = group, fill = PlotReg)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  geom_polygon(data = StudyMap.Pilb.df, 
               aes(x = long, y = lat, group = group, fill = PlotSubReg)) + 
  scale_fill_manual(values = c(MapPal)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 9, 
                                    margin = margin(0,10,0,0), 
                                    face ="bold"),
        axis.title.x = element_text(size = 9, 
                                    margin = margin(10,0,0,0), 
                                    face ="bold"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")


load(paste0(gl.Path, "/RelClean.gl.Nt.rdata"))
gl <- sp.gl
load(paste0(Out.Path, "/Tess.obj.FullDat.best100reps.K1_7.Nt_SNPs.R"))
  
# Retrieve tess3 Q matrix rep with the lowest cross-entropy criterion
Q.matrix.Tess <- qmatrix(Tess.obj, K = 2)
  
# Create palette (rearrange for different spp so clusters are coloured consistently by location)
Pal <- c("#594DA4", "#46C697")

# Reformat data for plotting
Q.matrix.Tess.df <- cbind(gl@ind.names, as.data.frame(Q.matrix.Tess[,]))
colnames(Q.matrix.Tess.df)[1] <- c("sampleID")
  
# Sort by cluster
  
# Assign individuals to clusters
Q.matrix.Tess.df$GenCluster <- colnames(Q.matrix.Tess.df[,2:3])[apply(Q.matrix.Tess.df[, 2:3],1,which.max)]
Qmat.K.sorted <-  Q.matrix.Tess.df[order(Q.matrix.Tess.df$GenCluster, Q.matrix.Tess.df$V2), -(which(colnames(Q.matrix.Tess.df) == "GenCluster"))]

# Add new dummy IDs that are in the correct order for plotting
Qmat.K.sorted$sampleID.ordered <- c(paste0("Dummy_00", 1:9), paste0("Dummy_0", 10:99), paste0("Dummy_", 100:length(gl@ind.names)))

Qmat.K.sorted.gather <- tidyr::gather(Qmat.K.sorted, "popGroup", "prob", -sampleID, -sampleID.ordered)
  
  
# Plot
TessBarPlot <- ggplot(Qmat.K.sorted.gather, aes(factor(sampleID.ordered), prob, fill = factor(popGroup))) +
    scale_fill_manual(values = Pal) +
    geom_col(color = "transparent", size = 0.1, width = 1) +
    labs(y = " ") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", 
                                      fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          strip.text = element_text(face = "plain", size = 6, angle = 90),
          strip.background = element_rect(linetype = "solid", size = 0.1, 
                                          colour = "black", fill = "lightgrey"),
          axis.text = element_blank(),
          axis.title.y = element_text(size = 9, margin = margin(0,24,0,0), face ="bold"),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
  
  
# Assign individuals to clusters
TessAssign <- Qmat.K.sorted[ ,-(which(colnames(Qmat.K.sorted) == "sampleID.ordered"))]
colnames(TessAssign)[2:3] <- paste0("Clust", 1:2)
Tess.Admixed <- apply(TessAssign[, 2:3],1,max) < 0.7
TessAssign$Admixed <- ifelse(Tess.Admixed == TRUE, "Hybrid", "Assigned")
TessAssign$GenCluster <- colnames(TessAssign[,2:3])[apply(TessAssign[, 2:3],1,which.max)]

  
######################
# Tess PIE CHART MAP #
######################
  
# Define xy position for north arrow in plot
xmin.n <- xlims[1] + 0.1
ymin.n <- ylims[1]
xmax.n <- xlims[2]
ymax.n <- ylims[2] - 0.1

# Define xy position for scale bar in plot
xmin.sb <- xlims[1] + 0.25
ymin.sb <- ylims[1] + 0.25
xmax.sb <- xlims[2]
ymax.sb <- ylims[2]

Fig.h <- 4.74
Fig.w <- 7.2
  
# Get mean of individual coordinates so we have a single set of coordinates for each pop
q.matrix.df <- cbind(Qmat.K.sorted, gl@other$ind.metrics[match(Qmat.K.sorted$sampleID, gl@other$ind.metrics$id), c("lon", "lat", "Buffer_5km")])
coord.pop.mean <- q.matrix.df %>% group_by(Buffer_5km) %>% summarise(Site.lon = mean(lon), Site.lat = mean(lat))
  
# Calculate population average ancestry proportions and create an array with population coordinates:
q.matrix.df.pop <- q.matrix.df %>%
  group_by(Buffer_5km) %>%
  summarise(Mean_AP1 = mean(V1), Mean_AP2 = mean(V2))    
  
# Reformat coordinates and q matrix for plotting
coord.pie <- left_join(coord.pop.mean, q.matrix.df.pop, by = "Buffer_5km")
coord.pie$Clust <- ifelse(coord.pie$Mean_AP1 >= 1/K[i], "Clust1", "Clust2")
coord.pie <- coord.pie[order(coord.pie$Clust), ]
coord.pie$Group <- c(paste0("Group.0", 1:9), paste0("Group.", 10:nrow(coord.pie)))

# Tess pie charts of ancestry proportions on map
Tess.Pie.Map <- Map.plot +
    scale_fill_manual(values = c(MapPal, Pal)) +
    scalebar(location = "bottomleft", dist_unit = "km", 
             transform = TRUE, model = 'WGS84',
             x.min = xmin.sb, x.max = xmax.sb,
             y.min = ymin.sb, y.max = ymax.sb,
             dist = 100, st.dist = 0.025, st.size = 2.5,
             height = 0.02, border.size = 0.15) +
    north(location = "topleft", scale = 0.06, symbol = 12,
          x.min = xmin.n, x.max = xmax.n, 
          y.min = ymin.n, y.max = ymax.n) +
    geom_scatterpie(data = coord.pie, aes(x=Site.lon, y=Site.lat, group = Group), 
                    cols = colnames(coord.pie[ , 4:5]),  
                    pie_scale = 1.65, size = 0.1)
  
Tess.Combined <- ggarrange(TessBarPlot, Tess.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))
  
# Check that figures are okay for colourblind
cvd_grid(Tess.Combined)

# Save as jpeg
ggsave(filename = "TessPieBar_K2_Nt.jpg",
         plot = Tess.Combined,
         device = "jpeg",
         path = Out.Path,
         scale = 1,
         width = Fig.w,
         height = Fig.h,
         units = "in",
         dpi = 300)
  
