#############################################
#       Isolation-By-Resistance (IBR)       #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/LandGen_AridMammals/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

library(dplyr)
library(data.table)
library(raster)
library(splitstackshape)
library(ecodist)
library(reshape2)
library(stringr)
library(ResistanceGA)
library(parallel)
library(doParallel)
library(rgdal)
library(pointdexter)
library(viridis)
library(ggplot2)
library(ggnewscale)
library(ggsn)

####################
# SET UP VARIABLES #
####################

# Set path to a directory for saving data/outputs
Out.path <- "IBR/IBR_outputs/"


##################
# Import samples #
##################

Samps <- read.csv("IBR/Data/SampleMetaData_Cleaned.csv", stringsAsFactors = FALSE)

# Remove island samples
Samps <- Samps[!(Samps$pop %like% "Island" | Samps$pop %like% "island"), ]

# Make separate dfs by species
Samps.Nt <- Samps[Samps$species == "Ningaui timealeyi", ]
Samps.Pc <- Samps[Samps$species == "Pseudomys chapmani", ]
Samps.Ph <- Samps[Samps$species == "Pseudomys hermannsburgensis", ]

# Remove any outlier samples and any NAs
Samps.Nt <- Samps.Nt[!(is.na(Samps.Nt$year)), ]
Samps.Pc <- Samps.Pc[!(is.na(Samps.Pc$year)), ]
Samps.Ph <- Samps.Ph[!(is.na(Samps.Ph$year)), ]


####################################
# Get pixel coordinate from raster #
####################################

# Read in one of the processed rasters so that I can aggregate the genetic data based on the pixel size of the rasters 
Raster <- raster("Rasters_Shapefiles/SOMO29.asc")

sp <- c("Nt", "Pc", "Ph")

for (i in 1:length(sp)) {
  
  Samps.sp <- get(paste0("Samps.", sp[i]))
  
  # Get xy coordinates of raster pixels that overlap with samples
  cell <- data.frame(id = Samps.sp[, "id"], 
                     cell = raster::extract(Raster, Samps.sp[, c("UTMX", "UTMY")], cellnumbers=T)[, 1], 
                     Rast = raster::extract(Raster, Samps.sp[, c("UTMX", "UTMY")], cellnumbers=T)[, 2])
  
  # Find any samples outside of raster extent (to remove later)
  rmNAs <- cell$id[is.na(cell$Rast)]
  print(paste0(sp[i], " points outside raster: ", length(rmNAs)))
  xy <- cbind(cell, xyFromCell(Raster, cell[,2]))
  
  # How many unique points
  print(paste0(sp[i], " unique points: ", length(unique(xy$cell))))
  
  # Create a new column called "pixel" to identify duplicated coords easily (i.e. all coordinates that fall within the same pixel)
  # Create a new data frame with only the unique coordinates
  if(length(unique(xy$cell)) >= 100) {
    pixel.group <- data.frame(cell = unique(xy$cell), 
                              pixel = c(paste0("P00", 1:9), 
                                        paste0("P0", 10:99), 
                                        paste0("P", 100:length(unique(xy$cell)))))
  } else if(length(unique(xy$cell)) < 100) {
    pixel.group <- data.frame(cell = unique(xy$cell), 
                              pixel = c(paste0("P00", 1:9), 
                                        paste0("P0", 10:length(unique(xy$cell)))))
  } else {
    print("More than 199 samples")
    stop()
  }
  
  
  # Match these to the main data set
  Samps.sp <- cbind(Samps.sp, pixel.group[match(xy$cell, pixel.group$cell), "pixel"])
  colnames(Samps.sp)[ncol(Samps.sp)] <- "pixel.grp"
  
  # Remove samples outside of raster extent
  Samps.sp <- Samps.sp[!Samps.sp$id %in% rmNAs, ]
  
  # Save to workspace
  assign(paste0("Samps.", sp[i]), Samps.sp)
  
  rm(Samps.sp, cell, rmNAs, pixel.group, xy)
  
}


###################################################
# READ IN FILTERED GENETC DATA - GENLIGHT FORMAT  #
###################################################

# Load in genlight saved at the end of SNP filtering script
gl.list <- list.files("SNP_Filtering/SampleClean_outputs", pattern = "gl.", full.names = TRUE)

for (i in 1:length(gl.list)) {
  load(gl.list[i])
  assign(paste0(sp[i], ".gl.rmRel"), sp.gl)
}

#######################################
# GENERATE GENETIC DISTANCE FOR SNPS  #
#######################################

for (i in 1:length(sp)) {
  
  # Subset genlight by samp df
  gl <- get(paste0(sp[i], ".gl.rmRel"))
  Samps.sp <- get(paste0("Samps.", sp[i]))
  gl <- gl[gl@ind.names %in% Samps.sp$id, ]
  
  # Make sure ind metrics match ids (0 if match)
  if (sum(!(gl@ind.names == gl@other$ind.metrics$id)) > 0) {
    print("genlight mismatch error")
    stop()
  }
  
  # Rename Ids as the Coordinate group so it's easier to get mean
  indNames(gl) <- Samps.sp$pixel.grp[match(gl@ind.names, Samps.sp$id)]
  
  # Some of the genetic distance methods can't deal with missing data
  # I'll interpolate missing data by using the mean
  gl.mat <- round(tab(gl, NA.method="mean"), 0)
  
  ##############################
  # Calculate genetic metrics  #
  ##############################
  
  # According to Shirk et al. 2017
  # https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12684?casa_token=kDgBnlvK8aMAAAAA:cAGe_TzoPPF-u4Qf0p-3WXiB0jhTLo0vUcOV6suW05lCUJifsq8w-nvWECqKeNuigFtiPssskhxAxcuI
  # The best metrics for landscape genetic analysis are PC eigenvectors (as long as more than one axis is included), Euclidean distance, Bray-Curtis distance and Proportion of shared alleles.
  # The PCA based metrics performed best when dispersal was high and sample sizes were low.
  # For this reason I'll use these metrics in my analysis (excluding proportion shared alleles). 
  # I'll try all three, because these were tested for multi-allelic markers (meaning I'm not quite sure how these metrics perform when using bi-allelic SNPs)
  
  
  #### EUCLIDEAN GENETIC DISTANCE #### 
  
  # Justification for taking the average:
  # Gallego-Garc?a et al. 2019 (Mol Ecol)
  # Calculated pairwise Euclidean genetic distances, as this metric does not assume random mating. Because many sample sites were represented by multiple individuals, whereas others contained few or just one, we estimated the genetic distance across all sites as the mean genetic distance between individuals among sites. For example, the genetic distance between a hypothetical site A with three individuals (1, 2, and 3) and site B with one individual (1) was estimated as the mean between pairs A1-B1, A2-B1, and A3-B1.
  
  # Calculate Euclidean genetic distance between samples
  EucGenDist <- distance(gl.mat, method = "euclidean")
  
  # Convert lower triangular to table
  # Have already tested to make sure Euc dist order is the same as ResistanceGA lower convience function
  # I just need to make sure I sort by Individual 1, then individual 2 in the code below
  EucGenDist.mat <- as.matrix(EucGenDist)
  EucGenDist.tab <- melt(EucGenDist.mat)[melt(lower.tri(EucGenDist.mat))$value, ]
  names(EucGenDist.tab) <- c("Ind2", "Ind1", "EucDist")
  
  # Remove comparisons between same pixel group
  EucGenDist.tab <- EucGenDist.tab[!EucGenDist.tab$Ind1 == EucGenDist.tab$Ind2,]
  
  # Create a new column to group by (and take the mean of)
  Euc.pixdf <- data.frame(Ind2 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind2)), Ind1 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind1)))
  
  # Make a new column with both pixel groups so I can get the mean across duplicates
  EucGenDist.tab$Pair <- paste(ifelse(apply(Euc.pixdf, 1, which.min) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), ifelse(apply(Euc.pixdf, 1, which.max) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), sep = "_")
  
  # Group by pair and then take the mean
  EucGenDist.Av <- EucGenDist.tab %>%
    group_by(Pair) %>%
    summarise(EucDistMean = mean(EucDist))
  
  # Check correct number of values
  if (!(length(unique(Samps.sp$pixel.grp))*length(unique(Samps.sp$pixel.grp))-length(unique(Samps.sp$pixel.grp)))/2 == nrow(EucGenDist.Av)) {
    print("Distance matrix conversion error")
    stop()
  }
  
  # Create final vector to run in ResistanceGA
  EucDist.Vect <- EucGenDist.Av$EucDistMean
  
  
  ########################################
  # WRITE DATA FOR INPUT TO RESISTANCEGA #
  ########################################
  
  # Write genetic data
  write.csv(data.frame(EucDist = EucDist.Vect), 
            paste0("IBR/Data/EucDist_SNPs", sp[i], ".csv"), 
            row.names = FALSE)
  
  # Write coordinates
  Coords <- Samps.sp[!duplicated(Samps.sp$pixel.grp), c("UTMX", "UTMY", "pixel.grp")]
  Coords <- Coords[order(Coords$pixel.grp), c("UTMX", "UTMY")]
  
  # Check numbers match up
  if (!(nrow(data.frame(EucDist = EucDist.Vect)) == ((nrow(Coords)*nrow(Coords))-nrow(Coords))/2)) {
    print("Genetic metric df and coordinate df don't match")
    stop()
  }
  
  # Save coords
  write.csv(Coords, paste0("IBR/Data/Coords_", sp[i], ".csv"), row.names = FALSE)
  
}

########################################
# GENERATE GENETIC DISTANCE FOR MSATs  #
########################################

str.list <- list.files("IBB/Data", pattern = ".str", full.names = TRUE)

# Read in genetic data (MSATS)
# Have already converted to structure format in genalex
gi.Nt <- read.structure(file = str.list[1], n.ind = 180, n.loc = 12, onerowperind = TRUE, col.lab = 1, row.marknames = 1, col.pop = 0, col.others = 0)

gi.Pc <- read.structure(file = str.list[2], n.ind = 80, n.loc = 14, onerowperind = TRUE, col.lab = 1, row.marknames = 1, col.pop = 0, col.others = 0)

gi.Ph <- read.structure(file = str.list[3], n.ind = 173, n.loc = 14, onerowperind = TRUE, col.lab = 1, row.marknames = 1, col.pop = 0, col.others = 0)


for (i in 1:length(sp)) {
  
  # Subset genlight by samp df
  gi <- get(paste0("gi.", sp[i]))
  Samps.sp <- get(paste0("Samps.", sp[i]))
  gi <- gi[indNames(gi) %in% Samps.sp$id, ]
  gi@other$ind.metrics <- Samps.sp[match(indNames(gi), Samps.sp$id), ]
  
  # Rename Ids as the pixel group so it's easier to get mean
  indNames(gi) <- Samps.sp$pixel.grp[match(indNames(gi), Samps.sp$id)]
  
  # Some of the genetic distance methods can't deal with missing data
  # I'll interpolate missing data by using the mean
  # Keep all the naming the same (even though gi, not gl)
  gl.mat <- round(tab(gi, NA.method="mean"), 0)
  
  ##############################
  # Calculate genetic metrics  #
  ##############################
  
  # According to Shirk et al. 2017
  # https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12684?casa_token=kDgBnlvK8aMAAAAA:cAGe_TzoPPF-u4Qf0p-3WXiB0jhTLo0vUcOV6suW05lCUJifsq8w-nvWECqKeNuigFtiPssskhxAxcuI
  # The best metrics for landscape genetic analysis are PC eigenvectors (as long as more than one axis is included), Euclidean distance, Bray-Curtis distance and Proportion of shared alleles.
  # The PCA based metrics performed best when dispersal was high and sample sizes were low.
  # For this reason I'll use these metrics in my analysis (excluding proportion shared alleles). 
  # I'll try all three, because these were tested for multi-allelic markers (meaning I'm not quite sure how these metrics perform when using bi-allelic SNPs)
  
  
  #### EUCLIDEAN GENETIC DISTANCE #### 
  
  # Justification for taking the average:
  # Gallego-Garc?a et al. 2019 (Mol Ecol)
  # Calculated pairwise Euclidean genetic distances, as this metric does not assume random mating. Because many sample sites were represented by multiple individuals, whereas others contained few or just one, we estimated the genetic distance across all sites as the mean genetic distance between individuals among sites. For example, the genetic distance between a hypothetical site A with three individuals (1, 2, and 3) and site B with one individual (1) was estimated as the mean between pairs A1-B1, A2-B1, and A3-B1.
  
  # Calculate Euclidean genetic distance between samples
  EucGenDist <- distance(gl.mat, method = "euclidean")
  
  # Convert lower triangular to table
  # Have already tested to make sure Euc dist order is the same as ResistanceGA lower convience function
  # I just need to make sure I sort by Individual 1, then individual 2 in the code below
  EucGenDist.mat <- as.matrix(EucGenDist)
  EucGenDist.tab <- melt(EucGenDist.mat)[melt(lower.tri(EucGenDist.mat))$value, ]
  names(EucGenDist.tab) <- c("Ind2", "Ind1", "EucDist")
  
  # Remove comparisons between same pixel group
  EucGenDist.tab <- EucGenDist.tab[!EucGenDist.tab$Ind1 == EucGenDist.tab$Ind2,]
  
  # Create a new column to group by (and take the mean of)
  Euc.pixdf <- data.frame(Ind2 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind2)), Ind1 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind1)))
  
  # Make a new column with both pixel groups so I can get the mean across duplicates
  EucGenDist.tab$Pair <- paste(ifelse(apply(Euc.pixdf, 1, which.min) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), ifelse(apply(Euc.pixdf, 1, which.max) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), sep = "_")
  
  # Group by pair and then take the mean
  EucGenDist.Av <- EucGenDist.tab %>%
    group_by(Pair) %>%
    summarise(EucDistMean = mean(EucDist))
  
  # Check correct number of values
  if (!(length(unique(Samps.sp$pixel.grp))*length(unique(Samps.sp$pixel.grp))-length(unique(Samps.sp$pixel.grp)))/2 == nrow(EucGenDist.Av)) {
    print("Distance matrix conversion error")
    stop()
  }
  
  # Create final vector to run in ResistanceGA
  EucDist.Vect <- EucGenDist.Av$EucDistMean
  
  
  ########################################
  # WRITE DATA FOR INPUT TO RESISTANCEGA #
  ########################################
  
  # Write genetic data
  write.csv(data.frame(EucDist = EucDist.Vect), 
            paste0("IBR/Data/EucDist_MSAT", sp[i], ".csv"), 
            row.names = FALSE)
  
}



####################
# Run ResistanceGA #
####################

# Set data directory
data.dir <- "IBR/IBR_outputs/"

# Set path/filename to Coordinate file
Coords <- "IBR/Data/Coords_Nt.csv"

# Set path/filename to genetic file 
GenMet.df.name <- "IBR/Data/EucDist_MSATNt.csv"

# Raster folder
Rasts <- "Rasters_Shapefiles/"

# Define Genetic Metric (i.e. column name)
GenMet <- "EucDist"

# Set seed for single surface
seed <- 1234

# Replicate (should replicate 2x)
rep <- 1

# Set number of cores to use when running parallel (i.e. number of VCPUs available to instance)
Para <- 16

###############
# Import data #
###############

# Load in sample coordinates
Samps <- as.matrix(read.csv(Coords, stringsAsFactors = FALSE))

# Load in genetic distance (in vector format)
# Already ordered so matches coords
GenMet.df <- read.csv(GenMet.df.name)
Gen_vect <- GenMet.df[, GenMet]


#####################################################
# Run Single Surface Optimisation with ResistanceGA #
#####################################################

# I'm doing an initial single surface optimisation of all rasters to find out which raster performs best when layers are correlated (I'll drop the others for multisurface optimisation so that I'm only including uncorrelated variables)
# Create raster stack (only use those with a 5km resolution))
Rstack <- stack(list.files(Rasts, pattern = ".asc", full.names = TRUE))
Rstack.names <- names(Rstack)

# NOTE: PRELIMINARY! RESULTS FROM SINGLE SURFACE OPTIM NOT INCLUDED IN REPOSITORY

# Run a loop to analyse each surface separately
for (i in 1:length(Rstack.names)) {
  
  # Create/Set results directory:
  Results.dir <- paste0("IBR/IBR_outputs/SS_Optim/", Rstack.names[i])
  ifelse(!dir.exists(file.path(Results.dir)), dir.create(Results.dir, recursive = TRUE), FALSE)
  
  # Run single surfaces
  GA.inputs <- GA.prep(method = "LL",
                       ASCII.dir = Rstack[[i]],
                       Results.dir = paste0(Results.dir, "/"),
                       seed = seed + rep,
                       parallel = Para)
  
  gdist.inputs <- gdist.prep(n.Pops = nrow(Samps),
                             samples = Samps,
                             response = Gen_vect,
                             method = 'commuteDistance')
  
  # First run all single surfaces, Multi-surface is response variable
  SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs, 
                         GA.inputs = GA.inputs, 
                         dist_mod = TRUE, 
                         null_mod = TRUE)
  
}


####################
# Run all combined #
####################

# After single surface optimisation, I have decided which (uncorrelated) rasters to include in the final analysis
# I will run the all_comb function for multi-surface optimisation


# Directory containing final rasters
Rasts.final <- "" # Need to move final raster set to a new folder


GA.inputs <- GA.prep(method = "LL",
                     ASCII.dir = Rasts.final,
                     Results.dir = "all.comb",
                     seed = 1234 + rep,
                     parallel = Para)

gdist.inputs <- gdist.prep(n.Pops = nrow(Samps),
                           samples = Samps,
                           response = Gen_vect,
                           method = 'commuteDistance')

AC_RESULTS <- all_comb(gdist.inputs = gdist.inputs,
                       GA.inputs = GA.inputs, 
                       results.dir = "all.comb",
                       max.combination = 4,
                       iters = 1000,
                       replicate = 1,
                       sample.prop = 0.75,
                       dist_mod = TRUE,
                       null_mod = TRUE)

