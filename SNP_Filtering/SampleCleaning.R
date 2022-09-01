#############################################
#                 STEP 2                    #  
#        Sample Cleaning: Relatives         #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/LandGen_AridMammals/")

############
# Packages #
############

library(stringr)
library(sp)
library(raster)
library(related)
library(dartR)
library(dplyr)
library(sf)
library(rgdal)
library(RColorBrewer)
library(leaflet)

###########################
# READ IN DATA - ADD UTMS #
###########################

# List species codes
sp <- c("Nt", "Pc", "Ph")

# Define output folder
Out.path <- "SNP_Filtering/SampleClean_outputs/"

# Read in meta data
Metadata <- read.csv("SNP_Filtering/Data/SampleMetaData.csv", stringsAsFactors = FALSE)

# list genlights
gl.paths <- list.files("SNP_Filtering/Filtering_outputs/", pattern = "rdata", full.names = TRUE)

# Run loop that:
# loads filtered genlight
# subsets individuals in meta-data to match those retained in genlight after filtering
# converts lon/lat to UTMs

for (i in 1:length(sp)) {
  
  load(gl.paths[i])
  IndMD <- Metadata[Metadata$id %in% gl.FinalFilt@ind.names, ]
  
  Spatial.P <- SpatialPoints(IndMD[ , c("lon", "lat")])
  crs(Spatial.P) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  Spatial.P.UTM <- spTransform(Spatial.P, CRS("+proj=utm +zone=50 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
  colnames(Spatial.P.UTM@coords) <- c("UTMX", "UTMY")
  IndMD <- cbind(IndMD, Spatial.P.UTM@coords)
  
  assign(x = paste0(sp[i], ".Metadata"), value = IndMD)
  assign(x = paste0(sp[i], ".gl"), value = gl.FinalFilt)
  
  rm(gl.FinalFilt, IndMD, Spatial.P, Spatial.P.UTM)
  
}

#######################################
# CREATE BUFFERS FOR GROUPING SAMPLES #
#######################################

# Create a 5km buffer to group samples
Buff_m <- 5000
names(Buff_m) <- "Buffer_5km"

# Run a loop that groups samples withing the specified buffer distance and 
# saves unique buffer ID as a column in meta data
for (i in 1:length(sp)) {
  IndMD <- get(paste0(sp[i], ".Metadata"))
  
  # Add column called CoordLocation 
  # That groups samples at exact same coordinates
  IndMD <- left_join(IndMD, 
                     cbind(unique(IndMD[, c("UTMX", "UTMY")]), 
                           CoordLocation = 1:nrow(unique(IndMD[, c("UTMX", "UTMY")]))), 
                     by = c("UTMX", "UTMY"))
  
  
  # Create spatial points, generate buffers of different distances in a loop
  Coords.utm <- SpatialPoints(IndMD[, c("UTMX", "UTMY")], proj4string = crs("+proj=utm +zone=50 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
  
  
  print(paste0("Generating ", names(Buff_m), " for ", sp[i]))
    
  # Create a buffer around samples
  # Convert to sf so that I can split polygon layer so each buffer has a unique ID
  UTM.Buff <- st_cast(as(buffer(x = Coords.utm, Buff_m, dissolve = TRUE), 
                           "sf"),"POLYGON")
    
  # Name with unique ID
  UTM.Buff$BuffID <- 1:dim(UTM.Buff)[1]
    
  # Convert back to sp object
  # Get buffer polygon info at each sample point
  IndMD <- cbind(IndMD, over(Coords.utm, as(UTM.Buff, "Spatial")))
  colnames(IndMD)[which(colnames(IndMD) == "BuffID")] <- names(Buff_m)
  
  assign(x = paste0(sp[i], ".Metadata"), value = IndMD)
  rm(IndMD)
  
}


#############################
# ADD METADATA TO GENLIGHTS #
#############################

for (i in 1:length(sp)) {
  
  # Assign metadata and genlight to variable
  IndMD <- get(paste0(sp[i], ".Metadata"))
  sp.gl <- get(paste0(sp[i], ".gl"))
  
  # Replace ind.metrics slot in genlight with updated metadat 
  # making sure sample order matches
  sp.gl@other$ind.metrics <- IndMD[match(sp.gl@ind.names, IndMD$id), ]
  
  # Double check sample order matches
  if (sum(!sp.gl@other$ind.metrics$id == sp.gl@ind.names) == 0) {
    print(paste("Datasets match: ", sp[i]))
  } else {
    print(paste("Error: dataset mis-match for ", sp[i]))
    stop()
  }
  
  # Save to workspace and output cleaned/updated genlight
  assign(x = paste0(sp[i], ".gl"), value = sp.gl)
  
  # Remove loop variables
  rm(IndMD, sp.gl)
}


################################
# CREATE "RELATED" INPUT FILES #
################################

# List to run in loop
gl.list <- list(Nt.gl,
                Pc.gl,
                Ph.gl)

names(gl.list) <- c("Nt",
                    "Pc",
                    "Ph")

# Write out relatedness input files
for (i in 1:length(gl.list)) {
  
  demerel <- gl2demerelate(gl.list[[i]])
  demerel <- demerel[, 3:ncol(demerel)]
  colnames(demerel) <- NULL
  
  # Export as a txt file
  write.table(demerel, na = "0", file = paste0(Out.path, names(gl.list)[[i]], ".related.txt"), sep = "\t", row.names = TRUE)
  
}


################################
# CALCULATE WANG'S RELATEDNESS #
################################

# List files
rel.files <- list.files("SNP_Filtering/SampleClean_outputs",pattern = ".txt", full.names = TRUE)
names(rel.files) <- str_remove(str_remove(rel.files, "SNP_Filtering/SampleClean_outputs/"), ".related.txt")

# Run through a loop that calculates relatedness and outputs pairs of related individuals 
# where Wang's relatedness estimate equals or exceeds 0.24 (i.e. ~half-sibs and above)
for (i in 1:length(rel.files)) {
  
  # Read in text file
  rel.tab <- readgenotypedata(rel.files[i])
  
  # Calculate relatedness (island pop - so allow inbreeding)
  rel.dat <- coancestry(rel.tab$gdata, lynchli = 0, lynchrd = 0, quellergt = 0, ritland = 0, wang = 1)
  
  gl <- get(paste0(names(rel.files)[i], ".gl"))
  
  rel.dat <- rel.dat$relatedness[rel.dat$relatedness$wang >= 0.24, c(2:3, 6)]
  rel.dat <- cbind(rel.dat,
                   gl@other$ind.metrics[match(rel.dat$ind1.id, gl@ind.names), c("sex", "pop", "year")],
                   gl@other$ind.metrics[match(rel.dat$ind2.id, gl@ind.names), c("sex", "pop", "year")])
  
  colnames(rel.dat)[4:9] <- c("Ind.1.sex", 
                              "Ind.1.pop", 
                              "Ind.1.year",
                              "Ind.2.sex", 
                              "Ind.2.pop", 
                              "Ind.2.year" )
  
  # Output results
  if (nrow(rel.dat) > 0) {
    write.csv(rel.dat, paste0(Out.path, names(rel.files)[i], "_wang0.24.csv"), row.names = FALSE)
  }
  
}

# Find out dataset info that will help me decide which individuals to delete
for (i in 1:length(sp)) {
  
  gl <- get(paste0(sp[i], ".gl"))
  
  print(paste0("Sex/Year breakdown for ", sp[i]))
  print(table(gl@other$ind.metrics$sex, useNA = "ifany"))
  print(table(gl@other$ind.metrics$year, useNA = "ifany"))
}


#########################
# FILTER ON RELATEDNESS #
#########################

# Remove one of the pair of related individuals identified in the spread sheet
RelSS <- read.csv(paste0(Out.path, "SamplesToDelete_Related.csv"))

for (i in 1:length(unique(RelSS$Spp))) {
  Spp <- unique(RelSS$Spp)[i]
  rmID <- RelSS$ID[RelSS$Spp == Spp]
  
  gl <- get(paste0(Spp, ".gl"))
  gl <- gl[!gl@ind.names %in% rmID, ]
  assign(x = paste0(Spp, ".gl"), value = gl)
  
  rm(Spp, rmID, gl)
}


##########################
# OUTPUT FINAL GENLIGHTS #
##########################

for (i in 1:length(sp)) {
  
  # Assign genlight to variable
  sp.gl <- get(paste0(sp[i], ".gl"))
  
  # Output cleaned genlight
  save(sp.gl, file = paste0(Out.path, "RelClean.gl.", sp[i], ".rdata"))
  
  # Remove loop variables
  rm(sp.gl)
}





