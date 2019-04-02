library(INLA)
library(raster)
library(readxl)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(fields)
library(dismo)
library(spdep)
library(rJava)
library(ggplot2)

n.threads <- 4

# Used to speed up the script when parts of it are already done
force.run <- F

# FOLDERS
data.dir <- paste0(root.dir, "data/")
output.dir <- paste0(root.dir, "outputs/")
hbm.output.dir <- paste0(output.dir, "hbm/")
maxent.output.dir <- paste0(output.dir, "maxent/")
manuscript.dir <- paste(root.dir, "manuscript/", sep="")

dir.create(hbm.output.dir, recursive=T, showWarnings = F)
dir.create(maxent.output.dir, showWarnings = F)
dir.create(paste0(manuscript.dir, "figures"), recursive=T, showWarnings = F)
dir.create(paste0(manuscript.dir, "tables"), showWarnings = F)

# DATA VARIABLES
genetic.clusters <- c("gc1", "gc2", "gc3", "gc4")
variables <- c("bio1", "bio2", "bio3", "bio4", "bio8", "bio12", "bio15", "bio18")

# DATA THRESHOLD FOR MAXENT MODELS
th <- 0.5
rsac.nr.simulations <- 10000

# MANUSCRIPT
mean.gc.colour.gradient <- c("Floral White", "Wheat",
                             "Dark Olive Green 2", "Dark Olive Green 4", "Firebrick")
sd.gc.colour.gradient <- c("Snow", "Light Goldenrod 1",
                           "Gold 3", "Gold 4", "Dark Orchid 4")

mean.sp.colour.gradient <- c("Linen", "Bisque", "Pale Turquoise", "Dark Turquoise", "Dark Blue")
sd.sp.colour.gradient <- c("Snow", "Light Goldenrod 1",
                           "Gold 3", "Gold 4", "Dark Orchid 4")

