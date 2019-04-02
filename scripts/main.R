# Execute this script to generate all results and manuscript tables and figures.
# Warning: Long execution time depending on computing resources used, may take up to 2-3 days.
# Warning: Results will occupy on the order of 12 Gb.

# Set root dir to the directory containing the data and scripts folders
root.dir <- "<YOUR PATH>" # path to the directory where data and scripts are located
setwd(root.dir)
dir.create(paste0(root.dir, "/manuscript"))
dir.create(paste0(root.dir, "/outputs"))

source("scripts/init.R")
source("scripts/functions.R")
source("scripts/maxent.R")
source("scripts/hbm.R")
source("scripts/manuscript.R")

# Run Maxent models
executeMaxentModels(force.run)

# Run HBM models
executeHBMModels(force.run)

# Check for residual spatial autocorrelation
assessModelsRSAC(force.run)

# Write a csv with all results combined
getCombinedMeanResultsAll(0.5, "best-maxent-model", force.run)

# Generate manuscript tables and figures
manuscriptTables()
manuscriptFigures()
suppInformationTables()
suppInformationFigures()
