#### DATA PREPARATION ##############################################################################
# Species genetic cluster data
data <- read.csv(paste0(data.dir, "sp/ath_accessions.csv"))

# Boundaries
ip.boundaries <- getStudyAreaBoundaries()
crs.proj <- crs(ip.boundaries)

# Mesh
mesh <- createMesh(ip.boundaries, crs=crs.proj, force.run)

# Spde
spde <- getSPDE(mesh, force.run)

# Projection grid at 5 km resolution
cat("\n### Preparing the projection grid\n")
proj.grid.mat <- inla.mesh.projector(
    mesh, xlim=bbox(ip.boundaries)[1,],
    ylim=bbox(ip.boundaries)[2,],
    dims=c(round((bbox(ip.boundaries)[1,"max"]-bbox(ip.boundaries)[1,"min"])/5000,0),
           round((bbox(ip.boundaries)[2,"max"]-bbox(ip.boundaries)[2,"min"])/5000,0))
)

# Prediction points inside the study area boundaries
cat("\n### Preparing prediction points\n")
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, ip.boundaries@proj4string), ip.boundaries)
i.map <- !is.na(ov)

# Index of points with actual data in predictors, we use one of the raster predictors to check
i.bio <- which(!is.na(
    raster::extract(
        raster(paste0(data.dir, "/climate/2000/bio1.tif")),
        proj.grid.mat$lattice$loc[i.map,])
)
)

prediction.points <- data.frame(proj.grid.mat$lattice$loc[i.map, ][i.bio, ])
names(prediction.points) <- c("x", "y")

# Save for further use
f.out <- paste(hbm.output.dir, "prediction-at-points.csv", sep="")
write.csv(prediction.points, f.out, row.names = FALSE, quote = FALSE)
cat("==>", f.out, "\n")

#### MODEL DEVELOPMENT #############################################################################
developHBMs <- function(gc, data, mesh, spde, n.threads, force.run){
    cat("\n\n")
    cat(paste0(rep("=", 80), collapse = ""))
    cat(paste("\nModel for", toupper(gc)), "\n")
    cat(paste0(rep("-", 80), collapse = ""), "\n\n")
    #--- INLA ESTIMATION STACK --------------------------------------------------------------------#
    # Spatial
    stk.est <- getINLAEstimationStack(gc, mesh, data, spde)

    #--- COMPARISON OF ALL POSSIBLE MODELS --------------------------------------------------------#
    modelling.variables <- variables
    # gc3 does not result in a correct spatial model (see article)
    if(gc != 'gc3'){
        modelling.variables <- c(variables, "f(s, model=spde)")
    }
    n.models <- sum(sapply(1:length(modelling.variables),
                           FUN=function(x) dim(combn(length(modelling.variables), x))[2])) + 1
    models <- getAllModels(gc, data[[gc]], modelling.variables,
                           stk.est, n.models, n.threads, T)

    #--- SELECTION OF BEST MODEL AND FINAL MODEL ESTIMATION ---------------------------------------#
    #... Non-Spatial ..............................................................................#
    formula.ng <- selectBestModel(gc, models[["models.by.lcpo"]], 5, F)
    mod.est.ng <- fitModel(gc, formula.ng, stk.est, n.threads, force.run)
    makeEstimationModelReport(gc, mod.est.ng, "non-spatial")

    #... Spatial ..................................................................................#
    mod.est.g <-  NULL
    formula.g <- NA
    if(gc != 'gc3'){ # No spatial model for GC3 (see article)
        formula.g <- selectBestModel(gc, models[["models.by.lcpo"]], 5, T)
        mod.est.g <- fitModel(gc, formula.g, stk.est, n.threads, force.run)
        makeEstimationModelReport(gc, mod.est.g, "spatial")
        makePosteriorSpatialParametersReport(gc, mod.est.g, spde)
        plotPosteriorSpatialParameters(gc, mod.est.g, spde)
        projectPosteriorLatentFieldMeanAndSD(gc, proj.grid.mat, mod.est.g, T)
    }
    #--- MODELS SUMMARY ---------------------------------------------------------------------------#
    csvModelPerformanceSummary(gc, mod.est.ng, mod.est.g)

    model <- list("gc"=gc, "formula.ng"=formula.ng, "formula.g"=formula.g,
                  "mod.est.ng"=mod.est.ng, "mod.est.g"= mod.est.g,
                  "stk.est"=stk.est)
    return(model)
}

#### MODEL PREDICTIONS #############################################################################
makeHBMPredictions <- function(model, gcc.model, rcp, yr, force.run=F){
    # Prediction stack
    stk.pred.g <- getINLAPredictionStack(gcc.model, rcp, yr, spde, force.run)

    stk <- inla.stack(model[["stk.est"]], stk.pred.g)

    gc <- model[["gc"]]

    # Non-spatial
    formula.ng <- model[["formula.ng"]]
    mod.est.ng <- model[["mod.est.ng"]]
    cat(paste("Making non-spatial prediction model for", gc, gcc.model, rcp, yr), "\n")
    mod.pred.ng <- makePredictionModel(gc, gcc.model, rcp, yr, formula.ng,
                                       mod.est.ng, stk, force.run)
    makePredictiveModelReport(gc, gcc.model, rcp, yr, mod.pred.ng, mod.est.ng,
                              stk, type="non-spatial")
    predictScenario(gc, gcc.model, rcp, yr, proj.grid.mat, i.map, i.bio, mod.pred.ng, stk, F, T)

    # Spatial
    if(gc != 'gc3'){ # No spatial model for GC3 (see article)
        formula.g <- model[["formula.g"]]
        mod.est.g <- model[["mod.est.g"]]
        cat(paste("Making spatial prediction model for", gc, gcc.model, rcp, yr), "\n")
        mod.pred.g <- makePredictionModel(gc, gcc.model, rcp, yr, formula.g, mod.est.g,
                                          stk, force.run)
        makePredictiveModelReport(gc, gcc.model, rcp, yr, mod.pred.g, mod.est.g,
                                  stk, type="spatial")
        predictScenario(gc, gcc.model, rcp, yr, proj.grid.mat, i.map, i.bio, mod.pred.g, stk, T, T)
    }
}

#### MAIN ##########################################################################################
executeHBMModels <- function(force.run){
    for(gc in genetic.clusters){
        model <- developHBMs(gc, data, mesh, spde, 4, force.run)
        cat("\nPredicting to present time (2000)\n")
        makeHBMPredictions(model, "current", "na", "00") # Prediction to present time
        for(mo in c("mean")){
            for(rcp in c("26", "85")){
                for(yr in "70"){
                    cat(paste0(rep("-", 80), collapse = ""))
                    cat("\nMaking predictions ...\n")
                    cat("\nMaking predictions for", model[["gc"]], mo, rcp, yr, "\n")
                    cat(paste0(rep("-", 80), collapse = ""), "\n\n")
                    makeHBMPredictions(model, mo, rcp, yr, force.run)
                }
            }
        }
    }
    writePredictionResults()
    writeModelsPredictiveMetrics(data)
}
