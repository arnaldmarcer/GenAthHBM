#### DATA PREPARATION ##############################################################################
prepareMaxentAccessionsData <- function(data){
    for(gc in genetic.clusters){
        f.out <- paste(maxent.output.dir, gc, "_0.5_modelling-presences-bg.csv", sep="")
        if(file.exists(f.out))
            next
        cat("\nPreparing modelling data for", toupper(gc))
        data$pa <- ifelse(data[[gc]] >= th, 1, 0)
        occ <- data.frame(data[data$pa==1, c("pa", "x", "y", variables)])
        tifs <- paste(data.dir, "climate/2000/", variables, ".tif", sep="")
        env <- stack(tifs)
        spdf.env <- data.frame(as(env, "SpatialPixelsDataFrame"))
        # Get rid of cells with presences before selecting background points
        spdf.env <- anti_join(spdf.env, occ, by=c("x", "y"))
        set.seed(1234)
        bg <- sample_n(spdf.env, 10000)
        bg$pa <- 0
        bg <- bg[, c("pa", "x", "y", variables)]
        occ <- occ[, c("pa", "x", "y", variables)]
        df.mx <- rbind(occ, bg)
        write.csv(df.mx, f.out, row.names=F, quote=F)
        cat("\n==>", f.out)
    }
}

prepareMaxentPredictorsData <- function(){
    f.out <- paste0(maxent.output.dir, "predictor_values_1km-resolution_current-na-00.csv")
    if(!file.exists(f.out)){
        cat("\nPreparing predictors data")
        tifs <- paste0(data.dir, "climate/2000/", variables, ".tif")
        st.env <- stack(tifs)
        spdf.env <- data.frame(as(st.env, "SpatialPixelsDataFrame"))
        write.csv(spdf.env, f.out, quote=F, row.names=F)
        cat("\n==>", f.out)
    }

    for(rcp in c("2.6", "8.5")){
        f.out <- paste(maxent.output.dir, "predictor_values_1km-resolution_mean-",
                       gsub("\\.", "", rcp), "-70.csv", sep="")
        if(!file.exists(f.out)){
            tifs <- paste(data.dir, "climate/2070/mean", gsub("\\.", "", rcp),
                          "bi70/", variables, ".tif", sep="")
            st.env <- stack(tifs)
            spdf.env <- data.frame(as(st.env, "SpatialPixelsDataFrame"))
            write.csv(spdf.env, f.out, quote=F, row.names=F)
            cat("\n==>", f.out)
        }
    }
}

#### MODEL DEVELOPMENT #############################################################################
developMaxentModels <- function(data){
    models <- getAllPossibleModels(variables)
    for(model in models){
        vars <- unlist(str_split(model, " \\+ "))
        for(gc in genetic.clusters){
            path <- paste(maxent.output.dir, gc, "_", th, "_5-fold-crossvalidation_",
                          paste(vars, collapse="-"), "/", sep="")
            if(file.exists(path)){
                cat("\n--> Model already exists:", path)
                next
            }

            cat(paste("\nGenerating", toupper(gc), "model:", model))

            df.mx <- read.csv(paste(maxent.output.dir, gc, "_", th,
                                    "_modelling-presences-bg.csv", sep=""))

            unlink(path, recursive=TRUE)
            dir.create(path, recursive = T)

            df.env <- df.mx %>% select(one_of(vars))

            mx <- dismo::maxent(data.frame(df.env), p=df.mx$pa, path=path,
                                args=c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                       "hinge=TRUE", "threshold=FALSE",
                                       "addSamplesToBackground=FALSE", "replicates=5"))

            cat(paste("\n==>", toupper(gc), ":", model))
        }
    }
}

generateFinalBestMaxentModels <- function(df, best.models, force.run=F){
    for(gc in genetic.clusters){
        vars <- best.models[best.models$gc == gc & best.models$th == th  &
                                best.models$best == '*', "model"]
        if(nrow(vars) == 0)
            next
        vars <- unlist(str_split(vars, "-"))
        df <- read.csv(paste0(maxent.output.dir, gc, "_", th, "_modelling-presences-bg.csv"))
        cat("\nGenerating final model for: ", toupper(gc), " (",
            paste(vars, collapse=", "), ")", sep="")
        m <- getMaxentModel(df, vars, gc, th, "best-maxent-model", force.run)
        if(is.null(m))
            cat("\n--> Not enough data for", gc, th)
    }
}

#### MODEL PREDICTIONS #############################################################################
maxentPredict <- function(gc, th, gcc, rcp, year, file.suffix, force.run=F){
    f.out <- paste(maxent.output.dir, gc, "_", th, "_", gcc, "_", rcp, "_", year, "_prediction_",
                   file.suffix, ".csv", sep="")
    if(file.exists(f.out) & !force.run)
        return()

    f.model <- paste(maxent.output.dir, gc, "_", th, "_final-", file.suffix, ".rds", sep="")
    if(!file.exists(f.model))
        return(NULL)

    m <- readRDS(f.model)
    f.in <- paste(maxent.output.dir, "predictor_values_1km-resolution_", gcc, "-",
                  gsub("\\.", "", rcp), "-", year, ".csv", sep="")
    env <- data.frame(read_csv(f.in, col_types = cols(), progress = FALSE))
    env <- env[complete.cases(env), ]
    pred <- dismo::predict(m, env, progress="text")
    df.pred <- data.frame(cbind(pred=pred, env[, c("x", "y")]))

    write_csv(df.pred, f.out)
    cat("\n==>", f.out)

    return(df.pred)
}

# For best models file.suffix = "best-maxent-model"
# For full models file.suffix = "full-maxent-model"
makeMaxentPredictions <- function(file.suffix, force.run){
    for(gc in genetic.clusters){
        # Current predictions
        cat("\nMaking prediction for", file.suffix, toupper(gc), th, "current", "na", "00")
        maxentPredict(gc, th, "current", "na", "00", file.suffix, force.run)

        # Future predictions
        for(gcc in c("mean")){
            for(rcp in c("2.6", "8.5")){
                for(yr in c("70")){
                    cat("\nMaking prediction for", file.suffix, toupper(gc), th, gcc, rcp, yr)
                    maxentPredict(gc, th, gcc, rcp, yr, file.suffix, force.run)
                }
            }
        }
    }
}

#### MAIN ##########################################################################################
executeMaxentModels <- function(force.run=F){
    # Data preparation
    data <- read.csv(paste0(data.dir, "sp/ath_accessions.csv"))
    prepareMaxentAccessionsData(data)
    prepareMaxentPredictorsData()

    # Model development
    developMaxentModels(data)
    bm <- selectBestMaxentModels(force.run)
    generateFinalBestMaxentModels(data, bm, force.run)
    maxentVariableContributions(force.run)

    # Model predictions
    makeMaxentPredictions("best-maxent-model", force.run)

    cat("\n==> Finished Maxent models and results")
}

