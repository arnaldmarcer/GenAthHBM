#### STUDY AREA BOUNDARIES, MESH AND SPDE ##########################################################
getStudyAreaBoundaries <- function(force.run=F){
    f.out <- paste(hbm.output.dir, "ip_boundaries.rds", sep="")
    if(file.exists(f.out) & !force.run)
        return(readRDS(f.out))
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    ip.pol <- Polygons(list(Polygon(ip.df)), ID="Península Ibérica")
    ip.pol <- SpatialPolygons(list(ip.pol), pO = 1:length(list(ip.pol)),
                              proj4string=CRS(as.character(NA)))
    proj4string(ip.pol) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000
                          +ellps=GRS80 +units=m +no_defs"
    saveRDS(ip.pol, f.out)
    cat("==>", f.out, "\n")
    return(ip.pol)
}

createMesh <- function(ip.boundaries, crs.proj, force.run=F){
    f.mesh <- paste(hbm.output.dir, "mesh.rds", sep="")
    if(file.exists(f.mesh) & !force.run){
        mesh <- readRDS(f.mesh)
        return(mesh)
    }
    cat("\n###Creating the mesh\n")
    ip.inla <- inla.sp2segment(ip.boundaries, crs=crs.proj)

    mesh <- inla.mesh.2d(boundary=ip.inla,
                         max.edge=c(30000, 60000),
                         cutoff=10000,
                         offset=c(10000, 300000), crs=crs.proj)
    saveRDS(mesh, file=f.mesh)
    cat("==>", f.mesh, "\n")
    return(mesh)
}

# SPDE - Stochastic Partial Differential Equation
# We define SPDE over the generated mesh. It will be used to project our data
# Prior distributions:
#  * range: We calculate the mesh diameter and divide it by 2 (log(range0)), just
#           as a starting value. The prior is given in the logarithmic scale,
#           hence log(range) ~ Normal(log(range0), precision=0.25)
#  * sigma: This is the standard deviation of the spatial effects, usually
#           assigned to 1 (sigma0). The prior distribution is also given in the
#           logarithmic scale, hence:
#               log(sigma) ~ Normal (log(sigma0), precision=0.25).
#
# SPDE is defined as a function of two parameters, kappa and tau, which in turn
# depend on range and standard deviaton of the spatial effect. These parameters
# values are internally used for a better computation:
#     range <- sqrt(8)/kappa
#     sigma^2 <- 1/(4pi kappa^2 tau^2)
#
# However, we first need to give reasonable values for the priors given the
# knowledge we have of the study area. The same priors for the hyperparameters of
# the spatial effect can be given to all genetic clusters since we consider them
# vague.
getSPDE <- function(mesh, force.run=F){
    f.out.spde <- paste(hbm.output.dir, "spde.rds", sep="")
    if(file.exists(f.out.spde) & !force.run){
        spde <- readRDS(f.out.spde)
        return(spde)
    }
    cat("\n### Preparing the SPDE\n")
    # Range total diameter
    size = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))

    # Mean range
    range0 <-  size / 2

    # Range in between 10 km and total size
    range0_in <-c (10000, size)

    # Range's prior distribution
    range0.mar <-
        inla.tmarginal(function(x) exp(x),
                       transform(data.frame(x = log(range0) + seq(-2, 2, length=2021)),
                                 y = dnorm(x, mean = log(range0),
                                           sd = 2)))

    kappa0 <- sqrt(8) / range0

    # Save as plot
    f.out.sp.range <- paste(hbm.output.dir, "spatial_effect_range_a_priori_distribution.png", sep="")
    png(f.out.sp.range, width=950, height=800)
    par(mfrow=c(1,1))
    par(mar=c(4,4,2,0))
    plot(range0.mar$x, range0.mar$y)
    cat("==>", f.out.sp.range, "\n")
    dev.off()

    # Spatial effect variance
    # In Beta regression, the dependent variable values are bounded by 0 and 1 (0 and 1 not included)
    # We use the logit transform for the modelling process
    y_beta <- seq(0.001, 0.999, 0.001)
    y_beta_logit <- log(y_beta/(1-y_beta))

    sigma0 <- 1
    # We know the spatial effect's standard deviation can range from 0 to abs(-6.907) + 6.907 ~= 14
    sigma0_in <- c(0.01, 14)

    # We set sd=2, that is, a precision of 1/2^2=0.25
    sigma0.mar <-
        inla.tmarginal(function(x) exp(x),
                       transform(data.frame(x = log(sigma0) + seq(-2, 2,
                                                                  length=2000)),
                                 y = dnorm(x, mean = log(sigma0),
                                           sd = 2)))

    # Save as plot
    f.out.sp.variance <- paste(hbm.output.dir, "spatial_effect_sigma_a_priori_distribution.png", sep="")
    png(f.out.sp.variance, width=950, height=800)
    par(mfrow=c(1,1))
    par(mar=c(4,4,2,0))
    plot(sigma0.mar$x, sigma0.mar$y)
    cat("==>", f.out.sp.variance, "\n")
    dev.off()

    # Reparametrize to parameters used by INLA
    kappa0 <- sqrt(8) / range0
    tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma0)

    # Define SPDE on the mesh and set the a priori distributions
    spde = inla.spde2.matern(mesh,
                             B.tau = cbind(log(tau0), -1, +1),
                             B.kappa = cbind(log(kappa0), 0, -1),
                             theta.prior.mean = c(0, 0),
                             theta.prior.prec = c(0.25, 0.25) )

    # Save object for further use
    saveRDS(spde, f.out.spde)
    cat("==>", f.out.spde, "\n")
    return(spde)
}

#### MAXENT MODEL FITTING ##########################################################################
getAllPossibleModels <- function(variables){
    models <- c()
    for(n in 1:length(variables)){
        mc <- t(combn(variables, n))
        for(i in 1:nrow(mc)){
            models <- c(models, paste(mc[i,], collapse = ' + '))
        }
    }
    return(models)
}

generateEvaluationReport <- function(force.run=F){
    f.out <- paste(maxent.output.dir, "maxent_models_evaluation.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read.csv(f.out))

    models <- getAllPossibleModels(variables)

    df.res <- data.frame()
    for(model in models){
        for(gc in genetic.clusters){
            vars <- unlist(str_split(model, " \\+ "))
            dir <- paste(maxent.output.dir, gc, "_", th, "_5-fold-crossvalidation_",
                         paste(vars, collapse="-"), "/", sep="")
            if(!file.exists(dir))
                next
            f <- paste(dir, "maxentResults.csv", sep="")
            df <- read.csv(f)[, c("Species", "X.Training.samples", "X.Test.samples", "Iterations",
                                  "Test.gain", "Test.AUC", "AUC.Standard.Deviation")]
            names(df) <- c("species", "n.train", "n.test", "fold", "gain", "auc.mean", "auc.sd")
            df <- df[df$species != "species (average)",]
            df$fold <- apply(df, 1, FUN=function(x) strsplit(x[["species"]], "\\_")[[1]][2])
            df$gc <- gc
            df$th <- th
            df$model <- paste(vars, collapse="-")
            df$dir <- dir
            df <- df[, c("gc", "th", "fold", "n.train", "n.test", "gain", "auc.mean", "auc.sd",
                         "model", "dir")]
            df.res <- rbind(df.res, df)
        }
    }
    write.csv(df.res, f.out)
    cat("\n==>", f.out)
    return(df.res)
}

# It orders the best five models for each gc and groups according to auc and number of predictors
# and then selects the best one (lowest number of predictors)
selectBestMaxentModels <- function(force.run=F){
    ev <- generateEvaluationReport(force.run)
    # Best five according to auc
    sev <- ev %>% mutate(n=n.train + n.test) %>% group_by(gc, th, n, model) %>%
        summarise(auc=mean(auc.mean), auc.sd=mean(auc.sd)) %>%
        group_by(gc, th, n) %>% arrange(gc, th, desc(auc)) %>% filter(row_number() <= 5)
    sev <- sev %>% group_by(gc, th, n, model) %>%
        mutate(l=length(unlist(strsplit(as.character(model), '-'))))
    sev <- sev %>% group_by(gc, th) %>% mutate(best=ifelse(l==min(l),"*","")) %>%
        select(gc, th, n, model, auc, auc.sd, best)
    # In case two get the '*' per group
    sev <- sev %>% arrange(gc, th, desc(best)) %>% mutate(best=ifelse(row_number() == 1, '*', ''))
    return(sev)
}

getMaxentModel <- function(df, vars, gc, th, file.suffix, force.run=F){
    f.out <- paste(maxent.output.dir, gc, "_", th, "_final-", file.suffix, ".rds", sep="")
    if(file.exists(f.out) & !force.run){
        return(readRDS(f.out))
    }
    path <- paste(maxent.output.dir, gc, "_", th, "_", file.suffix, "/", sep="")
    if(file.exists(path))
        unlink(path, recursive=TRUE)

    dir.create(path, recursive = T)

    df.env <- data.frame(df %>% select(one_of(vars)))
    mx <- dismo::maxent(df.env, p=df$pa, path=path,
                        args=c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                               "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE"))

    saveRDS(mx, f.out)
    cat("\n==>", f.out)
    return(mx)
}

maxentVariableContributions <- function(force.run){
    f.out <- paste0(maxent.output.dir, "variable_contributions.csv")
    if(file.exists(f.out) & !force.run)
        return(read.csv(f.out))
    df.res <- data.frame()
    for(gc in genetic.clusters){
        f.maxent <- paste0(maxent.output.dir, gc, "_0.5_", "best-maxent-model/maxentResults.csv")
        df <- read.csv(f.maxent)
        df <- df %>% select(one_of(names(df)[grep("contribution", names(df))])) %>% gather() %>%
            mutate(Cluster=gc, Model="Maxent", metric="percentage.contribution",
                   key=gsub("\\.contribution", "", key), metric='Percentage contribution') %>%
            select(Cluster, Model, "Predictor"=key, metric, value)

        df.res <- rbind(df.res, df)
    }

    write.csv(df.res, f.out, row.names = F, quote = F)
    cat("\n==>", f.out)
    return(df.res)
}

#### HBM MODEL FITTING ############################################################################
getINLAEstimationStack <- function(gc, mesh, data, spde){
    f.out <- paste(hbm.output.dir, gc, "_inla_estimation_stack.rds", sep="")
    A.est.g <- inla.spde.make.A(mesh, loc=as.matrix(data[,2:3]))
    # Spatial index to be used in the model
    ind <- inla.spde.make.index(name='s', n.spde=spde$n.spde)
    stk.est.g <- inla.stack(data=list(resp=data[[gc]]),
                            A=list(A.est.g, 1),
                            effects=list(c(ind, list(intercept=1)), data[, variables]),
                            tag='est')
    saveRDS(stk.est.g, f.out)
    cat("==>", f.out, "\n")
    return(stk.est.g)
}

getAllModels <- function(gc, y, variables, stk.est.g, n.models, n.threads, force.run=F){
    f.out <- paste(hbm.output.dir, gc, "_models.rds", sep="")
    if(file.exists(f.out) & !force.run){
        models <- readRDS(f.out)
        return(models)
    } else {
        models <- buildModels(
            gc,
            resp=y,
            variables=variables,
            data=inla.stack.data(stk.est.g), n=n.models,
            family="beta",
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk.est.g),link=1),
            control.fixed=list(mean=0, prec=10^(-4), mean.intercept=0, prec.intercept=10^(-4)),
            control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=n.threads,
            control.inla=list(strategy="simplified.laplace"),
            verbose=FALSE)
        saveRDS(models, f.out)
        cat(paste("\n==>", f.out))
    }
    return(models)
}

# BUILD ALL MODELS
# Given a formula returns the best model according to DIC
# * gc: cluster name
# * resp: Response variable
# * variables: Predictor variables in the models
# * data: Modelling data
# * n: Number of models to be listed in the output
# * family: distribution family of the response variable, binomial by default
# * ...: INLA function arguments
buildModels <- function(gc, resp, variables, data, n, family="binomial", ...){

    sel.terms <- switch('terms', terms=variables)

    # All combinations of m elements in v
    comb.terms <- function(m, v=sel.terms) {
        if(m==0) return('resp ~ -1 + intercept')
        else {
            combis <- apply(combn(v, m), 2, paste, collapse=' + ')
            return(paste('resp ~ -1 + intercept', combis, sep=' + '))
        }
    }

    # List of all possible models
    f.list <- unlist(sapply(0:length(sel.terms), comb.terms))

    # Try each model in f.list and keep the DIC value
    dic <- numeric()
    lcpo <- numeric()
    waic <- numeric()
    for(i in 1:length(f.list)){
        cat("\n")
        cat(paste0("--> Model ", str_pad(i, 3, pad='0'), "/", length(f.list)), ": ", f.list[i])
        res <- try(inla(formula = eval(parse(text=f.list[i])),
                        family=family,
                        data=data, ...),
                   silent=TRUE)
        if(inherits(res, 'try-error')){
            dic[i] <- NA
            lcpo[i] <- NA
            waic[i] <- NA
        }else{
            dic[i] <- res$dic$dic
            lcpo[i] <- -mean(log(res$cpo$cpo))
            waic[i] <- res$waic$waic

        }

        formula <- f.list[i]
        f.name <- trim(sub("-", "",
                           gsub(" \\+ ", "-",
                                gsub(")", "",
                                     gsub("f\\(s, model=", "",
                                          gsub("resp ~ -1 \\+ intercept", "int-", formula))))))
        individual.models.dir <- paste0(hbm.output.dir, "individual_models/", gc)
        dir.create(individual.models.dir, recursive=T, showWarnings = F)
        f.out <- paste0(individual.models.dir, "/m_", i, "_", f.name, ".rds")
        saveRDS(res, f.out)
        cat("\n  ==>", f.out)
    }
    if(any(is.na(dic))){
        cat("\nWARNING There are models with a non-sense spatial effect:\n")
        cat(f.list[which(is.na(dic))], sep="\n")
    }
    # List models ordered by DIC
    models.dic <- data.frame(f.list[order(dic)[1:n]],
                             dic[order(dic)[1:n]],
                             waic[order(dic)[1:n]],
                             lcpo[order(dic)[1:n]])
    colnames(models.dic) <- c("Models", "DIC", "WAIC", "LCPO")

    # List models ordered by WAIC
    models.waic <- data.frame(f.list[order(waic)[1:n]],
                              dic[order(waic)[1:n]],
                              waic[order(waic)[1:n]],
                              lcpo[order(waic)[1:n]])
    colnames(models.waic)<-c("Models", "DIC", "WAIC", "LCPO")

    # List models ordered by LCPO
    models.lcpo <- data.frame(f.list[order(lcpo)[1:n]],
                              dic[order(lcpo)[1:n]],
                              waic[order(lcpo)[1:n]],
                              lcpo[order(lcpo)[1:n]])
    colnames(models.lcpo)<-c("Models", "DIC", "WAIC", "LCPO")


    models <- list(models.dic, models.waic, models.lcpo)
    names(models) <- c("models.by.dic", "models.by.waic", "models.by.lcpo")
    return(models)
}

selectBestModel <- function(gc, models, top.n, spatial=T){
    selected.models <- getBestModels(gc, models, top.n, spatial)
    formula <- selected.models %>% filter(best=='*') %>% dplyr::select(Models) %>%
        first() %>% as.character()

    return(formula)
}

getBestModels <- function(gc, models, top.n, spatial=T){
    if(spatial){
        models <- models %>% filter(str_detect(Models, "spde"))
    }else{
        models <- models %>% filter(!str_detect(Models, "spde"))
    }
    selected.models <- models %>% mutate(gc=gc) %>%
        arrange(LCPO) %>% filter(row_number() <= top.n)
    selected.models <- selected.models %>% group_by(Models) %>%
        mutate(l=length(unlist(strsplit(as.character(Models), '\\+')))) %>%
        arrange(l) %>% ungroup() %>% mutate(best=ifelse(row_number()==1, '*', '')) %>%
        dplyr::select(gc, Models, DIC, WAIC, LCPO, best, -l) %>%
        arrange(gc, LCPO)
    return(selected.models)
}

fitModel <- function(gc, formula, stk.est.g, n.threads, force.run){
    if(str_detect(formula, "spde")){
        f.out <- paste(hbm.output.dir, gc, "_spatial_fit.rds", sep="")
    } else {
        f.out <- paste(hbm.output.dir, gc, "_non-spatial_fit.rds", sep="")
    }
    if(file.exists(f.out) & !force.run)
        return(readRDS(f.out))

    model <- inla(formula=as.formula(formula),
                  data=inla.stack.data(stk.est.g),
                  family="beta" ,
                  control.compute = list(config=TRUE, cpo=TRUE, dic=TRUE,
                                         return.marginals=TRUE, waic=TRUE),
                  control.inla=list(strategy='simplified.laplace'),
                  control.predictor=list(A=inla.stack.A(stk.est.g), compute=TRUE),
                  control.fixed=list(mean=0, prec=10^(-4),
                                     mean.intercept=0, prec.intercept=10^(-4)),
                  num.threads = n.threads,
                  verbose=FALSE)
    saveRDS(model, f.out)
    cat("\n==>", f.out, "\n")
    return(model)
}

makeEstimationModelReport <- function(gc, model, type="NA"){
    f.out <- paste(hbm.output.dir, gc, "_", type, "_fit_report.txt", sep="")

    capture.output(paste(rep('#', 80), collapse = ""), file=f.out)
    capture.output(print("summary(model)"), file=f.out, append = T)
    capture.output(summary(model), file=f.out, append = T)

    ### Model variable relevance
    capture.output(paste(rep('#', 80), collapse = ""), file=f.out, append = T)
    for(v in 2:length(names(model$marginals.fixed))){
        var <- names(model$marginals.fixed)[v]

        capture.output(print(paste("1-inla.pmarginal(0, model$marginals.fixed$", var, ")", sep="")),
                       file=f.out, append = T)
        capture.output(print(1-inla.pmarginal(0, model$marginals.fixed[[var]])),
                       file=f.out, append = T)
        capture.output(print(paste("inla.pmarginal(0, model$marginals.fixed$", var, ")", sep="")),
                       file=f.out, append = T)
        capture.output(print(inla.pmarginal(0, model$marginals.fixed[[var]])),
                       file=f.out, append = T)
    }

    ### Model comparison metrics
    capture.output(paste(rep('#', 80), collapse = ""), file=f.out, append = T)
    capture.output(print("model$dic$dic"), file=f.out, append = T)
    capture.output(model$dic$dic, file=f.out, append = T)
    capture.output(print("model$waic$waic"), file=f.out, append = T)
    capture.output(model$waic$waic, file=f.out, append = T)
    capture.output(print("-mean(log(model$cpo$cpo))"), file=f.out, append = T)
    capture.output(-mean(log(model$cpo$cpo)), file=f.out, append = T)
    cat("\n==>", f.out, "\n")

    ### Keep verbose to know check INLA computation
    f.out <- paste(hbm.output.dir, gc, "_", type, "_fit_verbose.txt", sep="")
    write.table(model$logfile, file=f.out)
    cat("\n==>", f.out, "\n")
}

makePosteriorSpatialParametersReport <- function(gc, mod.est.g, spde){
    f.out <- paste(hbm.output.dir, gc, "_selected_model_posterior.txt", sep="")
    spde.result = inla.spde2.result(mod.est.g, "s", spde, do.transform=TRUE)

    # Spatial effects variance
    sigma2 <- inla.emarginal(function(x) x, spde.result$marginals.variance.nominal[[1]])

    # Spatial effects standard deviation
    sigma <- inla.emarginal(function(x) sqrt(x), spde.result$marginals.variance.nominal[[1]])

    # kappa
    kappa <- inla.emarginal(function(x) x, spde.result$marginals.kappa[[1]])

    # Spatial effects range
    range <- inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]])
    tau <- inla.emarginal(function(x) x, spde.result$marginals.tau[[1]]) #tau

    # Credibility intervals
    kappa_in <- inla.hpdmarginal(0.95, spde.result$marginals.kappa[[1]])
    sigma2_in <- inla.hpdmarginal(0.95, spde.result$marginals.variance.nominal[[1]])
    range_in <- inla.hpdmarginal(0.95, spde.result$marginals.range.nominal[[1]])
    tau_in <- inla.hpdmarginal(0.95, spde.result$marginals.tau[[1]])

    fileConn <- file(f.out)
    writeLines(c("Spatial effects variance", sigma2, "",
                 "Spatial effects standard deviation", sigma, "",
                 "Kappa", kappa, "",
                 "Spatial effects range", range, "",
                 "Tau", tau, ""), fileConn)
    close(fileConn)
    cat("==>", f.out, "\n")
}

plotPosteriorSpatialParameters <- function(gc, mod.est.g, spde){
    f.out <- paste(hbm.output.dir, gc, "_selected_model_posterior.png", sep="")
    png(f.out, width=950, height=950)
    par(mfrow=c(2,2), mar=c(3,3,1,0.5))
    spde.result = inla.spde2.result(mod.est.g, "s", spde, do.transform=TRUE)
    plot(spde.result$marginals.variance.nominal[[1]], type='l', main="marginals.variance.nominal")
    plot(spde.result$marginals.kappa[[1]], type='l', main="kappa")
    plot(spde.result$marginals.range.nominal[[1]], type='l', main="marginals.range.nominal")
    plot(spde.result$marginals.tau[[1]], type='l', main="tau")
    dev.off()
    cat("==>", f.out, "\n")
}

projectPosteriorLatentFieldMeanAndSD <- function(gc, proj.grid.mat, mod.est.g, makeImagePlot=T){
    xmedia <- inla.mesh.project(proj.grid.mat, mod.est.g$summary.random$s$mean)
    xsd <- inla.mesh.project(proj.grid.mat, mod.est.g$summary.random$s$sd)

    # Only inside study area
    xmedia[!i.map] <- xsd[!i.map] <- NA

    f.out.avg <- paste(hbm.output.dir, gc, "_latent_field_avg.rds", sep="")
    saveRDS(xmedia, f.out.avg)
    cat("==>", f.out.avg, "\n")

    f.out.sd <- paste(hbm.output.dir, gc, "_latent_field_sd.rds", sep="")
    saveRDS(xsd, f.out.sd)
    cat("==>", f.out.sd, "\n")

    if(makeImagePlot){
        par(mfrow=c(1,1))
        par(mar=c(0,0,0,0))
        f.out <- paste(hbm.output.dir, gc, "_latent_field_avg.png", sep="")
        png(f.out, width=950, height=800)
        image.plot(xmedia)
        dev.off()
        cat("==>", f.out, "\n")

        f.out <- paste(hbm.output.dir, gc, "_latent_field_sd.png", sep="")
        png(f.out, width=950, height=800)
        image.plot(xsd)
        dev.off()
        cat("==>", f.out, "\n")
    }
}

hbmVariableCoefficients <- function(force.run){
    f.out <- paste0(hbm.output.dir, "variable_coefficients.csv")
    if(file.exists(f.out) & !force.run)
        return(read.csv(f.out))
    df.res <- data.frame()
    for(gc in genetic.clusters){
        for(sp in c("Spatial", "Non-spatial")){
            f.inla <- paste0(hbm.output.dir, gc, "_", tolower(sp), "_fit.rds")
            if(!file.exists(f.inla))
                next
            df <- readRDS(f.inla)$summary.fixed
            df$Predictor <- row.names(df)
            df <- df %>% mutate(Cluster=toupper(gc), Model=paste(sp, "HBMs")) %>%
                select(Cluster, Model, Predictor, mean, sd, "0.025quant", "0.5quant", "0.975quant", "mode", "kld")
            df.res <- rbind(df.res, df)
        }
    }
    write.csv(df.res, f.out, row.names = F, quote = F)
    cat("\n==>", f.out)
}

csvModelPerformanceSummary <- function(gc, mod.est.ng, mod.est.g){
    if(!is.null(mod.est.ng)){
        dic.ng <- mod.est.ng$dic$dic
        waic.ng <- mod.est.ng$waic$waic
        lcpo.ng <- - mean(log(mod.est.ng$cpo$cpo))
    } else {
        formula.ng <- NA
        dic.ng <- NA
        waic.ng <- NA
        lcpo.ng <- NA
    }
    if(!is.null(mod.est.g)){
        dic.g <- mod.est.g$dic$dic
        waic.g <- mod.est.g$waic$waic
        lcpo.g <- - mean(log(mod.est.g$cpo$cpo))
    } else {
        formula.g <- NA
        dic.g <- NA
        waic.g <- NA
        lcpo.g <- NA
    }

    summary <- data.frame(MODEL=c("Spatial model", "Best non-spatial model"),
                          DIC=c(dic.g, dic.ng),
                          WAIC=c(waic.g, waic.ng),
                          LCPO=c(lcpo.g, lcpo.ng))
    f.out <- paste(hbm.output.dir, gc, "_models_summary.csv", sep="")
    write.csv(summary, f.out, quote=F, row.names = F)
    cat("==>", f.out, "\n")
}

#### MODEL PREDICTIONS #############################################################################
preparePredictorValuesAtPredictionPoints <- function(gcc.model, rcp, yr){
    if(yr=="00"){
        tifs <- paste0(data.dir, "climate/2000/", variables, ".tif")
    } else {
        tifs <- paste0(data.dir, "climate/2070/", gcc.model,
                       rcp, "bi", yr, "/", variables, ".tif")
    }

    st.bioclim <- stack(tifs)
    points <- read.csv(paste0(hbm.output.dir, "prediction-at-points.csv"))
    env <- cbind(points, raster::extract(st.bioclim, points))
    f.out <- paste0(hbm.output.dir, "predictor_values_1km-resolution_",
                    gcc.model, "-", rcp, "-", yr, ".csv")
    write_csv(env, f.out)
    cat("==>", f.out, "\n")
    return(env)
}

getPredictorValues <- function(gcc.model, rcp, year, force.run=FALSE){
    f.out <- paste(hbm.output.dir, "predictor_values_1km-resolution_",
                   gcc.model, "-", rcp, "-", year, ".csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols(), progress = FALSE))
    env <- preparePredictorValuesAtPredictionPoints(gcc.model, rcp, year)
    return(env)
}

getINLAPredictionStack <- function(gcc.model, rcp, yr, spde, force.run=F){
    f.out <- paste(hbm.output.dir, gcc.model, "_", rcp, "_", yr, "_inla_prediction_stack.rds", sep="")
    if(!file.exists(f.out) | force.run){
        # Data frame of values at prediction points plus intercept column
        df.pred <- getPredictorValues(gcc.model, rcp, yr)

        A.pred.g <- inla.spde.make.A(mesh, loc=as.matrix(df.pred[, c("x", "y")]))

        stk.pred.g <- inla.stack(
            data=list(resp=NA),
            A=list(A.pred.g, 1),
            effects=list(s=1:spde$n.spde,
                         cbind(intercept=rep(1, nrow(df.pred)), df.pred[, variables])),
            tag='pred')
        saveRDS(stk.pred.g, f.out)
        cat("==>", f.out, "\n")
    } else {
        stk.pred.g <- readRDS(f.out)
    }
    return(stk.pred.g)
}

makePredictionModel <- function(gc, gcc.model, rcp, yr, formula, mod.est, stk, force.run){
    if(str_detect(formula, "spde")){
        f.out.predictive.model <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                                         "_spatial_predictive_model.rds")
    } else {
        f.out.predictive.model <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                                         "_non-spatial_predictive_model.rds")
    }

    if(file.exists(f.out.predictive.model) & !force.run){
        mod.pred <- readRDS(f.out.predictive.model)
    } else {
        # Predictive model
        mod.pred <- inla(as.formula(formula),
                         data=inla.stack.data(stk), family="beta" ,
                         control.compute = list(config=TRUE, cpo=TRUE, dic=TRUE, waic=TRUE),
                         control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),
                         control.inla=list(strategy="simplified.laplace"),
                         control.mode=list(theta=mod.est$mode$theta, restart=TRUE),
                         control.fixed=list(mean=0, prec=10^(-4), mean.intercept=0,
                                            prec.intercept=10^(-4)),
                         control.results=list(return.marginals.random=FALSE,
                                              return.marginals.predictor=FALSE),
                         num.threads = n.threads,
                         verbose=FALSE)

        saveRDS(mod.pred, f.out.predictive.model)
        cat("==>", f.out.predictive.model, "\n")
    }
    return(mod.pred)
}

makePredictiveModelReport <- function(gc, gcc.model, rcp, yr, model.pred,
                                      model.est, stk.g, type="NA"){
    f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr, "_", type,
                    "_model_prediction_report.txt")
    fileConn <- file(f.out)
    writeLines(c("model.pred$dic$dic", model.pred$dic$dic, "",
                 "model.pred$waic$waic", model.pred$waic$waic, "",
                 "-mean(log(model.pred$cpo$cpo))", -mean(log(model.pred$cpo$cpo)), ""), fileConn)
    close(fileConn)

    capture.output(print(""), file=f.out, append = T)
    capture.output(paste(rep("-",80), collapse = ""), file=f.out, append = T)
    capture.output(print("summary(model.pred)"), file=f.out, append = T)
    capture.output(summary(model.pred), file=f.out, append = T)

    capture.output(print(""), file=f.out, append = T)
    capture.output(paste(rep("-",80), collapse = ""), file=f.out, append = T)
    capture.output(print("summary(model.est)"), file=f.out, append = T)
    capture.output(summary(model.est), file=f.out, append = T)

    capture.output(print(""), file=f.out, append = T)
    capture.output(paste(rep("-",80), collapse = ""), file=f.out, append = T)
    capture.output(print("Mu"), file=f.out, append = T)
    capture.output(
        print("summary(model.pred$summary.fitted.values$mean[inla.stack.index(stk.g,'pred')$data])"),
        file=f.out, append = T)
    capture.output(
        summary(model.pred$summary.fitted.values$mean[inla.stack.index(stk.g, 'pred')$data]),
        file=f.out, append = T)

    capture.output(print(""), file=f.out, append = T)
    capture.output(paste(rep("-",80), collapse = ""), file=f.out, append = T)
    capture.output(print("Sd"), file=f.out, append = T)
    capture.output(
        print("summary(model.pred$summary.fitted.values$sd[inla.stack.index(stk.g, 'pred')$data])"),
        file=f.out, append = T)
    capture.output(
        summary(model.pred$summary.fitted.values$sd[inla.stack.index(stk.g, 'pred')$data]),
        file=f.out, append = T)
    cat("==>", f.out, "\n")

    # Keep verbose
    f.out <- paste(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr, "_", type,
                   "_model_prediction_verbose.txt", sep="")
    write.table(model.pred$logfile, file=f.out)
    cat("==>", f.out, "\n")
}

predictScenario <- function(gc, gcc.model, rcp, yr, proj.grid.mat, i.map, i.bio,
                            mod.pred, stk, spatial=T, makeImagePlot=T){
    # Predictions over mean
    # Create a matrix to visualize
    mu.mean <- mu.sd <- matrix(NA, proj.grid.mat$lattice$dims[1], proj.grid.mat$lattice$dims[2])

    # Set matrix values
    # Predictive distribution mean for mu
    mu.mean[i.map][i.bio] <-
        c(mod.pred$summary.fitted.values$mean[inla.stack.index(stk, 'pred')$data])

    if(spatial){
        f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                        "_spatial_prediction_mu_mean_matrix.rds")
    } else {
        f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                        "_non-spatial_prediction_mu_mean_matrix.rds")
    }
    saveRDS(mu.mean, f.out)
    cat("==>", f.out, "\n")

    # Predictive distribution standard deviation for mu
    mu.sd[i.map][i.bio] <-
        c(mod.pred$summary.fitted.values$sd[inla.stack.index(stk, 'pred')$data])

    if(spatial){
        f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                        "_spatial_prediction_mu_sd_matrix.rds")
    } else {
        f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                        "_non-spatial_prediction_mu_sd_matrix.rds")
    }
    saveRDS(mu.sd, f.out)
    cat("==>", f.out, "\n")

    if(makeImagePlot){
        par(mfrow=c(1,1))
        par(mar=c(0,0,0,0))
        if(spatial){
            f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                            "_spatial_prediction_mu_mean.png")
        } else {
            f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                            "_non-spatial_prediction_mu_mean.png")
        }
        png(f.out, width=950, height=800)
        image.plot(mu.mean)
        dev.off()
        cat("==>", f.out, "\n")

        if(spatial){
            f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                            "_spatial_prediction_mu_sd.png")
        } else {
            f.out <- paste0(hbm.output.dir, gc, "_", gcc.model, "_", rcp, "_", yr,
                            "_non-spatial_prediction_mu_sd.png")
        }
        png(f.out, width=950, height=800)
        image.plot(mu.sd)
        dev.off()
        cat("==>", f.out, "\n")

    }
}

#### RESULTS PREPARATION ###########################################################################
# Generates a csv file by genetic cluster with current and future predictions
preparePredictionResultsHelpingFunction <- function(mu.mean, mu.sd){
    iberia <- getStudyAreaBoundaries()
    # We need to rotate the matrix counterclockwise, if not it is printed flipped
    mu.mean <- t(mu.mean)[,nrow(mu.mean):1]
    mu.mean <- t(mu.mean)[,nrow(mu.mean):1]
    mu.mean <- t(mu.mean)[,nrow(mu.mean):1]
    r.mean <- raster(mu.mean)
    extent(r.mean) <- bbox(iberia)
    crs(r.mean) <- "+init=epsg:3035"
    df.mean <- data.frame(as(r.mean, "SpatialPixelsDataFrame"))
    names(df.mean) <- c("v", "x", "y")
    df.mean$var <- "mean"

    mu.sd <- t(mu.sd)[,nrow(mu.sd):1]
    mu.sd <- t(mu.sd)[,nrow(mu.sd):1]
    mu.sd <- t(mu.sd)[,nrow(mu.sd):1]
    r.sd <- raster(mu.sd)
    extent(r.sd) <- bbox(iberia)
    crs(r.sd) <- "+init=epsg:3035"
    df.sd <- data.frame(as(r.sd, "SpatialPixelsDataFrame"))
    names(df.sd) <- c("v", "x", "y")
    df.sd$var <- "sd"

    df <- rbind(df.mean, df.sd)

    return(df)
}

preparePredictionResults <- function(gc, force.run=F){
    f.out <- paste(hbm.output.dir, gc, "_prediction_results.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols()))

    df.res <- data.frame()
    for(sp in c("spatial", "non-spatial")){
        f.mean <- paste(hbm.output.dir,
                        gc, "_current_na_00_", sp, "_prediction_mu_mean_matrix.rds", sep="")
        f.sd <- paste(hbm.output.dir,
                      gc, "_current_na_00_", sp, "_prediction_mu_sd_matrix.rds", sep="")
        if(file.exists(f.mean) & file.exists(f.sd)){
            mu.mean <- readRDS(f.mean)
            mu.sd <- readRDS(f.sd)
            df.current <- preparePredictionResultsHelpingFunction(mu.mean, mu.sd)
            df.current$gcc <- "Current"
            df.current$rcp <- "-"
            df.current$year <- "2000"
            df.current$type <- sp
            df.current$gc <- gc
            df.res <- rbind(df.res, df.current)
        }

        for(rcp in c("26", "85")){
            print(paste("Preparing mean model ", rcp, " 2070 results for", gc))
            f.mean <- paste(hbm.output.dir,
                            gc, "_mean_", rcp, "_70_", sp, "_prediction_mu_mean_matrix.rds", sep="")
            f.sd <- paste(hbm.output.dir,
                          gc, "_mean_", rcp, "_70_", sp, "_prediction_mu_sd_matrix.rds", sep="")
            if(!file.exists(f.mean) | !file.exists(f.sd)){
                print(paste("File", f.mean, "or", f.sd, "do not exist, skipping ..."))
                next
            }
            mu.mean <- readRDS(f.mean)
            mu.sd <- readRDS(f.sd)
            df <- preparePredictionResultsHelpingFunction(mu.mean, mu.sd)
            df$gcc <- "Mean model"
            df$rcp <- rcp
            df$year <- "2070"
            df$type <- sp
            df$gc <- gc
            df.res <- rbind(df.res, df)
        }
    }
    df.res$gcc <- as.factor(df.res$gcc)
    df.res$rcp <- as.factor(df.res$rcp)
    levels(df.res$rcp) <- c("-", "RCP 2.6", "RCP 8.5")
    write_csv(df.res, f.out)
    cat("\n==>", f.out)
    return(df.res)
}

writePredictionResults <- function(){
    cat("\n", paste0(rep("-",80), collapse=""))
    cat("\n PREPARING PREDICTION RESULTS")
    cat("\n", paste0(rep("-",80), collapse=""))
    df.gc <- data.frame()
    for(gc in genetic.clusters){
        cat("\n  --> ", gc)
        df <- preparePredictionResults(gc, force.run)
        df$gc <- gc
        df.gc <- bind_rows(df.gc, df)
    }
    write_csv(df.gc, f.out)
    cat(paste("\n==>", f.out), "\n")
}

mergePredictionResults <- function(force.run=F){
    f.out <- paste(hbm.output.dir, "prediction_results.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols(), progress=F))

    df.gc <- data.frame()
    for(gc in genetic.clusters){
        df <- preparePredictionResults(gc, force.run)
        df$gc <- gc
        df.gc <- bind_rows(df.gc, df)
    }

    write_csv(df.gc, f.out)
    cat("\n==>", f.out)
    return(df.gc)
}

prepareLatentFieldResults <- function(gc, force.run=F){
    f.out <- paste(hbm.output.dir, tolower(gc), "_latent_field_results.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols()))

    f.mean <- paste(hbm.output.dir, tolower(gc), "_latent_field_avg.rds", sep="")
    f.sd <- paste(hbm.output.dir, tolower(gc), "_latent_field_sd.rds", sep="")

    if(!file.exists(f.mean) | !file.exists(f.sd)){
        print(paste("File", f.mean, "or", f.sd, "do not exist, skipping ..."))
        return()
    }

    mu.mean <- readRDS(f.mean)
    mu.sd <- readRDS(f.sd)
    df <- preparePredictionResultsHelpingFunction(mu.mean, mu.sd)
    df$gc <- gc

    write_csv(df, f.out)
    cat("\n==>", f.out)
    return(df)
}

mergeLatentFieldResults <- function(force.run=F){
    f.out <- paste(hbm.output.dir, "latent_field_results.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols(), progress=F))

    df.gc <- data.frame()
    for(gc in genetic.clusters){
        df <- prepareLatentFieldResults(gc, force.run)
        if(is.null(df))
            next
        df$gc <- gc
        df.gc <- bind_rows(df.gc, df)
    }

    write_csv(df.gc, f.out)
    cat("\n==>", f.out)
    return(df.gc)
}

getMaxentPrediction <- function(gc, th, gcc, rcp, yr, file.suffix){
    f <- paste0(maxent.output.dir, gc, "_", th, "_", gcc, "_", rcp, "_", yr, "_prediction_", file.suffix, ".csv", sep="")
    if(!file.exists(f))
        return()
    df <- read_csv(f, col_types=cols(), progress = FALSE)
    return(df)
}

getMaxentThresholdedCurrentAndFutureMeanPredictions <- function(genetic.clusters, thresholds, gccs,
                                                                rcps, yrs, file.suffix){
    df.res <- tibble()
    for(gc in genetic.clusters){
        for(gcc in gccs){
            for(rcp in rcps){
                for(yr in yrs){
                    for(th in thresholds){
                        df <- getMaxentPrediction(gc, th, gcc, rcp, yr, file.suffix)
                        if(is.null(df))
                            next
                        df$gc <- gc
                        df$gcc <- gcc
                        df$rcp <- rcp
                        df$yr <- yr
                        df$th <- th
                        if(!is.null(df))
                            df.res <- bind_rows(df.res, df)
                    }
                }
            }
        }
    }
    return(df.res)
}
# Maxent plus INLA (spatial and non-spatial)
getCombinedMeanResultsAll <- function(th, file.suffix, force.run=F){
    f.out <- paste(output.dir, "hbm_maxent-th-", th, "_combined-mean-results_all.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read_csv(f.out, col_types = cols(), progress=F))

    df.inla <- mergePredictionResults()
    df.inla <- df.inla %>% mutate(gcc=ifelse(gcc == 'Mean model', 'Mean', gcc), method='INLA')
    df.inla <- df.inla %>% select(x, y, var, v, gc, gcc, rcp, year, method, type)

    df.maxent <- getMaxentThresholdedCurrentAndFutureMeanPredictions(genetic.clusters, th,
                                                                     c("current", "mean"),
                                                                     c("na", "2.6", "8.5"),
                                                                     c("00", "70"), file.suffix)
    df.maxent <- df.maxent %>% mutate(gcc=ifelse(gcc == 'current', 'Current', gcc),
                                      gcc=ifelse(gcc == 'mean', 'Mean', gcc),
                                      rcp=ifelse(rcp == 'na', '-', rcp),
                                      rcp=ifelse(rcp == '2.6', 'RCP 2.6', rcp),
                                      rcp=ifelse(rcp == '8.5', 'RCP 8.5', rcp),
                                      var="mean",
                                      yr=as.integer(paste("20", yr, sep="")),
                                      method='Maxent',
                                      type = 'non-spatial')
    df.maxent <- df.maxent %>% rename(year=yr, v=pred) %>% select(x, y, var, v, gc, gcc, rcp, year, method, type)
    df.all <- bind_rows(df.inla, df.maxent)
    write_csv(df.all, f.out)
    cat("\n==>", f.out)
    return(df.all)
}

#### MODELS PREDICTIVE METRICS #####################################################################
# mae and rmse functions
mae <- function(o, p) {round(sum(abs(p - o)) / length(o),3)}
rmse <- function(o, p) {round(sqrt(mean((p - o)^2)),3)}

writeModelsPredictiveMetrics <- function(data){
    df.res.full <- data.frame()
    cat("\n", paste0(rep("-",80), collapse=""))
    cat("\n ASSESSING PREDICTIVE PERFORMANCE FOR MODELS WITH FULL DATASET")
    cat("\n", paste0(rep("-",80), collapse=""))
    for(gc in c("gc1", "gc2", "gc3", "gc4")){
        cat("\n  --> ", gc)
        f.non.sp.est.model <- paste0(hbm.output.dir, gc, "_non-spatial_fit.rds", sep="")
        non.sp.est.model <- readRDS(f.non.sp.est.model)
        predictors <- non.sp.est.model$names.fixed[str_detect(non.sp.est.model$names.fixed, "bio")]
        formula <- paste("resp ~ -1 + intercept +", paste(predictors, collapse=' + '))
        fitted.non.sp <- non.sp.est.model$summary.fitted.values$mean[1:301]
        mae.full.non.sp <- mae(data[[gc]], fitted.non.sp)
        rmse.full.non.sp <- rmse(data[[gc]], fitted.non.sp)

        df.res.full <- rbind(df.res.full, cbind(gc=gc, model='non-spatial',
                                                formula=formula, metric='mae',
                                                value=mae.full.non.sp))
        df.res.full <- rbind(df.res.full, cbind(gc=gc, model='non-spatial',
                                                formula=formula, metric='rmse',
                                                value=rmse.full.non.sp))

        if(gc != "gc3"){
            f.sp.est.model <- paste0(hbm.output.dir, gc, "_spatial_fit.rds", sep="")
            sp.est.model <- readRDS(f.sp.est.model)
            predictors <- sp.est.model$names.fixed[str_detect(sp.est.model$names.fixed, "bio")]
            formula <- paste("resp ~ -1 + intercept +", paste(predictors, collapse=' + '), "+ f(s, model=spde)")
            fitted.sp <- sp.est.model$summary.fitted.values$mean[1:301]
            mae.full.sp <- mae(data[[gc]], fitted.sp)
            rmse.full.sp <- rmse(data[[gc]], fitted.sp)

            df.res.full <- rbind(df.res.full, cbind(gc=gc, model='spatial',
                                                    formula=formula, metric='mae',
                                                    value=mae.full.sp))
            df.res.full <- rbind(df.res.full, cbind(gc=gc, model='spatial',
                                                    formula=formula, metric='rmse',
                                                    value=rmse.full.sp))
        }
    }
    f.out <- paste0(hbm.output.dir, "models_predictive_performance.csv")
    write.csv(df.res.full, f.out, quote = T, row.names = F)
    cat("\n==>", f.out)
}

#### SPATIAL AUTOCORRELATION ###################################################
# gets a data frame of x,y,gc.value,pred,obs where gc.value is the original gc
# membership value and pred is the predicted probability of presence given a
# threshold of 0.5
getMaxentPredictedObserved <- function(gc){
    f.in <- paste(maxent.output.dir, tolower(gc), "_0.5_best-maxent-model/species_samplePredictions.csv", sep= "")
    mx.res <- read.csv(f.in)
    po <- cbind(pred=mx.res[["Logistic.prediction"]], obs=1)

    f.in <- paste(maxent.output.dir, gc, "_0.5_modelling-presences-bg.csv", sep="")
    mod.data <- read.csv(f.in) %>% filter(pa == 1)

    po <- data.frame(cbind(mod.data[, c("x", "y")], po))
    return(po)
}

assessModelsRSAC <- function(force.run=F){
    f.out <- paste0(output.dir, "RSAC_evaluation.csv", sep="")
    if(file.exists(f.out) & !force.run)
        return(read.csv(f.out))

    # Maxent residual spatial autocorrelation
    df.res <- data.frame()
    for(gc in genetic.clusters){
        cat("\nMaxent - checking RSAC for", toupper(gc))
        po <- getMaxentPredictedObserved(gc)
        residuals <- 1 - po$pred
        rsac <- checkRSAC(po$x, po$y, residuals, rsac.nr.simulations)
        rsac$gc <- gc
        rsac$model.type <- "Maxent"
        df.res <- bind_rows(df.res, rsac)
    }

    # Non-spatial INLA residual spatial autocorrelation
    df <- read.csv(paste0(data.dir, "sp/ath_accessions.csv"))
    for(gc in toupper(genetic.clusters)){
        for(sp in c(F, T)){
            f.in <- paste0(hbm.output.dir, tolower(gc), "_", ifelse(sp, "spatial", "non-spatial"), "_fit.rds")
            if(!file.exists(f.in))
                next
            cat("\nHBM -", ifelse(sp, "spatial", "non-spatial"), "- Checking RSAC for", toupper(gc))
            model <- readRDS(f.in)
            fitted <- model$summary.fitted.values$mean[1:301]
            # fitted <- getFittedValues(gc, sp)
            if(is.null(fitted))
                next
            po <- data.frame(cbind(pred=fitted, obs=df[[tolower(gc)]]))
            residuals <- po$obs-po$pred
            rsac <- checkRSAC(df$x, df$y, residuals, rsac.nr.simulations)
            rsac$gc <- gc
            rsac$model.type <- "INLA"
            rsac$type = ifelse(sp, "spatial", "non-spatial")
            df.res <- bind_rows(df.res, rsac)
        }
    }

    write.csv(df.res, f.out, quote=F, row.names=F)
    cat("\n==>", f.out)
    return(df.res)
}

checkRSAC <- function(x, y, residuals, nr.simulations){
    dists <- as.matrix(dist(cbind(x, y)))
    dists.inv <- 1/dists
    diag(dists.inv) <- 0
    m <- as.matrix(dists.inv)
    mw <- mat2listw(m)
    set.seed(1234)
    sim1 <- moran.mc(residuals, listw=mw, nsim=nr.simulations)
    # print(sim1)
    moran.test <- "failed"
    if(round(sim1$p.value,2) >= 0.05){
        moran.test <- "passed"
    }
    # print(moran.test)
    df <- data.frame(cbind(statistic=sim1$statistic, parameter=sim1$parameter, p.value=sim1$p.value,
                           alternative=sim1$alternative, method=sim1$method, nr.simulations=length(sim1$res) - 1))
    return(df)
}

