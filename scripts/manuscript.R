#### MANUSCRIPT TABLES #############################################################################
#--- Table 1 ---
table_1 <- function(){
    cat("\nGenerating manuscript table 1")
    df.inla <- read.csv(paste0(hbm.output.dir, "variable_coefficients.csv")) %>%
        mutate(metric="variable coefficient") %>%
        select(Cluster, Model, Predictor, metric, value=mean)
    df.maxent <- read.csv(paste0(maxent.output.dir, "variable_contributions.csv"))
    df <- rbind(df.maxent, df.inla) %>% filter(Predictor != 'intercept') %>%
        mutate(value=round(value,3)) %>%
        select(Cluster, Model, Predictor, value) %>%
        spread(Predictor, value) %>% select(Cluster, Model, variables) %>% replace(., is.na(.), "-") %>%
        mutate(Cluster=gsub('G', '', toupper(Cluster))) %>% arrange(Cluster, Model)
    f.out <- paste0(manuscript.dir, "tables/table_1.csv")
    write.csv(df, f.out, row.names = F, quote = F)
    cat("\n  ==>", f.out, "\n")
}

#--- Table 2 ---
table_2 <- function(){
    cat("\nGenerating manuscript table 2")
    df.metrics <- read.csv(paste0(hbm.output.dir, "models_predictive_performance.csv")) %>%
        spread(metric, value) %>%
        mutate(gc=gsub('G', '', toupper(gc)), mae=sprintf("%.3f", as.numeric(as.character(mae))),
               rmse=sprintf("%.3f", as.numeric(as.character(rmse)))) %>%
        rename(Cluster=gc, HBM=model, Model=formula, MAE=mae, RMSE=rmse)
    f.out <- paste0(manuscript.dir, "tables/table_2.csv")
    write.csv(df.metrics, f.out, row.names = F, quote = F)
    cat("\n  ==>", f.out, "\n")
}

#--- Table 3 ---
table_3 <- function(){
    cat("\nGenerating manuscript table 3")
    df <- read_csv(paste0(output.dir, "hbm_maxent-th-0.5_combined-mean-results_all.csv"),
                   col_types = cols(), progress = F)
    nr.cells <- df %>% filter(gc == 'gc1' & method == 'INLA' & type == 'non-spatial' &
                                  year == '2000' & var == 'mean') %>% nrow()

    df.res <- df %>% filter(var=='mean') %>% group_by(gc, type, method, rcp, year) %>%
        summarise(n=n(), sum=sum(v)) %>% group_by(gc, type, method) %>%
        mutate(cum.prob=round(nr.cells/n*sum,2),
               change=round((cum.prob/cum.prob[1] - 1) * 100, 2)) %>% ungroup() %>%
        mutate(gc=gsub('G','', toupper(gc)), GCC=ifelse(year==2000, year, rcp),
               type=ifelse(method == 'Maxent', '', type),
               method=ifelse(method == 'INLA', paste0(type, ' HBMs'), method)) %>%
        select(Cluster=gc, GCC, method, cum.prob, change) %>% arrange(Cluster, method, GCC)
    f.out <- paste0(manuscript.dir, "tables/table_3.csv")
    write.csv(df.res, f.out, row.names = F, quote = F)
    cat("\n  ==>", f.out, "\n")
}

#### MANUSCRIPT FIGURES ############################################################################
#--- Figure 1 ---
figure_map_genetic_cluster_data <- function(ge){
    df <- read.csv(paste0(data.dir, "sp/ath_accessions.csv"))
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    df <- df %>% select(x, y, gc1, gc2, gc3, gc4) %>% gather(gc, value, 3:6) %>%
        mutate(gc=toupper(gc)) %>% filter(gc == toupper(ge))
    p <- ggplot(df) + geom_path(data=ip.df, aes(x=x, y=y), size=1.5)
    p <- p + geom_point(aes(x=x, y=y, size=value), colour="black", fill="grey45", shape=21)
    p <- p + scale_size_continuous("",
                                   limits = c(0,1),
                                   range = c(0, 15))
    p <- p + coord_equal()
    p <- p + theme(panel.grid = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   axis.ticks = element_blank(), axis.line = element_blank(),
                   plot.margin = unit( c(0,0,0,0) , units = "cm" ),
                   panel.background = element_rect(fill = "white"),
                   strip.background=element_rect(fill="white"),
                   strip.text.x = element_blank(),
                   legend.text=element_text(family="Times New Roman", size=30),
                   legend.key=element_blank())
    p <- p + guides(fill=guide_colorbar(barwidth = unit(0.8, "cm"), barheight = unit(10, "cm")))
    return(p)
}

figure_map_maxent_occurrences_thresholded <- function(ge, threshold = 0.5){
    df <- read.csv(paste0(data.dir, "sp/ath_accessions.csv"))
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    df <- df %>% select(x, y, gc1, gc2, gc3, gc4) %>% gather(gc, value, 3:6) %>% filter(value >= 0.5) %>%
        mutate(gc=toupper(gc)) %>% filter(gc == toupper(ge))
    p <- ggplot(df) + geom_path(data=ip.df, aes(x=x, y=y), size = 1.5)
    p <- p + geom_point(aes(x=x, y=y), colour="black", fill="grey45", shape=21, size=10)
    p <- p + coord_equal()
    p <- p + theme(panel.grid = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   axis.ticks = element_blank(), axis.line = element_blank(),
                   plot.margin = unit( c(0,0,0,0) , units = "cm" ),
                   panel.background = element_rect(fill = "white"),
                   strip.background=element_rect(fill="white"),
                   strip.text.x = element_blank(),
                   legend.title = element_text(family="Times New Roman", size=11),
                   legend.text=element_text(family="Times New Roman", size=9),
                   legend.key=element_blank())
    return(p)
}

figure_1 <- function(){
    cat("\nGenerating manuscript figure 1")
    dir.create(paste0(manuscript.dir, "figures/figure_1"), showWarnings = F)
    for(gc in toupper(genetic.clusters)){
        # Maps, continuous
        f.out <- paste(manuscript.dir, "figures/figure_1/figure_1a_", gsub('G','', toupper(gc)), "_accessions-continuous.png", sep="")
        png(f.out, width=1200, height=936, units='px')
        p <- figure_map_genetic_cluster_data(gc)
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
        # Maps, binary
        f.out <- paste(manuscript.dir, "figures/figure_1/figure_1b_", gsub('G','', toupper(gc)), "_accessions-binary.png", sep="")
        png(f.out, width=1200, height=936, units='px')
        p <- figure_map_maxent_occurrences_thresholded(gc, 0.5)
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
    }
    cat("\n")
}

#--- Figure 2 ---
getMaxentCurrentPredictions <- function(gc, file.suffix){
    df.pred.mx <- data.frame()

    f <- paste(maxent.output.dir, gc, "_", th, "_current_na_00_prediction_", file.suffix, ".csv", sep="")

    if(!file.exists(f))
        next
    df <- read_csv(f, col_types = cols(), progress = FALSE)
    df$gc <- gc
    df$th <- th
    df.pred.mx <- rbind(df.pred.mx, df)

    return(df.pred.mx)
}

getDataForFigureCurrentPredictions <- function(){
    df <-
        # MAXENT
        bind_rows(getMaxentCurrentPredictions("gc1", "best-maxent-model"),
                  getMaxentCurrentPredictions("gc2", "best-maxent-model"),
                  getMaxentCurrentPredictions("gc3", "best-maxent-model"),
                  getMaxentCurrentPredictions("gc4", "best-maxent-model")) %>%
        filter(th==0.5) %>% rename(v=pred) %>% mutate(var='-', method="Maxent") %>% select(x, y, v, gc, method) %>%
        mutate(var="", gc=toupper(gc)) %>% select(x,y,var,v,gc,method) %>%
        # INLA NON SPATIAL
        bind_rows(
            mergePredictionResults() %>% filter(type == 'non-spatial', year == 2000) %>%
                mutate(method='INLA', gc=toupper(gc)) %>% select(x, y, var, v, gc, method)
        ) %>%
        # INLA SPATIAL
        bind_rows(
            mergePredictionResults() %>% filter(type == 'spatial', year == 2000) %>%
                bind_rows(mergePredictionResults() %>%
                              filter(gc == "gc1" & type == 'spatial' & year == 2000) %>%
                              mutate(gc = 'GC3', v = NA)) %>%
                mutate(method='spINLA', gc=toupper(gc)) %>% select(x, y, var, v, gc, method)
        ) %>% mutate(method=factor(method, c("Maxent", "INLA", "spINLA")))
    return(df)
}

figure_current_predictions <- function(df, ge, me, va, limits, with.legend=F){
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    dfcp <- df %>% filter(method == me & var == va & gc==toupper(ge))
    if(va != 'sd'){
        colour.gradient <- mean.gc.colour.gradient
    } else {
        colour.gradient <- sd.gc.colour.gradient
    }
    p <- ggplot(dfcp) + geom_raster(aes(x=x, y=y, fill=v)) + geom_path(data=ip.df, aes(x=x, y=y), size=0.3)

    p <- p + scale_fill_gradientn(NULL, colours = colour.gradient, limits=limits)
    p <- p + coord_equal()
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                   axis.title = element_blank(), panel.grid = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   legend.position = "none",
                   plot.margin = unit( c(0,0,0,0) , units = "cm" ),
                   strip.background=element_rect(fill="white"),
                   legend.text=element_text(family="Times New Roman", size=44))
    if(with.legend){
        p <- p + theme(legend.position = c(0.9,0.25))
        p <- p + guides(fill=guide_colorbar(barwidth = unit(0.75, "cm"), barheight = unit(13, "cm")))
    } else {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}

figure_2 <- function(with.legend=F){
    cat("\nGenerating manuscript figure 2")
    dir.create(paste0(manuscript.dir, "figures/figure_2"), showWarnings = F)
    df <- getDataForFigureCurrentPredictions()
    ## MEANS ------------------------------------------------------------------
    # Maxent
    for(gc in toupper(genetic.clusters)){
        f.out <- paste(manuscript.dir, "figures/figure_2/figure_2a_", gsub('G','', toupper(gc)), "-maxent.png", sep="")
        png(f.out, width=1200, height=936, units='px')
        p <- figure_current_predictions(df, gc, "Maxent", "", c(0,1), with.legend)
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
    }
    # INLA
    for(method in c("INLA", "spINLA")){
        for(gc in toupper(genetic.clusters)){
            f.out <- paste0(manuscript.dir, "figures/figure_2/figure_2a_", gsub('G','', toupper(gc)), "-",
                           ifelse(method=='INLA','non-spatial_HBM', 'spatial_HBM'), "-mean.png")
            png(f.out, width=1200, height=936, units='px')
            p <- figure_current_predictions(df, gc, method, "mean", c(0,1), with.legend)
            print(p)
            dev.off()
            cat("\n  ==>", f.out)
        }
    }
    ## SD ---------------------------------------------------------------------
    for(gc in toupper(genetic.clusters)){
        # INLA
        f.out <- paste0(manuscript.dir, "figures/figure_2/figure_2b_", gsub('G','', toupper(gc)),
                       "-non-spatial_HBM-sd.png")
        png(f.out, width=1200, height=936, units='px')
        p <- figure_current_predictions(df, gc, "INLA", "sd", c(0, 0.2), with.legend)
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
        # spINLA
        f.out <- paste0(manuscript.dir, "figures/figure_2/figure_2b_", gsub('G','', toupper(gc)),
                       "-spatial_HBM-sd.png")
        p <- figure_current_predictions(df, gc, "spINLA", "sd", c(0, 0.2), with.legend)
        png(f.out, width=1200, height=936, units='px')
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
    }
    cat("\n")
}

#--- Figure 3 ---
getDataForFigureCurrentPredictionsLatentField <- function(){
    df.inla.lf <- mergeLatentFieldResults() %>% filter(gc == 'gc1') %>% mutate(gc = 'gc3', v = NA) %>%
        bind_rows(mergeLatentFieldResults() %>% filter(gc != "gc3")) %>%
        mutate(method='spINLA', gc=toupper(gc)) %>%
        select(x, y, var, v, gc, method)
    return(df.inla.lf)
}

figure_inla_latent_field <- function(df, ge, va, with.legend = F){
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    dfcp <- df %>% filter(var == va & gc==toupper(ge))
    if(va != 'sd'){
        colour.gradient <- mean.sp.colour.gradient
    } else {
        colour.gradient <- sd.sp.colour.gradient
    }
    p <- ggplot(dfcp) + geom_raster(aes(x=x, y=y, fill=v))
    p <- p + geom_path(data=ip.df, aes(x=x, y=y), size=0.3)
    p <- p + scale_fill_gradientn(NULL, colours = colour.gradient,
                                  limits=c(floor(min(dfcp$v)*10)/10, ceiling(max(dfcp$v)*10)/10),
                                  breaks=c(floor(min(dfcp$v)*10)/10, ceiling(max(dfcp$v)*10)/10))
    p <- p + coord_equal()
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                   axis.title = element_blank(), panel.grid = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   plot.margin = unit( c(0,0,0,0) , units = "cm" ),
                   strip.background=element_rect(fill="white"),
                   strip.text = element_text(size = 12, face = 'bold'),
                   legend.text=element_text(family="Times New Roman", size=80))

    p <- p + guides(fill=guide_colorbar(barwidth = unit(15, "mm"), barheight = unit(150, "mm")))

    if(with.legend){
        p <- p + theme(legend.position = c(0.9,0.25))
        p <- p + guides(fill=guide_colorbar(barwidth = unit(15, "mm"), barheight = unit(150, "mm")))
    } else {
        p <- p + theme(legend.position = "none")
    }

    return(p)
}

figure_3 <- function(with.legend=F){
    cat("\nGenerating manuscript figure 3")
    dir.create(paste0(manuscript.dir, "figures/figure_3"), showWarnings = F)
    df <- getDataForFigureCurrentPredictionsLatentField()
    for(gc in toupper(genetic.clusters)){
        # spINLA latent field mean
        f.out <- paste(manuscript.dir, "figures/figure_3/figure_3_", gsub('G','', toupper(gc)), "-spatial_HBM-latent_field_mean.png", sep="")
        p <- figure_inla_latent_field(df, gc, "mean", with.legend)
        png(f.out, width=1200, height=936, units='px')
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
        # spINLA latent field sd
        f.out <- paste(manuscript.dir, "figures/figure_3/figure_3_", gsub('G','', toupper(gc)), "-spatial_HBM-latent_field_sd.png", sep="")
        p <- figure_inla_latent_field(df, gc, "sd", with.legend)
        png(f.out, width=1200, height=936, units='px')
        print(p)
        dev.off()
        cat("\n  ==>", f.out)
    }
    cat("\n")
}

#--- Figure 4 ---
figure_future_predictions_maps <- function(df, with.legend=F){
    ip.df <- read_csv(paste0(data.dir, "study_area/ip.csv"),
                      col_types = cols(x=col_double(), y=col_double()))
    p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=v))
    p <- p + geom_path(data=ip.df, aes(x=x, y=y), size=0.3)
    p <- p + scale_fill_gradientn(NULL, colours = mean.gc.colour.gradient, limits=c(0,1))
    p <- p + coord_equal()
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                   axis.title = element_blank(), panel.grid = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   legend.text=element_text(family="Times New Roman", size=16),
                   strip.background=element_rect(fill="white"),
                   strip.text = element_text(size = 12, face = 'bold'),
                   strip.text.x = element_blank())

    if(with.legend){
        p <- p + theme(legend.position = c(0.9,0.25))
        p <- p + guides(fill=guide_colorbar(barwidth = unit(0.75, "cm"), barheight = unit(13, "cm")))
    } else {
        p <- p + theme(legend.position = "none")
    }

    return(p)
}

figure_4 <- function(with.legend = F){
    cat("\nGenerating manuscript figure 4")
    dir.create(paste0(manuscript.dir, "figures/figure_4"), showWarnings = F)
    df <- getCombinedMeanResultsAll(0.5, "best-maxent-model")
    # GC3 spatial models are not valid, we set it to NA to still appear on figure but as grey map
    df[df$gc == 'gc3' & df$type == 'spatial', 'v'] <- NA
    dfsc <- df %>% mutate(gc=toupper(gc), method=factor(method, c("Maxent", "INLA"))) %>%
        mutate(rcp = case_when(.$rcp == '-' ~ '2000', .$rcp == 'RCP 2.6' ~ '2070 - RCP 2.6',
                               .$rcp == 'RCP 8.5' ~ '2070 - RCP 8.5'),
               method = case_when(.$method == 'INLA' & .$type == 'spatial' ~ 'spatial_HBM',
                                  .$method == 'INLA' & .$type != 'spatial' ~ 'non-spatial_HBM',
                                  .$method == 'Maxent' & .$type != 'spatial' ~ 'Maxent'))

    # Maxent ---------------------------------------------------------------------
    for(ge in toupper(genetic.clusters)){
        for(me in c("Maxent", "non-spatial_HBM", "spatial_HBM")){
            for(scenario in c("2000", "2070 - RCP 2.6", "2070 - RCP 8.5")){
                dfp <- dfsc %>% filter(gc == toupper(ge) & method == me & rcp == scenario & var == 'mean' & (v > 0 | is.na(v)))
                f.out <- paste(manuscript.dir, "figures/figure_4/figure_4_", gsub('G', '', toupper(ge)), "-", me, "-",
                               str_replace(str_replace_all(scenario, " ", ""), "\\.", ""), ".png", sep="")
                pf <- figure_future_predictions_maps(dfp, with.legend)
                png(f.out, width=1200, height=936, units="px")
                print(pf)
                dev.off()
                cat("\n  ==>", f.out)
            }
        }
    }
    cat("\n")
}

#### TABLES SUPPLEMENTARY INFORMATION ##############################################################
#--- Table S1 ---
table_S1 <- function(){
    cat("\nGenerating supplementary information table S1")
    bm <- selectBestMaxentModels(F)
    df <- bm %>% mutate(model=paste("Y ~ ", str_replace_all(toupper(model), "-", " + "), sep="")) %>%
        rename(Cluster=gc, N=n, Model=model, AUC=auc, SD=auc.sd, Best=best) %>%
        arrange(Cluster, desc(AUC)) %>%
        ungroup() %>% mutate(Cluster=toupper(Cluster)) %>%
        select(Cluster, N, Model, AUC, SD, Best)
    f.out <- paste0(manuscript.dir, "tables/table_S1.csv")
    write.csv(df, f.out, row.names=F, quote=T)
    cat("\n  ==>", f.out, "\n")
}

#--- Table S2 ---
table_S2 <- function(){
    cat("\nGenerating supplementary information table S2")
    df <- assessModelsRSAC(F) %>%
        select(gc, model.type, type, statistic, p.value) %>%
        mutate(moran.test=ifelse(round(as.numeric(p.value),5) >= 0.05, 'passed', 'failed')) %>%
        rename(GC=gc, Method=model.type, Spatial=type, Statistic=statistic, "P-value"=p.value, "Moran's I"=moran.test)
    f.out <- paste0(manuscript.dir, "tables/table_S2.csv")
    write.csv(df, f.out, row.names=F, quote=T)
    cat("\n  ==>", f.out, "\n")
}

#--- Table S3 ---
table_S3 <- function(){
    cat("\nGenerating supplementary information table S3")
    df.res <- data.frame()
    for(gc in genetic.clusters){
        df <- readRDS(paste(hbm.output.dir, gc, "_models.rds", sep=""))[[3]]
        dfsp <- getBestModels(gc, df, 5, T)
        df.res <- bind_rows(df.res, dfsp)
        dfnsp <- getBestModels(gc, df, 5, F)
        df.res <- bind_rows(df.res, dfnsp)
    }

    df.res <- df.res %>% rowwise() %>% mutate(gc=toupper(gc)) %>% rename("BEST"="best", "Cluster"="gc", "Model"="Models") %>%
        mutate(Type=ifelse(str_detect(Model, "spde"), "Spatial HBMs", "Non-spatial HBMs")) %>%
        select(Type, Cluster, Model, DIC, WAIC, LCPO, BEST) %>% arrange(Type)
    f.out <- paste0(manuscript.dir, "tables/table_S3.csv")
    write.csv(df.res, f.out, row.names=F, quote=T)
    cat("\n  ==>", f.out, "\n")
}

#--- Table S4 ---
table_S4 <- function(){
    cat("\nGenerating supplementary information table S4")
    for(gc in genetic.clusters){
        # A - Non-spatial fits
        ns.fit <- readRDS(paste0(hbm.output.dir, gc, "_non-spatial_fit.rds"))
        df.ns <- txtRound(round(ns.fit$summary.fixed, 4)[,1:6],4)
        f.out <- paste0(manuscript.dir, "tables/table_S4_A_", toupper(gsub('g','',gc)), ".csv")
        write.csv(df.ns, f.out, row.names=T, quote=T)
        cat("\n  ==>", f.out)
        # B Spatial fits
        if(gc != 'gc3'){
            sp.fit <- readRDS(paste0(hbm.output.dir, gc, "_spatial_fit.rds"))
            df.sp <- txtRound(round(sp.fit$summary.fixed, 4)[,1:6],4)
            f.out <- paste0(manuscript.dir, "tables/table_S4_B_", toupper(gsub('g','',gc)), ".csv")
            write.csv(df.sp, f.out, row.names=T, quote=T)
            cat("\n  ==>", f.out)
        }

    }
}

#--- Table S5 ---
table_S5 <- function(){
    cat("\nGenerating supplementary information table S5")
    non.sp.means <- c()
    sp.means <- c()
    for(gc in genetic.clusters){
        # Non-spatial
        ns.fit <- readRDS(paste0(hbm.output.dir, gc, "_non-spatial_fit.rds"))
        mean <- round(ns.fit$summary.hyperpar["precision parameter for the beta observations", 'mean'], 3)
        non.sp.means <- c(non.sp.means, mean)
        # Spatial
        if(gc == 'gc3'){
            mean <- NA
        } else {
            sp.fit <- readRDS(paste0(hbm.output.dir, gc, "_spatial_fit.rds"))
            mean <- round(sp.fit$summary.hyperpar["precision parameter for the beta observations", 'mean'], 3)
        }
        sp.means <- c(sp.means, mean)
    }
    df <- data.frame(rbind(non.sp.means, sp.means))
    names(df) <- toupper(gsub('g', '', genetic.clusters))
    row.names(df) <- c('Non-spatial HBMs', 'Spatial HBMs')
    f.out <- paste0(manuscript.dir, "tables/table_S5.csv")
    write.csv(df, f.out, row.names=T, quote=T)
    cat("\n  ==>", f.out)
}

#### FIGURES SUPPLEMENTARY INFORMATION #############################################################
#--- Figure S1 ---
figure_S1 <- function(){
    cat("\nGenerating supplementary information figure S1")
    dir.create(paste0(manuscript.dir, "figures/figure_S1"), showWarnings = F)
    mesh <- readRDS(paste0(hbm.output.dir, "mesh.rds"))
    f.out <- paste(manuscript.dir, "figures/figure_S1/figure_S1_mesh.png", sep="")
    png(f.out, width=3800, height=3800, units="px", res=600)
    plot(mesh, main="")
    dev.off()
    cat("\n  ==>", f.out, "\n")
}

#--- Figure S2 ---
figure_current_predictions_density_overlap <- function(df, ge, with.legend=F){
    dff <- df %>% mutate(gc = toupper(gc), method=ifelse(method == 'INLA', 'HBM', method)) %>%
        filter(var == 'mean' & year == 2000 & gc == toupper(ge))
    dff <- dff %>% mutate(method=ifelse(method == 'Maxent', method, paste0(type, " ", method))) %>%
        mutate(method=factor(method, c("Maxent", "non-spatial HBM", "spatial HBM")))
    qt <- dff %>% filter(v >= 0) %>% group_by(method) %>%
        summarise(q.0=quantile(v, 0),
                  q.5 = quantile(v, 0.5),
                  q.75 = quantile(v, 0.75),
                  q.90 = quantile(v, 0.90),
                  q.1 = quantile(v, 1)) %>%
        gather(key=quantile, value=value, c(q.0, q.5, q.75, q.90, q.1))
    qt$quantile = factor(qt$quantile, c("q.0", "q.5", "q.75", "q.90", "q.1"))

    method.colours <- c("Magenta4", "Dark Orange", "Saddle Brown")
    p <- ggplot(dff)
    p <- p + geom_density(aes(x=v, y=..scaled.., fill=method, colour=method), size=0.4, alpha=0.3)
    p <- p + geom_point(data=qt %>% filter(quantile %in% c("q.75")),
                        aes(x=value, y=-0.03, colour=method, fill=method), size=3, shape=24)
    p <- p + scale_fill_manual("", values=method.colours)
    p <- p + scale_colour_manual("", values=method.colours)

    p <- p + scale_x_continuous("", limits = c(0,1))
    p <- p + scale_y_continuous("")
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                   axis.title = element_blank(), panel.grid = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.5, size=18),
                   panel.background = element_rect(fill = "white"),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   plot.margin = unit(c(1, 1, 0.5, 1), "lines"),
                   axis.line = element_line(),
                   strip.background=element_rect(fill="white"))

    if(with.legend){
        p <- p + theme(legend.position = c(0.8,0.8))
    } else {
        p <- p + theme(legend.position = "none")
    }

    return(p)
}

figure_S2 <- function(with.legend = F){
    cat("\nGenerating supplementary information figure S2")
    dir.create(paste0(manuscript.dir, "figures/figure_S2"), showWarnings = F)
    # Spatial GC3 models are not valid, we filter them out
    df <- getCombinedMeanResultsAll(0.5, "best-maxent-model") %>% filter(!(gc == 'gc3' & type == 'spatial'))
    for(gc in toupper(genetic.clusters)){
        f.out <- paste(manuscript.dir, "figures/figure_S2/figure_S2_", gsub('G', '', toupper(gc)), "_current_density_distributions.png", sep="")
        p <- figure_current_predictions_density_overlap(df, gc, with.legend)
        png(f.out, width=950, height=1050, units="px", res=300)
        plot(p)
        dev.off()
        cat("\n ==>", f.out)
    }
    cat("\n")
}

#--- Figure S3 ---
figure_future_predictions_density_rcp_overlap <- function(dfp, ge, me, with.legend=F){
    dff <- dfp %>% filter(var == 'mean' & gc == toupper(ge) & method == me)

    qt <- dff %>% filter(v >= 0) %>% group_by(method, rcp) %>%
        summarise(q.0=quantile(v, 0),
                  q.5 = quantile(v, 0.5),
                  q.75 = quantile(v, 0.75),
                  q.90 = quantile(v, 0.90),
                  q.1 = quantile(v, 1)) %>%
        gather(key=quantile, value=value, c(q.0, q.5, q.75, q.90, q.1))
    qt$quantile = factor(qt$quantile, c("q.0", "q.5", "q.75", "q.90", "q.1"))

    rcp.colours <- c("#00BA38", "#619CFF", "#F8766D")
    p <- ggplot(dff)
    p <- p + geom_density(aes(x=v, y=..scaled.., fill=rcp, colour=rcp), size=0.4, alpha=0.3)
    p <- p + geom_point(data=qt %>% filter(quantile %in% c("q.75")),
                        aes(x=value, y=-0.03, colour=rcp, fill=rcp), size=3, shape=24)
    p <- p + scale_fill_manual("", values=rcp.colours)
    p <- p + scale_colour_manual("", values=rcp.colours)
    p <- p + scale_x_continuous("", limits = c(0,1))
    p <- p + scale_y_continuous("")
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                   axis.title = element_blank(), panel.grid = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.5, size=18),
                   panel.background = element_rect(fill = "white"),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   plot.margin = unit(c(1, 1, 0.5, 1), "lines"),
                   axis.line = element_line(),
                   strip.background=element_rect(fill="white"))

    if(with.legend){
        p <- p + theme(legend.position = c(0.8,0.8))
    } else {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}

figure_S3 <- function(with.legend=F){
    cat("\nGenerating supplementary information figure S3")
    dir.create(paste0(manuscript.dir, "figures/figure_S3"), showWarnings = F)
    # GC3 spatial models are not valid, we filter them out
    df <- getCombinedMeanResultsAll(0.5, "best-maxent-model") %>%
        filter(!(gc== 'gc3' & type == 'spatial'))
    dfp <- df %>% mutate(gc=toupper(gc), method=ifelse(method=='INLA', 'HBM', method)) %>%
        mutate(method=factor(method, c("Maxent", "HBM"))) %>%
        mutate(rcp = case_when(.$rcp == '-' ~ '2000', .$rcp == 'RCP 2.6' ~ '2070 - RCP 2.6',
                               .$rcp == 'RCP 8.5' ~ '2070 - RCP 8.5'),
               method = case_when(.$method == 'HBM' & .$type == 'spatial' ~ 'spatial_HBM',
                                  .$method == 'HBM' & .$type != 'spatial' ~ 'non-spatial_HBM',
                                  .$method == 'Maxent' & .$type != 'spatial' ~ 'maxent'))

    for(ge in toupper(genetic.clusters)){
        for(me in c("maxent", "non-spatial_HBM", "spatial_HBM")){
            if(ge == 'GC3' & me == 'spatial_HBM')
                next
            f.out <- paste0(manuscript.dir, "figures/figure_S3/figure_S3_",
                            gsub('G', '', toupper(ge)), "-", me, ".png")

            png(f.out, width=950, height=1050, units="px", res=300)
            p <- figure_future_predictions_density_rcp_overlap(dfp, ge, me, with.legend)
            print(p)
            dev.off()
            cat("\n  ==>", f.out)
        }
    }
    cat("\n")
}

#### MAIN ##########################################################################################
manuscriptTables <- function(){
    table_1()
    table_2()
    table_3()
}

manuscriptFigures <- function(){
    figure_1()
    figure_2(T)
    figure_3(T)
    figure_4(T)
}

suppInformationTables <- function(){
    table_S1()
    table_S2()
    table_S3()
    table_S4()
    table_S5()
}

suppInformationFigures <- function(){
    figure_S1()
    figure_S2(T)
    figure_S3(T)
}
