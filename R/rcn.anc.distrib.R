# Ancestors distribution
rcn.anc.distrib <- function(anc.niche,
                            pca.wrapper,
                            bioclim.past.data,
                            bioclim.ext='tif',
                            start.boost=.05, boost.step=.1
                            ) {

# in:
#           Phylogeny-based climatic niche
#           RCN PCA wrapper object
# params:
#           path to bioclim variables (past conditions)
#           extension of gis layers
#           start.boost
# out: ancestral distribution object
    if(!inherits(anc.niche,"Phylogeny-based climatic niche"))
        stop("anc.niche should be an object of class \"Phylogeny-based climatic niche\".")
    if(!inherits(pca.wrapper,"RCN PCA wrapper"))
        stop("pca.wrapper should be an object of class \"RCN PCA wrapper\".")

    .pca2 <- pca.wrapper$prcomp
    .pca <- pca.wrapper$dudi.pca
    .pca.dim <- ncol(.pca$li)
    .phylo.states <- anc.niche$states
    .predictors.boost <- start.boost
    .nodes <- anc.niche$nodes
    result <- list()

    message('Restoring predictor values for ancestor nodes.')
    .pca.variables <- factoextra::facto_summarize(.pca2, 'var', axes=1:.pca.dim)
    .pca.eig <- factoextra::get_eigenvalue(.pca2)[1:.pca.dim, 1]
    .pca.theo_contrib <- sum(100/length(.pca.variables$contrib) * .pca.eig)
    .paleo.use <- which(.pca.variables$contrib > .pca.theo_contrib)

    if (class(.phylo.states) == 'list') {
        .anc.restored <- list()
        for (i in 1:length(.phylo.states)) {
            .anc.restored[[i]] <- t(t(.phylo.states[[i]] %*% t(.pca2$rotation[,1:.pca.dim])) * .pca2$scale + .pca2$center)
        }
    } else {
        .anc.restored <- t(t(.phylo.states %*% t(.pca2$rotation[,1:.pca.dim])) * .pca2$scale + .pca2$center)
    }
    result$anc.presvals <- .anc.restored


    # Load ancient bioclim data
    .paleo.geoFiles <- list.files(
        path=bioclim.past.data,
        pattern=paste0('*.', bioclim.ext), full.names=T
        )

    .paleo.stack <- raster::stack(.paleo.geoFiles[.paleo.use])

    # Search for places matching restored predictors
    if (class(.anc.restored) == 'list') {
        result$distrib <- list()
        result$boost <- list()
        result$points <- list()
        for (i in 1:length(.anc.restored)) {
            repeat {
                .paleo.distrib <- GetDistrByPred(.paleo.stack,
                                              predictors=.anc.restored[[i]][,.paleo.use],
                                              boost=.predictors.boost,
                                              checkIfStop=T)
#                 .paleo.distrib <- OverlapLayers(.paleo.list, stopIfNA=T)
                if ( is.null(.paleo.distrib) || all(is.na(raster::values(.paleo.distrib))) ) {
                    message('Empty result. Trying again.')
                    .predictors.boost <- .predictors.boost + boost.step
                    message('Now boosting presvals by ',.predictors.boost)
                } else {
                    break
                }
            }
            result$distrib[i] <- .paleo.distrib
            result$boost[i] <- .predictors.boost
            # Extract points
            .paleo.distrib[.paleo.distrib == 0] <- NA
            result$points[[i]] <- raster::rasterToPoints(.paleo.distrib)
            .predictrs.boost <- start.boost
        }
        names(result$distrib) <- .nodes
        names(result$points)  <- .nodes
        names(result$anc.presvals) <- .nodes
    } else {
        repeat {
            .paleo.distrib <- GetDistrByPred(.paleo.stack,
                                          predictors=.anc.restored[,.paleo.use],
                                          boost=.predictors.boost,
                                          checkIfStop=T)
#             .paleo.distrib <- OverlapLayers(.paleo.list)
            if ( is.null(.paleo.distrib) || all(is.na(values(.paleo.distrib))) ) {
                message('Empty result. Trying again.')
                .predictors.boost <- .predictors.boost + boost.step
            } else {
                break
            }
        }
        # Extract points
        result$distrib <- .paleo.distrib
        result$boost <- .predictors.boost
        .paleo.distrib[.paleo.distrib == 0] <- NA
        result$points <- raster::rasterToPoints(.paleo.distrib)
    }

    # result$nodes <- .nodes

    class(result) <- 'Ancestral distribution data'
    return(result)
}






