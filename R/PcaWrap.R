as.pca.wrapper <- function(result.presvals, pca, pca2, pca.expl, pca.loadings){
    # Create an object of class "RCN PCA wrapper"
    result <- list();
    result$presvals <- result.presvals
    result$dudi.pca <- pca
    result$prcomp <- pca2
    result$explain <- pca.expl
    result$loadings <- pca.loadings
    class(result) <- 'pca.wrapper'
    return(result)
}

PcaWrap <- function(occurence.points, bioclim.data, bioclim.ext='tif'){
#   in: paths to occurence data, bioclim variables, extension of gis layers
#   class name: RCN PCA wrapper
    # Data processing
    .list.locations <- read.csv(
        occurence.points,
        stringsAsFactors=F, header=T, sep=",", dec="."
    )
    .list.geoFiles <- list.files(
        path=bioclim.data,
        pattern=paste0('*.', bioclim.ext), full.names=T
    )
    message('Extracting predictors values.')
    .list.presvals <- raster::extract( #predictors values
        raster::stack(.list.geoFiles),
        data.frame(.list.locations$Long, .list.locations$Lat)
    )
    .result.presvals <- cbind(
        cbind(
            .list.locations$ID,
            .list.locations$Type
        ),
        cbind(
            .list.locations$Clade,
            .list.presvals
        )
    )
    colnames(.result.presvals)[1] <- "ID"
    colnames(.result.presvals)[2] <- "Type"
    colnames(.result.presvals)[3] <- "Clade"
    .result.presvals <- .result.presvals[complete.cases(.result.presvals),]
    if (nrow(.result.presvals) == 0) stop('No predictor values for given coordinates!')

    # PCA
    message('Performing PCA.')
    .list.types <- unique(as.data.frame(.result.presvals)$Type)
    .pca.data <- .list.presvals[complete.cases(.list.presvals),]
    .pca <- ade4::dudi.pca( .pca.data, nf=3, scannf=F, scale=T  )
    .pca.expl <- 100*.pca$eig/sum(.pca$eig)
    .pca.loadings <- .pca$c1

    # alternative pca for data restoration:
    .pca2 <- prcomp(.pca.data, scale.=T)
    result <- as.pca.wrapper(.result.presvals, .pca, .pca2, .pca.expl, .pca.loadings)
    message('Done.')
    return(result);
}