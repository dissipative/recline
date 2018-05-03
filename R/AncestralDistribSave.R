#' Save ancestor distribution to CSV and GeoTiff files
#'
#' @param anc.distrib an object of "AncestralDistrib" class.
#' @param output.path a character string with path to save produced files.
#'
#' @return null
#' @export
#' @importFrom raster stack writeRaster
#' @importFrom utils write.csv
#' @examples
AncestralDistribSave <- function(anc.distrib, output.path){
# in: anc.distrib, output path

    if(!inherits(anc.distrib, 'AncestralDistrib'))
        stop('anc.distrib should be an object of class \"AncestralDistrib\".')

    if (is.null(output.path)) {
        output.path <- getwd()
    } else {
        output.path <- gsub('/$', '', output.path)
    }
    base.path = "AncestralDistrib"

	CheckDir(file.path(output.path, base.path))

    if (class(anc.distrib$distrib) == 'list') {
        ancestral_distrib <- raster::stack(anc.distrib$distrib)
    } else {
        ancestral_distrib <- anc.distrib$distrib
    }

    if (length(anc.distrib$nodes) > 1) {

        # multiple-node case
        for (i in 1:length(anc.distrib$points)) {

            # Save produced coordinates
            write.csv(anc.distrib$points[[i]],
                file.path(output.path, base.path,
                          paste0('coords.x',
                                 anc.distrib$boost[[i]],'.node',
                                 anc.distrib$nodes[i],'.csv')),
                quote = F
                )
        }

        # Save produced layers
        raster::writeRaster(
            ancestral_distrib,
            file.path(output.path, base.path, 'distrib.tif'),
            format = "GTiff", bylayer = T, overwrite = T,
            suffix = paste0('node', anc.distrib$nodes, 'x', anc.distrib$boost)
            )
    } else {
        # single-node case

        # Save produced coordinates
        write.csv(anc.distrib$points,
            file.path(output.path, base.path,
                      paste0('coords.x',
                             anc.distrib$boost,'.node',
                             anc.distrib$nodes,'.csv')),
            quote = F
            )

        # Save produced layer
        raster::writeRaster(
            ancestral_distrib,
            file.path(output.path, base.path,
                      paste0('distrib.x',
                             anc.distrib$boost, '.node',
                             anc.distrib$nodes, '.tif')),
            format = "GTiff", overwrite = T)
    }
}
