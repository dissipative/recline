# Ancestors distribution - save to files
AncDistrib.Save <- function(anc.distrib, output.path){
# in: anc.distrib, output path

    if(!inherits(anc.distrib, 'anc.distrib'))
        stop('anc.distrib should be an object of class \"anc.distrib\".')

    if (is.null(output.path)) {
        output.path <- getwd()
    } else {
        output.path <- gsub('/$', '', output.path)
    }

	CheckDir(file.path(output.path, 'anc.distrib'))

    if (class(anc.distrib$distrib) == 'list') {
        ancestral_distrib <- raster::stack(anc.distrib$distrib)
    } else {
        ancestral_distrib <- anc.distrib$distrib
    }

    if (length(anc.distrib$nodes) > 1) { # multiple-node case
        for (i in 1:length(anc.distrib$points)) {


            # Save produced coordinates
            write.csv(anc.distrib$points[[i]],
                file=paste0(output.path, '/anc.distrib/coords.x', anc.distrib$boost[[i]],'.node', anc.distrib$nodes[i],'.csv'),
                quote=F
                )
        }

        # Save produced layers
        raster::writeRaster(
            ancestral_distrib,
            file.path(output.path, '/anc.distrib/distrib.tif'),
            format="GTiff", bylayer=T, overwrite=T,
            suffix=paste0('node', anc.distrib$nodes, 'x', anc.distrib$boost)
            )
    } else { # single-node case

        # Save produced coordinates
        write.csv(anc.distrib$points,
            file=paste0(output.path, '/anc.distrib/coords.x', anc.distrib$boost,'.node', anc.distrib$nodes,'.csv'),
            quote=F
            )

        # Save produced layer
        raster::writeRaster(ancestral_distrib,
            paste0(output.path, '/anc.distrib/distrib.x', anc.distrib$boost,'.node', anc.distrib$nodes,'.tif'),
            format="GTiff", overwrite=T)
    }
}