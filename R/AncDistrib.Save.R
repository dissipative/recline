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
    .anc.restored <- anc.distrib$anc.presvals
    .points <- anc.distrib$points
    if (class(anc.distrib$distrib) == 'list') {
        .paleo.distrib <- raster::stack(anc.distrib$distrib)
    } else {
        .paleo.distrib <- anc.distrib$distrib
    }
    .phylo.ancNode <- anc.distrib$nodes
    .boost <- anc.distrib$boost

    if (length(.phylo.ancNode) > 1) { # multiple-node case
        for (i in 1:length(.anc.restored)) {
            write.csv(.anc.restored[[i]],
                file=paste0(output.path, '/presvals.node', .phylo.ancNode[i],'.csv'),
                quote=F
                )

            # Save produced coordinates
            write.csv(.points[[i]],
                file=paste0(output.path, '/coords.x', .boost[[i]],'.node', .phylo.ancNode[i],'.csv'),
                quote=F
                )
        }

        # Save produced layers
        raster::writeRaster(
            .paleo.distrib,
            file.path(output.path, '/distrib.tif'),
            format="GTiff", bylayer=T, overwrite=T,
            suffix=paste0('node', .phylo.ancNode, 'x', .boost)
            )
    } else { # single-node case
        write.csv(.anc.restored,
            file=paste0(output.path, '/presvals.node', .phylo.ancNode,'.csv'),
            quote=F
            )

        # Save produced coordinates
        write.csv(.points,
            file=paste0(output.path, '/coords.x', .boost,'.node', .phylo.ancNode,'.csv'),
            quote=F
            )

        # Save produced layer
        raster::writeRaster(.paleo.distrib,
            paste0(output.path, '/distrib.x', .boost,'.node', .phylo.ancNode,'.tif'),
            format="GTiff", overwrite=T)
    }
}