# Ancestors distribution - save to files

    if (is.null(output.path)) {
        output.path <- getwd()
    } else {
        output.path <- gsub('/$', '', output.path)
    }
    CheckDir(paste0(config.output.path, '/ancestors'))

rcn.saveCsv(.phylo.states,
			paste0(config.output.path, '/ancestors/pcaEstim.node', .phylo.ancNode, '.csv')
);

rcn.saveCsv(.anc.restored[,.paleo.use],
			paste0(config.output.path, '/ancestors/presvals.node', .phylo.ancNode,'.csv')
			);

	# Save produced coordinates
rcn.saveCsv(.paleo.pts, filename=paste0(config.output.path, '/ancestors/coords.x', config.predictors.boost,'.node', .phylo.ancNode,'.csv'))

# Save produced layer
writeRaster(.paleo.distrib,
			paste0(config.output.path, '/ancestors/distrib.x', config.predictors.boost,'.node', .phylo.ancNode,'.tif'),
			format="GTiff", overwrite=T)