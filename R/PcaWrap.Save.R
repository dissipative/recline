PcaWrap.Save <- function(pca.wrapper, output.path) {
    # in: output path, dudi.pca, prcomp, result.presvals, pca.loadings, pca.expl
    if(!inherits(pca.wrapper, 'pca.wrapper'))
        stop('pca.wrapper should be an object of class \"pca.wrapper\".')

    output.path <- gsub('/$', '', output.path)
    CheckDir(file.path(output.path, 'pca.plots'))
    CheckDir(file.path(output.path, 'pca.data'))

    .pca <- pca.wrapper$dudi.pca
    .pca2 <- pca.wrapper$prcomp
    .result.presvals <- as.data.frame(pca.wrapper$presvals)
    .pca.loadings <- pca.wrapper$loadings
    .pca.expl <-pca.wrapper$explain

    # Save graphs
    svg(file.path(output.path, 'pca.plots/pca.evaluation.svg'),
        width=10, height=10)
    EvPlot(.pca$eig)
    dev.off()

    .pca.dim <- ncol(.pca$li)
    .pca.summ <- factoextra::facto_summarize(.pca2, element="var", axes=1:.pca.dim)

    svg(file.path(output.path, 'pca.plots/ind.biplot.svg'), width=10, height=10)
    factoextra::fviz_pca_ind(.pca2, geom = "point",
                 habillage=.result.presvals$Clade, addEllipses=TRUE,
                 ellipse.level= 0.95) +
        ggplot2::scale_color_brewer(palette="Set1") +
        ggplot2::theme_minimal()
    dev.off()

    svg(file.path(output.path, 'pca.plots/var.biplot.svg'), width=10, height=10)
    factoextra::fviz_pca_var(.pca2, col.var="contrib") +
        ggplot2::scale_color_gradient2(low="white", mid="blue",
                              high="red", midpoint = min(.pca.summ$contrib)) +
        ggplot2::theme_minimal()
    dev.off()

    .pca.plotList <- list()
    for (i in 1:.pca.dim) .pca.plotList[[i]] <- factoextra::fviz_contrib(.pca, choice='var', axes=i)
    for (i in 1:.pca.dim) {
        svg(paste0(output.path, '/pca.plots/contrib.dim', i ,'.svg'), width=10, height=7)
        print(.pca.plotList[[i]])
        dev.off()
    }
    svg(file.path(output.path, 'pca.plots/contrib.all.svg'), width=10, height=7)
    print(factoextra::fviz_contrib(.pca2, choice='var', axes=1:.pca.dim))
    dev.off()

    # Save PCA data
    .result.pca <- cbind(.result.presvals[,1:3], .pca2$x[,1:.pca.dim])
    write.csv(.result.presvals,
                file.path(output.path, '/pca.data/input.presvals.csv'), quote=F, row.names=F)
    write.csv(.result.pca,
                file.path(output.path, '/pca.data/result.global.csv'), quote=F, row.names=F)

    .list.types <- unique(.result.presvals$Type)
    for (i in 1:length(.list.types)) {
        condition = ( .result.presvals$Type == .list.types[i] )
        filename = paste0(output.path, '/pca.data/result.', .list.types[i],'.csv')
        write.csv( .result.pca[condition,], filename, quote=F, row.names=F )
    }

    # PCA loadings table:
    if ( exists('.pca.loadings') ) {
        colnames(.pca.loadings) <- paste0('PC', 1:ncol(.pca.loadings))
        .result.pcaLoad <- rbind(.pca.loadings, .pca.expl[1:ncol(.pca.loadings)]
        )
        row.names(.result.pcaLoad)[nrow(.result.pcaLoad)] <- 'Var explainaition %'
        write.csv(.result.pcaLoad, file.path(output.path, 'pca.data/loadings.csv'), quote=F)
    }
}