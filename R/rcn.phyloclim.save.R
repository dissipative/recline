rcn.plyloclim.save <- function(phyloclim, pca.wrapper, output.path){
# in: phylogeny based climatic niche, output path

    if(!inherits(phyloclim,"Phylogeny-based climatic niche"))
        stop("phyloclim should be an object of class \"phylogeny based climatic niche\".")
    if(!inherits(pca.wrapper,"RCN PCA wrapper"))
        stop("pca.wrapper should be an object of class \"RCN PCA wrapper\".")

    output.path <- gsub('/$', '', output.path)
    CheckDir(paste0(output.path, '/phylo.plots'))

    .phylo.pca <- phyloclim$phylo.pca
    .phylo.states <- phyloclim$states
    .phylo.ancNode <- phyloclim$nodes
    .pca2 <- pca.wrapper$prcomp
    .result.presvals <- as.data.frame(pca.wrapper$presvals)

    .pca.dim <- ncol(.phylo.pca)
    # Plot of continual traits on the best tree
    for (i in 1:.pca.dim) {
        factor <- .phylo.pca[,i]
        names(factor) <- rownames(.phylo.pca)
        svg(paste0(output.path, '/phylo.plots/contMapPC', i ,'.svg'), width=10, height=7)
        phytools::contMap(phyloclim$bestTree, factor, fsize=0.2, lwd=1)
        dev.off()
    }

    # Save models info
    if ('models' %in% names(phyloclim)) {
        CheckDir(paste0(output.path, '/phylo.data'))
        SaveCsv(as.data.frame(phyloclim$models), paste0(output.path, '/phylo.data/models.csv'), row.names=F)
    }

    # Save plot of PCA with ancestor space
    if (class(.phylo.states) == 'list') {
        states <- reshape2::melt(.phylo.states)
        states <- data.frame(states$value[states$Var2==1], states$value[states$Var2==2], states$L1[states$Var2==1])
        colnames(states) <- c('X', 'Y', 'Node')
        colours <- grDevices::rainbow(length(.phylo.states), v=.3)
        col <- match(states$Node, unique(states$Node))
        for (i in 1:length(colours)) {
            col[col == i] <- colours[i]
        }
    } else {
        states <- as.data.frame(.phylo.states[,1:2])
        colnames(states) <- c('X', 'Y')
        col <- rainbow(1, v=.3)
    }
    .plot <- factoextra::fviz_pca_ind(.pca2, geom = "point",
                         habillage=.result.presvals$Clade, addEllipses=F) +
        ggplot2::geom_point(data=states,
                                  mapping=ggplot2::aes_string(x=states$X, y=states$Y),
                                  color=col,
                                  ) +
        ggplot2::theme_minimal()
    ggplot2::ggsave(filename=(paste0(output.path, '/phylo.plots/ancPCA.svg')),
                plot=.plot, width=10, height=10)
}