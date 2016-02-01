AncNiche.Save <- function(anc.niche, pca.wrapper, output.path){
# in: phylogeny based climatic niche, output path

    if(!inherits(anc.niche, 'anc.niche'))
        stop("anc.niche should be an object of class \"anc.niche\".")
    if(!inherits(pca.wrapper,'pca.wrapper'))
        stop("pca.wrapper should be an object of class \"pca.wrapper\".")

    output.path <- gsub('/$', '', output.path)
    CheckDir(file.path(output.path, 'phylo.plots'))
    CheckDir(file.path(output.path, 'phylo.data'))

    .phylo.pca <- anc.niche$phylo.pca
    .phylo.states <- anc.niche$states
    .phylo.ancNode <- anc.niche$nodes
    .pca2 <- pca.wrapper$prcomp
    .result.presvals <- as.data.frame(pca.wrapper$presvals)

    .pca.dim <- ncol(.phylo.pca)
    # Plot of continual traits on the best tree
    for (i in 1:.pca.dim) {
        factor <- .phylo.pca[,i]
        names(factor) <- rownames(.phylo.pca)
        svg(paste0(output.path, '/phylo.plots/contMapPC', i ,'.svg'), width=10, height=7)
        phytools::contMap(anc.niche$bestTree, factor, fsize=0.2, lwd=1)
        dev.off()
    }

    # Save models info
    if ('models' %in% names(anc.niche)) {
        write.csv(anc.niche$models, file.path(output.path, 'phylo.data/models.csv'), quote=F, row.names=F)
    }
    # # Save states
    if (class(.phylo.states) == 'list') {
        for (i in 1:length(.phylo.states)) {
            write.csv(.phylo.states[i],
                file.path(output.path, paste0('phylo.data/pcaEstim.node', .phylo.ancNode[i], '.csv')),
                quote=F
                )
        }
    } else {
        write.csv(.phylo.states,
                file.path(output.path, paste0('phylo.data/pcaEstim.node', .phylo.ancNode, '.csv')),
                quote=F
                )
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
    ggplot2::ggsave(filename=(file.path(output.path, 'phylo.plots/ancPCA.svg')),
                plot=.plot, width=10, height=10)
}