rcn.phyloclim <- function(trees,
                          best.tree,
                          pca.wrapper,
                          outgroup,
                          nodes,
                          trees.sample=100,
                          trees.burnin=20,
                          do.lambda=T,
                          do.models=F) {
#     in:
#         trees,
#         best tree,
#         RCN PCA wrapper object
#     params:
#         lambda y/n,
#         model est. y/n
#         trees sample fraction
#         outgroup tips (comma separated string)
#         internal nodes to estimate states (combined characters or numbers)
    if(!inherits(pca.wrapper,"RCN PCA wrapper"))
        stop("pca.wrapper should be an object of class \"RCN PCA wrapper\".")
    if(!inherits(trees,"multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    if(!inherits(best.tree,"phylo"))
        stop("best.tree should be an object of class \"phylo\".")
    .pca <- pca.wrapper$dudi.pca
    .pca2 <- pca.wrapper$prcomp
    .result.presvals <- as.data.frame(pca.wrapper$presvals)
    result <- list()

    # Process outgroup data - remove spaces and convert to combined characters
    if (inherits(outgroup, 'character') && length(outgroup) == 1)
        outgroup <- strsplit(gsub(" ", "", outgroup), ",")[[1]]

    # Process trees
    .bestTree <- ape::ladderize(best.tree)
    .bestTree <- DropTip(.bestTree, outgroup)
    .newtrees <- GetFraction(trees, burnin=trees.burnin, sample=trees.sample)
    .newtrees <- DropTip(.newtrees, outgroup)
    result$bestTree <- .bestTree

    # Select pca data with ID = tree tips:
    .pca.dim <- ncol(.pca$li)
    .pcaSelect <- match(.bestTree$tip.label, .result.presvals$ID)
    .phylo.pca <- EmptyDF(colnames=1:.pca.dim)
    for (i in 1:length(.pcaSelect)) {
        .phylo.pca[i,] <- .pca2$x[.pcaSelect[i], 1:.pca.dim]
    }
    rownames(.phylo.pca) <- .bestTree$tip.label
    # return .phylo.pca
    result$phylo.pca <- .phylo.pca

    # Check phylogenetic signal (Lambda) from selected PCA data
    if (do.lambda) {
        for (j in 1:.pca.dim) {
            factor <- .phylo.pca[,j]
            names(factor) <- rownames(.phylo.pca)
            .temp.lambda <- MultiplePhylosig(.newtrees, factor)
            if (j > 1) {
                .lambda <- cbind(.lambda, .temp.lambda)
            }    else {
                .lambda <- .temp.lambda
            }
        }
        k <- 1
        for (i in 1:(.pca.dim*2)) {
            if (i%%2 == 0) next
            colnames(.lambda)[i] <- paste0('PC', k, ' lambda')
            colnames(.lambda)[i+1] <- paste0('PC', k, ' P')
            k <- k+1
        }
        # return .lambda
        result$phylo.signal <- .lambda
    }

    # Estimate ancestral states
    if (is.null(nodes))
        nodes <- 'root'
    result$nodes <- nodes
    if (do.models) {
        for (j in 1:.pca.dim) {
            message('PC', j, ' started..')
            factor <- .phylo.pca[,j]
            names(factor) <- rownames(.phylo.pca)
            .temp <- GetNodeAncSlow(.newtrees, factor, node=nodes)
            if (j > 1) {
                .states <- cbind(.states, .temp$states)
                .models <- cbind(.models, .temp$selected)

            } else {
                .states <- .temp$states
                .models <- .temp$selected
            }
        }
        result$models <- .models
    } else {
        .states <- GetNodeAncFast(.newtrees, node=nodes, factors=.phylo.pca)
    }
    # return .states, .models
    result$states <- .states
    class(result) <- 'Phylogeny-based climatic niche'
    return(result)
}