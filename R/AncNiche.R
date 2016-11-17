AncNiche <- function(trees,
                     best.tree,
                     pca.wrapper,
                     predictors,
                     outgroup = F,
                     nodes,
                     trees.sample = 100,
                     trees.burnin = 20,
                     do.lambda = T,
                     do.models = F) {
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

    if (!length(predictors))
        stop('Predictors list is empty')
    if (!inherits(pca.wrapper, 'pca.wrapper'))
        stop('pca.wrapper should be an object of class \"pca.wrapper\".')
    if (!inherits(trees, 'multiPhylo'))
        stop('trees should be an object of class \"multiPhylo\".')
    if (!inherits(best.tree, 'phylo'))
        stop('best.tree should be an object of class \"phylo\".')

    result <- list()
    result$factors <- SelectPresvals(predictors, pca.wrapper)
    .bestTree <- ladderize(best.tree)
    .newtrees <-
        GetFraction(trees, burnin = trees.burnin, sample = trees.sample)
    if (is.list(outgroup) || outgroup != F) {
        # Process outgroup data - remove spaces and convert to combined characters
        if (inherits(outgroup, 'character') &&
            length(outgroup) == 1)
            outgroup <- strsplit(gsub(" ", "", outgroup), ",")[[1]]
        .bestTree <- DropTip(.bestTree, outgroup)
        .newtrees <- DropTip(.newtrees, outgroup)
    }
    result$bestTree <- .bestTree

    if (do.lambda) {
        # check signal

        for (i in 1:ncol(result$factors)) {
            message('Estimating phylogenetic signal for ',
                    colnames(result$factors)[i])
            factor <-
                structure(result$factors[, i], names = pca.wrapper$presvals[, 1])
            .temp.lambda <- MultiplePhylosig(.newtrees, factor)
            if (i > 1) {
                .lambda <- cbind(.lambda, .temp.lambda)
            }    else {
                .lambda <- .temp.lambda
            }
        }

        k <- 1

        for (i in 1:(ncol(result$factors) * 2)) {
            if (i %% 2 == 0)
                next
            colnames(.lambda)[i] <-
                paste(colnames(result$factors)[k], 'lambda')
            colnames(.lambda)[i + 1] <-
                paste(colnames(result$factors)[k], 'P')
            k <- k + 1
        }

        result$phylo.signal <- .lambda

        # check false signal
        .rem.trees <- vector()

        for (i in 1:(ncol(result$factors) * 2)) {
            if (i %% 2 == 0)
                .rem.trees = c(as.numeric(rownames(.lambda[.lambda[, i] > 0.05,])), .rem.trees)
        }
        .newtrees <- .newtrees[-unique(.rem.trees)]
    }

    if (!exists('nodes') || is.null(nodes))
        nodes <- 'root'
    result$nodes <- nodes

    .states <- list()

    if (do.models) {
        for (j in 1:ncol(result$factors)) {
            message(colnames(result$factors)[j], ' started')
            factor <- structure(result$factors[, j], names = pca.wrapper$presvals[, 1])
            .temp <- GetNodeAncSlow(.newtrees, factor, node = nodes)
            if (length(nodes) > 1) {
                # add states for selected nodes
                for (i in 1:length(nodes)) {
                    if (j > 1) {
                        .states[[i]] <- cbind(.states[[i]], .temp$states[[i]])
                    } else {
                        .states[[paste0('node', as.character(nodes[i]))]] <-
                            .temp$states[[i]]
                    }
                }
            } else {
                if (j > 1) {
                    .states <- cbind(.states, .temp$states)
                } else {
                    .states <- .temp$states
                }
            } # then append models information
            if (j > 1) {
                .models <- cbind(.models, .temp$selected)

            } else {
                .models <- .temp$selected
            }
        }
        if (length(nodes) > 1) {
            for (k in 1:length(nodes)) {
                colnames(.states[[k]]) <- colnames(factors)
            }
        } else {
            colnames(.states) <- colnames(factors)
        }
        result$models <- .models
    } else {
        .factors <- result$factors
        rownames(.factors) <- pca.wrapper$presvals[, 1]
        .states <- GetNodeAncFast(.newtrees, node = nodes, factors = .factors)
    }

    result$states <- .states
    class(result) <- 'anc.niche'
    return(result)

}