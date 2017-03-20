AncestralNiche <- function(trees,
                           best.tree,
                           predictors.values,
                           outgroup = F,
                           nodes,
                           trees.sample = 100,
                           trees.burnin = 20,
                           do.lambda = T,
                           do.models = F) {

    #     in:
    #         trees,
    #         best tree,
    #         predictors values in matrix
    #     params:
    #         lambda y/n,
    #         model est. y/n
    #         trees sample fraction
    #         outgroup tips (comma separated string)
    #         internal nodes to estimate states (combined characters or numbers)

    if (!inherits(trees, 'multiPhylo'))
        stop('trees should be an object of class \"multiPhylo\".')
    if (!inherits(best.tree, 'phylo'))
        stop('best.tree should be an object of class \"phylo\".')

    result <- list()
    bestTree <- ladderize(best.tree)
    newtrees <- GetFraction(trees, burnin = trees.burnin, sample = trees.sample)
    if (is.list(outgroup) || outgroup != F) {
        # Process outgroup data - remove spaces and convert to combined characters
        if (inherits(outgroup, 'character') &&
            length(outgroup) == 1)
            outgroup <- strsplit(gsub(" ", "", outgroup), ",")[[1]]
        bestTree <- DropTip(bestTree, outgroup)
        newtrees <- DropTip(newtrees, outgroup)
    }
    result$bestTree <- bestTree

    if (do.lambda) {

        # check signal
        lambda <- matrix(nrow = 100, ncol = 0)

        for (i in 1:ncol(predictors.values)) {
            message('Estimating phylogenetic signal for ',
                    colnames(predictors.values)[i])

            lambda <- cbind(lambda, MultiplePhylosig(newtrees, predictors.values[, i]))
        }

        k <- 1
        for (i in 1:(ncol(predictors.values) * 2)) {
            if (i %% 2 == 0)
                next
            colnames(lambda)[i] <-
                paste(colnames(predictors.values)[k], 'lambda')
            colnames(lambda)[i + 1] <-
                paste(colnames(predictors.values)[k], 'P')
            k <- k + 1
        }

        result$phylo.signal <- lambda

        # check false signal
        .rem.trees <- vector()

        for (i in 1:(ncol(predictors.values) * 2)) {
            if (i %% 2 == 0)
                .rem.trees = c(as.numeric(rownames(lambda[lambda[, i] > 0.05,])), .rem.trees)
        }
        newtrees <- newtrees[-unique(.rem.trees)]
    }

    if (!exists('nodes') || is.null(nodes))
        nodes <- 'root'
    result$nodes <- nodes

    states <- list()

    if (do.models) {
        for (j in 1:ncol(predictors.values)) {
            message(colnames(predictors.values)[j], ' started')
            ancestors <- GetNodeAncSlow(newtrees, predictors.values[, j], node = nodes)
            if (length(nodes) > 1) {
                # add states for selected nodes
                for (i in 1:length(nodes)) {
                    if (j > 1) {
                        states[[i]] <- cbind(states[[i]], ancestors$states[[i]])
                    } else {
                        states[[paste0('node', as.character(nodes[i]))]] <-
                            ancestors$states[[i]]
                    }
                }
            } else {
                if (j > 1) {
                    states <- cbind(states, ancestors$states)
                } else {
                    states <- ancestors$states
                }
            } # then append models information
            if (j > 1) {
                models <- cbind(models, ancestors$selected)

            } else {
                models <- ancestors$selected
            }
        }
        if (length(nodes) > 1) {
            for (k in 1:length(nodes)) {
                colnames(states[[k]]) <- colnames(factors)
            }
        } else {
            colnames(states) <- colnames(factors)
        }
        result$models <- models
    } else {
        states <- GetNodeAncFast(newtrees, node = nodes, factors = predictors.values)
    }

    result$states <- states
    class(result) <- 'AncestalNiche'
    return(result)

}