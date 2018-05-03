# internal functions (non-imported)
requireNamespace('ape')
requireNamespace('raster')
requireNamespace('phytools')
requireNamespace('phangorn')
requireNamespace('geiger')

CheckDir <- function(path) {
    # Check if directory exist
    if (!dir.exists(path)) {
        dir.create(path)
        message('Created directory \"', path, '\"')
    }
}

EmptyDF <- function(colnames = c('a', 'b', 'c')) {
    # Create empty data.frame with fixed column names
    newDF <-
        data.frame(matrix(NA, nrow = 0, ncol = length(colnames)), stringsAsFactors =
                       FALSE)
    colnames(newDF) <- colnames
    return(newDF)
}

#' @importFrom graphics par legend barplot abline
EvPlot <- function(ev) {
    # Check PC significance levels
    # use graphics
    # Args:
    #   ev: a vector of eigenvalues of PCA object
    # source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot")
    # Borcard, D., Gillet, F. & Legendre, P. 2011. Numerical Ecology with R.
    # Springer. Supplementary material
    # Broken stick model (MacArthur 1957)
    n <- length(ev)
    bsm <- data.frame(j = seq(1:n), p = 0)
    bsm$p[1] <- 1 / n
    for (i in 2:n)
        bsm$p[i] <- bsm$p[i - 1] + (1 / (n + 1 - i))
    bsm$p <- 100 * bsm$p / n
    # Plot eigenvalues and % of variation for each axis
    op <- par(mfrow = c(2, 1))
    barplot(ev,
            main = "Eigenvalues",
            col = "bisque",
            las = 2)
    abline(h = mean(ev), col = "red")
    legend(
        "topright",
        "Average eigenvalue",
        lwd = 1,
        col = 2,
        bty = "n"
    )
    barplot(
        t(cbind(100 * ev / sum(ev), bsm$p[n:1])),
        beside = TRUE,
        main = "% variation",
        col = c("bisque", 2),
        las = 2
    )
    legend(
        "topright",
        c("% eigenvalue", "Broken stick model"),
        pch = 15,
        col = c("bisque", 2),
        bty = "n"
    )
    par(op)
}

#' @export
#' @title Select predictor values
#'
#' @description Extract predictors values only for selected variables
#' after selection done by \code{\link{SuggestVariables}}.
#'
#' @param vars a character vector with names of selected variables.
#' @param predvals predictor values obtained with
#' \code{\link{GetPredictorsValues}}.
#'
#' @return A named numeric vector.
SelectPredVals <- function(vars, predvals) {
    # Extract predictors values only for selected variables from pca.wrapper
    # Args:
    #   vars: names of variables
    #   predvals: predictors values
    v <- vector()
    for (i in 1:length(vars))
        v[i] <- grep(vars[i], colnames(predvals))
    result <- predvals[, v]
    rownames(result) <- rownames(predvals)
    class(result) <- 'numeric'
    if (!is.null(as.data.frame(predvals)$ID))
        rownames(result) <- as.data.frame(predvals)$ID
    return (result)
}

#' @importFrom raster extract stack
#' @importFrom utils read.csv
#' @export
#' @title Get predictor values for each occurence point.
#'
#' @description Predictor values is the values of each geospatial layer
#' (BIOCLIM) for each occurence point of the species of interest.
#'
#' @param occurence.points a character string with the path to CSV file
#' containing occurence points. Necessary columns in the CSV file are:
#' Type, ID, Long, Lat, Clade.
#' @param bioclim.data a character string with the path to geospatial layers
#' (BIOCLIM or similar variables). The format of layers must be compatible with
#' \code{stack} and \code{extract} methods of \code{raster} package.
#' @param bioclim.ext a character string with the file extention of geospatial
#' layers given in bioclim.data argument. By default is "tif" (geotif files).
#'
#' @return list with two collections,
#' \code{$raw} (extracted BIOCLIM data for each occurence point) and
#' \code{$values} (including columns ID, Type, and Clade of occurence.points +
#' data from \code{$raw} collection).
GetPredictorsValues <-
    function(occurence.points,
             bioclim.data,
             bioclim.ext = 'tif') {
        locations <- read.csv(
            file = occurence.points,
            stringsAsFactors = F,
            header = T,
            sep = ",",
            dec = "."
        )

        spatial.layers <- list.files(
            path = bioclim.data,
            pattern = paste0('*.', bioclim.ext),
            full.names = T
        )

        result <- list()

        result$raw <- raster::extract(raster::stack(spatial.layers),
                                      data.frame(locations$Long, locations$Lat))
        result$values <- cbind(cbind(locations$ID,
                                     locations$Type),
                               cbind(locations$Clade,
                                     result$raw))
        colnames(result$values)[1:3] <- c('ID', 'Type', 'Clade')

        result$values <- result$values[complete.cases(result$values), ]
        if (nrow(result$values) == 0)
            stop('No predictor values for given coordinates!')

        return(result)
    }

GetLowestContribVar <- function(pca, pc_to_use = c(1, 2, 3)) {
    if (!inherits(pca, 'prcomp'))
        stop('pca must be the object of the class "prcomp"')

    sdev <- pca$sdev
    loadings <- abs(pca$rotation[, 1:length(pc_to_use)])

    pca.abs <-
        matrix(NA,
               ncol = length(pc_to_use),
               nrow = nrow(loadings))
    contrib_sum <- vector(length = nrow(pca.abs))

    for (i in 1:length(pc_to_use))
        pca.abs[, i] <- loadings[, i] * sdev[i]

    for (j in 1:nrow(pca.abs))
        contrib_sum[j] <- sum(pca.abs[j, ] ^ 2)

    names(contrib_sum) <- row.names(loadings)
    factor_name <- names(sort(contrib_sum)[1])

    return(factor_name)
}

#' Suggest predictors variables
#'
#' @description Suggest predictors variables for prediction of ancestral
#' niche & distibution by PCA and Broken Stick model. Suggested variables
#' contain the maximum information about studied species.
#'
#' @param predictors.values numeric matrix or data.frame of predictor values
#' for used occurence points.
#'
#' @return a character vector with the names of suggested variables.
#'
#' @export
#' @importFrom stats complete.cases prcomp
#'
#' @examples
#' predictors <- GetPredictorsValues(occurence.points, bioclim.present)
#' predictors.num <- predictors$values[1:3112,4:ncol(predictors$values)]
#' predictors.use <- SuggestVariables(predictors.num)
#' predictors.selected <- SelectPredVals(predictors.use, predictors$values)
SuggestVariables <- function(predictors.values) {

    message('Performing PCA...')

    pca.data <-
        predictors.values[complete.cases(predictors.values), ]
    pca.data <- apply(t(pca.data), 1, as.numeric)
    pca <- prcomp(pca.data, scale. = T)

    # get eigenvalues
    ev <- pca$sdev ^ 2

    # broken stick model:
    ev.n <- length(ev)
    bsm <- data.frame(j = seq(1:ev.n), p = 0)
    bsm$p[1] <- 1 / ev.n
    for (i in 2:ev.n)
        bsm$p[i] <- bsm$p[i - 1] + (1 / (ev.n  + 1 - i))
    bsm$p <- 100 * bsm$p / ev.n
    bsm <- cbind(100 * ev / sum(ev), bsm$p[ev.n:1])

    # choose the main PCs
    pc_to_use <- which(bsm[, 1] > bsm[, 2])
    main_loadings <- abs(pca$rotation[, pc_to_use])

    suggested_variables <- list()
    for (i in 1:ncol(main_loadings))
        suggested_variables[i] <-
        names(which(main_loadings[, i] == max(main_loadings[, i])))

    # additional variable!
    suggested_variables[ncol(main_loadings) + 1] <-
        GetLowestContribVar(pca, pc_to_use)

    suggested_variables <- unlist(suggested_variables)
    message('Suggested variables are: ',
            paste(suggested_variables, collapse = ', '))

    return(suggested_variables)
}

#' Get fraction of phylogenetic trees
#'
#' @description Perform burn-in and get random sapmle inside set of trees.
#'
#' @param trees a "multiPhylo" object with multiple trees.
#' @param burnin the percentage of starting trees to discard.
#' @param sample a number of trees to sample randomly after burn-in.
#'
#' @return a "multiPhylo" object with new trees selection.
#' @export
#'
#' @examples
#' # get random 100 trees after discarding first 20%
#' forest.sample <- GetFraction(forest)
GetFraction <- function(trees,
                        burnin = 20,
                        sample = 100) {

    # Perform burn-in and get random sapmle inside set of trees
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   burnin: the percentage of starting trees to discard
    #   sample: a number of trees to sample randomly after burn-in

    if (!inherits(trees, "multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    if (burnin > 0) {
        burnin <- burnin / 100
        newtrees <- trees[(burnin * length(trees)):length(trees)]
    } else {
        newtrees <- trees
    }
    if (sample != 0)
        newtrees <- sample(newtrees, size = sample)
    return(newtrees)
}

#' @importFrom ape drop.tip
DropTip <- function(phylo, tipsToDrop = c(1, 2)) {
    # Drop tip from tree or set of trees
    # Wrapper function for ape drop.tip
    # Args:
    #   phylo: a "phylo" or  "multiPhylo" object with phylogenetic tree or trees
    #   tipsToDrop: a vector of tree tips to drop from the tree
    if (!inherits(phylo, "phylo") && !inherits(phylo, "multiPhylo"))
        stop("phylo should be an object of class \"phylo\" or \"multiPhylo\".")
    if (class(phylo) == 'phylo')
        return(drop.tip(phylo, tipsToDrop))
    if (class(phylo) == 'multiPhylo') {
        newtrees <- lapply(phylo, drop.tip, tip = tipsToDrop)
        class(newtrees) <- "multiPhylo"
        return(newtrees)
    }
}

PhyloRescale <- function (trees, times = 1000) {
    # Function to rescale brach lengths on phylogenetic trees
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   times: a number, in which times trees branch lenghts should be scaled
    #   (e.g. 1000 for kya to mya scaling)
    if (!inherits(trees, "multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    newtrees <- list()
    for (x in 1:length(trees)) {
        trees[[x]]$edge.length <- trees[[x]]$edge.length * times
        newtrees[[x]] <- trees[[x]]
    }
    class(newtrees) <- "multiPhylo"
    return(newtrees)
}


#' Multiple phylogenetic signal test
#'
#' @param trees a "multiPhylo" object with multiple trees.
#' @param testdata a vector with continuous trait with names matching
#' trees tips
#'
#' @return data.frame.
#'
#' @importFrom phytools phylosig
#' @importFrom stats sd
MultiplePhylosig <- function (trees, testdata) {

    # Wrapper function for phylogenetic signal test (lambda estimation)
    # use phytools
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   testdata: a vector with continuous trait with names matching trees tips

    phylotest <- as.data.frame(t(
        sapply(
            trees,
            phytools::phylosig,
            x = testdata,
            test = T,
            method = "lambda"
        )
    ))
    phylotest <- as.data.frame(cbind(unlist(phylotest$lambda),
                                     unlist(phylotest$P)))
    colnames(phylotest) <- c('lambda', 'P')
    message("mean Lambda = ",
            round(mean(phylotest$lambda), 2),
            " +/- ", round(sd(phylotest$lambda), 2))
    return(phylotest)
}

#' @importFrom phytools fastAnc
#' @importFrom phangorn getRoot
GetNodeAncFast <- function(trees, node = 'root', factors) {
    # Estimation of factor values for a given node(s) in a set of trees
    # use phangorn, phytools
    # Args:
    #   trees: a "multiPhylo" object w. multiple trees
    #   node: internal node(s) for values estimation
    #   factors: a set of continuous factors in data.frame or matrix
    if (!inherits(trees, "multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    message('Estimating ancestor characters with BM model.')
    # remember if there are several nodes:
    nodes <- (length(node) > 1)
    # now work with each tree
    for (i in 1:length(trees)) {
        cat('.')
        # if default node is root
        if (length(node) == 1 && node == 'root')
            node <- getRoot(trees[[i]])
        # result class differs in case of node numbers
        if (nodes) {
            temp.anc <- list()
        } else {
            temp.anc <- vector()
        }
        # obtain ancestor values for each factor dimension
        for (j in 1:ncol(factors)) {
            factor <- structure(factors[, j], names = rownames(factors))
            factor <-
                factor[names(factor) %in% trees[[i]]$tip.label]
            temp.fastAnc <- fastAnc(trees[[i]], factor)
            if (nodes) {
                # collect characters for each selected node in a list
                for (k in 1:length(node)) {
                    temp.val <- temp.fastAnc[names(temp.fastAnc) == node[k]]
                    if (j > 1) {
                        temp.anc[[k]][j] <- temp.val
                    } else {
                        temp.anc[paste0('node', as.character(node[k]))] <-
                            temp.val
                    }
                }
            } else {
                # if selected the only node - combine all estimations in vector
                temp.fastAnc <-
                    temp.fastAnc[names(temp.fastAnc) == node]
                temp.anc[j]  <- temp.fastAnc
            }
        }
        if (i > 1) {
            if (nodes) {
                for (k in 1:length(node)) {
                    states[[k]] <- rbind(states[[k]], temp.anc[[k]])
                }
            } else {
                states <- rbind(states, temp.anc)
            }
        } else {
            if (nodes) {
                states <- list()

                for (k in 1:length(node)) {
                    states[[k]] <- temp.anc[[k]]
                    names(states)[k] <-
                        paste0('node', as.character(node[k]))
                }
            } else {
                states <- temp.anc
            }

        }
    }
    if (nodes) {
        for (k in 1:length(node)) {
            rownames(states[[k]]) <- c(1:length(trees))
            colnames(states[[k]]) <- colnames(factors)
        }
    } else {
        rownames(states) <- c(1:length(trees))
        colnames(states) <- colnames(factors)
    }
    message('Done')
    return(states)
}

#' @importFrom geiger fitContinuous rescale
#' @importFrom phangorn getRoot
#' @importFrom phytools fastAnc
GetNodeAncSlow <-
    function (trees,
              factor,
              node = 'root',
              ncores = NULL) {
        # Estimation of factor values for a given node(s) in a set of trees
        # with fitting different models of char evolution (BM, EB and OU)
        # use phangorn, phytools, geiger, parallel
        # Args:
        #   trees: a "multiPhylo" object w. multiple trees
        #   node: internal node(s) for values estimation
        #   factor: a set of contionuous factors in data.frame
        if (!inherits(trees, "multiPhylo"))
            stop("trees should be an object of class \"multiPhylo\".")
        message('Estimating ancestor characters with best-fit model')
        cont.model <- character()
        if (length(node) > 1) {
            cont.states <- list()
        } else {
            cont.states <- vector()
        }
        # rescale data if need
        if (min(factor) < 1 || max(factor) > 1) {
            abs.max <- max(abs(factor))
            factor <- factor / abs.max
        } else {
            abs.max <- 1
        }
        # work with each tree in set
        message('Fitting models and obtaining ancestral values.')
        for (i in 1:length(trees)) {
            this.factor <- factor[names(factor) %in% trees[[i]]$tip.label]
            # fit models
            bm <- fitContinuous(trees[[i]],
                                this.factor,
                                model = "BM",
                                ncores = ncores)
            ou <- fitContinuous(trees[[i]],
                                this.factor,
                                model = "OU",
                                ncores = ncores)
            eb <- fitContinuous(trees[[i]],
                                this.factor,
                                model = "EB",
                                ncores = ncores)
            # sort them by AICc
            aicc <-
                structure(c(bm$opt$aicc, ou$opt$aicc, eb$opt$aicc),
                          names = c("BM", "OU", "EB"))
            selected <- names(sort(aicc)[1])
            # if diff too small - change selected
            if (abs((sort(aicc)[1] - sort(aicc)[2])) < 4 &&
                names(sort(aicc)[2]) == 'BM' ||
                abs((sort(aicc)[1] - sort(aicc)[2])) < 4 &&
                abs((sort(aicc)[2] - sort(aicc)[3])) < 4)
                selected <- 'BM'
            cont.model[i] <- selected
            # message("Tree no.", i, " preferred model: ", selected);
            # message( 'BM:', round(aicc[1], 2), ' OU:', round(aicc[2], 2),
            # ' EB:', round(aicc[3], 2) );
            # now fit tree to model
            cat('.')
            if (selected == 'OU') {
                tree <- rescale(trees[[i]],
                                model = 'OU',
                                alpha = ou$opt$alpha)
            } else if (selected == 'EB') {
                tree <-
                    rescale(
                        trees[[i]],
                        model = 'EB',
                        a = eb$opt$a,
                        sigsq = bm$opt$sigsq
                    )
            } else {
                tree <- rescale(trees[[i]],
                                model = 'BM',
                                sigsq = bm$opt$sigsq)
            }
            if (length(node) == 1 && node == 'root')
                node <- phangorn::getRoot(tree)
            temp.fastAnc <- fastAnc(tree, this.factor) * abs.max
            if (length(node) > 1) {
                # collect characters for each selected node in a list
                for (k in 1:length(node)) {
                    temp.val <- temp.fastAnc[names(temp.fastAnc) == node[k]]
                    if (i > 1) {
                        cont.states[[k]][i] <- temp.val
                    } else {
                        cont.states[paste0('node', as.character(node[k]))] <-
                            temp.val
                    }
                }
            } else {
                # if selected the only node - combine all estimations in vector
                temp.fastAnc   <-
                    temp.fastAnc[names(temp.fastAnc) == node]
                cont.states[[i]] <- temp.fastAnc
            }
        }
        message('Done')
        result <- list()
        result$selected <- cont.model
        result$states <- cont.states
        return(result)
    }

GetDistrByPred <- function(stack,
                           predictors,
                           boost = .05,
                           checkIfStop = F) {

    # Get approximate distribution for each layer in stack by predictors values
    # use raster
    # Args:
    #   stack: a RasterStack object
    #   predictors: a data.frame with predictors values
    #   boost: percentage of each predictor value increment (+/-, %)
    #   checkIfStop: check for useful information in overlapped layers

    if (!inherits(stack, "RasterStack"))
        stop("stack should be an object of class \"RasterStack\".")
    layer.list <- list()
    message('Approximating distibution in each layer')
    for (i in 1:length(stack[1])) {
        cat('.')
        layerVal.list <- unique(round(predictors[, i]))
        for (j in 1:length(layerVal.list)) {
            layerVal <- layerVal.list[j]
            if (boost > 0) {
                if (layerVal > 0) {
                    layerVal.min <- layerVal - layerVal * boost
                    layerVal.max <- layerVal + layerVal * boost
                } else {
                    layerVal.min <- layerVal + layerVal * boost
                    layerVal.max <- layerVal - layerVal * boost
                }
                layer.temp.min <- stack[[i]] > layerVal.min
                layer.temp.max <- stack[[i]] < layerVal.max
                layer.temp <- layer.temp.min * layer.temp.max
                layer.temp[layer.temp > 0] <- 1
            } else {
                layer.temp <- (stack[[i]] == layerVal)
            }
            if (j > 1) {
                layer.merged <- layer.merged + layer.temp
            } else {
                layer.merged <- layer.temp
            }
        }
        layer.merged    <- layer.merged / max(values(layer.merged), na.rm = T)
        layer.list[[i]]   <- layer.merged

        if (checkIfStop) {
            layer.result <- OverlapLayers(layer.list)
            exit         <- all(is.na(values(layer.result)))
        } else {
            exit <- F
        }
        if (exit) {
            message('\nWeak signal')
            break
        }
    }
    if (!exit) message('\nDone')
    if (checkIfStop) return (layer.result)
    return (layer.list)
}

OverlapLayers <- function(list, stopIfNA = F) {
    # Overlaps list of Raster layers with same size to single layer
    # use raster
    # Args:
    #   list: list of Raster objects produced by GetDistrByPred function
    #   stopIfNA: stop if all values in resulting layer is NA
    for (i in 1:length(list)) {
        if (i > 1) {
            layer.result <- layer.result * list[[i]]
        } else {
            layer.result <- list[[i]]
        }
        if (all(is.na(values(layer.result))) && stopIfNA) {
            message('Breaking merge at ', i, ' layer')
            break
        }
    }
    layer.result <- layer.result / max(values(layer.result), na.rm = T)
    return(layer.result)
}