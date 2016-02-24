# internal functions (non-imported)

CheckDir <- function(path){
    # Check if directory exist
    if (!dir.exists(path)) {
        dir.create(path)
        message('Created directory \"', path, '\"')
    }
}

EmptyDF <- function(colnames=c('a', 'b', 'c')){
    # Create empty data.frame with fixed column names
    newDF <- data.frame(matrix(NA,nrow=0,ncol=length(colnames)), stringsAsFactors=FALSE)
    colnames(newDF) <- colnames
    return(newDF)
}

EvPlot <- function(ev) {
    # Check PC significance levels
    # use graphics
    # Args:
    #   ev: a vector of eigenvalues of PCA object
    # source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
    # Broken stick model (MacArthur 1957)
    n <- length(ev)
    bsm <- data.frame(j=seq(1:n), p=0)
    bsm$p[1] <- 1/n
    for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
    bsm$p <- 100*bsm$p/n
    # Plot eigenvalues and % of variation for each axis
    op <- par(mfrow=c(2,1))
    barplot(ev, main="Eigenvalues", col="bisque", las=2)
    abline(h=mean(ev), col="red")
    legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
    barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,
            main="% variation", col=c("bisque",2), las=2)
    legend("topright", c("% eigenvalue", "Broken stick model"),
           pch=15, col=c("bisque",2), bty="n")
    par(op)
}

GetFraction <- function(trees, burnin=20, sample=100) {
    # Perform burn-in and get random sapmle inside set of trees
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   burnin: the percentage of starting trees to discard
    #   sample: a number of trees to sample randomly after burn-in
    if(!inherits(trees,"multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    if (burnin > 0) {
        burnin <- burnin/100
        newtrees <- trees[(burnin*length(trees)):length(trees)]
    } else {
        newtrees <- trees
    }
    if (sample != 0) newtrees <- sample(newtrees, size=sample)
    return(newtrees)
}

DropTip <- function(phylo, tipsToDrop=c(1,2)) {
    # Drop tip from tree or set of trees
    # Wrapper function for ape drop.tip
    # Args:
    #   phylo: a "phylo" or  "multiPhylo" object with phylogenetic tree or trees
    #   tipsToDrop: a vector of tree tips to drop from the tree
    if(!inherits(phylo, "phylo") && !inherits(phylo, "multiPhylo"))
        stop("phylo should be an object of class \"phylo\" or \"multiPhylo\".")
    if (class(phylo) == 'phylo')
        return( drop.tip(phylo, tipsToDrop) )
    if (class(phylo) == 'multiPhylo'){
        newtrees <- lapply(phylo, drop.tip, tip=tipsToDrop)
        class(newtrees) <- "multiPhylo"
        return(newtrees)
    }
}

PhyloRescale <- function (trees, times=1000) {
    # Function to rescale brach lengths on phylogenetic trees
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   times: a number, in which times trees branch lenghts should be scaled
    #   (e.g. 1000 for kya to mya scaling)
    if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"multiPhylo\".")
    newtrees <- list()
    for ( x in 1:length(trees) ) {
        trees[[x]]$edge.length <- trees[[x]]$edge.length*times
        newtrees[[x]] <- trees[[x]]
    }
    class(newtrees) <- "multiPhylo"
    return(newtrees)
}

MultiplePhylosig <- function (trees, testdata) {
    # Wrapper function for phylogenetic signal test (lambda estimation)
    # use phytools
    # Args:
    #   trees: a "multiPhylo" object with multiple trees
    #   testdata: a vector with continuous trait with names matching trees tips
    phylotest <- as.data.frame(t(sapply(trees, phylosig, x=testdata, test=T, method="lambda")))
    phylotest <- as.data.frame(cbind(unlist(phylotest$lambda), unlist(phylotest$P)))
    colnames(phylotest) <- c('lambda', 'P')
    message("mean Lambda = ",
            round(mean(phylotest$lambda), 2),
            " +/- ", round(sd(phylotest$lambda), 2)
            )
    return(phylotest)
}

GetNodeAncFast <- function(trees, node='root', factors) {
    # Estimation of factor values for a given node(s) in a set of trees
    # use phangorn, phytools
    # Args:
    #   trees: a "multiPhylo" object w. multiple trees
    #   node: internal node(s) for values estimation
    #   factors: a set of contionuous factors in data.frame or matrix
    if(!inherits(trees,"multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    cat('Estimating ancestor characters with BM model.')
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
        for (j in 1:length(factors)) {
            factor <- structure(factors[,j], names=rownames(factors))
            factor <- factor[names(factor) %in% trees[[i]]$tip.label]
            temp.fastAnc <- fastAnc(trees[[i]], factor)
            if (nodes) {
                # collect characters for each selected node in a list
                for (k in 1:length(node)) {
                    temp.val <- temp.fastAnc[names(temp.fastAnc) == node[k]]
                    if (j > 1) {
                        temp.anc[[k]][j] <- temp.val
                    } else {
                        temp.anc[paste0('node', as.character(node[k]))] <- temp.val
                    }
                }
            } else {
                # if selected the only node - combine all estimations in vector
                temp.fastAnc <- temp.fastAnc[names(temp.fastAnc) == node]
                temp.anc[j]  <- temp.fastAnc
            }
        }
        if (i > 1) {
            if (nodes) {
                for (k in 1:length(node)) {
                    ancNiche[[k]] <- rbind(ancNiche[[k]], temp.anc[[k]])
                }
            } else {
                ancNiche <- rbind(ancNiche, temp.anc)
            }
        } else {
            if (nodes) {
                ancNiche <- list();
                for (k in 1:length(node)) {
                    ancNiche[[k]] <- temp.anc[[k]]
                    names(ancNiche)[k] <- paste0('node', as.character(node[k]))
                }
            } else {
                ancNiche <- temp.anc
            }

        }
    }
    if (nodes) {
        for (k in 1:length(node)){
            rownames(ancNiche[[k]]) <- c(1:length(trees))
        }
    } else {
        rownames(ancNiche) <- c(1:length(trees))
        colnames(ancNiche) <- c(1:length(factors))
    }
    cat('Done\n')
    return(ancNiche)
}

GetNodeAncSlow <- function (trees, factor, node='root', ncores=NULL) {
    # Estimation of factor values for a given node(s) in a set of trees
    # with fitting different models of char evolution (BM, EB and OU)
    # use phangorn, phytools, geiger, parallel
    # Args:
    #   trees: a "multiPhylo" object w. multiple trees
    #   node: internal node(s) for values estimation
    #   factor: a set of contionuous factors in data.frame
    if(!inherits(trees, "multiPhylo"))
        stop("trees should be an object of class \"multiPhylo\".")
    message('Estimating ancestor characters with best-fit model')
    cont.model <- character()
    if (length(node) > 1) {
        cont.states <- list()
    } else {
        cont.states <- vector()
    }
    # rescale tree if need
    if ( max(trees[[1]]$edge.length) < 1 )
        trees <- PhyloRescale(trees, times=100)
    # work with each tree in set
    cat('Fitting models and obtaining ancestral values.')
    for (i in 1:length(trees)) {
        this.factor <- factor[names(factor) %in% trees[[i]]$tip.label]
        # fit models
        bm <- fitContinuous(trees[[i]], this.factor, model="BM", ncores=ncores)
        ou <- fitContinuous(trees[[i]], this.factor, model="OU", ncores=ncores)
        eb <- fitContinuous(trees[[i]], this.factor, model="EB", ncores=ncores)
        # sort them by AICc
        aicc <- c(bm$opt$aicc, ou$opt$aicc, eb$opt$aicc)
        names(aicc) <- c("BM","OU","EB")
        selected <- names( sort(aicc)[1] )
        # if something wrong - change selected
        if ((sort(aicc)[1] - sort(aicc)[2]) < 4)
            selected <- names( sort(aicc)[2] )
        if ((selected == 'OU' && ou$opt$alpha == ou$bnd[1,2])
            && (selected == 'EB' && eb$opt$a == eb$bnd[1,2]))
            selected <- 'BM';
        if (selected == 'OU' && ou$opt$alpha == ou$bnd[1,2])
            selected <- names(sort(aicc[names(aicc)!='OU']))[1]
        if (selected == 'EB' && eb$opt$a == eb$bnd[1,2])
            selected <- names(sort(aicc[names(aicc)!='EB']))[1]
        cont.model[i] <- selected
        # message("Tree no.", i, "preferred model:", selected);
        # message( 'BM:', aicc[1], 'OU:', aicc[2], 'EB:', aicc[3] );
        # now fit tree to model
        cat('.')
        if (selected == 'OU') {
            tree <- rescale(trees[[i]], model='OU', alpha=ou$opt$alpha)
        } else if (selected == 'EB') {
            tree <- rescale(trees[[i]], model='EB', a=eb$opt$a, sigsq=bm$opt$sigsq)
        } else {
            tree <- rescale(trees[[i]], model='BM', sigsq=bm$opt$sigsq)
        }
        if (length(node) == 1 && node == 'root')
            node <- getRoot(tree)
        temp.fastAnc <- fastAnc(tree, this.factor)
        if (length(node) > 1) {
            # collect characters for each selected node in a list
            for (k in 1:length(node)) {
                temp.val <- temp.fastAnc[names(temp.fastAnc) == node[k]]
                if (i > 1) {
                    cont.states[[k]][i] <- temp.val
                } else {
                    cont.states[paste0('node', as.character(node[k]))] <- temp.val
                }
            }
        } else {
            # if selected the only node - combine all estimations in vector
            temp.fastAnc   <- temp.fastAnc[names(temp.fastAnc) == node]
            cont.states[[i]] <- temp.fastAnc
        }
    }
    cat('Done\n')
    result <- list()
    result$selected <- cont.model
    result$states <- cont.states
    return(result)
}

GetDistrByPred <- function(stack, predictors,
                           boost=.05, checkIfStop=F) {
    # Get approximate distribution for each layer in stack by predictors values
    # use raster
    # Args:
    #   stack: a RasterStack object
    #   predictors: a data.frame with predictors values
    #   boost: percentage of each predictor value extension (+/-, %)
    if(!inherits(stack, "RasterStack"))
        stop("stack should be an object of class \"RasterStack\".")
    layer.list <- list()
    cat('Approximating distibution in each layer.')
    for (i in 1:length(stack[1])) {
        cat('.')
        layerVal.list <- unique(round(predictors[,i]))
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
                layer.temp <- layer.temp.min*layer.temp.max
                layer.temp[layer.temp > 0] <- 1
            } else {
                layer.temp <- stack[[i]] == layerVal
            }
            if (j > 1) {
                layer.merged <- layer.merged+layer.temp
            } else {
                layer.merged <- layer.temp
            }
        }
        layer.merged <- layer.merged/max(values(layer.merged), na.rm=T)
        layer.list[i] <- layer.merged
        if (checkIfStop) {
            layer.result <- OverlapLayers(layer.list)
            exit <- all(is.na(values(layer.result)))
        } else {
            exit <- F
        }
        if ( exit ) {
            cat('weak signal\n')
            break
        }
    }

    if (!exit && checkIfStop)
        cat('done\n\n')
        return(layer.result)
    if (!exit)
        cat('done\n\n')
        return(layer.list)
}

OverlapLayers <- function(list, stopIfNA=F) {
    # Overlaps list of Raster layers with same size to ingle layer
    # use raster
    # Args:
    #   list: list of Raster objects produced by GetDistrByPred function
    #   stopIfNA: stop if all values in resulting layer is NA
    for (i in 1:length(list)) {
        if ( i > 1) {
            layer.result <- layer.result*list[[i]]
        } else {
            layer.result <- list[[i]]
        }
        if ( all(is.na(values(layer.result))) && stopIfNA ) {
            message('Breaking merge at ', i,' layer')
            break
        }
    }
    layer.result <- layer.result/max(values(layer.result), na.rm=T)
    return(layer.result)
}