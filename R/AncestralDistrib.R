#' Ancestors distribution
#'
#' @param anc.niche object of class "AncestralNiche"
#' @param predictors predictors to use. A vector obtained with
#' \code{\link{SuggestVariables}}.
#' @param bioclim.past.data a character string with the path to geospatial
#' layers in other timeframe (BIOCLIM or similar variables). The format of
#' layers must be compatible with \code{stack} and \code{extract} methods
#' of \code{raster} package.
#' @param bioclim.ext a character string with the file extention of geospatial
#' layers given in bioclim.past.data argument. By default is "tif" (geotif files).
#' @param start.boost a number to begin doing factor boosting if there are no
#' compatible conditions in bioclim.past.data.
#' @param boost.step a number to sequential increase of boosting level
#'
#' @return an object of class AncestralDistrib
#'
#' @export
#' @importFrom raster stack values rasterToPoints
#'
#' @examples
AncestralDistrib <- function(anc.niche,
                       predictors,
                       bioclim.past.data,
                       bioclim.ext = 'tif',
                       start.boost = 0,
                       boost.step = .1) {
    # in:
    #           Phylogeny-based climatic niche
    #           RCN PCA wrapper object
    # params:
    #           path to bioclim variables (past conditions)
    #           extension of gis layers
    #           boost
    # out: ancestral distribution object
    if (!inherits(anc.niche, 'AncestalNiche'))
        stop('anc.niche should be an object of class \"AncestalNiche\".')

    result <- list()
    result$nodes <- anc.niche$nodes

    states <- anc.niche$states
    boost <- start.boost

    # Load ancient bioclim data
    geo_files <- list.files(
        path = bioclim.past.data,
        pattern = paste0('*.', bioclim.ext),
        full.names = T
    )

    layers_to_use <- vector()

    for (i in 1:length(predictors)) {
        layers_to_use[i] <- grep(predictors[i], geo_files)
    }

    geo_stack <- raster::stack(geo_files[layers_to_use])

    # Search for places matching restored predictors
    if (class(states) == 'list') {
        result$distrib <- list()
        result$boost <- list()
        result$points <- list()
        for (i in 1:length(states)) {
            message('Node = ', result$nodes[i], '. Starting boost = ', boost)
            distrib <-
                InternalBoost(geo_stack, states[[i]], boost, boost.step, TRUE)

            result$distrib[[i]] <- distrib$ad
            result$boost[i] <- distrib$boost

            # Extract points
            values(distrib$ad)[values(distrib$ad == 0)] <- NA
            result$points[[i]] <- rasterToPoints(distrib$ad)
            boost <- start.boost
        }
        names(result$distrib) <- result$nodes
        names(result$points)  <- result$nodes
    } else {
        distrib <- InternalBoost(geo_stack, states, boost, boost.step, TRUE)

        # Extract points
        result$distrib <- distrib$ad
        result$boost <- distrib$boost

        values(distrib$ad)[values(distrib$ad == 0)] <- NA
        result$points <- rasterToPoints(distrib$ad)
    }
    class(result) <- 'AncestralDistrib'
    return(result)
}

InternalBoost <- function(geo_stack,
                              states,
                              boost,
                              boost.step,
                              checkIfStop = T) {
    result <- list()

    repeat {
        ancestal_distrib <- GetDistrByPred(
            geo_stack,
            predictors = states,
            boost = boost,
            checkIfStop = T
        )
        if (is.null(ancestal_distrib) ||
            all(is.na(raster::values(ancestal_distrib)))) {
            boost <- boost + boost.step
            message('Boosting predictors values range by ', boost)
        } else {
            break
        }
    }

    result$ad <- ancestal_distrib
    result$boost <- boost

    return(result)
}