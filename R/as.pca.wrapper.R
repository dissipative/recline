as.pca.wrapper <- function(result.presvals, pca, pca2, pca.expl, pca.loadings){
    # Create an object of class "RCN PCA wrapper"
    result <- list();
    result$presvals <- result.presvals
    result$dudi.pca <- pca
    result$prcomp <- pca2
    result$explain <- pca.expl
    result$loadings <- pca.loadings
    class(result) <- "RCN PCA wrapper"
    return(result)
}