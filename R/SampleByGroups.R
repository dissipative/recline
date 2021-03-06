SampleByGroups <- function(phylo,
                           grouping,
                           replicates = 10,
                           stay.max = 2) {
    # Sample tips in tree randomly by grouping data.frame
    # Necessary columns in grouping: tips, group
    # Args:
    #   phylo: a "phylo" or  "multiPhylo" object with phylogenetic tree or trees
    #   grouping: data.frame with selected groups to sample from
    #	replicates: number of random samples
    #   stay.max: maximum number of nodes to stay inside of group
    if (!inherits(phylo, "phylo") &&
        !inherits(phylo, "multiPhylo"))
        stop("phylo should be an object of class \"phylo\" or \"multiPhylo\".")
    gs <- unique(grouping$group)
    for (j in 1:replicates) {
        message('Sampling nodes inside groups. Replicate ', j)
        replicate.drop <- vector()
        for (i in 1:length(gs)) {
            this.tips <- as.vector(grouping[grouping$group == gs[i],]$tips)
            if (length(this.tips) > stay.max) {
                # cut any group to hold only pair of sequences,
                # if less - stay unchanged
                if (stay.max > 1) {
                    for (k in 1:(stay.max - 1)) {
                        this.tips <- this.tips[this.tips != sample(this.tips, 1)]
                    }
                }
                this.drop <-
                    this.tips[this.tips != sample(this.tips, 1)]
                replicate.drop <- c(replicate.drop, this.drop)
            }
        }
        this.phylo <- DropTip(phylo, replicate.drop)
        if (j > 1) {
            simple.forest <- c(this.phylo, simple.forest)
        } else {
            simple.forest <- this.phylo
        }
    }
    if (replicates > 1)
        class(simple.forest) <- "multiPhylo"
    return(simple.forest)
}
