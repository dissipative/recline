# test.final
# remove after package be ready

source('R/Recline.R')
source('R/AncestralNiche.R')

source('R/AncNiche.Save.R')
source('R/AncestralDistrib.R')
source('R/AncDistrib.Save.R')
source('R/SampleByGroups.R')

require(ape)
require(ade4)
require(phytools)
require(geiger)
require(factoextra)
require(raster)
require(phangorn)

# get cropped geo-layers
bioclim.present <- '/Users/artem/Cloud/current/gis-data/MAP_LAYERS/BIOCLIM/Present/2.5arcmin/mnem'
bioclim.past <- '/Users/artem/Cloud/current/gis-data/MAP_LAYERS/BIOCLIM/Past/LIG-120k-2.5-mnem/'

occurence.points = '../racoon/data/input.clades.csv'

# load calibrated phylogeny
mnem.tree <- ape::read.nexus('../../phylo/02_beast/2015-12-26/pmn.tree')
mnem.forest <- ape::read.nexus('../../phylo/02_beast/2015-12-26/pmn_combined_2015-12-26.trees')

# get random 100 trees after discarding first 20%
mnem.forest.sample <- GetFraction(mnem.forest)

# then normalize sampling across the tree tips
grouping <- read.csv('../racoon/data/grouping.csv')
simple.tree <- SampleByGroups(mnem.tree, grouping, replicates = 1)
simple.forest <- SampleByGroups(mnem.forest.sample, grouping, replicates = 10)

# visualize
plotTree(simple.tree, fsize=.5)
plotTree(simple.forest, node.numbers=T, fsize=.5)

# input
outgroup <- c('ariadne01', 'ariadne02')
nodes <- c(38, 58, 59)

# get predictors from occurence points + geo-layer
predictors <- GetPredictorsValues(occurence.points, bioclim.present)
predictors.num <- predictors$values[1:3112,4:ncol(predictors$values)] # only numeric values
predictors.use <- SuggestVariables(predictors.num)
predictors.selected <- SelectPredVals(predictors.use, predictors$values)

# restore occurences and climatic conditions by variables
# use: bio2, bio11, bio12, bio15
simple.niche <- AncestralNiche(
    simple.forest,
    simple.tree,
    predictors.selected,
    outgroup,
    nodes,
    trees.sample = 100,
    trees.burnin = 20,
    do.lambda = T,
    do.models = T
)

simple.distrib <- AncestralDistrib(simple.niche, predictors.use, bioclim.past)

save.image()
