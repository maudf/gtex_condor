#makeCondorMatrixPlot.R

### Load libraries
library(condor)
library(igraph)

### Load variables
load("code/variables_definition.R")

### Set parameters
community.file <- paste0("all_tissues_communities_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")

### Load data
load(paste0(cluster.dir, community.file))

### Plot networks as matrix
for(tissue.name in names(community.file)){
    plotfile = paste0(figure.dir, tissue.name,"_matrixPlot.pdf")
    cols = rep("dodgerblue",length(unique(communities[[`tissue.name`]]$red.memb$com)))
    
    pdf(plotfile)
    condor.plot.communities(condor.modularity,color_list=cols)
    dev.off()
}
