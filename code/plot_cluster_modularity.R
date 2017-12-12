### Maud Fagny
### 2015/12/21
### plot_cluster_modularity.R
### make modularity plots 
### _________________________________

### Set variables
load("variables_definition.R")
communities.file <- paste0('all_tissues_communities_fdr', FDRcis, FDRtrans, "_", window, "MB.Rdata")
modularity.file <- paste0("modularity.values_fdr",FDRcis, FDRtrans, "_", window, "MB.Rdata")
modularity.pdf <- paste0('alltissues_fdr',  FDRcis, FDRtrans, "_", window, 'MB_modularity.pdf')

### Load data

load(tissue.file)
load(cluster.dir, communities.file)

modularity <- unlist(lapply(communities, function(x){max(x$modularity)}))
save(modularity, file=paste0(cluster.dir, modularity.file))

### Plot modularity
mod <- sort(modularity, decreasing=T)
n <- as.character(Tissues[names(modularity)[order(modularity, decreasing=T)],2])
pdf(paste(figure.dir, modularity.pdf, sep=""),
    width=2, height=4)
par(las=1, mar=c(4, 3, 2, 1)+0.1, mfrow=c(1,4))

barplot(mod, horiz=T, main = 'Modularity',  xlab = 'Modularity',
        cex.names=0.8, names=n, xlim = c(0,1), col='darkorchid')
dev.off()
