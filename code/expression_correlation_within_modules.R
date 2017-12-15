### expression_correlation_within_modules.R

### Load variables
source("code/variables_definition.R")
library(RColorBrewer)

### Set parameters
community.file <- paste0( 'all_tissues_communities_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
genes.file <- paste0('all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
cluster.results <- paste0( "all_tissues_expr_cor_clusters_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
cluster.pdf <- paste0( "comp_condor_expr_cor_clusters_fdr", FDRcis, FDRtrans, "_", window, "MB.pdf")

### Load data
load(normalized.rnaseq)
load(paste0(cluster.dir, community.file))
load(paste0(cluster.dir, genes.file))
load(tissue.file)

### Compute gene expression correlation within each community
expr.cor <- list()
expression <- list()
for (f in 1:length(genes)){
    cat("Running tissue", names(genes)[f], "\n")
    expres <- read.delim(paste0("rna_dir/", names(genes)[f], "_norm.tsv"),
                         header=T, sep="\t", row.names=1)
    expression[[f]] <- expres
    expr.cor[[f]] <- list()
    for(com in 1:length(genes[[f]])){
        cat("Community", com, "\n")
        expres.tmp <- expres[genes[[f]][[com]],]
        expr.cor[[f]][[com]] <- cor(t(expres.tmp))
    }
}

save(expression,
     file=paste0(rna.dir, "all_tissues_expression_level.RData"))
save(expr.cor,
     file=paste0(cluster.dir, "all_tissues_correlation_communities.RData"))

cor.list <- lapply(expr.cor, function(x){
    lapply(x, function(y){
        as.numeric(y[upper.tri(y)])
    })})



###Plot gene expression correlation distribution

colors <- c(brewer.pal(n = 8, name = "Dark2"),  brewer.pal(n = (length(cor.list)-8), name = "Accent") )
pdf(paste0(general.path, figure.path, "communities_correlation_distribution_",
           names(cor.list)[x], "_", y, ".pdf"))
par( mar=c(4,4,1,1)+.1)
plot(c(-1, 1), c(0,10), type='n', ylab="density", xlab="Pearson r")
for(x in 1:length(cor.list)){
    print(as.character(Tissues[names(cor.list[x]),1]))
    for(y in 1:length(cor.list[[x]])){
        print(y)
        lines(density(cor.list[[x]][[y]]), col=colors[x])
    }
}
legend("topleft", legend=Tissues[names(cor.list),2],
               col=colors, lty=1, pch=-1, bty='n')
dev.off()

### Plot gene expression correlation median for each community as a function of community size
a <- unlist(lapply(cor.list, function(x){lapply(x, function(y){median(abs(y))})}))
b <- unlist(lapply(genes, function(x){lapply(x, length)}))

pdf(paste0(figure.dir, "correlation_vs_nbgenes_log.pdf"), width=8, height=11)
par(mfrow=c(5,4), mar=c(4,4,3,1)+.1)
plot(b, a, pch=20, log="x",  main="",
     ylab="Median pairwise |Pearson r|", xlab="Nb of genes")
dev.off()

