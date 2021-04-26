###extract_snp_genes_edges.R

### Load libraries
library(condor)

### Load variables
source("code/variables_definition.R")
load(paste0(tissue.file))

community.file <- paste0("all_tissues_communities_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
snps.file <- paste0('all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, "MB.Rdata")
genes.file <- paste0('all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, "MB.Rdata")
edges.file <- paste0('all_tissues_edges_fdr', FDRcis, FDRtrans, "_", window, "MB.Rdata")

### Load communities
list.communityfiles <- list.files(cluster.dir, pattern="_modularity.Rdata")

for(f in list.communityfiles){
    load(paste(cluster.dir, f))
    tissue <- gsub(paste0("fdr", FDRcis, FDRtrans, "_", window, "MB_edges_modularity.Rdata"), "", f)
    communities[[`tissue`]] <- condor.modularity
}
save(communities, file=paste(cluster.dir, community.file))

###Extract snps
snps <- lapply(communities, function(x){tapply(as.character(x$red.memb$red.names), x$red.memb$com, function(y){y})})

###Extract genes
genes <- lapply(communities, function(x){tapply(as.character(x$blue.memb$blue.names), x$blue.memb$com, function(y){y})})

###Extract edges
edges <- lapply(communities, function(x){
    g <- x$blue.memb$com
    names(g) <- as.character(x$blue.memb$blue.names)
    s <- x$red.memb$com
    names(s) <- as.character(x$red.memb$red.names)
    e <- x$edges
    e$snp.com <- s[e$red]
    e$gen.com <- g[e$blue]
    f <- e[e$snp.com == e$gen.com, ]
    tapply(paste(f$red, f$blue, sep="_"), f$snp.com, function(x){x})
})


save(snps, file=paste0(cluster.dir, snps.file))
save(genes, file=paste0(cluster.dir, genes.file))
save(edges, file=paste0(cluster.dir, edges.file))
