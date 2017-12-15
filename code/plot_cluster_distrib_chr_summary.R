### Maud Fagny
### 2015/12/21
### plot_cluster_distrib_chr_summary.R
### Plot nb of community with SNPs and genes from at least 2 chromosomes 
### ________________________________________________________


### Load data
source("code/variables_definition.R")
library(data.table)

### Set variables
eqtls.file <- paste0("all_tissues_eqtls_fdr",  FDRcis, FDRtrans, "_", window, "MB.Rdata")
edges.file <- paste0('all_tissues_edges_fdr',  FDRcis, FDRtrans, "_", window, 'MB.Rdata')
chr.data.file <- paste0('summary_cluster_chr_fdr',  FDRcis, FDRtrans, "_", window, 'MB.Rdata')
chr.snps.data.txt <- paste0('summary_cluster_snps_chr_fdr',  FDRcis, FDRtrans, "_", window, 'MB.txt')
chr.genes.data.txt <- paste0('summary_cluster_genes_chr_fdr',  FDRcis, FDRtrans, "_", window, 'MB.txt')
chr.data.txt <- paste0('summary_cluster_chr_fdr',  FDRcis, FDRtrans, "_", window, 'MB.txt')
chr.data.pdf <- paste0('summary_chr_bycommunities_fdr',  FDRcis, FDRtrans, "_", window, 'MB.pdf')
sup2.chr.data.pdf <- paste0('summary_chr_sup2_bycommunities_fdr',  FDRcis, FDRtrans, "_", window, 'MB.pdf')

###Load data
load(paste0(eqtl.dir, eqtls.file))
load(paste0(tissue.file))
load(paste0(samples.file))

### Functions
### For one community, calculate number of chromosomes to which SNPs/Genes/all nodes map
count.chr <- function(edg, qtl){
    a <- unique(data.frame(qtl[edges %in% edg,], stringsAsFactors=F)$Chr)
    b <- unique(data.frame(qtl[edges %in% edg], stringsAsFactors=F)$chr)
    c("SNP"=length(a), "Genes"=length(b), "All"=length(unique(c(a, b))))
}

### For each tissue, for each community, calculate number of chromosomes to which SNPs/Genes/all nodes map
nb.chr <- list()
for(tissue in names(eqtl)){
    print(tissue)
    qtl <- eqtl[[`tissue`]]
    qtl$edges <- paste(qtl$RSID, qtl$genes, sep='_')
    qtl <- data.table(qtl)
    setkey(qtl, edges)
    nb.chr[[`tissue`]] <- matrix(unlist(lapply(edges[[`tissue`]], count.chr, qtl)), ncol=3, byrow=T)
}
save(nb.chr, file=paste0(cluster.dir, chr.data.file))


### Plot distribution of number of chromosomes by community

###Create table summarizing proportion of communities with element in 1, 2, ..., 22 chromosomes for one tissue.
tab.chr.snps <- data.frame( matrix(0, ncol=22, nrow=length(nb.chr)) )
tab.chr.genes <- data.frame( matrix(0, ncol=22, nrow=length(nb.chr)) )
tab.chr.all <- data.frame( matrix(0, ncol=22, nrow=length(nb.chr)) )
colnames(tab.chr.snps ) <- colnames(tab.chr.genes ) <- colnames(tab.chr.all ) <- as.character(1:22)
rownames(tab.chr.snps ) <- rownames(tab.chr.genes ) <- rownames(tab.chr.all ) <- Tissues[names(nb.chr),2]

for(i in 1:length(nb.chr)){
    tmp.snps <- table(nb.chr[[i]][,1])
    tmp.genes <- table(nb.chr[[i]][,2])
    tmp.all <- table(nb.chr[[i]][,3])
    tab.chr.snps[i,names(tmp.snps)] <- tmp.snps/sum(tmp.snps)
    tab.chr.genes[i,names(tmp.genes)] <- tmp.genes/sum(tmp.genes)
    tab.chr.all[i,names(tmp.all)] <- tmp.all/sum(tmp.all)
}
write.table(tab.chr.all, file=paste0(general.path, cluster.path, chr.data.txt),
            quote=F, sep="\t")
write.table(tab.chr.snps, file=paste0(general.path, cluster.path, chr.snps.data.txt),
            quote=F, sep="\t")
write.table(tab.chr.genes, file=paste0(general.path, cluster.path, chr.genes.data.txt),
            quote=F, sep="\t")

### Plot proportion of communities with SNPs and genes from more than 2 chromosomes

b<-cbind(1-tab.chr.snps[,1], 1-tab.chr.genes[,1], 1-tab.chr.all[,1])
rownames(b) <- rownames(tab.chr.snps)
colnames(b) <- c("SNPs", "Genes", "Both")

pdf(paste0(general.path, figure.path, sup2.chr.data.pdf, sep=""), width=8, height=6)
par(las=1, mar=c(4,5,3,1)+.1)
barplot(t(b[,3:1]), horiz=T, beside=T, cex.names=0.8,
        col=c( "green3", "blue","red"), xlab="Proportion of community with SNPs and genes in >2chr",
        xlim=0:1)
legend("right", legend=colnames(b), fill=c("red", "blue", "green3"), bty='n')

dev.off()
