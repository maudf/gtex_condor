### find_tissue_specific_qtls_gcc.R.R

###Load variables
load("variables_definition.R")

### Set parameters
snps.file <- paste0( 'all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
genes.file <- paste0( 'all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
edges.file <- paste0( 'all_tissues_edges_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
eqtl.file <- paste0( 'all_tissues_eqtls_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
go.results.modif <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_final.txt')

snps.mult.file <- paste0( 'snps_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
genes.mult.file <- paste0( 'genes_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
edges.mult.file <- paste0( 'edges_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
ontology.mult.file <- paste0( 'ontology_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')


### Load data
load(tissue.file)
load(paste0(cluster.dir, snps.file))
load(paste0(cluster.dir, genes.file))
load(paste0(cluster.dir, edges.file))
go.tab <- read.delim(paste0(cluster.dir, go.results.modif), stringsAsFactors=F)

### Extract list of SNPs, Genes, Edges and GO Terms in GCC
snps.tmp <- lapply(snps, unlist)
genes.tmp <- lapply(genes, unlist)
edges.tmp <- lapply(edges, unlist)
tis <- unlist(lapply(go.tab$Tissue, function(x, T){rownames(T)[T[,1]==x]}, T = Tissues))
go.tmp <- tapply(go.tab$GOID, tis, function(x){unique(x)})

### Compute multiplicity of SNPs, Genes, Edges and GO Terms in GCC
cat("running multiplicity for \n")
cat("\t snps\n")
s <- table(unlist(snps.tmp))
cat("\t genes\n")
g <- table(unlist(genes.tmp))
cat("\t edges\n")
e <- table(unlist(edges.tmp))
cat("\t ontology\n")
o <- table(unlist(go.tmp))

### Save results
save(s, file=paste0(cluster.dir. snps.mult.file))
save(g, file=paste0(cluster.dir. genes.mult.file))
save(e, file=paste0(cluster.dir. edges.mult.file))
save(o, file=paste0(cluster.dir. ontology.mult.file))
