###identify_snps_neighbors.R

### Load variables
source("code/variables_definition.R")

### Load packages
library(data.table)

### Load data
load(anno.snps.file)
load(anno.genes.file)

### Set up data tables
anno.snps$start <- pmax(anno.snps$position-1000000, 0)
anno.snps$end <- anno.snps$position+1000000
snps <- data.table(anno.snps)
setkey(snps, chromosome_name, start, end)

anno.genes$start_position <- as.numeric(anno.genes$start_position)
anno.genes$end_position <- as.numeric(anno.genes$end_position)
anno.genes$chromosome_name <- as.numeric(anno.genes$chromosome_name)
genes <- data.table(anno.genes)
setkey(genes, chromosome_name, start_position, end_position)

### Cross data tables
ans = foverlaps(snps, genes, type="any", nomatch=0L)
nb.rs <- ans[, .N, by = rs_id]
nb.rs <- data.frame(nb.rs, stringsAsFactors=F)
null.rs <- anno.snps$rs_id[!(anno.snps$rs_id %in% nb.rs$rs_id)]
nb.rs <- rbind(nb.rs, data.frame("rs_id"=null.rs, "N"=rep(0, length(null.rs)), stringsAsFactors=F))

save(nb.rs, file=paste0(eqtl.dir, snp.neighbors.file))

