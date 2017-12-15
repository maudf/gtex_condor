###count_eqtls.R

### Load variables
source("code/variables_definition.R")

### Set variables
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")

### Load data
load(tissue.file)
load(paste0(eqtl.dir, eqtls.file))
count.unique.qtl <- function(qtl, type, fdr.th){
    rsid <- qtl$RSID[(qtl$cis.or.trans %in% type) & (qtl$FDR<=fdr.th)]
    length(unique(rsid))
}

### Count snps for each tissue
nb.snps=NULL
for(n in names(eqtl)){
    nb.snps <- c(nb.snps,
                 strsplit(system(paste0("wc -l ", geno.dir, n, ".pos"), intern=T), " ")[[1]][1])
}

### Compute proportion SNPs that are eQTLs for 4 FDR levels
prop.table <- NULL
for(f in c(0.05, 0.1, 0.15, 0.2)){
    for(n in names(eqtl)){
        nb.cis <- count.unique.qtl(eqtl[[`n`]], type="cis", fdr.th=f)
        prop.cis <- nb.cis/nb.snps[`n`]
        nb.trans <- lapply(eqtl, count.unique.qtl, type="trans")
        prop.trans <- nb.trans/nb.snps[`n`]
        prop.table <- rbind(prop.table, c(as.character(Tissues[n, 2]), f,
                                          nb.cis, prop.cis,
                                          nb.trans, prop.trans))
    }
}
colnames(prop.table) <- c("Tissue", "FDR", "#cis-eQTLs", "% cis-eQTLs", "#trans-eQTLs", "% trans-eQTLs")
write.table(prop.table, file=paste(eqtl.dir, "proportion_eqtls.csv"), sep=",", row.names=F, quote=F)

