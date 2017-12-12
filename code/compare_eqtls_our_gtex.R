###compare_eqtls_our_gtex.R

### Load variables
load("code/variables_definition.R")

### Set variables
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
comp.gtex.our <- paste0("comparison_with_gtex_fdr", FDRcis, FDRtrans, "_", window, "MB_fdr0.05.txt")

### Load data
load(tissue.file)
load(gtexqtl.file)
load(paste0(eqtl.dir , eqtls.file))
load(tissue.file)

### Functions
extract.qtls <- function(qtl,fdr){
    tab <- qtl
   if(length(grep("VariantID", colnames(qtl)))>0){
       tab <- tab[tab$cis.or.trans=="cis",]
       r <- paste(tab$VariantID, tab$genes, sep="_")
       r <- r[tab$FDR<=fdr]
    } else {
        r <- paste(tab$variant_id, tab$gene_id, sep="_")
    }    
    return(r)
}

### Extract eqtls from gtex data and our data

our.eqtl.list <- lapply(eqtl, extract.qtls, fdr=0.05)
gtex.eqtl.list  <- lapply(gtex.eqtl, extract.qtls, fdr=0.05)
gtex.eqtl.list$skin <- unique(c(gtex.eqtl.list$skin_not_sun_exposed_suprapubic,
                                gtex.eqtl.list$skin_sun_exposed_lower_leg))

res <- matrix(ncol=3, nrow=0)
for(n in names(our.eqtl.list)){
    res <- rbind(res, c(sum(gtex.eqtl.list[[`n`]] %in% our.eqtl.list[[`n`]]),
                        length( gtex.eqtl.list[[`n`]]),
                        length( our.eqtl.list[[`n`]]))
                 )
}
colnames(res) <- c("Common", "GTEx", "Our_analysis")
rownames(res) <- as.character(Tissues[ names(our.eqtl.list), 1])
write.table(res, file = paste0(general.path, eqtl.path, comp.gtex.our), quote=F, sep="\t")
