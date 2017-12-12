###summarize_eqtls.R

### Load variables
load("code/variables_definition.R")
load(paste0(tissue.file))

eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")

### Obtain list of eqtls files
list.files.eqtls <- list.files(eqtl.dir, pattern="*_edges.txt")

### Read eqtls and retrieve number of samples

eqtl <- list()
nb.samples <- c()

for(f in 1:length(list.files.eqtls)){

    tissue <- gsub(paste0("*fdr", FDRcis, FDRtrans, "_", window, "MB_edges.txt"), "", f)

    cat("Loading tissue ", tissue, "...\n", sep="")
    qtl.tmp <- read.table(paste0(eqtl.dir, list.files.eqtls[f]), header=T, stringsAsFactors=F, quote="")
    colnames(qtl.tmp)[2:3] <- colnames(qtl.tmp)[3:2]
    if(("RS_ID_dbSNP142_CHG37p13" %in% colnames(qtl.tmp)) & ("RS_ID_dbSNP135_original_VCF" %in% colnames(qtl.tmp))){
        qtl.tmp$RSID <- qtl.tmp$RS_ID_dbSNP142_CHG37p13
        qtl.tmp$RSID[qtl.tmp$RSID=="."] <- qtl.tmp$RS_ID_dbSNP135_original_VCF[qtl.tmp$RSID=="."]
    }
    eqtl[[`tissue`]] <- qtl.tmp
    rm(qtl.tmp)
    
    geno <- scan( pipe( paste0("head -1 ", general.path, geno.path, tissue, ".dosage")), what = character(0))
    rna <- scan( pipe( paste0("head -1 ", general.path, rna.path, tissue, "_norm.tsv")), what = character(0))
    geno <- gsub("(GTEX-[^-]*)-.*", "\\1", geno, perl=TRUE)
    rna <- gsub("(GTEX-[^-]*)-.*", "\\1", rna, perl=TRUE)
    nb.samples <- c(nb.samples, sum(geno %in% rna))
}
names(nb.samples) = names(eqtl)

cat("Save files...\n", sep="")
save(eqtl, file=paste0(eqtl.dir, eqtls.file))
save(nb.samples, file=paste0(eqtl.dir, samples.file))
