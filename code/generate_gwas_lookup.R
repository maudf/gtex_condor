### generate_gwas_lookup.R

### Load variables
source("code/variables_definition.R")

### Load data
gwas_data <- read.delim( gwas.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gwas_data <- gwas_data[gwas_data$PVALUE_MLOG>=8,]
gwas_data <- gwas_data[!is.na(gwas_data$PVALUE_MLOG),]

### Generate lookup table data
tab <- data.frame(gwas_data[, c("SNPS", "DISEASE.TRAIT")], stringsAsFactors=F)
tab <- unique(tab)

write.table(tab, file=paste0(gwas.dir, gwas.diseases.file),
            row.names=F, col.names=F, quote=F, sep="\t")
