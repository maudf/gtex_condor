### generate_circos.R

###Load variables
load("code/variables_definition.R")

### Set parameters
snps.file <- paste0( 'all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
genes.file <- paste0( 'all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
edges.file <- paste0( 'all_tissues_edges_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
eqtl.file <- paste0( 'all_tissues_eqtls_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
term.file <- paste0("data/gwas_metabolism_terms.txt")

snp.circos <- paste0( code.dir, 'circos_plot/snps.txt')
gene.circos <- paste0( code.dir, 'circos_plot/genes.txt')
edges.trans.circos <- paste0( code.dir, 'circos_plot/edges_trans.txt')
edges.cis.circos <- paste0( code.dir, 'circos_plot/edges_cis.txt')
snps.circos.description <- paste0(circos.dir, 'snp_description.txt')
genes.circos.description <- paste0(circos.dir, 'gene_description.txt')

tissue<- "heart_left_ventricle"
com <- 86

### Load data
load(paste0(cluster.dir, snps.file))
load(paste0(cluster.dir, genes.file))
load(paste0(cluster.dir, edges.file))
load( gene.anno.file)
load( snp.anno.file)
load( corres.file)
load( eqtl.dir, eqtl.file)
gwas_data <- read.delim( gwas.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gwas_data <- gwas_data[gwas_data$PVALUE_MLOG>=8,]
gwas_data <- gwas_data[!is.na(gwas_data$PVALUE_MLOG),]
gwas.diseases<-read.delim(paste0(gwas.dir, gwas.diseases.file), header=F, stringsAsFactors=F)
gwas.term <- scan(term.file, sep="\t", what=character(0))

# Extract genes involved in cellular respiration
genes.cell.resp <- scan( genes.hlv.86.file, sep="\t", what=character(0))
genes.cell.resp.ensg <- gsub("[^ENSG0-9][0-9]*$", "", unique(rownames(gene.tab)[gene.tab[,1] %in% genes.cell.resp]))

# Extract snps involved in metabolism
gwas.snps <- gwas_data$SNPS[gwas_data$SNPS %in% snps[[`tissue`]][[com]]]
diseases <- gwas.diseases[gwas.diseases[,1] %in% gwas.snps,]
gwas.snps.metabolism <- unique(diseases[diseases[,2] %in% gwas.term,1])

### Create snps table
s <- snps[[`tissue`]][[com]]
anno.s <- anno.snps[rownames(anno.snps) %in% s,]
tab.snp <- cbind(paste0("hs", anno.s$chromosome_name),
                    (anno.s$position-1), (anno.s$position+1),
                    rep("fill_color=vlgreen,stroke_color=vlgreen,stroke_thickness=5", length(s)))
tab.snp[rownames(anno.s) %in% gwas.snps.metabolism, 4] <- "fill_color=vdgreen,stroke_color=vdgreen,stroke_thickness=10"
tab.snp.final <- rbind(tab.snp[grep("vlgreen", tab.snp[,4]),], tab.snp[grep("vdgreen", tab.snp[,4]),])

### Create genes table

g <- gsub("[^ENSG0-9][0-9]*$", "", genes[[`tissue`]][[com]])
anno.g <- anno.genes[rownames(anno.genes) %in% g,]
tab.gene <- cbind(paste0("hs", anno.g$chromosome_name),
                  anno.g$start_position, anno.g$end_position,
                  rep("fill_color=vlblue,stroke_color=vlblue,stroke_thickness=5", length(g)) )
tab.gene[rownames(anno.g) %in% genes.cell.resp.ensg, 4] <- "fill_color=vdblue,stroke_color=vdblue,stroke_thickness=10"
tab.gene.final <- rbind(tab.gene[grep("vlblue", tab.gene[,4]),], tab.gene[grep("vdblue", tab.gene[,4]),])

### Create cis and trans edges table

qtl <- eqtl[[`tissue`]]
rownames(qtl) <- paste(qtl$RSID, gsub("[^A-Za-z0-9][0-9]*$", "", qtl$genes), sep="_")
e <- gsub("[^A-Za-z0-9][0-9]*$", "", edges[[`tissue`]][[com]])
e.m <- matrix(unlist(strsplit(e, "_")), ncol=2, byrow=T)
tab.edges <- cbind(paste0("hs", anno.s[e.m[,1], "chromosome_name"]),
                    (anno.s[e.m[,1], "position"]-1),
                    (anno.s[e.m[,1], "position"]+1),
                    paste0("hs", anno.g[e.m[,2], "chromosome_name"]),
                    anno.g[e.m[,2], "start_position"],
                    anno.g[e.m[,2], "end_position"],
                    rep( "color=lgrey,thickness=1", length(e)))

tab.edges[(e.m[,1] %in% gwas.snps.metabolism) | (e.m[,2] %in% genes.cell.resp.ensg), 7] <- "color=dred,thickness=3"

tab.edges.cis <- rbind(tab.edges[(qtl[rownames(qtl) %in% e, "cis.or.trans"]=="cis") & (1:nrow(tab.edges) %in% grep("lgrey",tab.edges[,7])),], tab.edges[(qtl[rownames(qtl) %in% e, "cis.or.trans"]=="cis") & (1:nrow(tab.edges) %in% grep("dred",tab.edges[,7])),])

tab.edges.trans <- rbind(tab.edges[(qtl[rownames(qtl) %in% e, "cis.or.trans"]=="trans") & (1:nrow(tab.edges) %in% grep("lgrey",tab.edges[,7])),], tab.edges[(qtl[rownames(qtl) %in% e, "cis.or.trans"]=="trans") & (1:nrow(tab.edges) %in% grep("dred",tab.edges[,7])),])

### Write circos input tables and circos annotation files
write.table(tab.snp.final, snp.circos,
            row.names=F, col.names=F, quote=F)
write.table(anno.s[rownames(anno.s) %in% gwas.snps.metabolism, ],
            snps.circos.description,
            row.names=F, col.names=F, quote=F)

anno.g <- anno.g[genes.cell.resp.ensg,]
anno.g$HGNC <- gene.tab[anno.g$ensembl_gene_id,1]
write.table(anno.g, genes.circos.description,
            row.names=F, col.names=F, quote=F)
write.table(tab.gene.final, gene.circos,
            row.names=F, col.names=F, quote=F)

write.table(tab.edges.trans, edges.trans.circos,
            row.names=F, col.names=F, quote=F)
write.table(tab.edges.cis, edges.cis.circos,
            row.names=F, col.names=F, quote=F)

