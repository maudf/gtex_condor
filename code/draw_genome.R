### draw_genome.R

###Load variables
load("code/variables_definition.R")

### Set variables
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
genes.file <- paste0("all_tissues_genes_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")


### Load data
heart.epi <- read.table(data.dir, "E095_15_coreMarks_mnemonics.bed", header=F, stringsAsFactors=F)
esophagus.epi <- read.table(data.dir, "E079_15_coreMarks_mnemonics.bed", header=F, stringsAsFactors=F)
colnames(heart.epi) <- colnames(esophagus.epi) <- c("chr", "start", "end", "region")
load(eqtl.dir, eqtls.file)
load(cluster.dir, genes.file)
load(anno.snps.file)
load(anno.genes.file)
load(genes.corres.file)

### Functions
draw.genome <- function(g, epi, pdffile){
    require(biomaRt)
    g$ENSG <- g$ENSG[g$annotation.genes$chromosome_name %in% g$annotation.snp$chromosome_name]
    g$HGNC <- g$HGNC[g$annotation.genes$chromosome_name %in% g$annotation.snp$chromosome_name]
    names(g$HGNC) <- gsub("\\.[0-9]*", "",g$ENSG)
    g$annotation.genes <- g$annotation.genes[g$annotation.genes$chromosome_name %in% g$annotation.snp$chromosome_name,]

    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    exons <- getBM(attributes=c('ensembl_gene_id', 'start_position', 'end_position', 'ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end'), 
                   filters = 'ensembl_gene_id', 
                   values = gsub("\\.[0-9]*", "",
                                 g$ENSG),
                   mart = ensembl)

    snpppos <-  g$annotation.snp$position
    min.pos <- min(as.numeric(c(exons$exon_chrom_start, exons$exon_chrom_end,
                     exons$start_position,  exons$end_position,
                     snpppos)))
    max.pos <- max(as.numeric(c(exons$exon_chrom_start, exons$exon_chrom_end,
                     exons$start_position,  exons$end_position,
                     snpppos)))
    d <- max.pos-min.pos
    lim <- c(min.pos-d/20, max.pos+d/20)

    epi.hg18=data.table(epi)
    epi.lim <- data.frame(epi.hg18[(chr == paste0("chr", g$annotation.snp$chromosome_name)) &
                   (((start >= lim[1]) & (start <= lim[2])) |
                    ((end >= lim[1]) & (end <= lim[2]))),], stringsAsFactors=F)
    epi.lim[,4] <- as.character(epi.lim[,4])
    color <- c("1_TssA"=rgb(255,0,0,maxColorValue=255),
               "2_TssAFlnk"=rgb(255,69,0,maxColorValue=255),
               "3_TxFlnk"=rgb(50,205,50,maxColorValue=255),
               "4_Tx"=rgb(0,128,0,maxColorValue=255),
               "5_TxWk"=rgb(0,100,0,maxColorValue=255),
               "6_EnhG"=rgb(194,225,5,maxColorValue=255),
               "7_Enh"=rgb(255,255,0,maxColorValue=255),
               "8_ZNF/Rpts"=rgb(102,205,170,maxColorValue=255),
               "9_Het"=rgb(138,145,208,maxColorValue=255),
               "10_TssBiv"=rgb(205,92,92,maxColorValue=255),
               "11_BivFlnk"=rgb(233,150,122,maxColorValue=255),
               "12_EnhBiv"=rgb(189,183,107,maxColorValue=255),
               "13_ReprPC"=rgb(128,128,128,maxColorValue=255),
               "14_ReprPCWk"=rgb(192,192,192,maxColorValue=255),
               "15_Quies"=rgb(255,255,255,maxColorValue=255))
    
    pdf(pdffile, height=4, width=8)
    par(mar=c(1,1,2,1)+.1)
    plot(lim, c(0,1.5),
         xaxt='n', yaxt='n', ylab='',type='n', frame.plot=F)
    div = ifelse(d<200000, 20000, ifelse(d<1000000, 50000, 100000))
    mult=div/1000000
    axis.lim.1 <- floor(lim[1]/div) - 1
    axis.lim.2 <- floor(lim[2]/div) + 1
    axis(side=3, at=seq(axis.lim.1*div, axis.lim.2*div, div),
         labels= format(seq(axis.lim.1*mult, axis.lim.2*mult, mult ), nsmall=2))
    mtext(side=3, line=2, text=paste0("chr", g$annotation.snp$chromosome_name, " (Mb)"))
    a <- unique(exons[,1:3])
    f <- rep(0, nrow(a))
    names(f) <- gsub("\\.[0-9]*", "", g$ENSG)
    for(i in 1:nrow(a)){
        
        if (i>1) {
            f[i] <- sum(apply(a[1:(i-1),], 1, function(x,y){
                return(((y[1,2] >= x[2] & y[1,2] <= x[3]) |
                        (y[1,3] >= x[2] & y[1,3] <= x[3])))
            }, y=a[i,]))
        }
        lines( c(a[i,2], a[i,3]), (rep(0.7,2)-f[i]/3), col="blue")
        lines(rep(a[i,2], 2), (c(0.6,0.8)-f[i]/3), col="blue")
        lines(rep(a[i,3], 2), (c(0.6,0.8)-f[i]/3), col="blue")
        par(xpd=T)
        text(a[i,2], (0.7-f[i]/3), labels=g$HGNC[a[i,1]], font=2, pos=2, col="blue")
        par(xpd=NA)
    }
    
    for(i in 1:nrow(exons)){
        polygon(c(exons[i,5], exons[i,5], exons[i,6], exons[i,6]),
                (c(0.6, 0.8, 0.8, 0.6)-f[exons[i,1]]), col="blue", border="blue")
    }
    for(i in 1:nrow(epi.lim)){
        polygon(c(epi.lim[i,2], epi.lim[i,2], epi.lim[i,3], epi.lim[i,3]),
                c(1, 1.2, 1.2, 1), col=color[epi.lim[i,4]], border=color[epi.lim[i,4]])
        
    }
    lines(rep(snpppos, 2), c(0,1.25), col="black")
    lines(lim, rep(1.2, 2), col="black")
    lines(lim, rep(1, 2), col="black")
    points(snpppos, 0, col="black", pch=17)
    par(xpd=T)
    text(snpppos, 0, pos=1, labels=g$annotation.snp$rs_id, col="black")
    par(xpd=NA)
    polygon(c(lim, lim[2:1]), c(rep(1.5,2), rep(1.4,2)), col="gray90", border="gray90")
    lines(lim, rep(1.5, 2), col="black")
    lines(lim, rep(1.4, 2), col="black")
    text(sum(lim)/2, 1.45, col="black", paste0("chr", g$annotation.snp$chromosome_name))

    dev.off()
    
}
### Plot genome pictures for rs2546765 in heart left ventricle and rs4072037 in esophagus mucosa

draw.genome(g=find.info.snp(rsid="rs2546765", tis="heart_left_ventricle", eqtl=eqtl, anno.snps=anno.snps, anno.genes=anno.genes, gene.tab=gene.tab, genes=genes),
            epi=heart.epi, pdffile=paste0(figure.dir, "example_rs254765.pdf"))
draw.genome(g=find.info.snp(rsid="rs4072037", tis="esophagus_mucosa", eqtl=eqtl, anno.snps=anno.snps, anno.genes=anno.genes, gene.tab=gene.tab, genes=genes),
            epi=heart.epi, pdffile=paste0(figure.dir, "example_rs4072037.pdf"))
