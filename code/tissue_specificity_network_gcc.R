### tissue_specificity_network_gcc.R

###Load variables
source("code/variables_definition.R")
library(metap)

### Set parameters
snps.file <- paste0( 'all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
genes.file <- paste0( 'all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
edges.file <- paste0( 'all_tissues_edges_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
go.results.modif <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_final.txt')

snps.mult.file <- paste0( 'snps_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
genes.mult.file <- paste0( 'genes_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
edges.mult.file <- paste0( 'edges_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
ontology.mult.file <- paste0( 'ontology_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')

LDblock.rdata <- gsub("txt", "RData", LDblock.file)

snps.alltissues.mult.LD.file <- paste0( 'all_tissues_snps_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_LD_gcc.Rdata')
genes.alltissues.mult.file <- paste0( 'all_tissues_genes_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
edges.alltissues.mult.LD.file <- paste0( 'all_tissues_edges_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_LD_gcc.Rdata')
go.alltissues.mult.file <- paste0( 'all_tissues_ontology_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
summary.stats.mult <- "which_go_unique_common.RData"

### Load data
load(paste0(cluster.dir, snps.file))
load(paste0(cluster.dir, genes.file))
load(paste0(cluster.dir, edges.file))
go.tab <- read.delim(paste0(cluster.dir, go.results.modif), stringsAsFactors=F)

load(paste0(cluster.dir. snps.mult.file))
load(paste0(cluster.dir. genes.mult.file))
load(paste0(cluster.dir. edges.mult.file))
load(paste0(cluster.dir. ontology.mult.file))
load(paste0(geno.dir, LDblock.rdata))

### Functions
## Find proportion of SNPs/edges/Genes/Ontology for each multiplicity (1-12) in each community.
find.bycluster.multiplicity <- function(elt, mult){
    lapply(elt, function(x, m){
        return(m[names(m) %in% x])
    }, m=mult)
}
## Correcting for LD
find.bycluster.multiplicity.LD <- function(elt, mult, LDblocks){
    require(data.table)
    lapply(elt, function(x, m, LD){
        l <- m[names(m) %in% x]
        if(length(grep("ENSG", names(l)))>0){
            snpid <- gsub("_ENSG.*", "", names(l))
            geneid <- gsub(".*_ENSG", "ENSG", names(l))
            tmp <- data.table(data.frame("SNP"=snpid, "geneid"=geneid, "edges"=names(l), "multiplicity"=l, stringsAsFactors=F))
            setkey(tmp, SNP)
            data.merge <- data.frame(merge(tmp, LD, by='SNP'), stringsAsFactors=F)[,1:5]
            re <- as.character(tapply(data.merge$edges, paste(data.merge$block.id, data.merge$multiplicity, data.merge$geneid, sep="-"), function(x){sample(x, 1)}))
        } else {
            tmp <- data.table(data.frame("SNP"=names(l), "multiplicity"=l, stringsAsFactors=F))
            setkey(tmp, SNP)
            data.merge <- data.frame(merge(tmp, LD, by='SNP'), stringsAsFactors=F)[,1:3]
            re <- as.character(tapply(data.merge$SNP, paste(data.merge$block.id, data.merge$multiplicity, sep="-"), function(x){sample(x, 1)}))
        }
        l <- l[names(l) %in% re]
        return(l)
    }, m=mult, LD=LDblocks)
}

## Find unique elements in each community
get.unique <- function(elt, th){
    l <- lapply(elt, function(x, a){
        unlist(lapply(x, function(y, b){
            sum(y<=b)
        }, b=a))
    }, a=th)
    return(l)
}

## Find shared elements in each community
get.com <- function(elt, th){
    l <- lapply(elt, function(x, a){
        unlist(lapply(x, function(y, b){
            sum(y>=b)
        }, b=a))
    }, a=th)
    return(l)
}

## Obtain length of each element of a list
get.length <- function(elt){
    l <- lapply(elt, function(x){
        unlist(lapply(x, length))
    })
    return(l)
}

### Summarize all multiplicity for each tissue
get.sumstat <- function(elt, th1, th2){
    u.com <- get.unique(elt=elt, th=th1)
    c.com <- get.com(elt=elt, th=th2)
    l.com <- get.length(elt=elt)
    p.u.com <- lapply(names(u.com), function(a,x,y){
        x[[`a`]]/y[[`a`]]
    }, x=u.com, y=l.com)
    names(p.u.com) <- names(elt)
    p.c.com <- lapply(names(c.com), function(a,x,y){
        x[[`a`]]/y[[`a`]]
    }, x=c.com, y=l.com)
    names(p.c.com) <- names(elt)
    u <- unlist(lapply(u.com, sum))
    c <- unlist(lapply(c.com, sum))
    l <- unlist(lapply(l.com, sum))
    p.u <- u/l
    p.c <- c/l
    d.u.com <- lapply(names(p.u.com), function(a,x,y){
        x[[`a`]]-y[a]
    }, x=p.u.com, y=p.u)
    names(d.u.com) <- names(elt)
    d.c.com <- lapply(names(p.c.com), function(a,x,y){
        x[[`a`]]-y[a]
    }, x=p.c.com, y=p.c)
    names(d.c.com) <- names(elt)
    return(list("unique.com"=u.com, "common.com"=c.com, "length.com"=l.com, "prop.unique.com"=p.u.com, "prop.common.com"=p.c.com,
                "unique"=u, "common"=c, "length"=l, "prop.unique"=p.u, "prop.common"=p.c,
                "dist.unique"=d.u.com, "dist.common"=d.c.com))
}
## Compute ratio of tissue specific element in communities enriched for tissue-specific vs shared ontology terms 
ratio.mean.elt.unique.byelt <- function(pdffile, a, b, c, w.u, w.c, Tis ){
    find.distrib <- function(elt, w.u, w.c){
        v <- lapply(names(elt$prop.unique.com), function(n, x, y, z){
            l <- list(x$prop.unique.com[[`n`]][y[[`n`]][!(y[[`n`]] %in% z[[`n`]])]],
                      x$prop.unique.com[[`n`]][z[[`n`]][!(z[[`n`]] %in% y[[`n`]])]])
            
            pval = wilcox.test(l[[1]], l[[2]], alternative = "greater")$p.value
            m1 = mean(l[[1]])
            m2 = mean(l[[2]])
            m <- m1/m2
            f <- ifelse(pval <= 0.001, "***",
                 ifelse(pval <= 0.01, "**",
                 ifelse(pval <= 0.05, "*",
                        "NS")))
            return(c(m, f, pval))
            
        }, x=elt, y=w.u, z=w.c)
        matrix(unlist(v), ncol=3, byrow=T)
    }
    a.sum <- find.distrib(elt=a, w.u=w.u, w.c=w.c)
    b.sum <- find.distrib(elt=b, w.u=w.u, w.c=w.c)
    c.sum <- find.distrib(elt=c, w.u=w.u, w.c=w.c)

    data <- t(cbind(as.numeric(a.sum[13:1,1]), as.numeric(c.sum[13:1,1]), as.numeric(b.sum[13:1,1])))
    pval <- cbind(a.sum[13:1,2], c.sum[13:1,2], b.sum[13:1,2])
                                        # Plot results
    pdf(pdffile, width=5, height=8)
    par(las=1)
    bar <- barplot(data-1, beside=T, horiz=T, col=2:4, xlab="Ratio of tissue-specific elements (tissue-specific / shared communities)",
                   names=as.character(Tissues[names(nb.samples)[nb.samples>=200],2])[13:1], cex.names=0.8, xaxt='n')
    axis(side=1, at=(floor(par("usr")[1])+1):(floor(par("usr")[2])), labels=((floor(par("usr")[1])+1):(floor(par("usr")[2])))+1)
    legend("top", bty='n', legend=c("Genes", "Edges", "SNPs"), fill=4:2)
    par(xpd=T)
    for(i in 1:nrow(pval)){
        text(rep(par("usr")[2], 3), bar[,i], labels=as.character(pval[i,]), cex=0.6, pos=4)
    }
     
   # dev.off()
    return(list("snps"=a.sum, "genes"=b.sum, "edges"=c.sum))
}

### reshape GO results in an object similar to snps, genes and edges
length.com <- unlist(lapply(snps, length))
go.tab$Tissue <- unlist(lapply(go.tab$Tissue, function(x, T){rownames(T)[T[,1]==x]}, T = Tissues))
go <- tapply(paste(go.tab$GOID, go.tab$comm.id, length.com[go.tab$Tissue], sep="_"), go.tab$Tissue,
             function(x){
                 com <- as.numeric(unlist(lapply(x, function(y){strsplit(y, "_")[[1]][2]})) )
                 GOID <- unlist(lapply(x, function(y){strsplit(y, "_")[[1]][1]}))
                 com <- factor(com, levels=1:as.numeric(strsplit(x[1], "_")[[1]][3]))
                 return(tapply(GOID, com, function(z){z}))
             }
             )

### Find distribution of SNPs/Genes/Edges/Ontology modularity among communities
snps.mult.LD <- list()
genes.mult <- list()
edges.mult.LD <- list()
go.mult <- list()
for(n in names(nb.samples)[nb.samples>=200]){
    cat("running multiplicity for tissue", n, "\n")
    cat("\t snps\n")
    snps.mult.LD[[`n`]] <- find.bycluster.multiplicity.LD(elt=snps[[`n`]], mult=s, LDblocks=block.info)
    cat("\t genes\n")
    genes.mult[[`n`]] <- find.bycluster.multiplicity(elt=genes[[`n`]], mult=g)
    cat("\t edges\n")
    edges.mult.LD[[`n`]] <- find.bycluster.multiplicity.LD(elt=edges[[`n`]], mult=e, LDblocks=block.info)
    cat("\t ontology\n")
    go.mult[[`n`]] <- find.bycluster.multiplicity(elt=go[[`n`]], mult=o) 
}

save(snps.mult.LD, file=paste0(cluster.dir, snps.alltissues.mult.LD.file))
save(genes.mult, file=paste0(cluster.dir, genes.alltissues.mult.file))
save(edges.mult.LD, file=paste0(cluster.dir, edges.alltissues.mult.LD.file))
save(go.mult, file=paste0(cluster.dir, ontology.alltissues.mult.file))


### Obtain summary statistics about multiplicity for SNPs/Genes/Edges/Ontology at different thresholds
## Threshold 2: unique = present in 1 or 2 tissues, shared = present in at least 12 of the 13 tissues
snps.sumstats.2.LD <- get.sumstat(elt=snps.mult.LD, th1=2, th2=12)
genes.sumstats.2 <- get.sumstat(elt=genes.mult, th1=2, th2=12)
edges.sumstats.2.LD <- get.sumstat(elt=edges.mult.LD, th1=2, th2=12)
go.sumstats.2 <- get.sumstat(elt=go.mult, th1=2, th2=12)

### For each threshold, find the communities that are enriched for unique and shared ontology terms
## Threshold 2
which.go.unique.2 <- lapply(go.sumstats.2$dist.unique, function(x){which(x>0 & !is.na(x))})
which.go.common.2 <- lapply(go.sumstats.2$dist.common, function(x){which(x>0 & !is.na(x))})

save(snps.sumstats.2.LD, edges.sumstats.2.LD,
     genes.sumstats.2, go.sumstats.2,
     which.go.unique.2, which.go.common.2,
     file=paste0(cluster.dir, summary.stats.mult))

### Plot ratio of tissue specific element in communities enriched for tissue-specific vs shared ontology terms 
ratio.mean.elt.2.LD <- ratio.mean.elt.unique.byelt(pdffile=paste0(figure.dir, "ratio_mean_unique_2_LD_gcc.pdf"),
                                                   a=snps.sumstats.2.LD, b=genes.sumstats.2,
                                                   c=edges.sumstats.2.LD, w.u=which.go.unique.2, w.c=which.go.common.2, Tis=Tissues)

### Compute combined p-values
pval.combined.2 <- unlist(lapply(ratio.mean.elt.2.LD, function(x){ sumlog(as.numeric(x[,3]))$p }))
tab.pval.combined <- data.frame( "th2"=pval.combined.2)
write.table(tab.pval.combined, file=paste0(cluster.dir, "ratio_mean_unique_pval_combined.txt"), quote=F, sep="\t")

