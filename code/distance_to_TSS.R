### distance_to_TSS.R

### Load variables
source("code/variables_definition.R")

### Set parameters
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
disttssfile <- "dist_snp_tss.RData"
pdfeqtlcisfile <- paste0("dist_tss_cis_fdr",FDRcis, FDRtrans, "_", window, "MB.pdf")
pdfeqtltransfile <- paste0("dist_tss_trans_fdr",FDRcis, FDRtrans, "_", window, "MB.pdf")
quantilecistrans <- paste0("dist_tss_summary_cis_trans_fdr",FDRcis, FDRtrans, "_", window, "MB.txt")

### Load data
load(anno.snps.file)
load(anno.genes.file)
load(paste0(eqtl.dir, eqtls.file))
load(tissue.file)

###Functions

##Extract list of cis and trans SNPs
extract.cis.or.trans.snp <- function(x){
  cis <- unique(x$RSID[x$cis.or.trans=="cis"])
  trans <- unique(x$RSID[x$cis.or.trans=="trans"])
  return(list("cis"=cis, "trans"=trans))
}

## Compute distance between a snp and the nearest TSS. 
compute.distance <- function(snp, genes){
  d <- as.numeric(snp[2]) - as.numeric(genes$transcript_start)
  j <- which(abs(d)==min(abs(d)))[1]
  return(c(d[j], rownames(genes)[j]) )
}
## Plot distance to TSS of cis- and trans-eQTLs

plot.dist.tss <- function(s, dist.tss, step, main, xlim){
    par(mar=c(4,5,4,1)+0.1)
    d <- as.numeric(dist.tss$nearest.tss[ rownames(dist.tss) %in% s])
    h <- hist(d, breaks=seq(min(d)-step, max(d)+step, step), plot=F)
    h$density <- h$counts/sum(h$counts)
    plot(h, freq=F, xlab="Distance to TSS (kb)", xaxt='n', ylab="Frequency", main=main, col="dodgerblue", xlim=xlim)
    axis(side=1, at=seq(min(xlim), max(xlim), (xlim[2]-xlim[1])/4),
         labels=seq(min(y=xlim)/1000, max(xlim)/1000,  (xlim[2]-xlim[1])/4000))
}


### Extract distance to TSS for each snp
tss.genes <- anno.genes[,c("chromosome_name", "transcript_start")]
tss.genes$transcript_start[anno.genes$strand == -1] <- anno.genes$transcript_end[anno.genes$strand == -1]
tss.genes$transcript_start <- as.numeric(tss.genes$transcript_start)

dist.tss <- anno.snps[, c("chromosome_name", "position")]
dist.tss$nearest.tss <- rep(NA, nrow(dist.tss))
dist.tss$nearest.gene <- rep(NA, nrow(dist.tss))

for(chr in unique(dist.tss$chromosome_name)){
    cat("Running chromosome", chr, "\n")
    a <- t(apply(dist.tss[dist.tss$chromosome_name==chr,], 1, compute.distance,
               genes=tss.genes[tss.genes$chromosome_name==chr,]))
    dist.tss[dist.tss$chromosome_name==chr, c("nearest.tss", "nearest.gene")] <- a[,1:2]
}
save(dist.tss, file=paste0(eqtl.dir, disttssfile))

### Extract snps by tissue

snp <- lapply(eqtl, extract.cis.or.trans.snp)

### Extract 50% distance
eqtl <- lapply(eqtl, function(x){x[x$FDR<=0.05,]})
qtl.cis <- lapply(eqtl, function(x){unique(x$RSID[x$cis.or.trans=="cis"])})
qtl.trans <- lapply(eqtl, function(x){unique(x$RSID[x$cis.or.trans=="trans"])})
d.cis <- lapply(qtl.cis, function(x,d){d[rownames(d) %in% x,3]}, d=dist.tss)
d.trans <- lapply(qtl.trans, function(x,d){d[rownames(d) %in% x,3]}, d=dist.tss)
range(unlist(lapply(d.cis, function(x) quantile(abs(as.numeric(x)), 0.5)))[nb.samples>200])
range(unlist(lapply(d.trans, function(x) quantile(abs(as.numeric(x)), 0.5)))[nb.samples>200])

### Get quantile distribition of distance between cis- and trans-eQTLs and nearest TSS
q.cis <-matrix(ncol=9, nrow=0)
q.trans <-matrix(ncol=9, nrow=0)
for(i in 1:length(snp)){
    s <- snp[[i]]$cis
    s2 <- snp[[i]]$trans
    d <- as.numeric(dist.tss$nearest.tss[ rownames(dist.tss) %in% s])
    d2 <- as.numeric(dist.tss$nearest.tss[ rownames(dist.tss) %in% s2])
    
    q.cis <- rbind(q.cis, quantile(d, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)) )
    q.trans <- rbind(q.trans, quantile(d2, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)) )
}
rownames(q.cis) <- rownames(q.trans) <- names(eqtl)
write.table(rbind(q.cis, q.trans), file=paste0( eqtl.dir, quantilecistrans), quote=F, sep="\t")

### Plot distance to TSS of cis-eQTLs

pdf(paste( figure.dir, pdfeqtlcisfile, sep=""),width=8, height=11)
par(mfrow=c(3,2))
for(i in 1:length(snp)){
    plot.dist.tss(snp[[i]]$cis, dist.tss, step=2000,
                  main=Tissues[names(snp)[i],1], xlim=c(-100000, 100000))
}
dev.off()

### Plot distance to TSS of trans-eQTLs
pdf(paste(figure.dir, pdfeqtltransfile, sep=""),width=8, height=11)
par(mfrow=c(3,2))
for(i in 1:length(snp)){
    plot.dist.tss(snp[[i]]$trans, dist.tss, 2000,
                  main=Tissues[names(snp)[i],1], xlim=c(-100000, 100000))
}
dev.off()
