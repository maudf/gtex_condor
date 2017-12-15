### plot_OR_epigenomic_states_across_tissues.R

### Load variables
source("code/variables_definition.R")
library(gtools)
library(data.table) 

### Set Parameters
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
topo.file <- paste0("network_topology_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")

threshold.topo.file <- paste0("all_tissues_epi_topo_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")
chromatin.nbgenes.pdf <- paste0("chromatin_alltissues_nbgenes_Qi_com_fdr", FDRcis, FDRtrans, "_", window, "MB.pdf")
chromstate.txt.file <- paste0("chromatinstates_OR_alltissues_nbgenes_Qi_com_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
chromstate.res.file <- paste0("chromatinstates_OR_alltissues_Qi_com_fdr", FDRcis, FDRtrans, "_", window, "MBRData")
chromatin.trans.file <- paste0("chromatintranseQTL_OR_alltissues_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
cistrans.res.file <- paste0("chromatintranseQTL_OR_alltissues_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")
cistrans.chromstate.res.file <- paste0("clogit_chromatintranseQTL_OR_alltissues_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")

###Load data
load(paste0(tissue.file))
tmp <- load(paste0(eqtl.dir, topo.file))
topo <- get(tmp)
tmp <- load(paste0(eqtl.dir, eqtls.file))
qtls <- get(tmp)
tmp <- load(paste0(eqtl.dir, snp.neighbors.file)))
neighbors <- get(tmp)
tmp <- load(paste0(data.dir, epi.file))
chromstates <- get(tmp)
tmp <- load( tissue.file)
Tissues <- get(tmp)
rm(tmp)

### Functions.

###Function building count table from a data table for 2 categories and a threshold.
make.table <- function(data, categ="A", colnm="new.reg.score"){
    a <- numeric(length=nrow(data))
    a[data[,colnm] == categ] <- 1
    return(a)
}

###Function plotting odds ratio and confidence intervals.
plot.or <- function(x, xlab, ylab, main, col, l, xlim){
    if(nrow(x)==15){
        x=x[c(1, 8:15, 2:7),]
    }
    par(las=1)
    plot(xlim, c(0.5, (nrow(x)+0.5)), type='n', xlab=xlab, ylab="", main=main, yaxt='n', col=col)
    axis(side=2, at=1:nrow(x), labels=x[nrow(x):1,1], cex.axis=0.9)
    par(las=0)
    mtext(ylab, side=2, line=l)
    lines(c(1,1), c(0,(nrow(x)+3)), col='red', lwd=2)
    for(i in 1:nrow(x)){
        lines(c(x[(nrow(x)+1-i),5], x[(nrow(x)+1-i),6]), c(i,i), lwd=2, col=col)
        lines(c(x[(nrow(x)+1-i),5], x[(nrow(x)+1-i),5]), c(i-0.3,i+0.3), lwd=2, col=col[1])
        lines(c(x[(nrow(x)+1-i),6], x[(nrow(x)+1-i),6]), c(i-0.3,i+0.3), lwd=2, col=col[1])
    }
    points(x[,3], nrow(x):1, col="black", bg=col, pch=21, cex=1 )
}

###high degree and core-snps search
tis <- names(topo)
data <- NULL
for(i in tis){
    cat ("Running", i, ":\n")
    tmp <- merge(topo[[`i`]], regdb, by.x="red.names", by.y="RsID")
    if(i %in% names(chromstates)){
        tmp <- merge(tmp, chromstates[[`i`]][,4:5], by.x="red.names", by.y="snp")
    } else {tmp$region <- rep(NA, nrow(tmp))}
    tmp <- merge(tmp, neighbors, by.x="red.names", by.y="rs_id")

    tmp$cs.thres <- rep(0, nrow(tmp))
    thres.Qi <- quantile(tmp$Qi, 0.75)
    tmp$cs.thres[tmp$Qi>=thres.Qi] <- 1
    tmp$dg.thres <- rep(0, nrow(tmp))
    tmp$dg.thres[tmp$degree>=10] <- 1
    tmp$tis <- rep(i, nrow(tmp))
    data <- rbind(data, tmp)
}
save(data, file=paste0(epi.dir, threshold.topo.file))

### Enrichment in high core-scores and high degrees among chromatin states

library(survival)
cs.epi <- list()
dg.epi <- list()

for(ca in sort(unique(data$region))){
    cat("Running chromstats for category", ca, "\n")
    data.tmp.epi  = data[!is.na(data$region),]
    epi <- make.table(data = data.tmp.epi,
                      categ=ca, colnm="region")
    dat <- data.frame("chrom"=epi, "cen.cs"=data.tmp.epi$cs.thres, "cen.dg"=data.tmp.epi$dg.thres,
                      "com"=as.factor(data.tmp.epi$com), "nbgenes"=data.tmp.epi$N, "tissue"=data.tmp.epi$tis)
    cs.epi[[`ca`]] <- clogit(chrom ~ cen.cs + com + nbgenes + strata(tissue),
                             data=dat, method="approximate")
    dg.epi[[`ca`]] <- clogit(chrom ~ cen.dg + nbgenes + strata(tissue),
                             data=dat, method="approximate")
}
save(cs.epi, dg.epi, file=paste0(epi.dir, chromstate.res.file))

###Summarize high core-score/high degree enrichment results
CI.cs.epi <- list()
CI.dg.epi <- list()
res.epi <- NULL
for(ca in names(cs.epi)){
    s.cs.epi<- summary(cs.epi[[`ca`]])$coefficient
    s.dg.epi<- summary(dg.epi[[`ca`]])$coefficient
    CI.cs.epi[[`ca`]] <- exp(confint(cs.epi[[`ca`]]))
    CI.dg.epi[[`ca`]] <- exp(confint(dg.epi[[`ca`]]))
    res.epi <- rbind(res.epi, rbind(c(ca, "core-score", s.cs.epi[1,2], s.cs.epi[1,5],
                                      CI.cs.epi[[`ca`]][1,1], CI.cs.epi[[`ca`]][1,2]),
                                    c(ca, "degree", s.dg.epi[1,2], s.dg.epi[1,5],
                                      CI.dg.epi[[`ca`]][1,1], CI.dg.epi[[`ca`]][1,2])))
}
colnames(res.epi) <- c("categ", "centrality", "OR", "pval", "CI.2.5", "CI.97.5")
write.table(res.epi, file=paste0(epi.dir, chromstate.txt.file, quote=F, sep="\t"))


### Set up green and blue color scales
greens <- colorRampPalette(c("forestgreen", "white"))
blues <- colorRampPalette(c("royalblue4", "white"))

### Plot enrichment scores

pdf(paste0(figure.dir, chromatin.nbgenes.pdf), width = 8, height = 5)
par(mfrow=c(1,2))
par(mar=c(2,7,0,1)+.1)
plot.or(x=res.epi[res.epi[,2]=="core-score",], ylab = "Chromatin State", xlab="", main = "Core scores",
         col = greens(nrow(res.epi[res.epi[,2]=="core-score",])), xlim=c(0,3), l=5)
plot.or(x=res.epi[res.epi[,2]=="degree",], ylab = "Chromatin State", xlab="", main = "Degrees",
         col = blues(nrow(res.epi[res.epi[,2]=="degree",])), xlim=c(0,3), l=5)
dev.off()

###trans-eQTLs search
tis <- names(chromstates)
categories <- mixedsort(unique(chromstates[[`1`]]$region))
res.trans <- NULL
for(i in tis){
    cat ("Running", i, ":\n")
    eqtl <- qtls[[`i`]]
    eqtl$cis.or.trans <- as.character(eqtl$cis.or.trans)
    tab1 <- data.table(chromstates[[`i`]])
    tab3 <- data.table(unique(eqtl[!is.na(eqtl$cis.or.trans) & eqtl$cis.or.trans=="trans",
                                   c("RSID", "cis.or.trans")]))
    setkey(tab1, snp)
    setkey(tab3, RSID)

    data.trans <- data.frame(merge(tab1, tab3, by.x="snp",
                                   by.y="RSID", all.x=TRUE), stringsAsFactors=F)
    data.trans$cis.or.trans[is.na(data.trans$cis.or.trans)] <- "0"
    data.trans$cis.or.trans[ data.trans$cis.or.trans=="trans"] <- "1"
    data.trans$Tissue <- rep(i, nrow(data.trans))
    res.trans <- rbind(res.trans,  data.trans)
}
save(res.trans, file=paste0(epi.dir, cistrans.res.file))

### Enrichment in trans-eQTLs among chromatin states
trans.epi <- list()
for(ca in sort(unique(res.trans$region))){
    cat("Running chromstats for category", ca, "\n")
    data.tmp.epi  = res.trans[!is.na(res.trans$region),]
    epi <- make.table(data = data.tmp.epi,
                      categ=ca, colnm="region")
    dat <- data.frame("chrom"=epi, "trans"=as.numeric(data.tmp.epi$cis.or.trans), "tissue"=data.tmp.epi$Tissue)
    trans.epi[[`ca`]] <- clogit(chrom ~ trans + strata(tissue),
                                data=dat, method="approximate")
}
save(trans.epi, file=paste0(epi.dir, cistrans.chromstate.res.file))


###Summarize trans-eQTL enrichment results
CI.trans.epi <- list()
trans.res.epi <- NULL
for(ca in names(trans.epi)){
    s.trans.epi<- summary(trans.epi[[`ca`]])$coefficient
    CI.trans.epi[[`ca`]] <- exp(confint(trans.epi[[`ca`]]))
    trans.res.epi <- rbind(trans.res.epi, rbind(c(ca, "core-score", s.trans.epi[1,2], s.trans.epi[1,5],
                                                  CI.trans.epi[[`ca`]][1,1], CI.trans.epi[[`ca`]][1,2])))
}
colnames(trans.res.epi) <- c("categ", "centrality", "OR", "pval", "CI.2.5", "CI.97.5")
write.table(trans.res.epi, file=paste0(epi.dir, chromatin.trans.file), quote=F, sep="\t")

