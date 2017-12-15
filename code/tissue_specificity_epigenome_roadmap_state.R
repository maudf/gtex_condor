### tissue_specificity_epigenome_roadmap_state.R

###Load variables
source("code/variables_definition.R")
library(metap)

### Set parameters
snps.file <- paste0( 'all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
snps.mult.file <- paste0( 'snps_multiplicity_fdr', FDRcis, FDRtrans, "_", window, 'MB_gcc.Rdata')
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

load(paste0(data.dir, epi.file))
load(paste0(cluster.dir, summary.stats.mult))

### Functions


## Create multiplicity info
find.epigenomic.multiplicity <- function(x, m, epi){
    l <- m[names(m) %in% x]
    re <- data.frame(as.numeric(l), row.names=names(l))
    
    for(k in names(epi)){
        print(k)
        epi.tmp <- epi[[`k`]][rownames(epi[[`k`]]) %in% names(l),"region"]
        names(epi.tmp) <- rownames(epi[[`k`]])[rownames(epi[[`k`]]) %in% names(l)]
        re <- data.frame(re, epi.tmp[rownames(re)], stringsAsFactors=F)
    }
    colnames(re)=c("mult", names(epi))
    return(re)
}

## Compute epigenetic states transition tables
make.transition.table <- function(mult, LD=NULL){
    res <- lapply(names(mult), function(n, epi, LD){
        if(!is.null(LD)){
            require(data.table)
            tmp <- data.table(data.frame("SNP"=rownames(epi[[`n`]]), epi[[`n`]], stringsAsFactors=F))
            setkey(tmp, SNP)
            data.merge <- data.frame(merge(tmp, LD, by='SNP'), stringsAsFactors=F)[,1:11]
            re <- as.character(tapply(data.merge$SNP, apply(data.merge[,c(2,which(colnames(data.merge)==n),11)], 1, paste, collapse="-"), function(x){sample(x, 1)}))
        } else {
            re <- rownames(epi[[`n`]])
        }
        epi[[`n`]] <- epi[[`n`]][rownames(epi[[`n`]]) %in% re,]
        a <- t(apply(epi[[`n`]], 1, function(y, x, z){
            i <- which(z %in% x)
            c(y[1], paste(y[(2:length(y))[!(2:length(y) %in% i)]], y[i], sep="->"))
        }, x=n, z=colnames(epi[[`n`]])))
        colnames(a) <- colnames(epi[[`n`]])[!(colnames(epi[[`n`]]) %in% n)]
        return(a)
    }, epi=mult, LD=LD)
    names(res) <- names(mult)
    return(res)
}

make.transition.table.com <- function(mult.transition, elt, LD=NULL){
    res <- lapply(names(mult.transition), function(x, snps, epi, LD){
        lapply(snps[[`x`]], function(y, x, z, LD){
            a <- z[ (rownames(z) %in% y),]
            if(!is.null(LD)){
                require(data.table)
                tmp <- data.table(data.frame("SNP"=rownames(a), a, stringsAsFactors=F))
                setkey(tmp, SNP)
                data.merge <- data.frame(merge(tmp, LD, by='SNP'), stringsAsFactors=F)[,1:10]
                re <- as.character(tapply(data.merge$SNP,
                                          apply(cbind(data.merge[,c(2,10)],
                                                      gsub("->.*", "", data.merge[,3])), 1, paste, collapse="-"),
                                          function(x){sample(x, 1)}))
                a <- a[ (rownames(a) %in% re),]
            }
            return(a)
        }, x=x, z=epi[[`x`]], LD=LD)
    }, snps=elt, epi=mult.transition, LD=LD)
    names(res) <- names(mult.transition)
    return(res)
}

## Summarize epigenetic states transition tables
summarize.transition <- function(u, identical, activated, repressed, similar){
    lapply(u, function(x, identical, activated, repressed, similar){
        a <- list()
        for(n in colnames(x[,2:ncol(x)])){
            f <- factor(gsub(" ", "", x[,1]), levels=as.character(1:13))
            b <- matrix(unlist(lapply(tapply(x[,n], f, function(z,  identical, activated, repressed, similar){
                res <- c(sum(z %in% identical),
                         sum(z %in% similar),
                         sum(z %in% activated),
                         sum(z %in% repressed))
                return(res)
            }, identical=identical, activated=activated, repressed=repressed, similar=similar
            ), function(y){
                if(is.null(y)){
                    y <- rep(0,4)
                }
                return(y)
            })), ncol=4, byrow=T)
            colnames(b) <- c("id", "si", "ac", "re")
            
            w <- which(table(f)==0)
            
            a[[`n`]] <- b
        }
        return(a)
    }, identical=identical, activated=activated, repressed=repressed, similar=similar)
}

summarize.transition.bytype <- function(summary.transition, w){
    res <- lapply(names(summary.transition), function(n, u, w){
        l <- NULL
        for(i in w[[`n`]]){
            if(is.null(l)){
                l <- u[[`n`]][[i]]
            } else {
                r <- lapply(names(l), function(x, v, z, j){
                    z[[`x`]]+v[[j]][[`x`]]
                }, v=u[[`n`]], z=l, j=i)
                names(r) <- names(l)
                l <- r
            }
        }
        return(l)
    }, u=summary.transition, w=w)
    names(res) <- names(summary.transition)
    return(res)
}

sum.transition.acrosstissues <- function(transition){
    l <- transition[[1]]
    for(i in 2:length(transition)){
        l <- l + transition[[i]]
    }
    return(l)
}

## Extract contingency table for Odd ratio computation
extract.table <- function(x, y, u){
    data <- c(sum(x[u,3]),
              sum(x[u, c(1:2,4)]),
              sum(y[u,3]),
              sum(y[u, c(1:2,4)]))
    tab <- matrix(data, ncol=2)
    return(tab)
}

### Extract Epigenomic roadmap tissues corresponding to available tissues
epigenomic.roadmap <- list()
for(i in names(epigenomic.rodamap)){
    if(i %in% names(genes)){
        epigenomic.roadmap[[`i`]] <- epigenomic.rodamap[[`i`]]
    }
}

### Find list of epigenetic states
list.states <- sort(unique(epigenomic.roadmap$adipose_subcutaneous$region))
active <- list.states[c(1,8:14)]
silent <- list.states[c(5:7,15)]
bivalent <- list.states[c(2:4)]

### Find list of epigenetic states transition
comb.list <- t(combn(list.states, m=2))
comb.list <- rbind(comb.list, cbind(comb.list[,2], comb.list[,1]))
list.transition <- paste(comb.list[,1], comb.list[,2], sep="->")
identical <- paste(list.states, list.states, sep="->")
activated <- list.transition[which(((comb.list[,1] %in% silent) | (comb.list[,1] %in% bivalent)) &
                                   (comb.list[,2] %in% active))]
activated <- c(activated, list.transition[which((comb.list[,1] %in% silent) &
                                                (comb.list[,2] %in% bivalent))])
repressed <- list.transition[which(((comb.list[,1] %in% active) | (comb.list[,1] %in% bivalent)) &
                                   (comb.list[,2] %in% silent))]
repressed <- c(repressed, list.transition[which((comb.list[,1] %in% active) &
                                                (comb.list[,2] %in% bivalent))])
similar <- c(list.transition[which((comb.list[,1] %in% silent) & (comb.list[,2] %in% silent))],
             list.transition[which((comb.list[,1] %in% bivalent) & (comb.list[,2] %in% bivalent))],
             list.transition[which((comb.list[,1] %in% active) & (comb.list[,2] %in% active))])


### Find epigenetic state of SNPs in each tissue
epi.mult<-list()
for(n in names(epigenomic.roadmap)){
    print(n)
    epi.mult[[`n`]] <- find.epigenomic.multiplicity(x=unlist(snps[[`n`]]), m=s, epi=epigenomic.roadmap)
}
save(epi.mult, file=paste0(eqtl.dir, "location_snps_multiplicity_gcc.RData"))

### Find transitions of epigenetic state of SNPs between pairs of tissue

## Across all SNPs
epi.mult.transition <- make.transition.table(epi.mult)
summary.transition.all <- summarize.transition(u=epi.mult.transition, identical=identical, activated=activated, repressed=repressed, similar=similar)

epi.mult.transition.LD <- make.transition.table(epi.mult, LD=block.info)
summary.transition.all.LD <- summarize.transition(u=epi.mult.transition.LD, identical=identical, activated=activated, repressed=repressed, similar=similar)

## Split result by community
epi.mult.transition.LD.bycom <- make.transition.table.com(epi.mult.transition, elt=snps, LD=block.info)
summary.transition.com.LD <- lapply(epi.mult.transition.LD.bycom,
                                    function(u, identical, activated, repressed, similar){
                                        summarize.transition(u=u, identical=identical, activated=activated,
                                                             repressed=repressed, similar=similar)
                                    }, identical=identical, activated=activated,
                                    repressed=repressed, similar=similar)

## Summarize transitions of epigenetic states in communities enriched for unique or shared ontology terms
summary.transition.unique.LD <- summarize.transition.bytype(summary.transition=summary.transition.com.LD, w=which.go.unique.2)
summary.transition.common.LD <- summarize.transition.bytype(summary.transition=summary.transition.com.LD, w=which.go.common.2)

## Summarize transitions of epigenetic states across tissues for all, enriched in unique and shared ontology terms
summary.transition.all.acrosstissues.LD <- lapply(summary.transition.all.LD, sum.transition.acrosstissues)
summary.transition.unique.acrosstissues.LD <- lapply(summary.transition.unique.LD, sum.transition.acrosstissues)
summary.transition.common.acrosstissues.LD <- lapply(summary.transition.common.LD, sum.transition.acrosstissues)


save(summary.transition.all.acrosstissues.LD, summary.transition.unique.acrosstissues.LD, summary.transition.common.acrosstissues.LD,
     file=paste0(cluster.dir, "transition_state_all_common_unique.RData"))

### Compute Odds Ratios for enrichement in activation among TS SNPs in TS communities compared to shared communities.
odds.ratio.all.LD <- NULL
odds.ratio.common.LD <- NULL
contingency.table.common <- NULL
contingency.table.all <- NULL
for(i in 1:length(summary.transition.unique.acrosstissues.LD)){
    t1 <- extract.table(x=summary.transition.unique.acrosstissues.LD[[i]], y=summary.transition.common.acrosstissues.LD[[i]], u=1)
    f1 <- fisher.test(t1)
    odds.ratio.common.LD <- rbind(odds.ratio.common.LD, c("OR"=f1$estimate, "pval"=f1$p.value))
    t2 <- extract.table(x=summary.transition.unique.acrosstissues.LD[[i]], y=summary.transition.all.acrosstissues.LD[[i]], u=1)
    t2 <- cbind(t2[,1], t2[,2]-t2[,1])
    f2 <- fisher.test(t2)
    odds.ratio.all.LD <- rbind(odds.ratio.all.LD, c("OR"=f2$estimate, "pval"=f2$p.value))
    
    contingency.table.common <- c(contingency.table.common, as.numeric(t1))
    contingency.table.all <- c(contingency.table.all, as.numeric(t2))

}

### Compute global Odd ratio across tissues using the cochran-mantel-haenszel test
common.glob.contingency<-array(contingency.table.common,
                               dim=c(2, 2, length(summary.transition.unique.acrosstissues.LD)),
                               dimnames=list("activation"=c("Y", "N"), "status"=c("unique", "common"), "tissue"=names(summary.transition.unique.acrosstissues.LD)))
all.glob.contingency<-array(contingency.table.all,
                               dim=c(2, 2, length(summary.transition.unique.acrosstissues.LD)),
                               dimnames=list("activation"=c("Y", "N"), "status"=c("unique", "common"), "tissue"=names(summary.transition.unique.acrosstissues.LD)))
cochran.common <- mantelhaen.test(common.glob.contingency)
cochran.all <- mantelhaen.test(all.glob.contingency)

### Write results
rownames(odds.ratio.all.LD) <- rownames(odds.ratio.common.LD) <- names(summary.transition.unique.acrosstissues.LD)
write.table(rbind(cbind(odds.ratio.all.LD, odds.ratio.common.LD), c(cochran.all$estimate, cochran.all$p.value, cochran.common$estimate, cochran.common$p.value)),
            file=paste0(cluster.dir, "epigenetic_state_unique_SNPs_odds_ratio_LD_gcc.txt"),
            quote=F, sep="\t")

### Plot Odds Ratios for enrichement in activation among TS SNPs in TS communisties.
pdf(paste0(figure.dir, "epigenetic_state_unique_SNPs_gcc_barplot_LD.pdf"), width=5, height=8)
par(las=1)
b <- barplot(odds.ratio.common.LD[8:1,1]-1, horiz=T, xlim=c(-0.5, 2.1), xlab="Odds Ratio",
             col="lightseagreen", names=as.character(Tissues[rownames(odds.ratio.common.LD),2])[8:1], xaxt='n')
axis(side=1, at=seq(-0.5, 1.6, 0.5), labels=seq(-0.5, 1.6, 0.5)+1)
pval <- odds.ratio.common.LD[8:1,2]
 f <- ifelse(pval <= 0.001, "***",
                 ifelse(pval <= 0.01, "**",
                 ifelse(pval <= 0.05, "*",
                        "NS")))
par(xpd=T)
text(rep(par("usr")[2], 3), b[,1], labels=f, cex=0.6, pos=4)
dev.off()
