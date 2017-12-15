### plot_resampling_shared_GO_pvalues.R

### Load variables
source("code/variables_definition.R")
args = commandArgs(trailingOnly=TRUE)

### Set variables
genes.file <- paste0("all_tissues_genes_fdr", FDRcis, FDRtrans,
                           "_", window, "MB.Rdata")
go.results.txt <- paste0('alltissues_fdr', FDRcis, FDRtrans,
                         "_", window, 'MB_go_results_bp.txt')
 
### Load data
load(paste0(cluster.dir, genes.file))
load(tissue.file)

all.tis <- names(genes)
tislist<-rownames(Tissues)
names(tislist)<-as.character(Tissues[,1])
tislist <- tislist[tislist %in% all.tis]

go.bp <- read.table(paste0(cluster.dir, go.results.txt), header=T, sep="\t", stringsAsFactors=F)
a <- go.bp[go.bp$Tissue %in% names(tislist)  &
         go.bp$p.adjust<=0.05 & go.bp$OddsRatio>1 &
         go.bp$NumberTissues>=12,]
shared.com <- tapply(a$com, a$Tissue, unique)
names(shared.com) <- tislist[names(shared.com)]
a$Tissue <- tislist[a$Tissue]


### Functions

compute.resample.pval <- function(x, resample){
    tis <- x[1]
    com <- x[2]
    goid <- x[3]
    pval.real <- as.numeric(x[5])
    iter <- as.numeric(x[13])
    pval.res <- paste0("<", 1/iter)
    if (goid %in% names(resample[[`tis`]][[`com`]])){
        tmp <- resample[[`tis`]][[`com`]][[`goid`]]
        if (length(tmp)<iter){
            tmp <- c(tmp, rep(1, (iter-length(tmp))))
        }
        s <- sum(tmp<=pval.real)
        if(s>0){pval.res <- s/iter}
    }
    return(pval.res)
}

### read resampling files and store p-values.
resample <- list()
for(i in 1:length(tislist)){
    tis <- tislist[i]
    print(tis)
    resample.files <- list.files(resampling.GO.dir,
                                 pattern=paste0(tislist[i], "_resampling_GO_shared_results_*"))
    tmp <- NULL
    for(f in resample.files){
        print(f)
        tmp <- rbind(tmp, read.table(paste0(resampling.GO.dir, f),
                                     header=T, sep="\t", stringsAsFactors=F, quote=""))
    }
    tmp$comm.id <- shared.com[[`tis`]][tmp$comm.id]

    resample[[`tis`]] <- tapply(1:nrow(tmp), tmp$comm.id, function(x, tab){
        res <- tapply(tab$Pvalue[x], tab$GOID[x], function(y){y})
    }, tab=tmp)
    rm(tmp)
}

### Compute real p-value rank
a$iter <- 1000
resample.pval <- apply(a, 1, compute.resample.pval, resample=resample)
go.bp$resample.pval <- NA
go.bp$resample.pval[which(go.bp$Tissue %in% names(tislist)  &
         go.bp$p.adjust<=0.05 & go.bp$OddsRatio>1 &
         go.bp$NumberTissues>=12)] <- resample.pval
write.table(go.bp, file=paste0(cluster.dir, go.results.txt), row.names=F, sep="\t", quote=F)

### Plot distribution of resampled p-values for term GO:0010468
terms <- "GO:0010468"
b<- a[grep(terms, a$GOID),]
d<-tapply(b$Pvalue, b$Tissue, function(x){which(x==min(x))})
b.tmp <- NULL
for(i in 1:length(d)){
    b.tmp <- rbind(b.tmp, b[b$Tissue==names(d)[i],][d[i],])
}

pdf(paste0(figure.dir, "hist_GO_terms_resampling_", terms, ".pdf"),
    width=8, height=11.5)
par(mfrow=c(4,4), mar=c(4,5,3,0)+.1)
for(i in 1:nrow(b.tmp)){
    tis <- b.tmp[i,1]
    com <- as.character(b.tmp[i,2])
    goid <- b.tmp[i,3]
    iter <- as.numeric(b.tmp[i,13])
    tmp <- resample[[`tis`]][[`com`]][[`goid`]]
    if (length(tmp)<iter){
        tmp <- c(tmp, rep(1, (iter-length(tmp))))
    }
    h <- hist(tmp, plot=F, breaks=seq(0,1,0.01))
    h$density <- h$counts/sum(h$counts)
    plot(h, col="dodgerblue", xlab="P-values", ylab="Density", main=as.character(Tissues[tis,1]), freq=F)
    par(xpd=T)
    text(0.5,(max(h$density)-0)*21/20, label=paste0("P = ", sprintf("%.2e", b$Pvalue[i])), font=2)
}
dev.off()
