### plot_GO_results.R

### Load variables
source("code/variables_definition.R")

### Set parameters
snps.file <- paste0( 'all_tissues_snps_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
genes.file <- paste0( 'all_tissues_genes_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
edges.file <- paste0( 'all_tissues_edges_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')
go.results.txt <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_gccbkg.txt')
go.results.modif <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_final.txt')
go.results.bp <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_bp.txt')
go.results.mf <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_mf.txt')
go.results.cc <- paste0('alltissues_fdr', FDRcis, FDRtrans, "_", window, 'MB_go_results_cc.txt')
heatmap.pdffile <- paste0("alltissues_fdr", FDRcis, FDRtrans, "_", window, "MB_go_heatmap_gccbkg.pdf")
heatmapsub.pdffile <- paste0("alltissues_fdr", FDRcis, FDRtrans, "_", window, "_go_subheatmap_gccbkg.pdf")
clusters.tissues <- paste0( "alltissues_fdr", FDRcis, FDRtrans, "_", window, "_go_subheatmap_gccbkg_tissue_clusters.txt")
clusters.GO <- paste0("alltissues_fdr", FDRcis, FDRtrans, "_", window, "_go_subheatmap_gccbkg_GO_clusters.txt")
uniqueGO.file <-paste0("alltissues_fdr", FDRcis, FDRtrans, "_", window, "_go_gccbkg_uniqueGO_proportion.txt")
allGO.file <-paste0("alltissues_fdr", FDRcis, FDRtrans, "_", window, "_go_gccbkg_allGO_proportion.txt")
bubble.plot.unique.top10 <- paste0("Gene_Ontology_", FDRcis, FDRtrans, "_", window, "_top10_")
bubble.plot.all.top10 <- paste0("Gene_Ontology_", FDRcis, FDRtrans, "_", window, "_top10_all_")

### Load libraries
library(ggplot2)
library(gplots)
library(plyr)
library(RColorBrewer)
library(dendextend)

### Load data
go <- read.table(paste0(cluster.dir, go.results.txt), header=T, stringsAsFactors=F, sep="\t", quote="")
load( tissue.file)
load(paste0(cluster.dir, snps.file))
load(paste0(cluster.dir, genes.file))
load(paste0(cluster.dir, edges.file))

### Functions
### Count number of groups in the dendrogram.
count.group <- function(dendro){
  n <- 1
  flag <- 0
  while(flag==0){
    if( length(dendro)>1 ){
      dendro <- dendro[[2]]
      n <- n+1
    } else { flag <- 1 }
  }
  return(n)
}

### Calculate the log P-values of GO results and filter them using a P-value threshold
calc.logpval <- function(go){
    corres.term <- unique(go[,c("GOID", "Term")])
    list.com <- unique(paste(go$Tissue, go$comm.id, sep="_"))
    goterm <- unique(go$GOID)
    
    mat <- matrix(0, ncol=length(list.com), nrow=length(goterm))
    colnames(mat) <- list.com
    rownames(mat) <- goterm
    for (n in 1:ncol(mat)){
        print(n)
        temp <- go$p.adjust[paste(go$Tissue, go$comm.id, sep="_") == colnames(mat)[n]]
        names(temp) <- go$GOID[paste(go$Tissue, go$comm.id, sep="_") == colnames(mat)[n]]
        mat[names(temp),n] <- (-log10(temp))
    }
    l <- apply(mat, 1, function(x){sum(x>(-log10(.05)))})
    m <- apply(mat, 2, function(x){sum(x>(-log10(.05)))})

    #mat[mat<0] <- 0

    mat=mat[l>4,m>0]
    return(list("data"=mat, "corres.key"=corres.term))
}


### Separate GO analysis results in 3 tables for  BP/CC/MF 

go$Tissue <- Tissues[go$Tissue, 1]
go <- go[go$p.adjust<=0.05 & go$OddsRatio>1,]
go <- go[,c(1,2,4,5,11,3,6,8:10,7)]
go.bp <- go[go$Category=="BP",]
go.mf <- go[go$Category=="MF",]
go.cc <- go[go$Category=="CC",]
write.table(go, file=paste0(cluster.dir, go.results.modif), quote=F, sep="\t", row.names=F)
write.table(go.bp, file=paste0(cluster.dir, go.results.bp), quote=F, sep="\t", row.names=F)
write.table(go.mf, file=paste0(cluster.dir, go.results.mf), quote=F, sep="\t", row.names=F)
write.table(go.cc, file=paste0(cluster.dir, go.results.cc), quote=F, sep="\t", row.names=F)

### Search for unique and shared terms in go.bp
goid.tmp <- tapply(go.bp$GOID, go.bp$Tissue, function(x){unique(x)})
tmp <- lapply(names(goid.tmp), function(x,y){cbind(rep(x, length(y[[`x`]])), as.character(y[[`x`]]))}, goid.tmp)
goid <- NULL
for(i in 1:length(tmp)){
    goid <- rbind(goid, tmp[[i]])
}
list.go <- table(goid[,2])
GO.unique <- list.go[list.go==1]
GO.all <- list.go[list.go>=12]
go.bp.tmp <- go.bp
go.bp.tmp$go.counts <- list.go[go.bp$GOID]
go.bp.tmp$class <- paste0(go.bp.tmp$Tissue, go.bp.tmp$comm.id)
summary.go.unique <- tapply(go.bp.tmp$go.counts, go.bp.tmp$class, function(x){sum(x==1)/length(x)})
summary.go.all <- tapply(go.bp.tmp$go.counts, go.bp.tmp$class, function(x){sum(x==13)/length(x)})
go.bp.tmp$proportion.unique <- summary.go.unique[go.bp.tmp$class]
go.bp.tmp$proportion.all <- summary.go.all[go.bp.tmp$class]
res.unique <- go.bp.tmp[go.bp.tmp$go.counts==1,c(1:10,13)]
res.all <- go.bp.tmp[go.bp.tmp$go.counts==13,c(1:10,13)]
res.all$names <- paste(res.all$Tissue, res.all$comm.id, sep="_")
res.unique$names <- paste(res.unique$Tissue, res.unique$comm.id, sep="_")
s.unique <- tapply(rep(1, nrow(res.unique)), res.unique$names, sum)
s.all <- tapply(rep(1, nrow(res.all)), res.all$names, sum)
res.all$counts <- s.all[res.all$names]
res.unique$counts <- s.unique[res.unique$names]
write.table(res.unique, file=paste0(cluster.dir, uniqueGO.file), row.names=F, col.names=T, quote=F, sep="\t")
write.table(res.all, file=paste0(cluster.dir, allGO.file), row.names=F, col.names=T, quote=F, sep="\t")
write(names(GO.all), file=paste0(cluster.dir, "list_common_GOTerms.txt"), ncol=1)

### Make bubble plots for unique terms
res <- res.unique
res$logp <- -log10(res$p.adjust)
res$Term <- factor(res$Term, rev(as.character(res$Term)))
res$class <- paste(res$Tissue, res$comm.id)
res$OddsRatio[res$OddsRatio==Inf] <- max(res$OddsRatio[res$OddsRatio!=Inf])
res.unique$OddsRatio[res.unique$OddsRatio==Inf] <- max(res.unique$OddsRatio[res.unique$OddsRatio!=Inf])
### For each community separately
for(cl in c("Heart left ventricle 86", "Heart left ventricle 30", "Muscle skeletal 75", "Esophagus muscularis 83")){
    cat("Running bubble plot for", cl, "\n")
    ## top10 terms
    pdf(paste0(figure.dir, bubble.plot.unique.top10,
                  gsub(" ", "_", cl), ".pdf"), h=11, w=8)
    res.tmp <- res.tmp[order(res.tmp$p.adjust),]
    res.tmp <- res.tmp[1:min(c(10, nrow(res.tmp))),]
    g<-ggplot(res.tmp, aes(x=logp, y=Term, size=OddsRatio))
    print(g+geom_point( shape=21, colour="blue", fill="dodgerblue3")+
    scale_size(range = c(min(res.tmp$OddsRatio)*500/max(res.unique$OddsRatio), max(res.tmp$OddsRatio)*500/max(res.unique$OddsRatio)))+
    labs(
    x = "-log10(FDR qvalue)",
    y = "",
    size = "Odds Ratio"
    )+
    theme_bw()+
    xlim(0, max(res$logp))+
    theme(axis.text.y = element_text(colour="black", size = rel(0.8)))+
    theme(plot.margin=unit(c(1, 0, 2, 3), unit="pt"), legend.position="none") # test
    )
    dev.off()
}

### For all community together to get the legend
pdf(paste0(figure.dir, bubble.plot.all,".pdf"), h=11, w=8)
cl="Heart left ventricle 86"
res.tmp <- res[res$class==cl,]
g<-ggplot(res.tmp, aes(x=logp, y=Term, size=OddsRatio))
print(g+geom_point( shape=21, colour="blue", fill="dodgerblue3")+
guides(fill=guide_legend(title.theme = element_text(size = 15)))+
scale_size_continuous(range = c(min(res.tmp$OddsRatio)*500/max(res$OddsRatio),
max(res.tmp$OddsRatio)*500/max(res$OddsRatio)),
breaks = c(2, seq( 10, 40, 10)))+
labs(
x = "-log10(FDR qvalue)",
y = "",
size = "Odds Ratio"
)+
theme_bw()+
xlim(0, max(res$logp))+
theme(axis.text.y = element_text(colour="black", size = rel(0.5)))+
theme(plot.margin=unit(c(1, 0, 2, 0), unit="pt"),
legend.key.width=unit(2, unit="pt"),
legend.key.height=unit(8, unit="pt"),
legend.key = element_blank()) # test
)
dev.off()

### Plot Gene Ontology heatmap

mycol2 <- colorRampPalette(c("lightblue1", "blue4"))(20)
data.bp <- calc.logpval(go.bp[go.bp$GOID %in% names(list.go)[list.go>11],])
data.bp$data[data.bp$data > 10] <- 20
pdf(paste0(figure.dir, heatmap.pdffile), paper='a4')
H.bp <- heatmap.2(data.bp$data,col=mycol2, cexRow=0.1, cexCol=0.35, trace='none',
                  breaks=seq(0,10,0.5), key.title=NA, key.xlab="-log10(p)",
                    key.par=list("mar"=c(4,1,4,0)+0.1), key.ylab=NA, keysize=1,  density.info='none',
                    main="Biological Processes")
dev.off()
