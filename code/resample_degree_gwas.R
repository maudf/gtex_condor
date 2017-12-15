### resample_degree_gwas.R

###Load variables
source("code/variables_definition.R")

### Set parameters
nsam <- 10000
breaks <- c(1:19, 24, 29, 39)
mids <- c(1:19, 22, 27, 34.5, 40)

topo.file <- paste0("network_topology_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")
snps.neigh.file <- "nb_snps_neighbors.RData"
community.file <- paste0( 'all_tissues_communities_fdr', FDRcis, FDRtrans, "_", window, 'MB.Rdata')

gwas.degree.pvalues <-  paste0("gwas_degree_pvalues_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
gwas.degree.pvalues.txt <-  paste0("gwas_degree_pvalues_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
gwas.degree.or.txt <-  paste0("gwas_degree_oddsratio_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
obs.exp.degree.pdf <- paste0("ObsVsExp_GWAS_SNPs_degree_plot_fdr", FDRcis, FDRtrans, "_", window, "MB.pdf")

term.file <- paste0("gwas_", term, "_terms.txt")
gwas.term.degree.pvalues <-  paste0("gwas_degree_whole_blood_pvalues_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")
gwas.term.degree.pvalues.txt <-  paste0("gwas_degree_whole_blood_pvalues_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
obs.exp.degree.term.pdf <- paste0("ObsVsExp_GWAS_SNPs_degree_plot_whole_blood_gwasautoimmune_fdr", FDRcis, FDRtrans, "_", window, "MB.pdf")

### Load data

load(paste0(cluster.dir, topo.file))
load(paste0(eqtl.dir, snps.neigh.file))
load(paste0( cluster.dir, community.file))
load(tissue.file)

gwas_data <- read.delim( gwas.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gwas_data <- gwas_data[gwas_data$PVALUE_MLOG>=8,]
gwas_data <- gwas_data[!is.na(gwas_data$PVALUE_MLOG),]
gwas.term <- read.delim(paste0(data.dir, term.file), header=F, stringsAsFactors=F)
gwas.diseases<-read.delim(paste0(gwas.dir, gwas.diseases.file), header=F, stringsAsFactors=F)
gwas.spe <- unique(gwas.diseases[(gwas.diseases[,2] %in% gwas.brain.term[,1]) & (gwas.diseases[,1] %in% gwas_data$SNPS), 1 ])

### Functions

plot.ticks <-  function(x, ax){
    y1 <- floor(log10(range(x)))
    pow <- seq(y1[1],y1[2]+1)
    labels <- sapply(pow,function(i){ ifelse(i>=0, 10^i, paste("10", i, sep="")) } )
    axis(ax, 10^pow, labels=labels, tck=0.04, hadj=0.5)
    y.minor.ticks <- c()
    for ( i in pow ) { y.minor.ticks <- c(y.minor.ticks, 1:9*10^(i)) }
    axis(ax, y.minor.ticks, labels=FALSE, tck=0.02, padj=-0.5)
}

## Resample gwas distribution taking into account number of neighboring genes
resample.distrib <- function(gwas, all,  distrib, nsam, central, breaks){
    
    ## Find probability law of number of SNPs in 1Mb window for GWAS SNPs
    tmp.nb.rs <- distrib[( distrib$rs_id %in% unique(all) ),]
    network.distrib <- sort(unique(tmp.nb.rs[,2]))
    gwas.distrib <- table(distrib[( distrib$rs_id %in% gwas),2])
    law.distrib <- rep(0, length(network.distrib))
    names(law.distrib) <- network.distrib
    law.distrib[names(gwas.distrib)] <- as.numeric(gwas.distrib)
    law.distrib[law.distrib==0] <- 0.001

    ## Resample list of GWAS SNPs taking into account number of neighboring SNPs 
    s <- lapply(1:nsam, function(x, n, gwas, law){
        return(sample(n, gwas, prob=law, replace=T))},
                n=as.numeric(network.distrib), gwas=sum(gwas.distrib), law=law.distrib)

    sampling.pool <- tapply(tmp.nb.rs$rs_id, tmp.nb.rs$N, function(x){x})
    all.samples <- lapply(1:length(s), function(index, nb.sample, list.snp){
        tab <- table(nb.sample[[index]]) 
        resampling <- lapply(names(tab), function(x, y, z){sample(z[[`x`]], y[x])}, y=tab, z=list.snp)
        names(resampling) <- names(tab)
        return(unlist(resampling))
    }, nb.sample=s, list.snp=sampling.pool)
    deg <- lapply(all.samples, function(sam, top){top[rownames(top) %in% sam, "degree"]}, central)

    ## Obtain distribution of centrality scores for all
    all.deg <- central[rownames(central) %in% unique(all), "degree"]
    b <- c(0, breaks,  max(all.deg))
    h.all <- hist(all.deg, breaks=b, plot=F)
    distrib.all <- h.all$counts/sum(h.all$counts)

    ## Obtain distribution of centrality scores for gwas
    gwas.deg <- central[rownames(central) %in% gwas, "degree"]
    h.gwas <- hist(gwas.deg, breaks=b, plot=F)
    distrib.gwas <- h.gwas$counts/sum(h.gwas$counts)

    ## Obtain Ob/Exp for real and resampled GWAS.
    obs.exp.true <- distrib.gwas/distrib.all 

    obs.exp.sample <- matrix(unlist(lapply(deg, function(x, distrib.all, breaks){
        h.tmp <- hist(x, breaks=breaks, plot=F)
        distrib.tmp <- h.tmp$counts/sum(h.tmp$counts)
        return(distrib.tmp/distrib.all)
    }, distrib.all=distrib.all, breaks=b)
                                                      ), nrow=length(deg), byrow=T)
                                                      
    perc <- unlist(
        lapply(1:ncol(obs.exp.sample),
               function(x, obs, sam){
                   return(sum(sam[,x]>obs[x])/length(sam[,x]))
               },
               obs=obs.exp.true, sam=obs.exp.sample)
    )
    return(list("deg"=deg, "obs.exp.true"=obs.exp.true, "obs.exp.sample"=obs.exp.sample, "perc"=perc))
}

## Get P-values for observed/expected gwas in each degree category

get.pval <- function(communities, tissue_name, gwas.list, topo, nb.rs, nsam, breaks){
    colnames(communities[[`tissue_name`]]$edges) <- c("red", "blue")
    gwas.tissue <- unique(gwas.list[gwas.list %in% unique(communities[[`tissue_name`]]$edges$red)])
    res <- resample.distrib(gwas=gwas.tissue, all=communities[[`tissue_name`]]$edges$red,
                            distrib=nb.rs, nsam=nsam, central=topo[[`tissue_name`]], breaks=breaks)
    return(res)
}



###Obtain Expected/Observed nb of GWAS SNPs/degree.
### Shuffle GWAS SNPs and keep distribution of number of neighbors to obtain resampled p-values

## For all tissues and all gwas
degree.sample <- list()
obs.exp <- list()
obs.exp.resampled <- list()
pval <- list()
for(tissue_name in names(communities)){
    print(tissue_name)
    res.tmp <- get.pval(communities=communities, tissue_name=tissue_name,
                        gwas.list=gwas_data$SNPS, topo=topo, nb.rs=nb.rs, nsam=nsam, breaks=breaks)
    degree.sample[[`tissue_name`]] <- res.tmp$deg
    obs.exp[[`tissue_name`]] <- res.tmp$obs.exp.true
    obs.exp.resampled[[`tissue_name`]] <- res.tmp$obs.exp.sample
    pval[[`tissue_name`]] <- res.tmp$perc
}
save(degree.sample, obs.exp, obs.exp.resampled, pval,   file=paste0(gwas.dir, gwas.degree.pvalues))

### Plot results
colors.scale <- colorRampPalette(c("blue", "white", "red"))
colors.all <- colors.scale(7)

pdf(paste0(general.path, figure.path, obs.exp.degree.pdf), width=8, height=11)
par(mfrow=c(5,3))
par(mar=c(4,5,3,1)+0.1)
for(tissue_name in names(obs.exp)[nb.samples>=200]){
    print(tissue_name)
    obs.exp.true <- obs.exp[[`tissue_name`]]
    perc <- pval[[`tissue_name`]]
    colors.tmp <- rep(colors.all[4], length(perc))
    colors.tmp[!is.na(perc) & perc >= 0.95 & obs.exp.true < 1] <- colors.all[3]
    colors.tmp[!is.na(perc) & perc >= 0.99 & obs.exp.true < 1] <- colors.all[2]
    colors.tmp[!is.na(perc) & perc >= 0.999 & obs.exp.true < 1] <- colors.all[1]
    colors.tmp[!is.na(perc) & perc <= 0.05 & obs.exp.true > 1] <- colors.all[5]
    colors.tmp[!is.na(perc) & perc <= 0.01 & obs.exp.true > 1] <- colors.all[6]
    colors.tmp[!is.na(perc) & perc <= 0.001 & obs.exp.true > 1] <- colors.all[7]

    plot(mids, obs.exp.true,
         xlim=c(1,45), ylim=c(0, 12), log="x", xaxt="n",
         main=Tissues[tissue_name,1], cex.main=1,
         xlab="Degree", ylab="GWAS-SNPs (Obs/Exp)",
         col="black", bg= colors.tmp, pch=23)
    plot.ticks(mids, 1)
}
par(xpd=NA)
plot(c(0,10), c(0,12), type='n', frame.plot=F, axes=F, xlab="", ylab="")
step = 12/4
y0 = 12
labels <- c( "< 0.001", "< 0.01", "< 0.05", "")
for(i in 1:4){
    polygon(c(6, 7, 7, 6), c(y0,y0,y0-step,y0-step), col=colors.all[i], border=NA)
    polygon(c(4, 5, 5, 4), c(y0,y0,y0-step,y0-step), col=colors.all[(8-i)], border=NA)
    text(7, y0-step/2, labels=labels[i], pos=4)
    y0 = y0-step
}
text(7, 13, labels="P-values", pos=4)
polygon(c(6, 7, 7, 6), c(0,0,12,12))
polygon(c(4, 5, 5, 4), c(0,0,12,12))
text(2, 0, labels="Obs/Exp", pos=1)
text(4.5, 0, labels="> 1", pos=1)
text(6.5, 0, labels="< 1", pos=1)
par(xpd=F)
dev.off()

### Write table with enrichments values
obs.exp.all <- matrix(unlist(obs.exp), ncol=length(obs.exp[[1]]), byrow=T)
pval.all <- matrix(unlist(pval), ncol=length(pval[[1]]), byrow=T)
rownames(pval.all) <- rownames(obs.exp.all) <- Tissues[names(topo),1]
colnames(pval.all) <- colnames(obs.exp.all) <- c(as.character(1:19), "20-24", "25-29", "30-40", "> 40")
write.table(pval.all, file=paste0(gwas.dir, gwas.degree.pvalues.txt), quote=F, sep="\t")
write.table(obs.exp.all, file=paste0(gwas.dir, gwas.degree.or.txt), quote=F, sep="\t")

## For term-related gwas in a specific tissue

res.term <- get.pval(communities=communities, tissue_name=tissue.spe, gwas.list=gwas.term,
                      topo=topo, nb.rs=nb.rs, nsam=nsam, breaks=breaks)
save(res.term,  file=paste0(gwas.dir, gwas.term.degree.pvalues))

### Plot results
colors.scale <- colorRampPalette(c("blue", "white", "red"))
colors.all <- colors.scale(7)

tissue_name=tissue.spe
pdf(paste0(figure.dir, obs.exp.degree.term.pdf), width=8, height=5)
par(mar=c(4,5,3,1)+0.1)
obs.exp.true <- res.term$obs.exp.true
perc <- res.term$perc
#mids <- res.brain$mids.coord
colors.tmp <- rep(colors.all[4], length(perc))
colors.tmp[!is.na(perc) & perc >= 0.95 & obs.exp.true < 1] <- colors.all[3]
colors.tmp[!is.na(perc) & perc >= 0.99 & obs.exp.true < 1] <- colors.all[2]
colors.tmp[!is.na(perc) & perc >= 0.999 & obs.exp.true < 1] <- colors.all[1]
colors.tmp[!is.na(perc) & perc <= 0.05 & obs.exp.true > 1] <- colors.all[5]
colors.tmp[!is.na(perc) & perc <= 0.01 & obs.exp.true > 1] <- colors.all[6]
colors.tmp[!is.na(perc) & perc <= 0.001 & obs.exp.true > 1] <- colors.all[7]

plot(mids, obs.exp.true,
     xlim=c(1,45), ylim=c(0, 12), log="x", xaxt="n",
     main=Tissues[tissue_name,1], cex.main=1, cex=2,
     xlab="Degree", ylab="GWAS-SNPs (Obs/Exp)",
     col="black", bg= colors.tmp, pch=23)
plot.ticks(mids, 1)
dev.off()

### Write list of p-values for autoimmune diseases in whole blood
pval.term <- matrix(res.term$perc, nrow=1, byrow=T)
rownames(pval.term) <- "Whole Blood"
colnames(pval.term) <- c(as.character(1:19), "20-24", "25-29", "30-40", "40")
write.table(pval.term, file=paste0(gwas.dir, gwas.term.pvalues.txt), quote=F, sep="\t")
