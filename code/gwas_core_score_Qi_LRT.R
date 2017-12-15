### gwas_core_score_Qi_LRT.R

### Load data
source("code/variables_definition.R")

### Set parameters
term.file <- paste0("gwas_", term, "_terms.txt")
topo.file <- paste0("network_topology_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")
community.file <- paste0("all_tissues_communities_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")

LDblock.rdata <- gsub("txt", "RData", LDblock.file)
pdf.effect.size <- "gwas_core-score_effect_size.pdf"
gwas.corescore.LD.pdf <- paste0("cores_score_GWAS_filteredLD_all_fdr", FDRcis, FDRtrans, "_", window, "MB.pdf")
gwas.corescore.term.ld.pdf <- paste0("cores_score_GWAS_", tissue.spe, "_", term,
                                     "_fdr", FDRcis, FDRtrans, "_", window, "MB_LD.pdf")
core.score.data <- paste0("cores_score_GWAS_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
ld.lrt.file.txt <- paste0("cores_score_GWAS_LD_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")
median.ratio.file <- paste0("cores_score_GWAS_median_ratio_fdr", FDRcis, FDRtrans, "_", window, "MB.txt")

### Load data
library(condor)
library(igraph)
library(data.table)
library(broom)

load(tissue.file)
load(paste0(cluster.dir, topo.file))
data <- load(paste0(cluster.dir, community.file))
communities <- get(data)
rm(data)

gwas_data <- read.delim(gwas.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gwas_data <- gwas_data[gwas_data$PVALUE_MLOG>=8,]
gwas_data <- gwas_data[!is.na(gwas_data$PVALUE_MLOG),]

gwas.term <- read.delim(paste0(data.dir, term.file), header=F, stringsAsFactors=F)
gwas.diseases<-read.delim(paste0(gwas.dir, gwas.diseases.file), header=F, stringsAsFactors=F)
snp.gwas.term <- unique(gwas.diseases[gwas.diseases[,2] %in% gwas.term[,1], 1])
snp.gwas.term <- snp.gwas.term[snp.gwas.term %in% gwas_data$SNPS]

block.info <- data.table(read.delim(paste0(geno.dir, LDblock.file), header=TRUE, sep="\t", stringsAsFactors=FALSE))
setkey(block.info, SNP)
save(block.info, file=paste0(geno.dir, LDblock.rdata))

### Compute p-values using LRT for GWAS enrichment in high core-score correcting for LD
LRT.Qi.LD <- NULL
Qi.data.final.LD <- list()
for(tissue in names(topo)[nb.samples>=200]){
    cat('Running LRT for', tissue, '...\n')
    
    Qi.data <- topo[[`tissue`]]
   
    Qi.data$GWAS <- rep(0, nrow(Qi.data))
    Qi.data$GWAS[rownames(Qi.data) %in% gwas_data$SNPS] <- 1
    Qi.data$GWAS <- as.factor(Qi.data$GWAS)
    Qi.data$com <- as.factor(Qi.data$com)
    Qi.data <- data.table(Qi.data)
    setkey(Qi.data, red.names)
    Qi.data.merge <- data.frame(merge(Qi.data, block.info, by.x='red.names', by.y='SNP'), stringsAsFactors=F)[,1:7]
    Qi.GWAS <- Qi.data.merge[Qi.data.merge$GWAS==1,]
    Qi.non.GWAS <- Qi.data.merge[Qi.data.merge$GWAS==0,]
    Qi.GWAS.values <- tapply(Qi.GWAS$Qi, paste(Qi.GWAS$block.id, Qi.GWAS$com, sep="_"), median)
    Qi.non.GWAS.values <- tapply(Qi.non.GWAS$Qi, paste(Qi.non.GWAS$block.id, Qi.non.GWAS$com, sep="_"), median)

    Qi.data.final.LD[[`tissue`]] <- rbind(data.frame("Qi"=Qi.GWAS.values, "com"=gsub("[0-9]*_", "", names(Qi.GWAS.values)),
                                   GWAS=rep(1, length(Qi.GWAS.values))),
                        data.frame("Qi"=Qi.non.GWAS.values, "com"=gsub("[0-9]*_", "", names(Qi.non.GWAS.values)),
                                   GWAS=rep(0, length(Qi.non.GWAS.values)))
                        )
    l3 <- lm(Qi ~ com, Qi.data.final.LD[[`tissue`]])
    l4 <- lm(Qi ~ GWAS + com, Qi.data.final.LD[[`tissue`]])
    LRT.Qi.LD <- rbind(LRT.Qi.LD, c(l4$coeff[2], summary(l4)$coefficients[2,2],
                                    tidy(anova(l3, l4, test="LRT"))$Pr..Chi.[2]))
}
colnames(LRT.Qi.LD) <- c("Qi.coef", "Qi.se", "Qi.LRT.pval")

### Compute combined p-values across tissues
library(metap)
combined.pval <- sumlog(LRT.Qi.LD$Qi.LRT.pval/2)$p
combined.pval <- sumlog(LRT.Qi.LD$Qi.LRT.pval/2)$p
LRT.Qi.LD <- rbind(LRT.Qi.LD, c(NA, NA, combined.pval))
write.table(cbind(c(as.character(Tissues[names(topo),1]), "Combined"), LRT.Qi.LD),
            file=paste0(general.path, gwas.path, ld.lrt.file.txt),
            quote=F, row.names=F, sep="\t")

### For the specified term, compute p-values using LRT for GWAS enrichment in high core-score correcting for LD
Qi.data.term <- topo[[`tissue.spe`]]
Qi.data.term$GWAS <- rep(0, nrow(Qi.data.term))
Qi.data.term$GWAS[rownames(Qi.data.term) %in% snp.gwas.term] <- 1
Qi.data.term$GWAS[rownames(Qi.data.term) %in% gwas_data$SNPS & Qi.data.term$GWAS != 1] <- NA
Qi.data.term <- Qi.data.term[!is.na(Qi.data.term$GWAS),]
Qi.data.term$GWAS <- as.factor(Qi.data.term$GWAS)
Qi.data.term$com <- as.factor(Qi.data.term$com)

Qi.data.term.ld <- data.table(Qi.data.term)
setkey(Qi.data.term.ld, red.names)
Qi.data.term.ld.merge <- data.frame(merge(Qi.data.term.ld, block.info,
                                            by.x='red.names', by.y='SNP'), stringsAsFactors=F)[,1:7]
Qi.GWAS.term.ld <- Qi.data.term.ld.merge[Qi.data.term.ld.merge$GWAS==1,]
Qi.non.GWAS.term.ld <- Qi.data.term.ld.merge[Qi.data.term.ld.merge$GWAS==0,]
Qi.GWAS.values.term.ld <- tapply(Qi.GWAS.term.ld$Qi, paste(Qi.GWAS.term.ld$block.id, Qi.GWAS.term.ld$com, sep="_"), median)
Qi.non.GWAS.values.term.ld <- tapply(Qi.non.GWAS.term.ld$Qi,
                                        paste(Qi.non.GWAS.term.ld$block.id, Qi.non.GWAS.term.ld$com, sep="_"),
                                    median)
    
Qi.data.final.term.ld <- rbind(data.frame("Qi"=Qi.GWAS.values.term.ld,
                                "com"=gsub("[0-9]*_", "", names(Qi.GWAS.values.term.ld)),
                                GWAS=rep(1, length(Qi.GWAS.values.term.ld))),
                                data.frame("Qi"=Qi.non.GWAS.values.term.ld,
                                "com"=gsub("[0-9]*_", "", names(Qi.non.GWAS.values.term.ld)),
                                GWAS=rep(0, length(Qi.non.GWAS.values.term.ld)))
                                )
l3 <- lm(Qi ~ com, Qi.data.final.term.ld)
l4 <- lm(Qi ~ GWAS + com, Qi.data.final.term.ld)
LRT.Qi.term.LD <- c(l4$coeff[2], summary(l4)$coefficients[2,2],
                    tidy(anova(l3, l4, test="LRT"))$Pr..Chi.[2])
save(Qi.data.final.LD, Qi.data.final.term.ld, file=paste0( gwas.dir, core.score.data))


median.ratio <- function(x){
    return(median(x$Qi[x$GWAS==1])/median(x$Qi[x$GWAS==0]))
}

mr.term <- median.ratio(Qi.data.term)
mr.term.LD <- median.ratio(Qi.data.final.term.ld)
mr.all <- unlist(lapply(Qi.data.final, median.ratio))
mr.all.LD <- unlist(lapply(Qi.data.final.LD, median.ratio))

d <- data.frame("LDcorrected"=c(mr.all.LD, mr.term.LD))
rownames(d) <- c(as.character(Tissues[names(mr.all),1]), "Whole blood (Autoimmune diseases )")
write.table(d, paste0(gwas.dir, median.ratio.file), sep="\t", quote=F)


### Plot enrichment GWAS filtered for LD SNPs in high core-scores
tissue_names=names(Qi.data.final.LD)

pdf(paste0(figure.dir, gwas.corescore.LD.pdf), width=8, height=11)
par(mfrow=c(4,4))
par(mar=c(2,5,3,1)+0.1)
for(i in 1:length(tissue_names)){
    tissue_name=tissue_names[i]
    print(tissue_name)
    b <- boxplot(Qi.data.final.LD[[`tissue_name`]]$Qi ~ Qi.data.final.LD[[`tissue_name`]]$GWAS, plot=F)
    boxplot(Qi.data.final.LD[[`tissue_name`]]$Qi ~ Qi.data.final.LD[[`tissue_name`]]$GWAS,
    col=c("grey", "dodgerblue"), outline=F,
    notch=F, names=c("n-GWAS", "GWAS"),
    ylab="Qi", main=Tissues[tissue_name,1], ylim=c(0, max(b$stats)*1.15), cex.main=0.9, cex.names=0.7)
}
dev.off()

### Plot GWAS/core-score effect size.
pdf(paste0( figure.dir, pdf.effect.size), width=8, height=4)
par(las=1, mar=c(4, 4, 1, 6)+.1)
plot(c(0, max(LRT.Qi.LD[,1]+LRT.Qi.LD[,2])), c(0.5,(nrow(LRT.Qi.LD)+0.5)), xlab="GWAS effect size (x10-7)", ylab="", yaxt='n', xaxt='n', type='n')
axis(side=2, at=1:nrow(LRT.Qi.LD), labels=Tissues[names(topo)[nb.samples>=200],2][nrow(LRT.Qi.LD):1])
axis(side=1, at=c(0,5,10,15, 20)*1e-07, labels=c(0,5,10,15, 20))
points(LRT.Qi.LD[,1], nrow(LRT.Qi.LD):1, pch=16, col="dodgerblue1")
for(i in 1:nrow(LRT.Qi.LD)){
lines(rep(LRT.Qi.LD[i,1]+LRT.Qi.LD[i,2],2), nrow(LRT.Qi.LD)+1-i+c(0.1, -0.1), lwd=2, col="dodgerblue1")
lines(rep(LRT.Qi.LD[i,1]-LRT.Qi.LD[i,2],2), nrow(LRT.Qi.LD)+1-i+c(0.1, -0.1), lwd=2, col="dodgerblue1")
lines(LRT.Qi.LD[i,1]+LRT.Qi.LD[i,2]*c(1,-1), rep(nrow(LRT.Qi.LD)+1-i, 2), lwd=2, col="dodgerblue1")
}
par(xpd=T)
text(rep(par("usr")[2], nrow(LRT.Qi.LD)), nrow(LRT.Qi.LD):1, labels=sprintf("p = %.2e", LRT.Qi.LD[,3]), pos=4)
par(xpd=NA)
dev.off()

### Plot enrichment tissue-related GWAS SNPs in high core-scores for specified term


pdf(paste0(figure.dir, gwas.corescore.term.ld.pdf ), width=5, height=8)
par(mar=c(4,5,3,1)+0.1)
b <- boxplot(Qi.data.final.term.ld$Qi ~ Qi.data.final.term.ld$GWAS, outline=F, plot=F)

boxplot(Qi.data.final.term.ld$Qi ~ Qi.data.final.term.ld$GWAS, xlab="", ylab="Core-scores",
names=c("non-GWAS", "GWAS"), col=c("grey", "dodgerblue"),
outline=F, notch=T, ylim=c(0, max(b$stats)*1.15) )
lines(1:2, rep(max(b$stats)*1.05,2))

text(1.5, max(b$stats)*1.05, labels=sprintf("P = %.2e", LRT.Qi.term.LD[3]), pos=3)
dev.off()

