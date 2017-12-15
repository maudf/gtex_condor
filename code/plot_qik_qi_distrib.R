### plot_qik_qi_distrib.R

### Load variables
source("code/variables_definition.R")

### Set variables
communities.file <- paste0("all_tissues_communities_fdr", FDRcis, FDRtrans,
                           "_", window, "MB.Rdata")
### Load data
load(cluster.dir, communities.file)
load(tissue.file)

### Plot Qi=f(Qik)
pdf(figure.dir, "Qik_Qi_distrib.pdf", height=11.5, width=8)
par(mfrow=c(13,2))
for(tis in names(communities)){
    print(tis)
    qik <- communities[[`tis`]]$qscores$red.qscore$Qik
    com <- as.character(communities[[`tis`]]$qscores$red.qscore$com)
    sum.com <- names(sort(table(com)))
        
    length.com=1:length(sum.com)
    names(length.com)=sum.com
    l = length.com[com]
    qk.list <- communities[[`tis`]]$Qcom[,1]
    names(qk.list) <- as.character(communities[[`tis`]]$Qcom[,2])
    qk <- qk.list[com]
    qi <- qik*qk
    
    se <- function(x){sd(x, na.rm=TRUE) /  
                          sqrt(length(x[!is.na(x)])) }
    
        
    mqi <- tapply(qi, com, mean)[sum.com]
    sqi <- tapply(qi, com, se)[sum.com]
    mqik <- tapply(qik, com, mean)[sum.com]
    sqik <- tapply(qik, com, se)[sum.com]
    
    par(mar=c(2,5,1,0)+.1)
    plot(1:length(sum.com), log10(mqi), xlab="",
         ylab="Qi",
         pch=20, col="black", main="")	
    mtext(as.character(Tissues[tis,2]), side=2, line=4)
    for(i in 1:length(sum.com)){
        lines(rep(i,2), log10(mqi[i]+sqi[i]*c(-1,1)), lwd=2, col="black")
    }
    par(mar=c(2,4,1,0)+.1)        
    plot(1:length(sum.com), log10(mqik), xlab="",
         ylab="Qik", pch=20, col="red")
    for(i in 1:length(sum.com)){
        lines(rep(i,2), log10(mqik[i]+sqik[i]*c(-1,1)), lwd=2, col="red")
    }
}
dev.off()
    
