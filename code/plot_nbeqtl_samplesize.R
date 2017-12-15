###summarize_eqtls.R

### Load variables
source("code/variables_definition.R")

### Set parameters
eqtls.file <- paste0("all_tissues_eqtls_fdr", FDRcis, FDRtrans, "_", window, "MB.Rdata")

### Load data
load(paste0(tissue.file))
load(samples.file)
load(paste0(eqtl.dir, eqtls.file))

### Find number of cis and trans eQTLs
size.com <- lapply(eqtl, function(x){
    c("cis"=sum(x$cis.or.trans=="cis" & x$FDR<=0.05),
    "trans"=sum(x$cis.or.trans=="trans" & x$FDR<=0.05))
})
nb.cis <- unlist(lapply(size.com, function(x){x['cis']}))
nb.trans <- unlist(lapply(size.com, function(x){x['trans']}))
names(nb.cis) <- gsub("\\.cis", "", names(nb.cis))
names(nb.trans) <- gsub("\\.trans", "", names(nb.trans))

### Plot number of cis- and trans-eQTLs as a function of tissue sample size
pdf(paste0(figure.dir, "eqtls_nbsamples_log.pdf"),
    height=5, width=8)
par(mar=c(4,5,1,1)+0.1)
plot(nb.samples[names(nb.cis)], (nb.cis), pch=16, col="red", log='y',
     ylim=c(1000,700000), xlab="Nb samples", ylab="Nb eQTLs")
points(nb.samples[names(nb.trans)], (nb.trans), pch=17, col="blue")
legend("topleft", bty='n', legend=c("cis-eQTLs", "trans-eQTLs"), pch=16:17, col=c("red", "blue"))
dev.off()
