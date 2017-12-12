### create_topo_file.R

### Load data
load("variables_definition.R")
load(paste0(tissue.file))

### Set parameters
topo.file <- paste0("network_topology_fdr", FDRcis, FDRtrans, "_", window, "MB.RData")

topo <- lapply(communities, function(x){
    cs <- x$qscores$red.qscore
    tmp <- x$Qcoms
    rownames(tmp) <- tmp[,2]
    Qk <- tmp[x$qscores$red.qscore$com,1]
    rm(tmp)
    qi <- x$qscores$red.qscore$Qik*Qk
    cs$Qi <- qi
    rownames(cs) <- as.character(cs$red.names)
    s <- as.character(x$edges[,grep("red", colnames(x$edges))])
    tab <- table(s)
    cs$degree <- as.numeric(tab[rownames(cs)])
    return(cs)})
save(topo, file=paste0(cluster.dir, topo.file))
