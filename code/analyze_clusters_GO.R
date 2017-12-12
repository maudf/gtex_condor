### Maud Fagny
### 2015/12/21
### analyze_cluster_GO.R
### make modularity plots and gene ontology analysis
### ________________________________________________________

### Load variables
load("code/variables_definition.R")

### Set variables
communities.file <- paste0("all_tissues_communities_fdr", FDRcis, FDRtrans,
                           "_", window, "MB.Rdata")
go.results.file <- paste0('alltissues_fdr', FDRcis, FDRtrans,
                          "_", window, 'MB_go_results_', gene.bg, 'bkg.Rdata')
go.results.txt <- paste0('alltissues_fdr', FDRcis, FDRtrans,
                         "_", window, 'MB_go_results_', gene.bg, 'bkg.txt')
 
### Load data
source(paste0(general.path, scripts.path, 'condor.GO.R'))
load(paste0(cluster.dir, communities.file))

### Enrichment in GO terms among communities
go.result <- list()
tissues <- names(communities)
for (tissue in tissues){
    cat("Computing GO and KEGG enrichment analyses for", tissue, "...\n")
    condor.modularity <- communities[[`tissue`]]
    condor.modularity$qscores$blue.qscore$blue.names <- gsub("\\.[0-9][0-9]*", "",
                                                             condor.modularity$qscores$blue.qscore$blue.names)
    condor.modularity$qscores$red.qscore$red.names <- gsub("\\.[0-9][0-9]*", "",
                                                           condor.modularity$qscores$red.qscore$red.names)
    condor.modularity$blue.memb$blue.names <-gsub("\\.[0-9][0-9]*", "",
                                                           condor.modularity$blue.memb$blue.names)
    condor.modularity$red.memb$red.names <-gsub("\\.[0-9][0-9]*", "",
                                                  condor.modularity$red.memb$red.names)

    if (gene.bg == "all"){
        genes <- read.table(paste(rna.dir, tissue, '_genes.tsv', sep = ''), header = T, sep = '\t', stringsAsFactors = F)
        all.gene <- gsub("\\.[0-9][0-9]*", "", genes$genes)
        rm(genes)
        } else {
        all.gene <- condor.modularity$qscores$blue.qscore$blue.names
    }
    go.result[[`tissue`]] <- condor.GO(condor.object = condor.modularity,
                                       go.ontology = c('BP', 'MF', 'CC'), unadj.p.cut = 0.1,
                                       min.overlap=5, symbol.universe = all.gene)
    rm(condor.modularity)
    save(go.result, file=paste0(cluster.dir, go.results.file))
}

ti=names(go.result)[1]
tab.GO <- data.frame("Tissue"=rep(ti, nrow(go.result[[`ti`]])), go.result[[`ti`]], stringsAsFactor=F)
for(ti in names(go.result)[2:length(names(go.result))]){
    tab.GO <- rbind(tab.GO, data.frame("Tissue"=rep(ti, nrow(go.result[[`ti`]])), go.result[[`ti`]], stringsAsFactor=F))
}

write.table(tab.GO, paste0( cluster.dir, go.results.txt),
            row.names=F, sep="\t", quote=F)
