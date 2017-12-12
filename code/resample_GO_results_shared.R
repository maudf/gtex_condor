### resample_GO_results_shared.R

### Load variables
load("code/variables_definition.R")
args = commandArgs(trailingOnly=TRUE)
tis <- args[1] # Tissue name
ind <- args[2] # index
sam <- args[3] # number of resamplings

### Set variables

genes.file <- paste0("all_tissues_genes_fdr", FDRcis, FDRtrans,
                           "_", window, "MB.Rdata")
go.results.txt <- paste0('alltissues_fdr', FDRcis, FDRtrans,
                         "_", window, 'MB_go_results_bp.txt')
 
### Load data
load(paste0(cluster.dir, genes.file))
load(tissue.file)
tislist<-rownames(Tissues)
names(tislist)<-as.character(Tissues[,1])

go.bp <- read.table(paste0(cluster.dir, go.results.txt), header=T, sep="\t", stringsAsFactors=F)
a <- go.bp[go.bp$p.adjust<=0.05 & go.bp$OddsRatio>1 & go.bp$NumberTissues>=12,]

### Extract communities with shared biological functions
shared.com <- tapply(a$com, a$Tissue, unique)
names(shared.com) <- tislist[names(shared.com)]

### Functions
resample.GO <- function(genes, tis, sam, shared.com, unadj.p.cut, min.overlap=4, ind){
#for other annotations, run ls("package:org.Hs.eg.db")
    require(GOstats)
    require(org.Hs.eg.db)
    x <- org.Hs.egENSEMBL
    mapped_genes <- mappedkeys(x)
    xx <- as.list(x[mapped_genes])
    base.sym <- unlist(xx)
    base.entrez <- unlist(lapply(seq_along(xx), function(i, y, n) { rep(n[i], length(y[[i]])) }, y=xx, n=names(xx)))

    all.genes <- gsub("\\.[0-9][0-9]*", "", unlist(genes[[`tis`]]))                           
    if(length(all.genes) != length(unique(all.genes))){ 
        stop("There are duplicate gene names in gene list!
                This is probably the result of an error in constructing the genes object")
    }
    
    all.entrez <- unique(base.entrez[base.sym %in% all.genes])

    stg.mc.entrez <- list()
    find.length <- function(x, coms){
        l<-NULL
        for(i in coms){
            l<-c(l, length(x[[i]]))
        }
        return(l)
    }
    length.com <- find.length(genes[[`tis`]], coms=shared.com[[`tis`]])
    samp <- lapply(1:sam, function(x, n, l){
        a <- sample(n, sum(l), replace = FALSE)
        res <- list()
        count=1
        for(k in 1:length(l)){
            res[[k]] <- a[count:(count+l[k]-1)]
            count=count+length.com[k]
        }
        return(res)
    }, n=all.genes, l=length.com)
    
    for(i in 1:length(samp)){
        cat("Sampling number", i, "\nCommunity number ")
        all.tabs <- data.frame()
        for(j in 1:length(samp[[i]])){
            cat( j, " ")
            this.comm <- samp[[i]][[j]]
            this.entrez <- unique(base.entrez[base.sym %in% this.comm])
    
            params<- new("GOHyperGParams", geneIds = this.entrez,
                         universeGeneIds = all.entrez,
                         ontology= 'BP',
                         pvalueCutoff = unadj.p.cut,
                         conditional = FALSE,
                         testDirection="over",
                         annotation="org.Hs.eg.db")
            
            tab.true <- NULL
            
            if (length(this.entrez)>min.overlap){
                p=params
                p@geneIds <- unique(p@geneIds[p@geneIds %in% p@universeGeneIds])
                cat2Entrez.list <- categoryToEntrezBuilder(p)
                if( length(cat2Entrez.list)>0 )
                {
                    go.res <- hyperGTest(params)
                    this.df <- summary(go.res)
                    this.df.p <- p.adjust(this.df[,2], method="BH")
                    
                    if (sum(this.df$Pvalue < unadj.p.cut)>0)
                    {
                        this.tab <- this.df[this.df$Pvalue < unadj.p.cut,]
                        tab.true <- rbind(tab.true,
                                          cbind(rep(tis, nrow(this.tab)), rep(i, nrow(this.tab)),
                                                this.df.p[this.df$Pvalue < unadj.p.cut],
                                                rep(j, nrow(this.tab)),
                                                this.tab))
                    }
                }
                all.tabs <- rbind(all.tabs,tab.true)
            }
            
        }
        cat("\n")
        names(all.tabs)[1:5] <- c("tissue", "sampling", "p.adjust","comm.id", "GOID")
        filename=paste0(resampling.GO.dir,
                        tis, "_resampling_GO_shared_results_",ind,".txt")
        if(i==1){
            write.table(all.tabs[all.tabs$Count > min.overlap,], file=filename, sep="\t", quote=F, row.names=F)
        } else {
            write.table(all.tabs[all.tabs$Count > min.overlap,], file=filename, sep="\t", quote=F,
                        append=T, row.names=F, col.names=F)
        }
    }
    

    
### A GOstats dependency,"graph", has conflicting functions with igraph, so I remove both after 
### running GOStats 
    if(sum(search() == "package:GOstats")){
        detach("package:GOstats"); detach("package:graph") }
}

resample.GO(genes=genes, tis=tis, sam=sam, shared.com=shared.com, unadj.p.cut=1, min.overlap=4, ind=ind)
