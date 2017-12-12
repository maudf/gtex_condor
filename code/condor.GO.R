#' Find Gene Ontology Term enrichment in condor communities
#' 
#' This function iterates through the genes in each community located in 
#' \code{condor.object$blue.memb} and tests the gene set using a one-sided
#' Fisher's exact test as described in \code{\link[GOstats]{GOHyperGParams}}
#' @param condor.object output of either \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @param go.ontology specify which gene ontology classes to include: 
#' \code{MF},\code{CC}, \code{BP}
#' @param unadj.p.cut unadjusted p-value threshold for recording results. 
#' All tests with an unadjusted p-value less than \code{unadj.p.cut} will 
#' not be output.
#' @param symbol.universe vector of genes used to define the set of 
#' possible genes to include in Fisher's Exact Test. 
#'
#' @param min.overlap The minimum overlap between GO gene group and genes
#' in community required to run the Fisher's Exact Test. Default is 4.
#' @return a \code{data.frame} containing the Community ID, GO term, number
#' of genes in the overlap, unadjusted p-value, 
#' @references \url{http://tools.medialab.sciences-po.fr/iwanthue/} for 
#'  a nice color generator at 
#' @note For the condor paper \url{http://arxiv.org/abs/1509.02816}, I used
#'   35 colors from the "Tarnish" palette with "hard" clustering
#' @examples
#' 
#' @import data.table
#' @import GOstats
#' @export

condor.GO <- function(condor.object,go.ontology,unadj.p.cut,min.overlap=4,symbol.universe=NULL){
#for other annotations, run ls("package:org.Hs.eg.db")
    require(GOstats)
    require(org.Hs.eg.db)
    require(data.table)
    x <- org.Hs.egENSEMBL
    mapped_genes <- mappedkeys(x)
    xx <- as.list(x[mapped_genes])
    base.sym <- unlist(xx)
    base.entrez <- unlist(lapply(seq_along(xx), function(i, y, n) { rep(n[i], length(y[[i]])) }, y=xx, n=names(xx)))
    if(is.null(symbol.universe)){
        print("Warning! Using all entrez genes as symbol universe, as none was provided")
        symbol.universe <- unique(base.sym)
    }
    all.entrez <- base.entrez[base.sym %in% symbol.universe]
    gene_com <- data.table(condor.object$blue.memb)
    setkey(gene_com,com)
    coms <- sort(unique(condor.object$blue.memb$com))

    stg.mc.entrez <- list()
    all.tabs <- data.frame()
                                
    for(j in go.ontology){
        print(j)
        for(i in coms){
        
            if(length(gene_com[J(i)]$blue.name) != length(unique(gene_com[J(i)]$blue.name))){ 
                stop("There are duplicate gene names in each community!
                This is probably the result of an error in constructing the condor object")
            }
            
            this.comm <- gene_com[J(i)]$blue.name
            this.entrez <- base.entrez[base.sym %in% this.comm]
            stg.mc.entrez[[i]] <- this.entrez  
        }
        this.comm <- stg.mc.entrez[[1]]
    
        params<- new("GOHyperGParams",geneIds = this.comm,
                     universeGeneIds = all.entrez,
                     ontology= j,
                     pvalueCutoff = 1,
                     conditional = FALSE,
                     testDirection="over",
                     annotation="org.Hs.eg.db")
    
        tab.true <- NULL
    
        for (i in 1:length(stg.mc.entrez))
            {
                print(i)
                this.comm <- stg.mc.entrez[[i]]
                if (length(this.comm)>min.overlap){
                    geneIds(params) <- this.comm
                                        #conditional(params) <- TRUE
                    p=params
                    p@geneIds <- unique(p@geneIds[p@geneIds %in% p@universeGeneIds])
                    cat2Entrez.list <- categoryToEntrezBuilder(p)
                    if( length(cat2Entrez.list)>0 )
                        {
                            go.res <- hyperGTest(params)
                            this.df <- summary(go.res)
                            this.df.p <- p.adjust(this.df[,2],method="BH")
                            
                            if (sum(this.df$Pvalue < unadj.p.cut)>0)
                                {
                                    this.tab <- this.df[this.df$Pvalue < unadj.p.cut,]
                                    tab.true <- rbind(tab.true,cbind(rep(j,nrow(this.tab)), this.df.p[this.df$Pvalue < unadj.p.cut],rep(i,nrow(this.tab)),this.tab))
                                }
                        }
                }
            }
        names(tab.true)[4] <- "GOID"
        all.tabs <- rbind(all.tabs,tab.true)
    }
    names(all.tabs)[1:3] <- c("Category", "p.adjust","comm.id")
    
#A GOstats dependency,"graph", has conflicting functions with igraph, so I remove both after 
                                        #running GOStats 
    if(sum(search() == "package:GOstats")){
        detach("package:GOstats"); detach("package:graph") }
    return(all.tabs[all.tabs$Count > min.overlap,])
}
