#eqtl_network_clustering.R

### Load libraries
library(condor)#, lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
library(igraph)

### Load variables
source("code/variables_definition.R")

### Set parameters
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
	cat("eqtl_network_clustering.R needs arguments: [eQTL file] [plot]\n")
	cat("e.g.: Rscript eqtl_network_clustering.R output/eqtls/skinfdr0.20.2_1MB_edges.txt TRUE\n")
	stop("Exit...\n")
}


eqtl.file <- args[1] # File containing eQTLs
plot <- args[2] # Whether to plot the graph [TRUE/FALSE]
output.folder <- cluster.dir  # Path to the folder to store output files
basename <- gsub(".*/([^/]*).txt","\\1", eqtl.file) 
print(basename)

if(!file.exists(eqtl.file)) {stop(paste(eqtl.file, "not found. Exit...\n"))}
if(!file.exists(output.folder)) {dir.create(output.folder)}



### Load network data 
eqtls <- read.table(eqtl.file, header = T, stringsAsFactors = F, quote = "")
if(("RS_ID_dbSNP142_CHG37p13" %in% colnames(eqtls)) & ("RS_ID_dbSNP135_original_VCF" %in% colnames(eqtls))){
    eqtls$RSID <- eqtls$RS_ID_dbSNP142_CHG37p13
    eqtls$RSID[eqtls$RS_ID_dbSNP142_CHG37p13=='.'] <- eqtls$RS_ID_dbSNP135_original_VCF[eqtls$RS_ID_dbSNP142_CHG37p13=='.']
}

elist<- data.frame("red" = eqtls$RSID, "blue" = eqtls$genes)
condor.object <- create.condor.object(elist)

### Compute clusters and node modularity
condor.result <- condor.cluster(condor.object)
print(paste(output.folder, basename, "_clusters.Rdata", sep = ""))
save(condor.result, file=paste(output.folder, basename, "_clusters.Rdata", sep = ""))
condor.modularity <- condor.qscore(condor.result)
save(condor.modularity, file=paste(output.folder, basename, "_modularity.Rdata", sep = ""))

### plot graph

if(plot){
	gtoy = graph.edgelist(as.matrix(elist), directed = F)
	set.graph.attribute(gtoy, "layout", layout.kamada.kawai(gtoy))
	V(gtoy)[c( unique(eqtls$genes), unique(eqtls$RS_ID_dbSNP142_CHG37p13) )]$color <- c(rep("red", 
				length(unique(eqtls$genes))), rep("blue", length(unique(eqtls$RS_ID_dbSNP142_CHG37p13))))
	
	pdf(paste(figure.dir, "/", basename, "_graph.pdf", sep = ""), family = "Helvetica") 
	plot(gtoy, vertex.label.dist = 2)
	dev.off()
}


