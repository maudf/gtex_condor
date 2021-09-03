#!/usr/bin/Rscript
### Maud Fagny
### 2021-08-13
### neutrality_stats.R
### Compute different statistics on phased VCF files
###-----------------------------------------------------

### Load packages
require("vcfR")
require("PopGenome")
require("getopt")

### Set variables
spec = matrix(c('help','h', 0, "logical",
        'directory','d', 1, "character",
        'chromosome', 'c', 1, "character",
        'simulation','s', 0, "logical",
        'ancestral','a', 1, "character",
        'populations','p', 1, "character",
        'width','w', 2, "integer", 
        'loci', 'l', 2, "integer"),
        byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

if ( is.null(opt$directory    ) ) {opt$directory    = './'     }
if ( is.null(opt$chromosome    ) ) {opt$chromosome    = '1'     }
if ( is.null(opt$simulation    ) ) {opt$simulation    = FALSE     }
if ( is.null(opt$populations    ) ) {opt$populations    = ''     }
if ( is.null(opt$ancestral    ) ) {opt$ancestral    = 'no'     }
if ( is.null(opt$width   ) ) { opt$width   = 100000    }
if ( is.null(opt$loci ) ) { opt$loci = 30 }

cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n",
    "chromosome =", opt$chromosome, "\n",
    "simulation =", opt$simulation, "\n",
    "ancestral =", opt$ancestral, "\n",
    "populations =", opt$populations, "\n",
    "chromosome length =", opt$width, "\n",
    "chromosome number =", opt$loci, "\n")

### Load sequences and compute statistics

vcfs <- list.files(path=opt$directory, pattern="*.vcf.gz$")
filenames=gsub("Rep_.*_", "Rep_*_", vcfs[1])

cat("Computing statistics for files", filenames, "in folder", opt$directory, "\n")

res <- list()
for (vcf in paste0(opt$directory, vcfs) ){

    ## Get replicate ID if simulation, file ID if not
    repl <- gsub(".*/", "", vcf)
    repl <- gsub("Rep_", "", repl)
    repl <- gsub("_.*.vcf.gz", "", repl)
    
    cat("\n File ID: ", repl, "\n")

    ## Load VCF file
    data <- readVCF(vcf, numcols=10000, tid=opt$chromosome, frompos=1, topos=opt$width*opt$loci*2) # Read chr 1
    data1 <- data
    
    ## Create ancestral individuals if none exists (for simulations only)
    pop=as.numeric(unlist(strsplit(opt$populations, split=",")))
    if(opt$ancestral=="no" & opt$simulation){
        a=data@region.data@biallelic.matrix[[1]][1:2,]
        rownames(a)=c("anc1", "anc1.2")
        a[a==1]=0
        data@region.data@biallelic.matrix[[1]]=as.ff(rbind(data@region.data@biallelic.matrix[[1]][,1:ncol(data@region.data@biallelic.matrix[[1]])], a))
        data@region.data@populations[[1]][[1]]=c(data@region.data@populations[[1]][[1]],
                                                 pop[length(pop)]*2+1, pop[length(pop)]*2+2)
        opt$ancestral=c("anc1")
        pop=c(pop, pop[length(pop)]+1)
    } else if (opt$ancestral=="no") {
        cat("\n!!! Warning!!!\nThese data do not have outgroups. Statistics such as Fay and Wu's H and Zeng's E and iHS will not be computed.\n!!! Warning!!!\n\n")
    }
    
    ## Set populations
    pop.list <- list()
    for(i in 1:length(pop)){
        pop.list[[i]] <- get.individuals(data)[[1]][seq(ifelse(i==1, 1, pop[(i-1)]*2+1), pop[i]*2, 2)]
    }
    data <- set.populations(data, pop.list , diploid = TRUE)
    pop.list1 <- pop.list
    pop.list1[[length(pop.list1)]] <- NULL
    data1 <- set.populations(data1, pop.list1, diploid = TRUE)

    ## Set outgroup
    if(opt$ancestral != 'no'){
        opt$ancestral <- unlist(strsplit(opt$ancestral, split=","))
        data=set.outgroup(data, new.outgroup = opt$ancestral, diploid = TRUE)
    }
    
    ## Compute neutrality stats
    data.chr <- sliding.window.transform(data, width=2500, jump=1250, type=2, whole.data=TRUE)
    data.chr <- neutrality.stats(data.chr)

    ## Compute stats by populations except for outgroup
    
    nbpop=ifelse(opt$ancestral=='no', length(pop), length(pop)-1) # Remove outgroup
    for(i in 1:nbpop){
        n <- paste0("neutrality.stat.P",i)
        stats.chr <- get.neutrality(data.chr)[[1]] # get neutrality statistics
        stats.chr=data.frame(stats.chr, stringsAsFactors=F)
        # Remove unusefull data from neutrality stat table
        stats.chr$Rozas.R_2=NULL
        stats.chr$Strobeck.S=NULL
        stats.chr$Fu.F_S=NULL
        # If not outgroup, remove Zheng's E and Fay and Wu's H
        if(opt$ancestral=='no'){
            stats.chr$Fay.Wu.H=NULL
            stats.chr$Zeng.E=NULL
        }
        # Extract window start and end
        stats.chr$start=as.numeric(gsub(" -.*", "", rownames(stats.chr)))
        stats.chr$end=as.numeric(gsub(" :.*", "", gsub(".* - ", "", rownames(stats.chr))))
        # If analysing SLIM simulation, compute chromosome
        if(opt$simulation){
            stats.chr$Chr=ifelse(((stats.chr$start) %/% ((opt$width*2)+1)) == ((stats.chr$end-2) %/% ((opt$width*2)+1)), ((stats.chr$start) %/% ((opt$width*2)+1))+1, NA)
            rownames(stats.chr)=NULL
        }
        res[[`n`]] <- rbind(res[[`n`]], data.frame(Replicate=rep(repl, nrow(stats.chr)), 
                                                      stats.chr, stringsAsFactors=F))
    }       

    ### Compute FST by nucleotide between pop1 and 2 # modify if more than 2 population of interest
    
    data <- detail.stats(data, site.FST=TRUE)
    FST <- data@region.stats@site.FST[[1]]
    res$FST <- rbind(res$FST, data.frame(Replicate=rep(repl, length(FST)),
                                         "chr"=ifelse(opt$simulation, (as.numeric(names(FST))-1) %/% ((opt$width*2*opt$loci) +1), opt$chromosome),
                                         "position"=(as.numeric(names(FST))), 
                                         "FST"=FST))
}

### Remove lines without data
if(!is.null(res$neutrality.stat.P1)){
    na.P1 <- which(is.na(res$neutrality.stat.P1$Chr))
    res$neutrality.stat.P1 <- res$neutrality.stat.P1[-na.P1,]
}

if(!is.null(res$neutrality.stat.P2)){
    na.P2 <- which(is.na(res$neutrality.stat.P2$Chr))
    res$neutrality.stat.P2 <- res$neutrality.stat.P2[-na.P2,]
}

if(!is.null(res$FST)){
    nan.fst <- which(is.nan(res$FST$FST))
    res$FST <- res$FST[-nan.fst,]
}
cat("\nSave neutrality stats and FST in", paste0(opt$directory, "summary_stats.rds"), "...\n")
saveRDS(res, file=paste0(opt$directory, "summary_stats.rds"))


