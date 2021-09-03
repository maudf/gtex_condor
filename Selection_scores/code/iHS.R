#!/usr/bin/Rscript
### Maud Fagny
### 2021-08-13
### neutrality_stats.R
### Compute different statistics on phased VCF files
###-----------------------------------------------------

### Load packages
require("rehh")
require("getopt")
require("this.path")
script.path=gsub("/[^/]*$", "/", this.path())
source(paste0(script.path, "ehh.R"))

### Set variables
spec = matrix(c('help','h', 0, "logical",
        'directory','d', 1, "character",
        'simulation','s', 0, "logical",
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
if ( is.null(opt$simulation    ) ) {opt$simulation    = FALSE     }
if ( is.null(opt$populations    ) ) {opt$populations    = ''     }
if ( is.null(opt$width   ) ) { opt$width   = 100000    }
if ( is.null(opt$loci ) ) { opt$loci = 30 }

cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n",
    "simulation =", opt$simulation, "\n",
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
    repl <- gsub(ifelse(opt$simulation, ".*Rep_", ".*/"), "", vcf)
    repl <- gsub("_.*.vcf.gz", "", repl)
    
    cat("\n File ID: ", repl, "\n")

     ### Read VCF file
    hh <- data2haplohh(hap_file = vcf,
                       vcf_reader = "data.table",
                       remove_multiple_markers=TRUE,
                       polarize_vcf=FALSE,
                       verbose = FALSE)
    
    ## Compute EHH and iHH
    for(i in 1:nbpop){
        n=paste0("iHS.P",i)
        sub <- rownames(hh@haplo)[ifelse(i==1, 1, (pop[(i-1)]*2+2)):(pop[i]*2)]
        hh.subset <- subset.haplohh(hh, select.hap =sub)
        mrk.opt <- seq(opt$width, opt$width*opt$loci*2, opt$width*2)
        
        for(j in 1:length(mrk.opt)){
            cat("Compute iHS for the", j, "th marker...\n")
            mrk <- which(abs(hh.subset@positions-mrk.opt[j])==min(abs(hh.subset@positions-mrk.opt[j])))[1]
            ehh <- calc_ehh(hh.subset, mrk = mrk, include_nhaplo = TRUE, 
                            discard_integration_at_border = FALSE)
            res[[`n`]] <- rbind(res[[`n`]], data.frame( Replicate=repl, 
                                                                        Chr=ifelse(opt$simulation, floor(mrk.opt[j]/(opt$width*2))+1, hh@chr.name), 
                                                                        Freq=ehh[[2]][2], 
                                                                        position=hh.subset@positions[mrk], 
                                                                        ihhA=ehh$ihh[1], ihhD=ehh$ihh[2], 
                                                                        ihs=log(ehh$ihh[1]/ehh$ihh[2])))
            rownames(res[[`n`]])=NULL
        }
    }
}

cat("\nSave iHS results  in", paste0(opt$directory, "iHS.rds"), "...\n")
saveRDS(res, file=paste0(opt$directory, "iHS.rds"))


