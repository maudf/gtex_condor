#!/usr/bin/Rscript
### Compute pcadapt statistics from VCFs

### Load packages
require("getopt")
require("pcadapt")
require("qvalue")

### Set variables
spec = matrix(c('help','h', 0, "logical",
                'directory','d', 1, "character"),
              byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$directory    ) ) {opt$directory    = './'     }
cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n")


### Load sequences and compute statistics

beds <- list.files(path=opt$directory, pattern="*.bed$")
filenames=gsub("Rep_.*_", "Rep_*_", beds[1])
cat("Computing pcadapt for files", filenames, "in folder", opt$directory, "\n")

res <- list(pcadapt=NULL, nb.outliers=NULL)
for (bed in paste0(opt$directory, beds) ){
  
  ## Get replicate ID if simulation, file ID if not
  repl <- gsub(".*/", "", bed)
  repl <- gsub("Rep_", "", repl)
  repl <- gsub("_.*.bed", "", repl)
  
  cat("\n File ID: ", repl, "\n")
  
  ## Load BED file
  cat("Loading", bed, "...\n")

  data <- read.pcadapt(input=bed, type='bed', type.out="matrix")
  
  cat("Computing pcadapt ...\n")
  
  data.m=bed2matrix(bedfile = bed)
  res$pcadapt[[repl]] <- pcadapt(  input=data,
            K = 7,
            min.maf = 0.01,
            ploidy = 2,
            LD.clumping = list(size = 500, thr = 0.1),
            pca.only = FALSE,
            tol = 1e-04
  )
  
  qval <- qvalue(res$pcadapt[[repl]]$pvalues)$qvalues
  alpha <- 0.05
  outliers <- which(qval < alpha)
  res$nb.outliers[[repl]]  <- length(outliers)
}
cat("Save PCadapt results  in", paste0(opt$directory, "pcadapt.rds"), "...\n")
saveRDS(res, file=paste0(opt$directory, "pcadapt.rds"))


