### create_matrix_eqtl_files.R

###Load variables
load("code/variables_definition.R")

### Load libraries
library(Biobase)
rm(list=ls())

### Load data
load(normalized.rnaseq)

###Set up environment variables
Sys.setenv(dosage_dir=dosage.dir)
Sys.setenv(output_dir=QC.dir)
Sys.setenv(expression_dir=rna.dir)
Sys.setenv(plink2=plink2)
Sys.setenv(vcfgenotype=vcfgenotypes)
set.seed(1)


obj = obj4
out = c()

###Functions
sample.to.subject <- function(s){
    unlist(lapply(strsplit(as.character(s),"-"),
                  function(x){
                      paste(x[1],x[2],sep="-")
                  }))
}

### Grab cell column indices
tissue_tmp = unique(pData(obj)$our_subtypes)

### Drop testis from consideration (the only sex organ with more than 150 samples)
tissues <- tissue_tmp[!(tissue_tmp == "testis")]
print(tissues)

### Sample ID to Subject ID lookup table for genotype data
geno_sample_id = as.character(read.delim(paste0(QC.dir, "geno_samples.tsv"),
                                         header=FALSE, stringsAsFactors=FALSE))
geno_subjid = sample.to.subject(geno_sample_id)

### Read in genetic PCs
all_pcs = read.delim(pc.data, stringsAsFactors = FALSE)
pc_subjid = sample.to.subject(all_pcs$FID)
all_pcs$SUBJID <- pc_subjid
genetic_pcs <- all_pcs[,c("SUBJID", paste0("C",1:number_pcs))]

for(tissue_id in tissues){
    tissue_indx <- grep(tissue_id,pData(obj)$our_subtypes)
    # define tissue object
    tissue  = obj[,c(tissue_indx)]
    n_matched = sum(geno_subjid %in% unique(pData(tissue)$SUBJID))
    
    if(n_matched > sample_threshold){
        
        ## Remove technical replicate RNA samples if there are any
        colnum = 1:ncol(tissue)
        names(colnum) <- pData(tissue)$id
        colSample <- sample(colnum)
        removeSamples = colSample[which(duplicated(names(colSample)))]
        
        if(length(removeSamples) > 0){
            ## Check to make sure we sampled correctly
            if(any(!(unique(pData(tissue)$id) %in% unique(pData(tissue[,-removeSamples])$id)))){
                stop("sampling error in RNA-seq duplicates")
            }
            tissue  = tissue[,-removeSamples]
            
            
            ## Number of matched gene expression and genotype files for tissue_id
            n_check1 = sum(pData(tissue)$SUBJID %in% geno_subjid)
            n_check2 <- n_matched
            print(n_matched)
            
            if(n_check1 != n_check2){
                print(paste("n_check1 is",n_check1))
                print(paste("n_check2 is",n_check2))
                stop("genotype sample number doesn't equal RNA-seq sample number")
            }
        }
        print(n_matched)        
        print(tissue_id)
        out <- rbind(out,c(tissue_id,n_matched))
        
        #Convert vcf to dosage matrix, subset on tissue
        Sys.setenv(tissue_id=tissue_id)
        write.table(data.frame(FID=0, IID=geno_sample_id[geno_subjid %in% pData(tissue)$SUBJID]),
                    file=paste0(QC.dir, "tmp_samples.tsv"), row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
        
        system("$plink2 --vcf $vcfgenotype --const-fid --keep $output_dir/tmp_samples.tsv --extract $output_dir/kept_snps_MAF05_maxmissing09.snplist --recode A-transpose --out $dosage_dir$tissue_id")
        
        
        
        ## Remove the FID (0) that PLINK adds to the sample IDs
        system("awk 'NR==1{gsub(/0_/,\"\")};{print}' $dosage_dir$tissue_id'.traw'> $dosage_dir$tissue_id'.tmp' ")
        system("head -1 $dosage_dir$tissue_id'.tmp' | cut -f 7- > $dosage_dir$tissue_id'.sample'")
        system("awk '{ print $2\"\t\"$1\"\t\"$4 }' $dosage_dir$tissue_id'.tmp' > $dosage_dir$tissue_id'.pos'")
        system("cut -f 2,7- $dosage_dir$tissue_id'.tmp' > $dosage_dir$tissue_id'.dosage'")
        dosage.indv = read.delim(paste0(dosage.dir,tissue_id,".sample"),
                                 header=FALSE,colClasses="character", stringsAsFactors=FALSE, sep="\t")
        dosage.header = sample.to.subject(dosage.indv)
        
        ## Make expression files
        tissue = tissue[,pData(tissue)$SUBJID %in% geno_subjid]
        ematrix = exprs(tissue)
	counts.samples = apply(ematrix, 1, function(x) sum(x>=reads_threshold))
	removeGenes = names(counts.samples[counts.samples<expr_sample])
	
	## Extract normalized gene counts
	nmatrix = assayData(tissue)[["qsmooth"]]
        nmatrix = nmatrix[!(rownames(nmatrix) %in% removeGenes),]
	
	## Correct for batches
        if(batch_correct == TRUE){
            if(!grepl("cells",tissue_id)){
                library(limma)
                batch  = factor(as.character(pData(tissue)$SMNABTCHT))
                nmatrix = removeBatchEffect(nmatrix,batch=batch)
            }
        }
        colnames(nmatrix) <- sample.to.subject(colnames(nmatrix))
        
        ## Expression, genotype, and covariate files must have the same column order
        nmatrix <- nmatrix[,dosage.header]
        pd = pData(tissue)
        rownames(pd) <- sample.to.subject(rownames(pd))
        fd = fData(tissue)
        
        names = c("genes-samples",colnames(nmatrix))
        nmatrix=cbind(rownames(nmatrix),nmatrix)
        nmatrix=rbind(names,nmatrix)
        
        pd = cbind(SUBJID=rownames(pd),pd[,c("gender","AGE","RACE")])
        fd = cbind(genes = rownames(fd),fd)
        
        ## Construct covariate matrix
        eqtl_covariates = merge(pd,genetic_pcs,by="SUBJID")
        rownames(eqtl_covariates) <- eqtl_covariates$SUBJID
        eqtl_covariates$SUBJID <- NULL

        ## Recode Male/Female as a number
        eqtl_covariates$gender <- as.numeric(factor(eqtl_covariates$gender))

        ## Make sure rows are ordered correctly
        eqtl_covariates <- eqtl_covariates[dosage.header,]
        write.table(t(eqtl_covariates), file=paste0(rna.dir, tissue_id, "_covariates.txt"),
                    row.names=TRUE, quote=FALSE, sep="\t")
        
        ## Write MatrixEQTL input files 
        write(t(nmatrix), file=paste0(rna.dir, tissue_id, "_norm.tsv"),
              sep="\t", ncolumns=ncol(nmatrix))
        write.table(pd, file=paste0(rna.dir, tissue_id, "_pheno.tsv"),
                    sep="\t", quote=FALSE, row.names=FALSE)
        write.table(fd[,c(1,5,3,4)], file=paste0(rna.dir, tissue_id, "_genes.tsv"),
                    sep="\t", quote=FALSE, row.names=FALSE)
    }
    
}


