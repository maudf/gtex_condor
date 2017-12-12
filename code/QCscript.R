### QCscript.R

### Load libraries
library(Biobase)
rm(list=ls())

### Load data
load("code/variables_definition.R")
load(normalized.rnaseq)

###Set up environment variables
Sys.setenv(plink2=plink2)
Sys.setenv(vcfgenotype=vcfgenotypes)
Sys.setenv(output_dir=QC.dir)

print(is.numeric(sample_threshold))
obj = obj4
out = c()

### Grab cell column indices
tissue_tmp = unique(pData(obj)$our_subtypes)

### Drop testis from consideration (the only sex organ with more than 150 samples)
tissues <- tissue_tmp[!(tissue_tmp == "testis")]

system("sed '18q;d' $vcfgenotype | cut -f 10- > $output_dir/geno_samples.tsv")
geno_sample_id = as.character(read.delim(paste0(QC.dir, "geno_samples.tsv"), header=FALSE,stringsAsFactors=FALSE))
geno_subjid = unlist(lapply(strsplit(geno_sample_id,"-"),function(x){paste(x[1],x[2],sep="-")}))

### Initialize samples to keep
tissue_old <- "all_samples"
Sys.setenv(tissue_old=tissue_old)

### Initial QC for MAF, missingness and indel removal for all 450 samples
system("$plink2 --vcf $vcfgenotype --snps-only --geno 0.1 --maf 0.05 --write-snplist --out $output_dir/$tissue_old")

### Remove duplicated SNPs, there are two SNP IDs that each show up twice
system("uniq -d $output_dir/$tissue_old'.snplist' > $output_dir/dups.snplist")
system("mv $output_dir/$tissue_old'.snplist' $output_dir/withdups.snplist")
system("grep -f $output_dir/dups.snplist -v  $output_dir/withdups.snplist > $output_dir/$tissue_old'.snplist'")

for(tissue_id in tissues){
    tissue_indx <- grep(tissue_id, pData(obj)$our_subtypes)
    ## Define tissue object
    tissue  <- obj[,c(tissue_indx)]
    
    removeSamples = which(duplicated(pData(tissue)$id))
    
    if(length(removeSamples) > 0){
        tissue  = tissue[,-removeSamples]
    }
    
    n_matched = sum(pData(tissue)$SUBJID %in% geno_subjid)
    
    out <- rbind(out,c(tissue_id,n_matched))
    if(n_matched > sample_threshold){
        ## Extract info
        tissue = tissue[,pData(tissue)$SUBJID %in% geno_subjid]
        print(tissue_id)
        ## Filter SNPs by MAF using tissue matched samples

        Sys.setenv(tissue_now=tissue_id)
        write.table(data.frame(FID=0, IID=geno_sample_id[geno_subjid %in% pData(tissue)$SUBJID]),
                      file=paste0(QC.dir, "tmp_samples.tsv"), row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
        system("$plink2 --vcf $vcfgenotype --const-fid --keep $output_dir/tmp_samples.tsv --extract $output_dir/$tissue_old'.snplist' --geno 0.1 --maf 0.05 --write-snplist --out $output_dir/$tissue_now")
        Sys.setenv(tissue_old=tissue_id)
        
        print(paste("Number of samples with expression + genotype in ",tissue_id,"is",n_matched))
        currentNumberOfSNPs = system(paste0("wc -l ",tissue_id,".snplist"),intern=TRUE)
        print(paste("Current number of SNPs",currentNumberOfSNPs))
    }
    
}

### Write out the final set of SNPs to keep:
system("$plink2 --vcf $vcfgenotype --const-fid --keep $output_dir/tmp_samples.tsv --extract $output_dir/$tissue_old'.snplist' --geno 0.1 --maf 0.05 --write-snplist --out $output_dir/kept_snps_MAF05_maxmissing09")
