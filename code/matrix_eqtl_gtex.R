### matrix_eqtl_gtex.R
### Take tissue name as input 

### Load libraries
library(MatrixEQTL)

### Load variables
load("code/variables_definition.R")

###Set up environment variables
output.dir = eqtl.dir
args <- commandArgs(trailingOnly=TRUE)
tissue_id <- args[1]
print(tissue_id)

### Read input files
SNP_file_name = paste0(geno.dir, tissue_id,".dosage");
snps_location_file_name = paste0(geno.dir, tissue_id,".pos");
expression_file_name = paste0(rna.dir, tissue_id,"_norm.tsv");
gene_location_file_name = paste0(rna.dir, tissue_id,"_genes.tsv");
covariates_file_name = paste0(rna.dir, stissue_id,"_covariates.txt");

### Set up output file names
output_file_name_cis = paste0(output.dir, tissue_id,"_cis.txt")
output_file_name_tra = paste0(output.dir, tissue_id,"_trans.txt")
errorCovariance = numeric();

### Load genotype data

#note: SlicedData is unique to MatrixEQTL
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;      #column labels
snps$fileSkipColumns = 1;   #row labels, SNP position info
snps$fileSliceSize = 2000;  #read in 2000 rows at a time
snps$LoadFile(SNP_file_name);
snps$LoadFile(SNP_file_name);

### Load gene expression

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 100000;
gene$LoadFile(expression_file_name);

### Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;
cvrt$fileSkipColumns = 1;
if(length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name);
}

### Load snps and genes positions

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

### Run eQTL mapping

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);



## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected local eQTLs:', '\n');


