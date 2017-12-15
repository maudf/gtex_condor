### Path to softwares

plink2 <- "/usr/local/bin/plink_v1.90beta/plink"

### Path to data

## Files
pca.data <- "data/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_PostImput_20genotPCs.txt"
vcfgenotypes <- "data/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf"
normalized.rnaseq <- "data/gtex_sub_noxymt.rdata" 
lookup <- "data/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt" # SNP ID/rsID lookup table
tissue.file <- "data/tissue_names.RData"
samples.file <- "data/nb_samples.RData"
gtexqtl.file <- "data/all_tissues_gtex_eqtls.RData"
anno.snps.file <- "data/annotations_snps.RData"
anno.genes.file <- "data/annotations_genes.RData"
genes.corres.file <- "data/genes_Ensembl_HGNC_corres.RData"
epi.file <- "data/epigenome_state_snps.RData"
gwas.file <- "data/gwas_catalog_curated.txt"
gwas.term.hlv.86.file <- "data/gwas_metabolism_terms.txt"
LDblock.file <- "data/LDblockInfo.txt"
gwas.degree.code.file <- "code/gwasDegreePlotFunctions.R"
genes.hlv.86.file <- "data/heart_left_ventricle_86_cellularrespiration_HGNC.txt"
snps.neigh.file <- "data/nb_snps_neighbors.RData"

## Folder
code.dir <- "code/"
data.dir <- "data/"
QC.dir <- "output/qc/"
rna.dir <- "output/expression/"
geno.dir <- "output/dosageMatrices/"
eqtl.dir <- "output/eqtls/"
cluster.dir <- "output/clusters/"
figure.dir <- "Figures/"
epi.dir <- "output/epigenome_roadmap/"
gwas.dir <- "output/gwas/"
resampling.GO.dir <- "output/clusters/GO_resampling/"
circos.dir <- "output/circos_annotations/"
### Variables

## QCscript.R
sample_threshold <- 200 

## Create_matrix_eqtl_files.R
number_pcs <- 3
batch_correct = TRUE
reads_threshold = 6
expr_sample = 10

## matrix_eqtl_gtex.R
pvOutputThreshold_cis = 1e-1; # P-values under which results are printed for cis-eQTLs
pvOutputThreshold_tra = 1e-4; # P-values under which results are printed for trans-eQTLs
cisDist = 1e6; # Maximum distance for local gene-SNP pairs for cis-eQTLs

## eqtl_annotation_merge_bytissue.sh
FDRcis <- 0.2 # FDR threshols for cis-eQTLs
FDRtrans <- 0.2 # FDR threshold for trans-eQTLs
window <- 1 # cis-window (in Mb)

## analyze_cluster_GO.R
gene.bg <- "gcc" # background set for Gene Ontology analysis (gcc=genes in eQTL network, all=all genes used in the eQTL analysis)

## gwas_core_score_Qi_LRT.R
term <- "autoimmmune" # traits or diseases for which core-score will be plotted
tissue.spe <- "whole_blood" # tissue corresponding to previous term 

## create_gwas_summary.R
gwas.diseases.file <- "gwas_snps_by_diseasestype.txt"
