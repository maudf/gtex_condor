#!/bin/bash
#Run this script to merge snp/gene annotation information. 
#FDR cut-off for cis and trans results
fdr_cis=$1
fdr_trans=$2
window=$3
tissuesfile=$4

gene_dir=$5
eqtl_dir=$6
lookup=$7

# Merge cis and trans results

while read tissues
do
  awk 'NR==1 {print $0"\t""cis.or.trans"}' $eqtl_dir$tissues'_'$window'MB_cis.txt' > $eqtl_dir$tissues$fdr_cis$fdr_trans'_'$window'MB.txt'
  awk -v afdr="$fdr_cis" '{if($6 < afdr) {print $0"\t""cis"}}' $eqtl_dir$tissues'_'$window'MB_cis.txt' >> $eqtl_dir$tissues$fdr_cis$fdr_trans'_'$window'MB.txt'
  awk -v afdr="$fdr_trans" '{if($6 < afdr) {print $0"\t""trans"}}' $eqtl_dir$tissues'_'$window'MB_trans.txt' >> $eqtl_dir$tissues$fdr_cis$fdr_trans'_'$window'MB.txt'
  #Join the snp and gene position information
  echo $tissues
  join -a2 -1 1 -2 2 <(tail -n +2 $gene_dir$tissues'_genes.tsv' | sort -k1,1 ) <(tail -n +2 $eqtl_dir$tissues$fdr_cis$fdr_trans'_'$window'MB.txt' | sort -k2,2) > $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp'
  join -a2 -1 3 -2 5 <(tail -n +2 $lookup | sort -k3,3 ) <(sort -k5,5 $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp') > $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp2' 
  cat <( paste -d"\t" <(head -1 $lookup | awk -v OFS="\t" -F"\t" '{print $3,$2,$1,$4,$5,$6,$7,$8}') <( head -1 $gene_dir$tissues'_genes.tsv' ) <( head -1 $eqtl_dir$tissues$fdr_cis$fdr_trans'_'$window'MB.txt' | cut -f3-)) $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp2' > $eqtl_dir$tissues'fdr'$fdr_cis$fdr_trans'_'$window'MB_edges.txt'
  rm $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp' $tissues$fdr_cis$fdr_trans'_'$window'MB.tmp2'
done < $tissuesfile
