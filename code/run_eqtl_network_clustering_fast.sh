#!/bin/bash
# Run this script from the /cccbstore-rc/projects/gtex/eqtl_networks/ directory

fdrcis=$1
fdrtrans=$2
window=$3
tissues_file=$4

code_dir=$5
eqtl_dir=$6

while read tissue
do
	echo "Running tissue $tissue"
	Rscript $code_dir/eqtl_network_clustering_fast.R "$eqtl_dir""$tissue"fdr'$fdr_cis$fdr_trans'_'$window'MB_edges.txt FALSE
done <$tissues_file
