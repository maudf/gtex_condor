#!/bin/bash

gene_dir=$1
code_dir=$2
tissuesfile=$3

if [ ! -e "$gene_dir" ]
then
    echo "Missing intput folder..."; exit
fi

if [ -e "$code_dir$tissuesfile" ]
then
    rm "$code_dir$tissuesfile"
fi

list_files=$(ls "$gene_dir"*norm.tsv)
while read line 
do
    basename $line | sed -e 's/_norm.tsv//g' >>"$code_dir$tissuesfile"
done <<<"$list_files"
