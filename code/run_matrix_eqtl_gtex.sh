#!/bin/bash

code_dir=$1
tissuesfile=$2

if [ ! -d $code_dir/log ]
then
   mkdir $code_dir/log
fi
   
var=$(cat $code_dir$tissuesfile)
while read line
do
    Rscript $code_dir/matrix_eqtl_gtex.sh $line >$code_dir/log/eqtl_calculation_"$line".log &
done <<<"var"
