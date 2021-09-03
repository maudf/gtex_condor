# Compute selection scores from VCF files

## Introduction
This document describes all the data files and scripts necessary to
replicates the selection scores computed in the PATTERNS project.

## General settings

### List of necessary R packages and softwares
Running the following scripts requires to have the following softwares and packages installed:

#### Softwares
Software               | Where to find it  
---------------------- | ---------------------------------------------------  
plink2 v1.90 or higher | [https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)  
R                      | [https://cran.r-project.org/](https://cran.r-project.org/)  
htslib v1.12 or higher | [http://www.htslib.org](http://www.htslib.org)  
SelectionHapStats      | [https://github.com/ngarud/SelectionHapStats/](https://github.com/ngarud/SelectionHapStats)  

#### R packages
Package        | Where to find it  
-------------- | ----------------------------------------------------------  
pcadapt     | [https://cran.r-project.org/web/packages/pcadapt/index.html](https://cran.r-project.org/web/packages/pcadapt/index.html)  
getopt         | [https://cran.r-project.org/web/packages/getopt/index.html](https://cran.r-project.org/web/packages/getopt/index.html)  
rehh        | [https://cran.r-project.org/web/packages/rehh/index.html](https://cran.r-project.org/web/packages/rehh/index.html) 
RColorBrewer   | [https://cran.r-project.org/web/packages/RColorBrewer/index.html](https://cran.r-project.org/web/packages/RColorBrewer/index.html) 
vcfR     | [https://cran.r-project.org/web/packages/vcfR/index.html](https://cran.r-project.org/web/packages/vcfR/index.html) 
PopGenome       | [https://cran.r-project.org/web/packages/PopGenome/index.html](https://cran.r-project.org/web/packages/PopGenome/index.html)
broom          | [https://cran.r-project.org/web/packages/broom/](https://cran.r-project.org/web/packages/broom/)  
metap          | [https://cran.r-project.org/web/packages/metap/](https://cran.r-project.org/web/packages/metap/)  
R Bioconductor | [http://bioconductor.org](http://bioconductor.org)  

#### R Bioconductor packages
Package        | Where to find it  
-------------- | -----------------------------------------------------------  
Biobase        | [http://bioconductor.org/packages/release/bioc/html/Biobase.html](http://bioconductor.org/packages/release/bioc/html/Biobase.html)  
limma          | [http://bioconductor.org/packages/release/bioc/html/limma.html](http://bioconductor.org/packages/release/bioc/html/limma.html)  


## List of required data files
Running these scripts requires to have the following data files:
* genotypes in phased VCF formats.

## Scripts pipeline
To compute the selection scores, please run the following pipeline.

### Prepare files
All files containing the sequences to be analyzed need to be stored in
a Results/ folder located in the Selection_scores/ folder.
For example, you can link the Results folder obtained from the
simulations:
```{bash qc, eval=FALSE}
ln -s $PWD/../Simulations/Results $PWD/
```
### Compute neutrality statistics and FST  for each VCF in a folder
This script computes some population genetics statistics, using the PopGenome R
packages on one or several VCF files:
* For each population: Tajima's <em>D</em>, Fu and Li's <em>D</em>
and <em>F</em>, Zheng's <em>E</em>, Fay and Wu's <em>H</em>
* For all pairs of populations: <em>F<sub>ST</sub></em>.

In this example, the statistics are computed on the Neutral and
Directional Opt=10 simulations. For each series of simulations, the
VCF files contains 100 samples from 2 populations. 30 independent loci
of 100kb were simulated, and SNPs of interest were located at position
50000 in each chromosomes. No outgroup is provided in the VCF file.
Results are saved in a RDS file in the directory where the VCF files
are stored.
**Warning:** These stats are computed by chromosome. Here are the
examples for chromosome 1

```{bash qc, eval=FALSE}
Rscript code/neutrality_stats.R -d ./Results/Neutral/ES_0.78_f_0.05/ -c 1 -s -a no -p 100,200 -w 50000 -l 30
Rscript code/neutrality_stats.R -d ./Results/Directional/Opt_10_ES_0.78_f_0.05/ -c 1 -s -a no -p 100,200 -w 50000 -l 30
```

### Compute iHS for each VCF in a folder
 This script computes <em>iHS</em> for each population and each SNP using
the rehh R package on one or several VCF files.

In this example, the <em>iHS</em>  is computed on the Neutral and
Directional Opt=10 simulations. For each series of simulations, the
VCF files contains 100 samples from 2 populations. 30 independent loci
of 100kb were simulated, and SNPs of interest were located at position
50000 in each chromosomes. No outgroup is provided in the VCF file.
Results are saved in a RDS file in the directory where the VCF files
are stored.

```{bash qc, eval=FALSE}
Rscript code/iHS.R -d ./Results/Neutral/ES_0.78_f_0.05/ -s -a no -p 100,200 -w 50000 -l 30
Rscript code/iHS.R -d ./Results/Directional/Opt_10_ES_0.78_f_0.05/ -c 1 -s TRUE -a no -p 100,200 -w 50000 -l 30
```

### Compute <em>PCadapt</em>
This script computes <em>PCadapt</em> scores for each population and each SNP using
the pcadapt R package on one or several BED/BIM/FAM files. For this,
you need to first convert the files using the convert\_VCF\_BED.sh
script (see ../Selection_scores/README.md).

```{bash qc, eval=FALSE}
Rscript code/neutrality_stats.R -d ./Results/Neutral/ES_0.78_f_0.05/ 
Rscript code/neutrality_stats.R -d ./Results/Directional/Opt_10_ES_0.78_f_0.05/ 
```

### Compute <em>H</em>12
This section explains how to compute H12 on converted files using
[https://cran.r-project.org/] (SelectionHapStats).
The H12 must be run for each position of interest using the following
command (here the example is for if SelectionHapStats installed in a
~/Documents/Softwares/ folder ):
```{bash qc, eval=FALSE}
SelectionHapStats=~/Documents/Softwares/SelectionHapStats/ # Change if necessary
python ${SelectionHapStats}/H12_H2H1.py Results/Directional/Opt_10_ES_0.78_f_0.05/Rep_1_OPT_10_ES_0.78_f_0.05_Directional.h12.txt 400 --singleWindow 50000 --window 200 >> Results/Directional/Opt_10_ES_0.78_f_0.05/Rep_1_OPT_10_ES_0.78_f_0.05_Directional.results.h12.txt 
```

It is possible to systematize it for all QTL for each replicates of each set of simulations using the run\_H12.sh script:
```{bash qc, eval=FALSE}
Rscript code/run_H12.sh -p ~/Documents/Softwares/SelectionHapStats/ -d ./Results/Neutral -s 0.78 -f 0.05 -d ./Results/Neutral/
```

