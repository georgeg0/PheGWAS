# PheGWAS: Three-dimensional approach to dynamically visualize GWAS across multiple phenotypes

# Introduction:
GWAS is Genome Wide Association Study which is done on many variants and one phenotype. The Manhattan plot is the popular way to visualize this( https://github.com/stephenturner/qqman ). On the other hand, PheWAS (https://rdrr.io/github/PheWAS/PheWAS/ ) is phenome wide association study which is done on a single variant and many phenotypes. PheGWAS is a new concept to do analysis on many variants and many phenotypes. This is a package to do interactive 3D visualization on many variants many phenotypes.

More details are found in the preprint: https://www.biorxiv.org/content/10.1101/694794v1. PheGWAS is an ongoing project from Department of Population Health and Genomics, University of Dundee. Efforts are being made to add new features that would make it more beneficial for comprehensive gene-disease association analysis.


# Installation:
## Prerequisites
Please Install plotly, devtools and Biomart before installing the PheGWAS package.

Install plotly
```
install.packages("plotly")
```

Install devtools
```
install.packages("devtools")
```

Install BioMart
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```
## Installation of PheGWAS package
```
install_github("georgeg0/PheGWAS")
```

Find vignette here
https://github.com/georgeg0/PheGWAS/blob/master/vignettes/PheGWAS-vignette.html

# Brief Documentation:
## Interactive PheGWAS landscape mode features:
•	X axis: Chromosome (level 1) and base pair positions (level 2). 

•	Z axis: The phenotypes name. 

•	Y axis: -log10(p-value) or other values such as the effect size can be specified. 

•	Info pop-up: Hovering over a point reveals information about the plotted SNP, the information includes Phenotype name, kbp position range, P value, Locus, SNP ID, Gene. The GENE name is mapped to the corresponding SNPID within the package by making use of the R package “BioMArt”.  

•	Interactive features like zoom in, zoom out, turntable rotation, rotation on an axis and pan scale features.  

•	The figure can be exported as static or interactive plot (which is recommended), in order to share with colleagues.  

## Exploring the GWAS summary data with PheGWAS:
### Preprocessing the dataframe’s:
The first stage of the exploration is to prepare the dataframe to pass to the PheGWAS exploration functions. The preparation of the dataframe is done using the function “processphegwas”. 

List of the GWAS summary stats files as dataframe’s are passed to this function, the GWAS summary stats file is to be of this format. 

If gene to rsid is not available. Gene will be mapped internally with the BioMart Module
```
CHR        BP       rsid <pvalue of the phenotype>
```
If gene to rsid is available, this could be passed as an input parameter.
```
CHR        BP       rsid gene  <pvalue of the phenotype> 
```
Here the p-value column name must be replaced by phenotype name, this column name is taken for the y axis of the landscape.

The PheGWAS exploration is undertaken in different levels (as of now its 2 different levels). The exploration is done using function “landscape”

### Level 1: Exploration on entire genome level:
The first level is the entire genome level. In this level the entire chromosomes are plotted in a 3D landscape view with each peak in the chromosome representing the highest pvalue in that chromosome. The X axis here is the chromosome number. Z axis and Y axis are phenotype name and -log10(p-value) respectively.

With this visualization the user gets an idea which chromosome contribute to multiple phenotypes. A chromosome contributing to multiple phenotypes has to be explored, as we are more interested to see if a certain region within the chromosome is contributing to multiple phenotypes. This is where the exploration is taken to the next level of the PheGWAS.

### Level 2: Exploration on individual chromosomes:
The second level is the single chromosome level. In this level the entire chromosome is divided into different groups based on bp units given by the user (100kbp is the default). The x axis here is the BP groups in the chromosome. Z axis and Y axis are phenotype name and -log10(p-value) respectively.

# Examples:

Checkout vignette here https://github.com/georgeg0/PheGWAS/blob/master/vignettes/PheGWAS-vignette.html
