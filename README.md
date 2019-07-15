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
If you downloaded the tar version of the package
```
install.packages("<path to the tar>",dependencies = TRUE, repos = NULL)
```
Eg:
```
install.packages("/Downloads/PheGWAS_0.1.0.tgz",dependencies = TRUE, repos = NULL)
```

Or if you have the unzipped version of the package then

```
devtools::install("<path to package folder name>", dependencies=TRUE)
```
Eg:
```
devtools::install("/Downloads/PheGWAS-master",dependencies = TRUE)
```

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
The 6 different data frames from 2 different consortium are made available in the package for getting used to the package and its usage.

## Following processed summary data are from the lipid consortium:
hdl - HDL GWAS summary dataset
ldl - LDL GWAS summary dataset
trigs - TRIG GWAS summary dataset
tchol - TOTAL CHOLESTROL GWAS summary dataset

Results are from H. Grallert et al., “Discovery and refinement of loci associated with lipid levels,” Nat. Genet., vol. 45, no. 11, pp. 1274–1283, 2013.

```
> head(hdl)
  CHR        BP       rsid     HDL                                                          gene
1   7  92383888       rs10 0.76460                                                          CDK6
2  12 126890980  rs1000000 0.33750 RP4-809F18.1;RP4-809F18.2;RP5-944M2.1;RP5-944M2.2;RP5-944M2.3
3   4  21618674 rs10000010 0.25460                            KCNIP4;RP11-556G22.1;RP11-556G22.3
4   4   1357325 rs10000012 0.26640                AC078852.1;AC078852.2;CRIPAK;MAEA;NKX1-1;UVSSA
5   4  37225069 rs10000013 0.63740                                              KIAA1239;MIR4801
6   4  84778125 rs10000017 0.05778                                                    RP11-8L2.1
```
```
> head(ldl)
  CHR        BP       rsid     LDL                                                          gene
1   7  92383888       rs10 0.03411                                                          CDK6
2  12 126890980  rs1000000 0.51210 RP4-809F18.1;RP4-809F18.2;RP5-944M2.1;RP5-944M2.2;RP5-944M2.3
3   4  21618674 rs10000010 0.13150                            KCNIP4;RP11-556G22.1;RP11-556G22.3
4   4   1357325 rs10000012 0.89260                AC078852.1;AC078852.2;CRIPAK;MAEA;NKX1-1;UVSSA
5   4  37225069 rs10000013 0.26410                                              KIAA1239;MIR4801
6   4  84778125 rs10000017 0.89340                                                    RP11-8L2.1
```
```
> head(trig)
  CHR        BP       rsid   TRIGS                                                          gene
1   7  92383888       rs10 0.39990                                                          CDK6
2  12 126890980  rs1000000 0.08781 RP4-809F18.1;RP4-809F18.2;RP5-944M2.1;RP5-944M2.2;RP5-944M2.3
3   4  21618674 rs10000010 0.26940                            KCNIP4;RP11-556G22.1;RP11-556G22.3
4   4   1357325 rs10000012 0.25140                AC078852.1;AC078852.2;CRIPAK;MAEA;NKX1-1;UVSSA
5   4  37225069 rs10000013 0.86350                                              KIAA1239;MIR4801
6   4  84778125 rs10000017 0.69120                                                    RP11-8L2.1
```
```
> head(tchol)
  CHR        BP       rsid TOTAL_CHOLESTROL                                                          gene
1   7  92383888       rs10          0.03129                                                          CDK6
2  12 126890980  rs1000000          0.99050 RP4-809F18.1;RP4-809F18.2;RP5-944M2.1;RP5-944M2.2;RP5-944M2.3
3   4  21618674 rs10000010          0.01953                            KCNIP4;RP11-556G22.1;RP11-556G22.3
4   4   1357325 rs10000012          0.94260                AC078852.1;AC078852.2;CRIPAK;MAEA;NKX1-1;UVSSA
5   4  37225069 rs10000013          0.24640                                              KIAA1239;MIR4801
6   4  84778125 rs10000017          0.98070                                                    RP11-8L2.1
```
## Following processed summary data are from the giant consortium:

bmimen - BMI for men 
bmiwomen - BMI for women

Results are from Randall JC, Winkler TW, Kutalik Z, Berndt SI, Jackson AU, Monda KL, Kilpeläinen TO, Esko T, Mägi R, Li S, et al. (2013). Sex-stratified genome-wide association studies including 270,000 individuals show sexual dimorphism in genetic loci for anthropometric traits. PLoS Genet 9: e1003500.
```
> head(bmimen)
  CHR        BP       rsid BMI_men
1  10   9960129  rs4747841  0.9400
2  10   9960259  rs4749917  0.9400
3  10 100012739   rs737656  0.7000
4  10 100012890   rs737657  0.7600
5  10 100013438 rs17524355  0.5700
6  10 100013563  rs7086391  0.0051
```
```
> head(bmiwomen)
  CHR        BP       rsid BMI_women
1  10   9960129  rs4747841      0.31
2  10   9960259  rs4749917      0.31
3  10 100012739   rs737656      0.28
4  10 100012890   rs737657      0.27
5  10 100013438 rs17524355      0.47
6  10 100013563  rs7086391      0.72
```

The dataframe’s are passed to processphegwas function as a list of dataframe’s.
```
library(PheGWAS)
x <- list(hdl,ldl,trig,tchol)
## This preprocess the dataframes’ for passing to the landscape function.
y <- processphegwas(x)

## pass the dataframe from the processphegwas
# 3D landscape visualization of all the phenotypes across the base pair #positions
landscape(y)

# 3D landscape visualization of all the phenotypes across the base pair #positions, showing only the crusts above a certain threshold
landscape(y,sliceval = 10)

# 3D landscape visualization of chromosome number 16
landscape(x,sliceval = 10,chromosome = 16)

# 3D landscape visualization of all the phenotypes across the base pair #positions, showing only the crusts that are within a certain region (between sliceval and upperlimit)
landscape(y,sliceval = 10 ,chromosome = 16,upperlimit = 50)

# 3D landscape visualization of chromosome number 16, with bp division on #100,000 (default is 1,000,000)
landscape(x ,chromosome = 16,bpdivision = 100,000)

# 3D landscape visualization of chromosome number 16, with bp division on #1,000,000 (default is 1,000,000), gene view active (refer paper for more details)

landscape(d, chromosome=16,geneview = TRUE)

## 3-D association table for many phenotypes
landscapetable(d, sliceval = 10, chromosome=3)
```
