# choros

choros is a software pipeline to estimate and correct the sequence biases 
present in ribosome profiling datasets. A preprint describing the method in 
greater detail can be found [on bioRxiv](https://doi.org/10.1101/2023.02.21.529452).

Briefly, choros performs the following:

1. Annotation of RPF alignments with: codon index corresponding to the ribosome 
A-site codon position; codons present in the ribosome A-, P-, and E-site 
positions; RPF 5' and 3' digest lengths (sequence length between RPF ends and 
A site); and RPF 5' and 3' bias sequences (nucleotides at RPF ends).

2. Generation of a training dataset of RPF counts on transcripts with high RPF
density for regression modeling.

3. Fitting of a negative binomial regression, modeling RPF counts on A-/P-/E-site 
codon identities, RPF 5' and 3' digest lengths, 5' and 3' bias sequences, and 
RPF GC-content.

4. Bias correction using scaling factors computed from regression modeling.

5. Computation and visualization of feature importance measurements per codon 
or nucleotide position, as performed in [iXnos](http://dx.doi.org/10.1038/s41594-018-0080-2)
to evaluate sequence biases before and after bias correction.

## Installation

1. Install [R](https://www.r-project.org/)

2. Install dependencies

Packages from CRAN:

```
install.packages(c("parallel", "doParallel", "foreach", "reshape", 
                   "ggplot2", "MASS", "patchwork", "RColorBrewer", 
                   "data.table", "devtools"))
```

Packages from Bioconductor:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsamtools")
BiocManager::install("Biostrings")
```

3. Install `choros`

```
devtools::install_github("lareaulab/choros", ref = "main",
                         build_vignettes = TRUE)
```

## Input files

- A .bam file of alignments of pre-processed ribosome protected fragments (RPFs)
- A FASTA file of transcript sequences (including 5' and 3' UTR regions)
- A transcriptome lengths file containing 5' UTR, CDS, and 3' UTR lengths
- An A-site offset rules file containing A-site offsets per read length and frame

## Running choros

A vignette showing usage and output of choros on a small sample dataset is included:
```
vignette("choros_vignette")
```

## Pre-processing

choros assumes that all read pre-processing steps have been performed by the 
user prior to regression modeling. These steps may include:

- adapter trimming
- filtering out of contaminant reads (i.e. ribosomal and noncoding RNAs)
- alignment to the transcriptome

In the case of reads with UMIs, the authors recommend additional UMI 
deduplication steps. See the [choros paper github repo](https://github.com/lareaulab/choros_paper)
for example scripts.

The authors additionally recommend computing posterior mapping weights for 
multi-mapping reads with [RSEM](https://deweylab.github.io/RSEM/). 
