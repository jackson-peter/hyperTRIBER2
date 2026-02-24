# hyperTRIBER2

R package for differential RNA editing analysis — an updated version of [hyperTRIBER](https://github.com/sarah-ku/hyperTRIBER).

## Summary

hyperTRIBER2 is an R package for detecting sites with significant differential RNA editing between conditions. It was originally developed for **hyperTRIBE** (targets of RNA binding proteins identified by editing), where an RBP is fused to a hyperactive ADAR domain, enabling detection of RBP-bound transcripts through A-to-I editing. The package is equally applicable to general differential RNA editing analyses.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DEXSeq", "GenomicRanges", "BiocParallel"))

library(devtools)
install_github("jackson-peter/hyperTRIBER2")
```

## Quick start

### 1. Generate base pileup

Before running the R package, generate base counts from your BAM files using `samtools mpileup` and the provided Perl script:

```bash
samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f reference_genome.fa \
  Samp1.bam Samp2.bam Samp3.bam | perl hyperTRIBE_mpileup2bases.pl > baseCounts.txt
```

### 2. Run the pipeline

```r
library(hyperTRIBER2)

data <- read.table("baseCounts.txt", header = FALSE)

samp.names <- c("Samp1", "Samp2", "Samp3", "Samp4", "Samp5", "Samp6")
design_vector <- c(Samp1="control", Samp2="control", Samp3="control",
                   Samp4="treat",   Samp5="treat",   Samp6="treat")

data_list <- extractCountData(data, samp.names, strand = FALSE)

# Filter sites
data_list <- restrict_data(data_list, design_vector, min_samp_treat = 2, min_count = 2, both_ways = FALSE)

# Fit model and get hits
generateCountFiles(data_list, out_dir = "./results/", design_vector = design_vector)
dxd.res <- make_test(out_dir = "./results/", design_vector = design_vector, ncores = 10)
posGR <- getHits(res = dxd.res, fdr = 0.1, design_vector = design_vector, data_list = data_list)

# Annotate
posGR <- addGenes(gtfGR = gtf, posGR = posGR, ncore = 10)
```


