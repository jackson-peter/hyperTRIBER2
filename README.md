# hyperTRIBER2

R package for differential RNA editing analysis — an updated version of [hyperTRIBER](https://github.com/sarah-ku/hyperTRIBER).

## Summary

hyperTRIBER2 is an R package for detecting sites with significant differential RNA editing between conditions. It was originally developed for **hyperTRIBE** (targets of RNA binding proteins identified by editing), where an RBP is fused to a hyperactive ADAR domain, enabling detection of RBP-bound transcripts through A-to-I editing. The package is equally applicable to general differential RNA editing analyses.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install dependencies
BiocManager::install(c(
  "DEXSeq",
  "GenomicRanges",
  "BiocParallel",
  "rtracklayer",
  "tximport",
  "dplyr",
  "purrr",
  "tidyverse"
))

# Install hyperTRIBER2 package
devtools::install_github("jackson-peter/hyperTRIBER2")

```

## Pipeline Overview

hyperTRIBER2 provides a complete pipeline for differential RNA editing analysis:
 - Data Loading: Import and preprocess data from Salmon quantification and mpileup files
 - Site Calling: Identify significant editing sites using DEXSeq
 - Filtering: Apply quality filters to editing sites
 - Annotation: Add gene and transcript annotations to significant sites

## Quick start

### 1. Prepare Data

The following files are required:
 - Salmon quantification files: Output from running Salmon on RNA data
 - Mpileup files: Base-level counts from sequencing data
 - Reference annotation: GTF/GFF files
 - Design file: A table describing samples and experimental conditions

### 2. Configuration file

Create an R file with settings like the example below. Replace placeholder paths with your actual file locations.

```R
# =============================================================================
# HyperTRIBER2 Pipeline Configuration
# =============================================================================

config <- list(

  # --- Run identity ---
  run_name = "hyperTRIBER2_example_run",   # used for logging/messages

  # --- Paths ---
  design_file  = "path/to/your/design_table.tsv",
  res_dir      = "path/to/your/results_directory",
  mpileup_dir  = "path/to/your/mpileup/output",
  salmon_dir   = "path/to/your/salmon/quantification",
  gtf_file     = "path/to/your/reference/annotation.gtf",
  gff_file     = "path/to/your/reference/annotation.gff",

  FvF = FALSE,

  # --- ADAR transgene name in the salmon quantification matrix ---
  adar_tx_name = "ADARclone",

  # --- restrict_data parameters ---
  min_samp_treat = 2,
  min_count      = 5,
  min_prop       = 0.0,
  both_ways      = TRUE,

  # --- Statistical thresholds ---
  fdr  = 0.1,

  # --- Parallelisation ---
  ncores_dexseq = 10,
  ncores_hits   = 40,

  # --- Edit types to consider ---
  # All 12 substitution types; filter to A>G / T>C happens downstream
  edits_of_interest = rbind(
    c("A","G"), c("G","A"),
    c("T","C"), c("C","T"),
    c("A","T"), c("T","A"),
    c("G","T"), c("T","G"),
    c("G","C"), c("C","G"),
    c("A","C"), c("C","A")
  ),

  # --- Final edit filter (ref targ pairs kept after getHits) ---
  # ADAR-specific: A>G on sense strand == T>C on antisense
  keep_edit_pairs = c("A G", "T C"),
  fc_thresh = 1,
  ep_thresh = 0.8,

  # --- GTF feature types to retain for annotation ---
  keep_gtf_types = c("3UTR", "5UTR", "CDS", "exon",
                     "start_codon", "stop_codon")
)

```


### 3. Run the pipeline

```bash
# Run the pipeline script
Rscript hyperTRIBER2_run.R hyperTRIBER2_config.R
```

## Key Functions

### Data loading
```R
# Load Salmon quantification data
txi <- load_salmon(
  salmon_dir = "path/to/your/salmon/quantification",
  design = design_df,
  sample_col = "sample"
)

# Load mpileup data
mp_data <- load_mpileup(mpileup_dir = "path/to/your/mpileup/output")

# Load annotations
annot <- load_annotations(
  gtf_file = "path/to/your/reference/annotation.gtf",
  gff_file = "path/to/your/reference/annotation.gff"
)


```

### Site Analysis
```
# Restrict data to sites of interest
analysis_data <- prepare_analysis(
  data_list_all = data_list,
  design_lists = design_lists,
  refBase = "A",
  min_samp_treat = 2,
  min_count = 5,
  min_prop = 0.0,
  both_ways = TRUE,
  edits_of_interest = c("A", "G")
)

# Call significant editing sites
posGR_list <- call_hits(
  dxd_list = dxd_list,
  design_vectors = design_vectors,
  data_restricted_lists = data_restricted_lists,
  locsGR = locsGR,
  edits_of_interest = c("A", "G"),
  fdr = 0.1,
  ncores = 40
)

# Filter editing sites
posGR_list_filtered <- filter_edits(
  posGR_list = posGR_list,
  keep_edit_pairs = c("A G", "T C"),
  ep_thresh = 0.8,
  fc_thresh = 1,
  FvF = FALSE
)

# Annotate sites with gene information
posGR_genes <- annotate_hits(
  posGR_list_filtered = posGR_list_filtered,
  design = design,
  txi = txi,
  gtfGR = annot$gtfGR,
  ids = annot$ids
)

```

### Utility Functions
```
# Save intermediate results
save_checkpoint(
  obj = design,
  res_dir = "path/to/your/results_directory/example_run",
  tag = "after_getHits"
)
```

### Example Usage
For a complete pipeline example, see the hyperTRIBER2_run.R script in the package repository. This script demonstrates the full workflow from data loading to final annotation.


### Configuration Parameters Explained
 - run_name: Name for this analysis run (used in output directories) &rarr; "hyperTRIBER2_example_run"
 - design_file: Path to TSV file describing samples and conditions &rarr; "path/to/your/design_table.tsv"
 - res_dir: Directory where results will be saved &rarr; "path/to/your/results_directory"
 - mpileup_dir: Directory containing mpileup output files &rarr; "path/to/your/mpileup/output"
 - salmon_dir: Directory containing Salmon quantification files &rarr; "path/to/your/salmon/quantification"
 - gtf_file: Path to GTF annotation file &rarr; "path/to/your/reference/annotation.gtf"
 - gff_file: Path to GFF annotation file &rarr; "path/to/your/reference/annotation.gff"
 - adar_tx_name: Name of ADAR transgene in Salmon quantification &rarr; "ADARclone"
 - min_samp_treat: Minimum number of samples with editing in treatment group &rarr; 2
 - min_count: Minimum total count for a site to be considered &rarr; 5
 - min_prop: Minimum proportion threshold for editing detection &rarr; 0.0
 - both_ways: Whether to test both directions of each edit type &rarr; TRUE
 - fdr: False discovery rate threshold for significance &rarr; 0.1
 - ncores_dexseq: Number of cores for DEXSeq analysis &rarr; 10
 - ncores_hits: Number of cores for hit calling &rarr; 40
 - edits_of_interest: All substitution types to consider initially &rarr; All 12 types
 - keep_edit_pairs: Edit types to keep after initial filtering &rarr; c("A G", "T C")
 - fc_thresh: Fold change threshold for filtering &rarr; 1
 - ep_thresh: Editing proportion threshold for filtering &rarr; 0.8
 - keep_gtf_types: GTF feature types to use for annotation &rarr; c("3UTR", "5UTR", "CDS", "exon", ...)

