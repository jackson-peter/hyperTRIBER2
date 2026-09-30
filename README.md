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
  "tximport"
))
install.packages(c("dplyr", "purrr", "tidyr", "readr", "glue", "doParallel", "foreach"))

# Install hyperTRIBER2 package
devtools::install_github("jackson-peter/hyperTRIBER2")

```

## Pipeline Overview

hyperTRIBER2 provides a complete pipeline for differential RNA editing analysis:
 - Data Loading: Import and preprocess data from Salmon quantification and mpileup files
 - Site Calling: Identify significant editing sites using DEXSeq
 - Filtering: Apply quality filters to editing sites
 - Annotation: Add gene and transcript annotations to significant sites. For unstranded
   data, A>G sites are matched only to + strand genes and T>C sites only to − strand genes
   (A-to-I on the transcript); sites with no gene on the matching strand get `gene = NA`.
   Gene symbols come from the GFF `symbol` attribute, else `Name`, else the gene ID.

## Quick start

### 1. Prepare Data

The following files are required:
 - Salmon quantification files: Output from running Salmon on RNA data
 - Mpileup files: Base-level counts from sequencing data
 - Reference annotation: GTF/GFF files
 - Design file: a TSV with one row per sample and experiment, with columns
   `sample` (matching the salmon `<sample>_quant/` folders and the BAM names in
   `mpileup_dir/bam_list.txt`), `experiment` (samples are compared within each
   experiment) and `condition` (`control` or `treat`)

Example design table:

```
sample	experiment	condition
R22	FvC	control
R23	FvC	control
R24	FvC	treat
R25	FvC	treat
```

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

  FvF      = FALSE,   # FALSE: keep only sites more edited in treat
  stranded = FALSE,   # TRUE only for 8-column-per-sample (strand-split) mpileup

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

The pipeline script is installed with the package; a template config is next to it
(`inst/scripts/config_example.R` in this repository).

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts", "hyperTRIBER2_run.R", package = "hyperTRIBER2"))')" hyperTRIBER2_config.R
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
```R
edits <- rbind(c("A", "G"), c("T", "C"))   # ref -> target matrix

# Per-sample base counts from the mpileup table
data_list_all <- extract_count_data(mp_data$mpileup_df, samp_names, all_samp_names)

# Build per-experiment design vectors and restrict to candidate sites
restricted <- build_design_and_restrict(
  design_df         = design_df,
  data_list_all     = data_list_all,
  refBase           = data_list_all[[1]]$ref,
  edits_of_interest = edits,
  min_samp_treat    = 2,
  min_count         = 5,
  min_prop          = 0.0,
  both_ways         = TRUE
)
design_vectors        <- restricted$design_vectors
data_restricted_lists <- restricted$data_restricted_lists

# DEXSeq test per experiment
write_count_files(data_restricted_lists, design_vectors, res_dir)
dxd_list <- run_dexseq(design_vectors, res_dir, ncores = 10)

# Call significant editing sites
posGR_list <- call_hits(
  dxd_list              = dxd_list,
  design_vectors        = design_vectors,
  data_restricted_lists = data_restricted_lists,
  locsGR                = mp_data$locsGR,
  edits_of_interest     = edits,
  fdr                   = 0.1,
  ncores                = 40
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
For a complete pipeline example, see `inst/scripts/hyperTRIBER2_run.R`. This script demonstrates the full workflow from data loading to final annotation.


### Configuration Parameters Explained
 - run_name: Name for this analysis run (used in output directories) &rarr; "hyperTRIBER2_example_run"
 - design_file: Path to TSV file describing samples and conditions &rarr; "path/to/your/design_table.tsv"
 - res_dir: Directory where results will be saved &rarr; "path/to/your/results_directory"
 - mpileup_dir: Directory containing mpileup output files &rarr; "path/to/your/mpileup/output"
 - salmon_dir: Directory containing Salmon quantification files &rarr; "path/to/your/salmon/quantification"
 - gtf_file: Path to GTF annotation file &rarr; "path/to/your/reference/annotation.gtf"
 - gff_file: Path to GFF annotation file &rarr; "path/to/your/reference/annotation.gff"
 - adar_tx_name: Name of ADAR transgene in Salmon quantification &rarr; "ADARclone"
 - FvF: If FALSE, keep only sites more edited in treat than control &rarr; FALSE
 - stranded: Whether the mpileup has strand-split counts (8 columns per sample) &rarr; FALSE
 - min_samp_treat: Minimum number of replicates with at least one edited read (applied per condition) &rarr; 2
 - min_count: Minimum total edited-base (target) read count across replicates &rarr; 5
 - min_prop: Minimum editing proportion, used when both_ways = TRUE &rarr; 0.0
 - both_ways: Keep sites passing the filters in either condition, not only in treat &rarr; TRUE
 - fdr: False discovery rate threshold for significance &rarr; 0.1
 - ncores_dexseq: Number of cores for DEXSeq analysis &rarr; 10
 - ncores_hits: Number of cores for hit calling &rarr; 40
 - edits_of_interest: All substitution types to consider initially &rarr; All 12 types
 - keep_edit_pairs: Edit types to keep after initial filtering &rarr; c("A G", "T C")
 - fc_thresh: Minimum absolute DEXSeq log2 fold change &rarr; 1
 - ep_thresh: Sites with treat editing proportion at or above this are removed (likely SNPs) &rarr; 0.8
 - keep_gtf_types: GTF feature types to use for annotation &rarr; c("3UTR", "5UTR", "CDS", "exon", ...)

