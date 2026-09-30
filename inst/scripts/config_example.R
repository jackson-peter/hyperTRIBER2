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
  
  # --- Design ---
  # Design table (TSV) columns: sample, experiment, condition ("control"/"treat")
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
  ep_thresh=0.8,
  
  # --- GTF feature types to retain for annotation ---
  keep_gtf_types = c("3UTR", "5UTR", "CDS", "exon",
                     "start_codon", "stop_codon")
)
