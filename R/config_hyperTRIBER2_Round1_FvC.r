# =============================================================================
# HyperTRIBER2 Pipeline Configuration
# =============================================================================

config <- list(
  
  # --- Run identity ---
  run_name = "hyperTRIBER2_FirstRound_F2C_TESTING",   # used for logging/messages
  
  # --- Paths ---
  design_file  = "/projects/renlab/people/qgr178/Projects/AGO1/Scripts/FirstRound/AGO1_design_table.tsv",
  res_dir      = "/maps/projects/renlab/people/qgr178/Projects/AGO1/FirstRound/Results",
  mpileup_dir  = "/projects/renlab/people/qgr178/Projects/AGO1/FirstRound/Results/variant/mpileup/",
  salmon_dir   = "/projects/renlab/people/qgr178/Projects/AGO1/FirstRound/Results/quantification/salmon/",
  gtf_file     = "/maps/projects/renlab/people/qgr178/Shared/Reference/A_thaliana/Araport11/Araport11_GTF_genes_transposons.20241001.gtf",
  gff_file     = "/maps/projects/renlab/people/qgr178/Shared/Reference/A_thaliana/Araport11/Araport11_GFF3_genes_transposons.20250813.gff",
  
  FvF=FALSE,

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
