# =============================================================================
# HyperTRIBE Pipeline
# Usage:
#   Rscript hyperTRIBER2_run.R path/to/hyperTRIBER2_config.R
# The script ships with the package; locate it with
#   system.file("scripts", "hyperTRIBER2_run.R", package = "hyperTRIBER2")
# =============================================================================

suppressPackageStartupMessages({
  library(hyperTRIBER2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- 0. Parse config path from command line  ----------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "hyperTRIBER2_config.R"
source(config_path)
stranded <- config$stranded %||% FALSE

message("=== HyperTRIBER pipeline: ", config$run_name, " ===")
out_dir <- file.path(config$res_dir, config$run_name, "hyperTRIBER2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(config_path, out_dir)

# ---- 1. Load data -----------------------------------------------------------
design_df <- read_tsv(config$design_file)


txi  <- load_salmon(config$salmon_dir, design_df)

adar_cat_vec <- categorise_adar(txi, config$adar_tx_name)

mp   <- load_mpileup(config$mpileup_dir)
bam_list_file <- file.path(config$mpileup_dir, "bam_list.txt")
mpileup_df <- mp$mpileup_df
locsGR     <- mp$locsGR

annot  <- load_annotations(config$gtf_file, config$gff_file,
                            config$keep_gtf_types)
gtfGR  <- annot$gtfGR
ids    <- annot$ids

# ---- 2.  count data -------------------------------------------------
message("ADAR expression category per sample:")
print(adar_cat_vec)
# Attach ADAR category to design
design_df <- design_df %>%
  left_join(
    tibble(sample = names(adar_cat_vec), ADAR_CAT = adar_cat_vec),
    by = "sample"
  )

unique_samp_names <- unique(design_df$sample)

all_samp_names <- readLines(bam_list_file) |>
  basename() |>
  stringr::str_remove("(_Aligned.sortedByCoord.out)?(_fwd|_rev)?\\.bam$")

data_list_all <- extract_count_data(mpileup_df, unique_samp_names, all_samp_names,
                                    stranded = stranded)
# Genome ref base per row, in the same orientation as the counts
refBase       <- data_list_all[[1]]$ref

out <- build_design_and_restrict(
  design_df          = design_df,
  data_list_all      = data_list_all,
  refBase            = refBase,
  edits_of_interest  = config$edits_of_interest,
  min_samp_treat     = config$min_samp_treat,
  min_count          = config$min_count,
  min_prop           = config$min_prop,
  both_ways          = config$both_ways
)

design_vectors        <- out$design_vectors
data_restricted_lists <- out$data_restricted_lists

# Wide-format design table for downstream steps
design <- design_df %>%
  select(-ADAR_CAT) %>%
  pivot_wider(
    id_cols    = experiment,
    names_from = condition,
    values_from = sample,
    values_fn  = list
  )
design$data_list <- data_restricted_lists[design$experiment]



# ---- 3. Generate DEXSeq count files ----------------------------------------

write_count_files(data_restricted_lists, design_vectors, out_dir,
                  stranded = stranded)

# ---- 4. Run DEXSeq ----------------------------------------------------------

dxd_list <- run_dexseq(design_vectors, out_dir,
                       ncores = config$ncores_dexseq, fdr = config$fdr)

design$dxd_list <- dxd_list[design$experiment]

obj <- design
save_checkpoint(obj, out_dir, "after_make_test")

# ---- 4bis. Run GLM per-site test -------------------------------------------

# glm_list <- lapply(setNames(design$experiment, design$experiment), function(exp) {
#   dl      <- design$data_list[design$experiment == exp][[1]]
#   dv      <- design_vectors[[exp]]
#   treat   <- names(dv[dv == "treat"])
#   control <- names(dv[dv == "control"])
#   # make_glm_test(dl, treat, control, family = config$glm_family %||% "quasibinomial")
#   make_glm_test(dl, treat, control, family = "binomial")
# })

# design$glm_res_binom <- glm_list[design$experiment]

# glm_list <- lapply(setNames(design$experiment, design$experiment), function(exp) {
#   dl      <- design$data_list[design$experiment == exp][[1]]
#   dv      <- design_vectors[[exp]]
#   treat   <- names(dv[dv == "treat"])
#   control <- names(dv[dv == "control"])
#   make_glm_test(dl, treat, control, family = "quasibinomial")
# })

# design$glm_res_quasi <- glm_list[design$experiment]

# obj <- design
# save_checkpoint(obj, out_dir, "after_glm_test")


# ---- 5. Call edit sites -----------------------------------------------------

posGR_list <- call_hits(
  dxd_list              = dxd_list,
  design_vectors        = design_vectors,
  data_restricted_lists = data_restricted_lists,
  locsGR                = locsGR,
  edits_of_interest     = config$edits_of_interest,
  fdr                   = config$fdr,
  symmetric             = config$symmetric %||% FALSE,
  ncores                = config$ncores_hits,
  stranded              = stranded
)

design$posGR <- posGR_list[design$experiment]
obj <- design
save_checkpoint(obj, out_dir, "after_getHits")

# ---- 6. Filter to ADAR edits ------------------------------------------------

posGR_list_filtered <- filter_edits(posGR_list, config$keep_edit_pairs, config$ep_thresh, config$fc_thresh, config$FvF)
design$posGR_filtered <- posGR_list_filtered[design$experiment]
obj <- design
save_checkpoint(obj, out_dir, "after_posGR_list_filtered")

# ---- 7. Annotate with gene information --------------------------------------

posGR_genes <- annotate_hits(
  posGR_list_filtered = posGR_list_filtered,
  design              = design,
  txi                 = txi,
  gtfGR               = gtfGR,
  ids                 = ids
)

design$posGR_genes <- posGR_genes[design$experiment]
obj <- design
save_checkpoint(obj, out_dir, "after_addGenes")

message("=== Pipeline complete: ", config$run_name, " ===")
