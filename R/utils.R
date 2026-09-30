#' Load and name salmon quantification files
#'
#' @param salmon_dir  Top-level directory containing per-sample *_quant/ folders
#' @param design A data frame with sample information (must contain a column named "sample")
#' @param sample_col Column name in design that contains sample names (default: "sample")
#' @param quant_suffix Suffix for quant files (default: "quant.sf")
#' @return tximport object with named columns
#' @export
load_salmon <- function(salmon_dir, design, sample_col = "sample", quant_suffix="quant.sf") {
  if (!sample_col %in% colnames(design)) stop("sample_col not found in design")

  valid_samples <- design[[sample_col]]

  # Get all quant.sf files with full paths
  quant_files <- list.files(
    salmon_dir,
    pattern = paste0( quant_suffix, "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(quant_files) == 0)
    stop("No quant.sf files found in: ", salmon_dir)

  # Extract sample names from directory structure
  dir_names <- basename(dirname(quant_files))
  samp_names <- sub("_quant$", "", dir_names)

  named_files <- stats::setNames(quant_files, samp_names)
  # Filter to only include samples in design
  named_files <- named_files[names(named_files) %in% valid_samples]

  # Identify missing samples
  missing_samples <- setdiff(valid_samples, names(named_files))
  if (length(missing_samples) > 0) {
    stop("Missing salmon files for samples in design: ", paste(missing_samples, collapse = ", "))
  }

  # Import in design order (a sample may appear in several experiments)
  named_files <- named_files[unique(valid_samples)]
  txi <- tximport::tximport(
    files = named_files,
    type = "salmon",
    txOut = TRUE
  )

  message("Successfully loaded ", length(named_files), " samples matching design file")
  txi
}



#' Categorise ADAR expression into low / high bins per sample
#'
#' @param txi tximport object
#' @param adar_tx_name Row name of the ADAR transgene transcript (default: "ADARclone")
#' @return Named character vector ("low"/"high") with sample names
#' @export
categorise_adar <- function(txi, adar_tx_name = "ADARclone") {

  if (!adar_tx_name %in% rownames(txi$abundance))
    stop("'", adar_tx_name, "' not found in txi abundance matrix.")

  adar      <- txi$abundance[adar_tx_name, ]
  adar_cat  <- cut(scale(log2(adar + 1)), 2)
  levels(adar_cat) <- c("low", "high")

  names(adar_cat)  <- names(adar)

  adar_cat
}

#' Load mpileup-derived text files and build matching GRanges
#'
#' @param mpileup_dir Directory containing output_*.txt files
#' @return list(mpileup_df, locsGR)
#' @export
load_mpileup <- function(mpileup_dir) {

  files <- list.files(mpileup_dir, pattern = "output_.*.txt$", full.names = TRUE)

  if (length(files) == 0)
    stop("No mpileup output files found in: ", mpileup_dir)

  df <- do.call(rbind, lapply(files, utils::read.table, header = FALSE,
                              colClasses = c(V1 = "character", V3 = "character")))
  df$V3 <- toupper(df$V3)
  rownames(df) <- paste(df[[1]], df[[2]], sep = "_")

  locsGR <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(df$V1),
    ranges   = IRanges::IRanges(df$V2, width = 1),
    ref      = df$V3,
    names    = paste(df$V1, "_", df$V2, sep = "")
  )
  names(locsGR) <- locsGR$names

  message("Mpileup: ", nrow(df), " positions loaded")

  list(mpileup_df = df, locsGR = locsGR)
}
#' Import and filter GTF + GFF annotations
#'
#' @param gtf_file Path to GTF file
#' @param gff_file Path to GFF3 file
#' @param keep_gtf_types Feature types to retain from GTF (default: c("3UTR","5UTR","CDS","exon","start_codon","stop_codon"))
#' @return list(gtfGR, genes, ids)
#' @export
load_annotations <- function(gtf_file, gff_file,
                             keep_gtf_types = c("3UTR","5UTR","CDS","exon",
                                                "start_codon","stop_codon")) {

  gtf    <- rtracklayer::import(gtf_file)
  gtfGR  <- gtf[S4Vectors::mcols(gtf)$type %in% keep_gtf_types]

  gff    <- rtracklayer::import(gff_file)
  genes  <- gff[gff$type == "gene"]

  # Gene symbol from `symbol` (Araport), else `Name`, else the gene ID
  gene_meta <- S4Vectors::mcols(genes)
  gene_id   <- as.character(gene_meta$ID)
  symbol_col <- intersect(c("symbol", "Name"), colnames(gene_meta))
  gene_symbol <- if (length(symbol_col) > 0) {
    as.character(gene_meta[[symbol_col[1]]])
  } else {
    gene_id
  }
  gene_symbol[is.na(gene_symbol)] <- gene_id[is.na(gene_symbol)]

  ids <- stats::setNames(gene_symbol, gene_id)

  message("Annotations: ", length(gtfGR), " GTF features | ",
          length(genes), " genes in GFF")

  list(gtfGR = gtfGR, genes = genes, ids = ids)
}

#' Build per-experiment design vectors and restrict count data
#'
#' @param design_df Tibble with columns: sample, experiment, condition
#' @param data_list_all Output of extract_count_data()
#' @param refBase Reference base vector aligned to data_list rows
#' @param edits_of_interest Edit matrix (12 x 2)
#' @param min_samp_treat Minimum treated samples with signal (default: 2)
#' @param min_count Minimum read count (default: 5)
#' @param min_prop Minimum edit proportion (default: 0.0)
#' @param both_ways Logical (default: TRUE)
#' @return list(design_vectors, data_restricted_lists)
#' @export
build_design_and_restrict <- function(design_df,
                                      data_list_all,
                                      refBase,
                                      edits_of_interest,
                                      min_samp_treat = 2,
                                      min_count      = 5,
                                      min_prop       = 0.0,
                                      both_ways      = TRUE) {

  bad_cond <- setdiff(unique(design_df$condition), c("control", "treat"))
  if (length(bad_cond) > 0) {
    stop("design condition must be 'control' or 'treat', found: ",
         paste(bad_cond, collapse = ", "))
  }

  design_lists <- design_df |>
    dplyr::group_by(experiment) |>
    dplyr::summarise(
      design = list(stats::setNames(condition, sample)),
      .groups = "drop"
    )

  design_vectors <- stats::setNames(design_lists$design, design_lists$experiment)

  data_restricted_lists <- purrr::imap(
    design_vectors,
    ~ restrict_data(
      data_list          = data_list_all[names(.x)],
      ref_base           = refBase,
      design_vector      = .x,
      min_samp_treat     = min_samp_treat,
      min_count          = min_count,
      min_prop           = min_prop,
      both_ways          = both_ways,
      edits_of_interest  = edits_of_interest
    )
  )

  list(design_vectors        = design_vectors,
       data_restricted_lists = data_restricted_lists)
}

#' Generate DEXSeq count files for all experiments
#'
#' @param data_restricted_lists Named list from restrict_data()
#' @param design_vectors Named list of named condition vectors
#' @param res_dir Top-level results directory
#' @param stranded Logical (default: FALSE)
#' @export
write_count_files <- function(data_restricted_lists,
                              design_vectors,
                              res_dir,
                              stranded = FALSE) {

  purrr::iwalk(
    data_restricted_lists,
    function(data_list, exper) {
      out <- file.path(res_dir, exper, "results")
      dir.create(out, recursive = TRUE, showWarnings = FALSE)

      generate_count_files(
        data_list     = data_list,
        design_vector = design_vectors[[exper]],
        out_dir       = out,
        stranded      = stranded
      )
      message("Count files written for: ", exper)
    }
  )
}
#' Run DEXSeq (make_test) for all experiments
#'
#' @param design_vectors Named list of named condition vectors
#' @param res_dir Top-level results directory
#' @param ncores Number of cores for DEXSeq (default: 10)
#' @param fdr FDR threshold used only for the progress message (default: 0.1)
#' @return Named list of DEXSeq result objects
#' @export
run_dexseq <- function(design_vectors, res_dir, ncores = 10, fdr = 0.1) {

  purrr::imap(
    design_vectors,
    function(design_vec, exper) {
      message("Running DEXSeq for: ", exper)

      dxd_res <- make_test(
        out_dir       = file.path(res_dir, exper, "results"),
        design_vector = design_vec,
        n_cores       = ncores
      )

      n_sig <- sum(dxd_res$padj < fdr, na.rm = TRUE)
      message("  -> ", n_sig, " sites with padj < ", fdr, " in ", exper)

      dxd_res
    }
  )
}
#' Call significant edit sites (getHits) for all experiments
#'
#' @param dxd_list Named list of DEXSeq results
#' @param design_vectors Named list of condition vectors
#' @param data_restricted_lists Named list from restrict_data()
#' @param locsGR GRanges of all mpileup positions
#' @param edits_of_interest Edit matrix
#' @param fdr FDR threshold (default: 0.1)
#' @param symmetric Logical (default: FALSE)
#' @param ncores Cores for getHits (default: 40)
#' @param stranded Logical (default: FALSE)
#' @return Named list of GRanges with edit sites
#' @export
call_hits <- function(dxd_list,
                      design_vectors,
                      data_restricted_lists,
                      locsGR,
                      edits_of_interest,
                      fdr = 0.1,
                      symmetric = FALSE,
                      ncores = 40,
                      stranded = FALSE) {

  purrr::imap(
    dxd_list,
    function(dxd_res, exper) {
      message("Running getHits for: ", exper)

      get_hits(
        dexseq_res        = dxd_res,
        stranded          = stranded,
        fdr               = fdr,
        n_cores           = ncores,
        add_meta          = TRUE,
        data_list         = data_restricted_lists[[exper]],
        edits_of_interest = edits_of_interest,
        design_vector     = design_vectors[[exper]],
        include_ref       = TRUE,
        ref_gr            = locsGR,
        symmetric         = symmetric
      )
    }
  )
}
#' Filter posGR list
#'
#' @param posGR_list Named list of GRanges from call_hits()
#' @param keep_edit_pairs Character vector of "ref targ" strings to keep (default: c("A G", "T C"))
#' @param ep_thresh Edit proportion threshold (default: 0.8)
#' @param fc_thresh Absolute log2 fold change threshold (default: 0)
#' @param FvF Logical; if FALSE keep only sites more edited in treat (default: FALSE)
#' @return Filtered named list of GRanges. Sites with NA proportion or fold
#'   change are dropped.
#' @export
filter_edits <- function(posGR_list,
                         keep_edit_pairs = c("A G", "T C"),
                         ep_thresh = 0.8,
                         fc_thresh = 0,
                         FvF = FALSE) {

  purrr::imap(
    posGR_list,
    function(posGR, exper) {
      posGR <- posGR[paste(posGR$ref, posGR$targ) %in% keep_edit_pairs]
      posGR <- posGR[which(posGR$prop < ep_thresh)]
      if (!FvF) {
        posGR <- posGR[which(posGR$fold_change > 0)]
      }
      if (fc_thresh != 0) {
        posGR <- posGR[which(abs(posGR$fold_change) > fc_thresh)]
      }

      message(exper, ": ", length(posGR), " sites after edit filter")
      posGR
    }
  )
}

#' Add gene annotations to filtered edit sites
#'
#' @param posGR_list_filtered Filtered named list of GRanges
#' @param design Wide-format design tibble (with control/treat cols)
#' @param txi tximport object
#' @param gtfGR Filtered GTF GRanges
#' @param ids Named vector gene_id -> gene_symbol
#' @return Named list of annotated GRanges. Unstranded A>G sites are only
#'   matched to + strand genes and T>C sites to - strand genes.
#' @export
annotate_hits <- function(posGR_list_filtered,
                          design,
                          txi,
                          gtfGR,
                          ids) {

  purrr::imap(
    posGR_list_filtered,
    function(posGR, exper_name) {
      if (is.null(posGR) || length(posGR) == 0) {
        message(exper_name, ": no sites to annotate, skipping")
        return(NULL)
      }

      i <- match(exper_name, design$experiment)
      ctrl_samples <- design$control[[i]]
      quant_vec    <- rowMeans(txi$abundance[, ctrl_samples, drop = FALSE], na.rm = TRUE)

      annotate_with_genes(
        pos_gr        = addStrandForHyperTRIBE(posGR),
        gtf_gr        = gtfGR,
        gene_ids      = ids,
        quant         = quant_vec,
        n_cores       = 1,
        assign_strand = TRUE,
        flank_distance = 1000
      )
    }
  )
}

#' Convenience wrapper: save a checkpoint .Rdat file
#'
#' @param obj Object to save
#' @param res_dir Directory
#' @param tag File name tag (e.g. "after_getHits")
#' @export
save_checkpoint <- function(obj, res_dir, tag) {
  path <- file.path(res_dir, paste0("design_", tag, ".Rdat"))
  save(obj, file = path)
  message("Checkpoint saved: ", path)
}
