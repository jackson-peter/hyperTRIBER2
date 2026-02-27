#' Extract count data from mpileup output
#'
#' Parses the output of the RNAeditR_mpileup2bases.pl/RNAeditR_mpileup2bases_stranded.pl
#' perl script into a tidy per-sample list of base count tables.
#'
#' @param dat Data frame read from the outputs of RNAeditR_mpileup2bases.pl (no header).
#'   Expected columns: chr, pos, ref, then per sample 4 columns (A, T, C, G)
#'   or 8 columns if stranded (A, T, C, G, a, t, c, g).
#' @param samp_names Character vector of sample names.
#' @param stranded Logical. If TRUE, expects 8 columns per sample. Default FALSE.
#'
#' @return A named list of tibbles (one per sample), each with columns:
#'   chr, pos, ref, A, T, C, G, site_id.
#'
#' @importFrom dplyr bind_cols mutate filter pull if_any all_of bind_rows
#' @importFrom purrr map set_names reduce
#' @importFrom glue glue
#' @importFrom rlang set_names
#'
#' @export
extract_count_data <- function(dat, samp_names, stranded = FALSE) {

  # -- Constants --
  n_meta_cols     <- 3
  base_cols_fwd   <- c("A", "T", "C", "G")
  base_cols_rev   <- c("a", "t", "c", "g")
  cols_per_sample <- ifelse(stranded, 8, 4)

  # -- If all_samp_names not provided, assume dat matches samp_names --
  if (is.null(all_samp_names)) {
    all_samp_names <- samp_names
  }

  n_all_samples <- length(all_samp_names)

  # -- Validate total columns against ALL samples --
  expected_ncol <- n_meta_cols + (n_all_samples * cols_per_sample)
  if (ncol(dat) != expected_ncol) {
    stop(glue::glue(
      "Column mismatch: expected {expected_ncol} ",
      "({n_meta_cols} meta + {n_all_samples} samples x {cols_per_sample}), ",
      "but got {ncol(dat)}."
    ))
  }

  # -- Validate requested samples exist --
  missing <- setdiff(samp_names, all_samp_names)
  if (length(missing) > 0) {
    stop(glue::glue(
      "Samples not found in all_samp_names: {paste(missing, collapse = ', ')}"
    ))
  }

  # -- Subset columns for requested samples only --
  # Find indices of requested samples within the full sample list
  samp_indices <- match(samp_names, all_samp_names)

  cols_to_keep <- c(
    1:n_meta_cols,
    unlist(lapply(samp_indices, function(i) {
      start <- n_meta_cols + (i - 1) * cols_per_sample + 1
      start:(start + cols_per_sample - 1)
    }))
  )

  dat <- dat[, cols_to_keep]


  # -- Input validation --
  expected_ncol <- n_meta_cols + (n_samples * cols_per_sample)
  if (ncol(dat) != expected_ncol) {
    stop(glue::glue(
      "Column mismatch: expected {expected_ncol} ",
      "({n_meta_cols} meta + {n_samples} samples x {cols_per_sample}), ",
      "but got {ncol(dat)}."
    ))
  }

  message(glue::glue(
    "Selected {length(samp_names)}/{n_all_samples} samples from mpileup file."
  ))

  n_samples <- length(samp_names)


  # -- Extract and name metadata --
  meta <- dat[, 1:n_meta_cols] |>
    rlang::set_names(c("chr", "pos", "ref")) |>
    dplyr::mutate(pos = as.integer(pos), chr = as.character(chr), ref = as.character(ref))

  # -- Per-sample extraction --
  data_list <- purrr::map(seq_along(samp_names), function(i) {

    start_col <- n_meta_cols + (i - 1) * cols_per_sample + 1
    end_col   <- n_meta_cols + i * cols_per_sample

    sample_counts <- dat[, start_col:end_col] |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))

    if (stranded) {
      colnames(sample_counts) <- c(base_cols_fwd, base_cols_rev)

      fwd_total <- sum(sample_counts[, base_cols_fwd], na.rm = TRUE)
      rev_total <- sum(sample_counts[, base_cols_rev], na.rm = TRUE)
      message(glue::glue(
        "{samp_names[i]}: fwd/rev ratio = {round(fwd_total / rev_total, 3)}"
      ))

      # Forward strand
      fwd <- dplyr::bind_cols(meta, sample_counts[, base_cols_fwd]) |>
        dplyr::mutate(strand = "+",
                      site_id = paste0(chr, "_", pos, ",+"))

      # Reverse strand: swap a<->t, c<->g
      rev_counts <- sample_counts[, c("t", "a", "g", "c")]
      colnames(rev_counts) <- base_cols_fwd

      rev <- dplyr::bind_cols(meta, rev_counts) |>
        dplyr::mutate(strand = "-",
                      site_id = paste0(chr, "_", pos, ",-"))

      dplyr::bind_rows(fwd, rev)

    } else {
      colnames(sample_counts) <- base_cols_fwd

      dplyr::bind_cols(meta, sample_counts) |>
        dplyr::mutate(site_id = paste0(chr, "_", pos))
    }
  }) |>
    purrr::set_names(samp_names)

  # -- Remove sites with NA in any sample --
  na_sites <- data_list |>
    purrr::map(~ .x |>
                 dplyr::filter(dplyr::if_any(dplyr::all_of(base_cols_fwd), is.na)) |>
                 dplyr::pull(site_id)
    ) |>
    purrr::reduce(union)

  if (length(na_sites) > 0) {
    message(glue::glue("Removing {length(na_sites)} sites with NA values."))
    data_list <- data_list |>
      purrr::map(~ .x |> dplyr::filter(!site_id %in% na_sites))
  } else {
    message("No NA sites found. All sites retained.")
  }

  # -- Summary --
  n_sites <- nrow(data_list[[1]])
  message(glue::glue(
    "Extracted {n_samples} samples, {n_sites} sites each. ",
    "Stranded: {stranded}."
  ))

  return(data_list)
}


#' Filter sites by minimum replicate support and editing evidence
#'
#' @param data_list Named list of tibbles from extract_count_data
#' @param ref_base Character vector of reference bases (one per site)
#' @param design_vector Character vector: "control" or "treat" for each sample
#' @param min_samp_treat Minimum number of replicates with editing > 0
#' @param min_count Minimum total edited-base count across replicates
#' @param min_prop Minimum editing proportion (only used when both_ways = TRUE)
#' @param both_ways If TRUE, keep sites passing filters in either condition
#' @param edits_of_interest Matrix of ref->target pairs to consider
#' @return Filtered data_list
#' @export
restrict_data <- function(data_list,
                          ref_base = NULL,
                          design_vector,
                          min_samp_treat = 3,
                          min_count = 2,
                          min_prop = 0.01,
                          both_ways = FALSE,
                          edits_of_interest = rbind(c("A", "G"), c("T", "C"))) {

  bases <- c("A", "T", "C", "G")

  # -- Extract ref_base from data if not supplied --
  if (is.null(ref_base)) {
    ref_base <- data_list[[1]]$ref
    message("Using 'ref' column from first sample as ref_base.")
  }

  # -- Build count matrices per base per condition --
  # Each element is a sites x replicates matrix
  build_count_matrix <- function(condition, base) {
    # subset data_list in either treat or control
    data_list[design_vector == condition] |>

      # subset by base and bind each replicate in a matrix (rep in columns, site in row)
      purrr::map(~ .x[[base]]) |>
      do.call(what = cbind)
  }

  # build Count matrix for each base for treat
  treat_counts <- purrr::set_names(bases) |>
    purrr::map(~ build_count_matrix("treat", .x))

  # build Count matrix for each base for control
  cont_counts <- purrr::set_names(bases) |>
    purrr::map(~ build_count_matrix("control", .x))

  # -- For each edit type, determine which sites to keep --
  # keep_per_edit is a list of vectors (one per edit type) with one element per site
  # specifying if there is editing and should be kept with our filtering options
  keep_per_edit <- seq_len(nrow(edits_of_interest)) |>
    purrr::map(function(i) {
      my_ref  <- edits_of_interest[i, 1]
      my_targ <- edits_of_interest[i, 2]

      # Editing proportions (target / (target + ref))
      treat_targ_total <- rowSums(treat_counts[[my_targ]])
      treat_ref_total  <- rowSums(treat_counts[[my_ref]])
      cont_targ_total  <- rowSums(cont_counts[[my_targ]])
      cont_ref_total   <- rowSums(cont_counts[[my_ref]])

      prop_treat   <- treat_targ_total / (treat_targ_total + treat_ref_total)
      prop_control <- cont_targ_total  / (cont_targ_total + cont_ref_total)
      prop_treat[is.nan(prop_treat)]     <- 0
      prop_control[is.nan(prop_control)] <- 0

      # Site cannot be its own target (ref != target)
      not_target_ref <- ref_base != my_targ

      if (!both_ways) {
        # WARNING HERE: original used >= 0 which is always TRUE for counts
        # Using > 0 here = replicates with at least 1 edited read
        reps_with_editing <- rowSums(treat_counts[[my_targ]] > 0)
        enough_reps  <- reps_with_editing >= min_samp_treat
        enough_count <- treat_targ_total >= min_count

        not_target_ref & enough_reps & enough_count

      } else {
        # Treatment
        treat_reps   <- rowSums(treat_counts[[my_targ]] > 0) >= min_samp_treat
        treat_pass   <- treat_reps & (treat_targ_total >= min_count) & (prop_treat > min_prop)

        # Control
        cont_reps    <- rowSums(cont_counts[[my_targ]] > 0) >= min_samp_treat
        cont_pass    <- cont_reps & (cont_targ_total >= min_count) & (prop_control > min_prop)

        not_target_ref & (treat_pass | cont_pass)
      }
    })

  # -- Keep site if ANY edit type passes --
  to_keep <- keep_per_edit |>
    purrr::reduce(`|`)

  n_before <- nrow(data_list[[1]])
  data_list <- data_list |>
    purrr::map(~ dplyr::slice(.x, which(to_keep)))

  n_after <- nrow(data_list[[1]])
  message(glue::glue(
    "restrict_data: {n_before} -> {n_after} sites ",
    "({n_before - n_after} removed). ",
    "both_ways={both_ways}, min_samp_treat={min_samp_treat}, ",
    "min_count={min_count}, min_prop={min_prop}."
  ))

  return(data_list)
}


#' Identify potential editing sites vs reference genome
#'
#' For each edit type (e.g. A->G), finds sites where the reference base
#' matches, enough samples show the ref and target bases, and few samples
#' show other (non-edit) bases.
#'
#' @param data_list List of count data frames (one per sample)
#' @param ref_gr GRanges with a `ref` column for reference bases
#' @param stranded Logical, whether data is strand-aware
#' @param edits_of_interest Matrix of ref->target base pairs
#' @param sample_type Which group to scan ("treat" or "control")
#' @param design_vector Named character vector of sample assignments
#' @param min_ref Minimum samples with ref base > 0
#' @param min_targ Minimum samples with target base > 0
#' @param max_other Maximum samples with non-ref/non-target bases > 0
#' @return Named list (per edit type) of site IDs passing filters
get_edits_vs_reference <- function(data_list,
                                   ref_gr,
                                   stranded = FALSE,
                                   edits_of_interest,
                                   sample_type = "treat",
                                   design_vector,
                                   min_ref = 3,
                                   min_targ = 3,
                                   max_other = 2) {

  dnames <- row.names(data_list[[1]])

  # -- 1. Parse positions and assign reference bases -------------------------
  pos_dat <- parse_site_positions(dnames, stranded)
  pos_dat$ref <- lookup_reference_base(pos_dat, ref_gr, stranded)

  # -- 2. Extract count matrices per base for target samples -----------------
  target_samples <- names(design_vector[design_vector == sample_type])

  base_counts <- c("A", "T", "C", "G") |>
    purrr::set_names() |>
    purrr::map(function(base) {
      mat <- do.call(cbind, lapply(data_list, function(x) x[, base]))
      row.names(mat) <- dnames
      mat[, target_samples, drop = FALSE]
    })

  # -- 3. For each edit type, find sites passing filters ---------------------
  poss_list <- seq_len(nrow(edits_of_interest)) |>
    purrr::set_names(apply(edits_of_interest, 1, paste, collapse = ":")) |>
    purrr::map(function(i) {
      a_ref  <- edits_of_interest[i, 1]
      a_targ <- edits_of_interest[i, 2]
      others <- setdiff(c("A", "T", "C", "G"), c(a_ref, a_targ))

      # Only look at sites where reference == expected ref base
      ref_sites <- which(pos_dat$ref == a_ref)

      count_nonzero_samples <- function(mat) {
        apply(mat, 1, function(x) sum(x > 0))
      }

      n_ref   <- count_nonzero_samples(base_counts[[a_ref]][ref_sites, , drop = FALSE])
      n_targ  <- count_nonzero_samples(base_counts[[a_targ]][ref_sites, , drop = FALSE])
      n_oth1  <- count_nonzero_samples(base_counts[[others[1]]][ref_sites, , drop = FALSE])
      n_oth2  <- count_nonzero_samples(base_counts[[others[2]]][ref_sites, , drop = FALSE])

      keep <- (n_ref >= min_ref) &
        (n_targ >= min_targ) &
        (n_oth1 <= max_other) &
        (n_oth2 <= max_other)

      names(which(keep))
    })

  return(poss_list)
}

#' Parse site names into chr, pos, (strand) data frame
#' @keywords internal
parse_site_positions <- function(dnames, stranded = FALSE) {
  if (stranded) {
    data.frame(
      chr    = gsub("(.+)_[0-9]+,[-+]", "\\1", dnames),
      pos    = as.numeric(gsub(".+_([0-9]+),[-+]", "\\1", dnames)),
      strand = gsub(".+_[0-9]+,([-+])", "\\1", dnames),
      row.names = dnames,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      chr = gsub("(.+)_[0-9]+", "\\1", dnames),
      pos = as.numeric(gsub(".+_([0-9]+)", "\\1", dnames)),
      row.names = dnames,
      stringsAsFactors = FALSE
    )
  }
}

#' Look up reference bases from GRanges, handling strand complementation
#' @keywords internal
lookup_reference_base <- function(pos_dat, ref_gr, stranded = FALSE) {

  names(ref_gr) <- paste(as.vector(GenomeInfoDb::seqnames(ref_gr)),
                         BiocGenerics::start(ref_gr), sep = "_")

  keys <- paste(pos_dat$chr, pos_dat$pos, sep = "_")
  ref_bases <- toupper(as.vector(ref_gr[keys]$ref))

  # Complement for minus strand
  if (stranded) {
    comp <- c(A = "T", T = "A", C = "G", G = "C")
    minus <- pos_dat$strand == "-"
    ref_bases[minus] <- comp[ref_bases[minus]]
  }

  ref_bases
}

#' Compute per-site editing proportions vs reference
#'
#' For each site, sums counts across all samples and returns the proportion
#' of reads supporting each base (A, T, C, G) plus total coverage.
#'
#' @param data_list List of count data frames (one per sample)
#' @param ref_gr GRanges with a `ref` column
#' @param stranded Logical
#' @return data.frame with chr, pos, (strand), ref, total, edit.A/T/C/G
editing_proportions_vs_reference <- function(data_list,
                                             ref_gr,
                                             stranded = FALSE) {

  dnames <- row.names(data_list[[1]])

  # -- 1. Parse positions and look up reference base -------------------------
  pos_dat <- parse_site_positions(dnames, stranded)
  pos_dat$ref <- lookup_reference_base(pos_dat, ref_gr, stranded)

  # -- 2. Sum counts across all samples per base -----------------------------
  bases <- c("A", "T", "C", "G")
  count_mat <- vapply(bases, function(b) {
    rowSums(do.call(cbind, lapply(data_list, function(x) x[, b])))
  }, numeric(length(dnames)))
  colnames(count_mat) <- bases

  # -- 3. Compute proportions ------------------------------------------------
  total <- rowSums(count_mat)
  props <- count_mat / total  # NaN where total == 0, that's fine

  pos_dat$total  <- total
  pos_dat$edit.A <- props[, "A"]
  pos_dat$edit.T <- props[, "T"]
  pos_dat$edit.C <- props[, "C"]
  pos_dat$edit.G <- props[, "G"]

  pos_dat
}

#' Extract significant editing hits from DEXSeq results
#'
#' Filters DEXSeq results by FDR, resolves multiple edit types per site,
#' and computes editing proportions in treat/control groups.
#'
#' @param dexseq_res DEXSeqResults data frame with padj, groupID, featureID, etc.
#' @param stranded Logical
#' @param fdr FDR threshold
#' @param n_cores Number of parallel cores
#' @param add_meta Logical, whether to compute per-site metadata
#' @param data_list List of per-sample count data frames
#' @param edits_of_interest Matrix of ref->target pairs
#' @param design_vector Named vector of sample group assignments
#' @param include_ref Logical, add reference base annotation
#' @param ref_gr GRanges with ref column
#' @param symmetric Logical, if TRUE use the "more edited" group for proportion
#' @return GRanges of significant hits with metadata columns
#' @export
get_hits <- function(dexseq_res,
                     stranded = FALSE,
                     fdr = 0.1,
                     n_cores = 20,
                     add_meta = TRUE,
                     data_list,
                     edits_of_interest = rbind(c("A", "G"), c("T", "C")),
                     design_vector,
                     include_ref = TRUE,
                     ref_gr = NULL,
                     symmetric = FALSE) {



  # -- 1. Filter results to edits of interest at given FDR -------------------
  target_bases <- unique(edits_of_interest[, 2])
  target_pattern <- paste0("E[", paste(target_bases, collapse = ""), "]")
  dexseq_res <- dexseq_res[grep(target_pattern, dexseq_res$featureID), ]
  sig_res <- dexseq_res[which(dexseq_res$padj <= fdr), ]
  hit_ids <- unique(sig_res$groupID)

  message(sprintf("Number of significant hits: %d", length(hit_ids)))

  # -- 2. Build GRanges of hit positions ------------------------------------
  pos_dat <- parse_site_positions(hit_ids, stranded)
  pos_gr <- GenomicRanges::GRanges(
    seqnames = pos_dat$chr,
    ranges   = IRanges::IRanges(start = pos_dat$pos, width = 1),
    strand   = if (stranded) pos_dat$strand else "*"
  )
  names(pos_gr) <- hit_ids

  # -- 3. Without metadata, just assign top edit base and return -------------
  if (!add_meta) {
    sig_res <- sig_res[order(sig_res$padj), ]
    top_feature <- tapply(sig_res$featureID, sig_res$groupID,
                          function(x) x[1])
    pos_gr$edit_base <- NA_character_
    pos_gr[names(top_feature)]$edit_base <- gsub("^E", "", top_feature)
    return(pos_gr)
  }

  # -- 4. Compute per-site metadata in parallel ------------------------------
  treat_idx  <- names(which(design_vector == "treat"))
  ctrl_idx   <- names(which(design_vector == "control"))

  if (n_cores > 1) {
    doParallel::registerDoParallel(cores = n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    meta_list <- foreach::foreach(i = seq_along(pos_gr)) %dopar% {
      compute_site_meta(
        site_id            = names(pos_gr)[i],
        sig_res            = sig_res,
        data_list          = data_list,
        edits_of_interest  = edits_of_interest,
        design_vector      = design_vector,
        treat_idx          = treat_idx,
        ctrl_idx           = ctrl_idx,
        symmetric          = symmetric
      )
    }
  } else {
    meta_list <- lapply(seq_along(pos_gr), function(i) {
      compute_site_meta(
        site_id            = names(pos_gr)[i],
        sig_res            = sig_res,
        data_list          = data_list,
        edits_of_interest  = edits_of_interest,
        design_vector      = design_vector,
        treat_idx          = treat_idx,
        ctrl_idx           = ctrl_idx,
        symmetric          = symmetric
      )
    })
  }

  meta <- do.call(rbind, meta_list)

  # -- 5. Optionally add reference base from GRanges ------------------------
  if (include_ref && !is.null(ref_gr)) {
    meta$base <- lookup_reference_base(pos_dat, ref_gr, stranded)
  }

  # -- 6. Attach metadata and filter self-edits -----------------------------
  if ("strand" %in% colnames(meta)) {
    colnames(meta)[colnames(meta) == "strand"] <- "gtf_strand"
  }
  GenomicRanges::mcols(pos_gr) <- meta
  pos_gr$ref  <- as.vector(pos_gr$ref)
  pos_gr$targ <- as.vector(pos_gr$targ)

  # Remove sites where ref == target (not real edits)
  pos_gr <- pos_gr[pos_gr$ref != pos_gr$targ]

  pos_gr
}


#' Compute metadata for a single hit site
#' @keywords internal
compute_site_meta <- function(site_id, sig_res, data_list,
                              edits_of_interest, design_vector,
                              treat_idx, ctrl_idx, symmetric) {

  bases <- c("A", "T", "C", "G")

  # Stack all samples for this site — index by site_id column, keep only base cols
  site_counts <- do.call(rbind, lapply(data_list, function(x) {
    row <- x[x$site_id == site_id, bases, drop = FALSE]
    if (nrow(row) == 0) return(setNames(rep(NA_integer_, 4), bases))
    row
  }))
  rownames(site_counts) <- names(data_list)

  # Determine reference as dominant base (original behavior)
  col_totals <- colSums(site_counts, na.rm = TRUE)
  ref <- names(which.max(col_totals))

  # Get significant hits for this site
  hit <- sig_res[sig_res$groupID == site_id, , drop = FALSE]

  # Resolve multiple edit types: keep only valid ref->target edits
  if (nrow(hit) > 1) {
    valid_targets <- edits_of_interest[edits_of_interest[, 1] == ref, 2]
    hit_targets <- gsub("^E", "", hit$featureID)
    keep <- hit_targets %in% valid_targets
    hit <- hit[keep, , drop = FALSE]

    if (nrow(hit) > 1) {
      message(sprintf("Multiple edits at %s, selecting most significant",
                      site_id))
      hit <- hit[which.min(hit$pvalue), , drop = FALSE]
    }
  }

  targ <- gsub("^E", "", hit$featureID)

  # Editing proportions: target / (ref + target)
  treat_ref  <- sum(site_counts[treat_idx, ref],  na.rm = TRUE)
  treat_targ <- sum(site_counts[treat_idx, targ], na.rm = TRUE)
  ctrl_ref   <- sum(site_counts[ctrl_idx, ref],   na.rm = TRUE)
  ctrl_targ  <- sum(site_counts[ctrl_idx, targ],  na.rm = TRUE)
  fc <- hit$log2fold_treat_control

  if (symmetric && fc < 0) {
    prop      <- ctrl_targ  / (ctrl_ref + ctrl_targ)
    prop_ctrl <- treat_targ / (ctrl_ref + ctrl_targ)
  } else {
    prop      <- treat_targ / (treat_ref + treat_targ)
    prop_ctrl <- ctrl_targ  / (ctrl_ref + ctrl_targ)
  }

  # Per-sample edit counts
  edit_counts <- vapply(data_list, function(x) {
    row <- x[x$site_id == site_id, targ, drop = TRUE]
    if (length(row) == 0) NA_integer_ else row
  }, numeric(1))

  data.frame(
    name         = site_id,
    ref          = ref,
    targ         = targ,
    prop         = prop,
    prop_ctrl    = prop_ctrl,
    padj         = hit$padj,
    pvalue       = hit$pvalue,
    control_par  = hit$control,
    treat_par    = hit$treat,
    fold_change  = fc,
    tags_treat   = sum(edit_counts[treat_idx], na.rm = TRUE),
    tags_control = sum(edit_counts[ctrl_idx],  na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}



addStrandForHyperTRIBE <- function(posGR)
{
  strand <- rep("*",length(posGR))
  strand[posGR$ref=="T" & posGR$targ=="C"] <- "-"
  strand[posGR$ref=="A" & posGR$targ=="G"] <- "+"
  strand(posGR) <- Rle(strand)
  return(posGR)
}

handleAnoType <- function(x) {
  x_uniq <- unique(x)
  # Remove redundant "exon" when a more specific type exists
  specific <- intersect(x_uniq, c("3UTR", "5UTR", "CDS"))
  if (length(specific) == 1) return(specific)
  if (length(x_uniq) == 1) return(x_uniq)
  stop("Ambiguous annotation types: ", paste(x_uniq, collapse = ", "))
}


#' Annotate editing sites with gene/transcript information from a GTF
#'
#' @param pos_gr GRanges of editing sites
#' @param gtf_gr GRanges from a parsed GTF (with $type, $gene_id, $transcript_id)
#' @param gene_ids Named vector mapping gene IDs to gene names
#' @param quant Named numeric vector of transcript-level expression (e.g. TPM),
#'        names should match transcript IDs in gtf_gr
#' @param assign_strand Logical; if TRUE, overwrite strand of pos_gr from GTF
#' @param n_cores Number of cores for parallel processing
#' @param flank_distance Distance to extend search if no direct overlap found
#' @return pos_gr with added metadata columns
#' @export
annotate_with_genes <- function(pos_gr,
                                gtf_gr,
                                gene_ids,
                                quant,
                                assign_strand = TRUE,
                                n_cores = 1,
                                flank_distance = 1000) {

  ## helper: resolve a single site ------------------------------------------
  annotate_single_site <- function(idx) {
    site <- pos_gr[idx]
    ols  <- gtf_gr[subjectHits(findOverlaps(site, gtf_gr, ignore.strand = FALSE))]
    out_of_range <- FALSE

    # If no direct overlap, try flanked window
    if (length(ols) == 0) {
      out_of_range <- TRUE
      ols <- gtf_gr[subjectHits(findOverlaps(site + flank_distance, gtf_gr,
                                             ignore.strand = FALSE))]
    }

    # No annotation at all
    if (length(ols) == 0) {
      return(list(
        info = c(gene = NA_character_, name = NA_character_, strand = "*",
                 transcripts = NA_character_, transcript_tpm = NA_character_,
                 transcript_types = NA_character_, out_of_range = as.character(out_of_range)),
        fprop      = NA_character_,
        new_strand = "*"
      ))
    }

    # --- Resolve gene (pick one if multiple) --------------------------------
    genes   <- unique(ols$gene_id)
    strands <- tapply(as.vector(strand(ols)), ols$gene_id, `[`, 1)

    if (length(genes) > 1) {
      genes <- resolve_gene_ambiguity(genes, ols, pos_gr, quant)
      strands <- strands[genes]
      ols <- ols[ols$gene_id %in% genes]
    }

    new_strand <- if (length(strands) > 0) strands[1] else "*"

    # --- Resolve transcript annotation types --------------------------------
    feature_type <- tapply(as.vector(ols$type), ols$transcript_id, handleAnoType)
    transcripts  <- names(feature_type)

    # --- Order transcripts by expression ------------------------------------
    tx_expr <- order_transcripts_by_expression(transcripts, quant)
    transcripts  <- names(tx_expr)
    feature_type <- feature_type[transcripts]

    # --- Compute relative position within feature ---------------------------
    fprop <- compute_feature_proportions(
      site, ols, feature_type, new_strand
    )

    gene_name <- gene_ids[genes]

    list(
      info = c(
        gene             = paste(genes, collapse = ","),
        name             = paste(gene_name, collapse = ","),
        strand           = paste(strands, collapse = ","),
        transcripts      = paste(transcripts, collapse = ","),
        transcript_tpm   = paste(round(tx_expr, 2), collapse = ","),
        transcript_types = paste(feature_type, collapse = ","),
        out_of_range     = out_of_range
      ),
      fprop      = paste(round(fprop, 4), collapse = ","),
      new_strand = new_strand
    )
  }

  ## run in parallel or serial -----------------------------------------------
  if (n_cores > 1) {
    doParallel::registerDoParallel(cores = n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    results <- foreach::foreach(
      i = seq_along(pos_gr),
    .packages = c("GenomicRanges", "IRanges", "S4Vectors")
    )%dopar% {
      annotate_single_site(i)
    }
  } else {
    results <- lapply(seq_along(pos_gr), annotate_single_site)
  }

  ## assemble results --------------------------------------------------------
  info_mat <- do.call(rbind, lapply(results, `[[`, "info"))

  pos_gr$gene             <- info_mat[, "gene"]
  pos_gr$name             <- info_mat[, "name"]
  pos_gr$gtf_strand       <- info_mat[, "strand"]
  pos_gr$transcripts      <- info_mat[, "transcripts"]
  pos_gr$transcript_tpm   <- info_mat[, "transcript_tpm"]
  pos_gr$transcript_types <- info_mat[, "transcript_types"]
  pos_gr$out_of_range     <- as.logical(info_mat[, "out_of_range"])
  pos_gr$feature_prop     <- vapply(results, `[[`, character(1), "fprop")

  if (assign_strand) {
    strand(pos_gr) <- vapply(results, `[[`, character(1), "new_strand")
  }

  pos_gr
}


#' Resolve gene ambiguity when a site overlaps multiple genes
#' @keywords internal
resolve_gene_ambiguity <- function(genes, ols, pos_gr, quant) {
  # Primary: pick gene with most overlapping features to all sites
  overlap_counts <- vapply(genes, function(g) {
    length(findOverlaps(ols[ols$gene_id == g], pos_gr))
  }, integer(1))

  if (sum(overlap_counts == max(overlap_counts)) == 1) {
    return(names(which.max(overlap_counts)))
  }

  # Tiebreak: highest total expression
  expr_sums <- vapply(genes, function(g) {
    matched <- grep(g, names(quant))
    if (length(matched) > 0) sum(quant[matched]) else 0
  }, numeric(1))

  genes[which.max(expr_sums)]
}


#' Order transcripts by expression, NAs last
#' @keywords internal
order_transcripts_by_expression <- function(transcripts, quant) {
  if (length(transcripts) == 0) {
    return(stats::setNames(NA_real_, ""))
  }

  tx_expr <- rep(NA_real_, length(transcripts))
  names(tx_expr) <- transcripts

  matched <- transcripts[transcripts %in% names(quant)]
  if (length(matched) > 0) {
    tx_expr[matched] <- quant[matched]
  }

  # Sort: expressed descending, then NAs
  tx_expr <- c(
    sort(tx_expr[!is.na(tx_expr)], decreasing = TRUE),
    tx_expr[is.na(tx_expr)]
  )
  tx_expr
}


#' Compute relative position of a site within its annotated feature
#' @keywords internal
compute_feature_proportions <- function(site, ols, feature_type, strand_val) {
  if (strand_val == "*" || length(ols) == 0) {
    fprop <- rep(NA_real_, length(feature_type))
    names(fprop) <- names(feature_type)
    return(fprop)
  }

  site_pos <- start(site)

  vapply(seq_along(feature_type), function(j) {
    tx_id <- names(feature_type)[j]
    ft    <- feature_type[j]
    feat  <- ols[ols$transcript_id == tx_id & ols$type == ft]

    if (length(feat) == 0) return(NA_real_)

    prop <- (site_pos - start(feat)) / (end(feat) - start(feat))
    if (strand_val == "-") prop <- 1 - prop
    mean(prop)
  }, numeric(1), USE.NAMES = TRUE) -> fprop

  names(fprop) <- names(feature_type)
  fprop
}



