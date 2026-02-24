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
  n_samples       <- length(samp_names)

  # -- Input validation --
  expected_ncol <- n_meta_cols + (n_samples * cols_per_sample)
  if (ncol(dat) != expected_ncol) {
    stop(glue::glue(
      "Column mismatch: expected {expected_ncol} ",
      "({n_meta_cols} meta + {n_samples} samples x {cols_per_sample}), ",
      "but got {ncol(dat)}."
    ))
  }

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


