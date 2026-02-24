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
