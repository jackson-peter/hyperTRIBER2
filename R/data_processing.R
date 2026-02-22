extractCountData <- function(dat, samp.names, stranded = FALSE) {
  library(tidyverse)

  # Calculate column indices for each sample
  base_cols <- if (stranded) 8 else 4
  col_indices <- map(1:length(samp.names), ~ {
    start <- 3 + ((.x - 1) * base_cols + 1)
    end <- start + base_cols - 1
    seq(start, end)
  })

  # Process each sample
  data.list <- imap(col_indices, ~ {
    samp_name <- samp.names[.y]

    # Extract data for current sample
    data_samp <- dat %>%
      select(all_of(1:2), all_of(.x)) %>%
      column_to_rownames("V1") %>%
      as.matrix() %>%
      setNames(if (stranded) {
        c("A", "T", "C", "G", "a", "t", "c", "g")
      } else {
        c("A", "T", "C", "G")
      })

    # Process stranded data
    if (stranded) {
      data_samp <- data_samp %>%
        as.data.frame() %>%
        mutate(
          strand = rep(c("+", "-"), each = nrow(dat)),
          position = rep(dat$V2, times = 2)
        ) %>%
        pivot_longer(
          cols = -c(strand, position),
          names_to = "base",
          values_to = "count"
        ) %>%
        pivot_wider(
          id_cols = c(position, strand),
          names_from = base,
          values_from = count
        ) %>%
        mutate(
          # Swap bases for negative strand
          `A-` = t,
          `T-` = a,
          `C-` = g,
          `G-` = c,
          across(c(A, T, C, G), ~ ifelse(strand == "-", NA, .x))
        ) %>%
        select(position, strand, A, T, C, G) %>%
        unite("id", position, strand, sep = ",") %>%
        column_to_rownames("id") %>%
        select(A, T, C, G) %>%
        as.matrix()
    }

    # Track NA rows to remove later
    na_rows <- which(is.na(data_samp), arr.ind = TRUE)[, 1]

    list(data = data_samp, na_rows = na_rows)
  })

  # Find all rows with NAs in any sample
  all_na_rows <- unique(unlist(map(data.list, ~ .x$na_rows)))

  # Remove NA rows from all samples
  data.list <- map(data.list, ~ {
    .x$data %>%
      `row.names<-`(.) %>%
      .[ -all_na_rows, ]
  }) %>%
    setNames(samp.names)

  return(data.list)
}


restrict_data <- function(data_list, ref_base = NULL, design_vector,
                          min_samp_treat = 3, min_count = 2, min_prop = 0.01,
                          both_ways = FALSE, edits_of_interest = rbind(c("A", "G"), c("T", "C"))) {

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (length(design_vector) != length(data_list)) stop("design_vector length must match data_list")
  if (!all(design_vector %in% c("control", "treat"))) stop("design_vector must contain only 'control' or 'treat'")

  # Split data into control and treatment groups
  split_data <- split(data_list, design_vector)
  cont_data <- split_data$control
  treat_data <- split_data$treat

  combine_counts <- function(df_list, base) {
    # Extract the specified base column from each dataframe
    base_columns <- lapply(df_list, function(df) df[, base])
    # Combine columns horizontally
    cbind_list <- do.call(cbind, base_columns)
    return(cbind_list)
  }

  cont_counts <- lapply(c("A", "T", "C", "G"), combine_counts, df_list = cont_data)
  treat_counts <- lapply(c("A", "T", "C", "G"), combine_counts, df_list = treat_data)
  names(cont_counts) <- names(treat_counts) <- c("A", "T", "C", "G")

  # Determine reference base (either provided or dominant from controls)
  dom_base <- if (is.null(ref_base)) {
    dom_df <- as.data.frame(lapply(cont_counts, rowSums))
    c("A","T","C","G")[max.col(dom_df, ties.method = "first")]
  } else {
    ref_base
  }

  # Initialize results list
  tokeep_list <- vector("list", nrow(edits_of_interest))

  for (i in seq_len(nrow(edits_of_interest))) {
    my_ref <- edits_of_interest[i, 1]
    my_targ <- edits_of_interest[i, 2]

    # Calculate proportions with protection against division by zero
    num_treat <- rowSums(treat_counts[[my_targ]])
    denom_treat <- num_treat + rowSums(treat_counts[[my_ref]])
    prop_treat <- ifelse(denom_treat == 0, 0, num_treat/denom_treat)

    num_control <- rowSums(cont_counts[[my_targ]])
    denom_control <- num_control + rowSums(cont_counts[[my_ref]])
    prop_control <- ifelse(denom_control == 0, 0, num_control/denom_control)

    # Count calculations
    count_treat <- rowSums(treat_counts[[my_targ]])
    count_control <- rowSums(cont_counts[[my_targ]])

    tokeep_treat <- rowSums(treat_counts[[my_targ]] > 0) >= min_samp_treat
    tokeep_control <- rowSums(cont_counts[[my_targ]] > 0) >= min_samp_treat

    # Filtering logic
    treat_ok <- tokeep_treat & (count_treat >= min_count) & (prop_treat > min_prop)
    control_ok <- tokeep_control & (count_control >= min_count) & (prop_control > min_prop)

    # Apply filtering based on both_ways parameter
    if (both_ways) {
      tokeep_list[[i]] <- (dom_base != my_targ) & (treat_ok | control_ok)
    } else {
      tokeep_list[[i]] <- (dom_base != my_targ) & treat_ok
    }
  }

  # Combine filtering results and apply to data
  use_table <- do.call(cbind, tokeep_list)
  to_keep <- rowSums(use_table) > 0

  return(lapply(data_list, function(x) x[to_keep, ]))
}
