#' Quick diagnostic plot/print of a single editing site
#'
#' @param data_list Named list of count matrices
#' @param design_vector Character vector of sample conditions
#' @param site_name Site ID to inspect (ignored if random = TRUE)
#' @param random Logical; pick a random site
#' @param ref_gr Optional GRanges with reference info
#' @param show_ref Logical; print reference info
plot_site_example <- function(data_list,
                              design_vector,
                              site_name = NULL,
                              random = TRUE,
                              ref_gr = NULL,
                              show_ref = TRUE) {

  if (random) {
    site_name <- sample(row.names(data_list[[1]]), 1)
  }
  stopifnot(!is.null(site_name))

  counts <- do.call(rbind, lapply(data_list, function(x) x[site_name, ]))
  cat("Raw counts:\n")
  print(counts)

  prop_by_group <- function(group) {
    idx <- which(design_vector == group)
    props <- t(apply(counts[idx, , drop = FALSE], 1, function(x) x / sum(x)))
    colMeans(props)
  }

  ctrl_prop  <- prop_by_group("control")
  treat_prop <- prop_by_group("treat")

  cat("\nNormalized proportions (control / treat):\n")
  print(round(rbind(control = ctrl_prop / sum(ctrl_prop),
                    treat   = treat_prop / sum(treat_prop)), 4))

  if (show_ref && !is.null(ref_gr)) {
    names(ref_gr) <- ref_gr$names
    cat("\nReference info:\n")
    print(ref_gr[site_name])
  }
}
