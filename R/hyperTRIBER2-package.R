#' @keywords internal
"_PACKAGE"

#' @importFrom BiocGenerics strand strand<- start end
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom foreach %dopar%
NULL

# Column names used in dplyr / tidyr pipelines
utils::globalVariables(c(
  "pos", "chr", "ref", "site_id", "base", "row_id", "count", "strand",
  "gene_id", "v9", "experiment", "condition", "sample", "i"
))
