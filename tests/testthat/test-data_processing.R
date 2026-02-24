library(testthat)

# Fix path: when running via test_file(), working dir is project root
# but when running via test(), it's tests/testthat/
# Use here::here() or testthat::test_path()
mpileup_files <- list.files(test_path("testdata"), pattern = "txt$", full.names = TRUE)
cat("Files found:", length(mpileup_files), "\n")

mpileup_df <- do.call(rbind, lapply(mpileup_files, read.table, header = FALSE))

cat("Dimensions:", dim(mpileup_df), "\n")

samp_names <- c("R22", "R23", "R24", "R25", "R26", "R27",
                "R28", "R29", "R30", "R31", "R32", "R33")

# ---- Test 1: Correct structure ----
test_that("extract_count_data returns correct structure with real data", {
  result <- extract_count_data(mpileup_df, samp_names, stranded = FALSE)

  # Result should be a names list by sample
  expect_type(result, "list")
  expect_length(result, 12) # number of samples
  expect_equal(nrow(result[["R22"]]), 10) # number of sites
  expect_named(result, samp_names)

  expect_true(all(c("A","T","C","G") %in% colnames(result[["R22"]])))
})

# ---- Test 2: Spot-check known values ----
test_that("extract_count_data correctly maps counts to samples", {
  result <- extract_count_data(mpileup_df, samp_names, stranded = FALSE)

  expect_equal(result[["R22"]]$A[1], 242)
  expect_equal(result[["R22"]]$T[1], 0)
  expect_equal(result[["R22"]]$C[1], 0)
  expect_equal(result[["R22"]]$G[1], 0)
  expect_equal(result[["R22"]]$C[3], 444)

  expect_equal(result[["R23"]]$A[1], 200)
  expect_equal(result[["R23"]]$A[5], 500)
  expect_equal(result[["R23"]]$T[5], 1)

  expect_equal(result[["R33"]]$A[1], 793)
  expect_equal(result[["R33"]]$G[1], 1)
  expect_equal(result[["R33"]]$G[8], 78)
})

# ---- Test 3: site_id is correct ----
test_that("site_id is correctly formed", {
  result <- extract_count_data(mpileup_df, samp_names, stranded = FALSE)

  expect_equal(result[["R22"]]$site_id[1], "ADARclone_9")
  expect_equal(result[["R22"]]$site_id[6], "Chr1_3879")
  expect_equal(result[["R22"]]$site_id[8], "Chr1_4061")
  expect_identical(result[["R22"]]$site_id, result[["R33"]]$site_id)
})
# ---- Test 4: Column count validation ----
test_that("extract_count_data rejects wrong number of columns", {
  bad_df <- mpileup_df[, 1:20]
  expect_error(extract_count_data(bad_df, samp_names, stranded = FALSE))
})

test_that("restrict_data keeps only sites with sufficient editing", {
  dl <- make_test_data()
  design <- c("treat", "treat", "treat", "control", "control")

  result <- restrict_data(
    data_list       = dl,
    design_vector   = design,
    min_samp_treat  = 3,
    min_count       = 2,
    min_prop        = 0.01
  )

  # Should keep sites 1 (A→G) and 3 (T→C)
  expect_equal(nrow(result[[1]]), 2)
  expect_equal(result[[1]]$pos, c(1L, 3L))
  # All samples should have same sites
  expect_equal(result[[1]]$site_id, result[[5]]$site_id)
})

test_that("restrict_data drops everything when no editing exists", {
  dl <- make_test_data()
  # Wipe all G and C counts to zero: no A>G or T>C editing
  dl <- purrr::map(dl, ~ dplyr::mutate(.x, G = 0L, C = 0L))
  design <- c("treat", "treat", "treat", "control", "control")

  result <- restrict_data(dl, design_vector = design, min_samp_treat = 3, min_count = 2)

  expect_equal(nrow(result[[1]]), 0)
})

test_that("min_samp_treat threshold works", {
  dl <- make_test_data()
  design <- c("treat", "treat", "treat", "control", "control")

  # Site 5 has editing in only 1 replicate (G=3,0,0)
  # With min_samp_treat=1, site 5 should now pass too
  result <- restrict_data(
    dl,
    design_vector  = design,
    min_samp_treat = 1,
    min_count      = 2
  )

  expect_true(5L %in% result[[1]]$pos)  # site 5 now included

  # With min_samp_treat=3, site 5 should be dropped
  result_strict <- restrict_data(
    dl,
    design_vector  = design,
    min_samp_treat = 3,
    min_count      = 2
  )

  expect_false(5L %in% result_strict[[1]]$pos)
})

test_that("ref == target sites are always excluded", {
  dl <- make_test_data()
  design <- c("treat", "treat", "treat", "control", "control")

  # Site 6: ref=G, which is the target for A→G → should never pass
  result <- restrict_data(
    dl,
    design_vector  = design,
    min_samp_treat = 1,
    min_count      = 1
  )

  expect_false(6L %in% result[[1]]$pos)
})

test_that("both_ways keeps sites with editing in control only", {
  dl <- make_test_data()
  # Kill all treat editing for site 1, but keep control editing
  dl[["T1"]]$G[1] <- 0L
  dl[["T2"]]$G[1] <- 0L
  dl[["T3"]]$G[1] <- 0L
  # Control still has G=2,1 at site 1

  design <- c("treat", "treat", "treat", "control", "control")

  # Without both_ways: site 1 should drop (no treat editing)
  result_oneway <- restrict_data(
    dl, design_vector = design,
    min_samp_treat = 2, min_count = 2, both_ways = FALSE
  )
  expect_false(1L %in% result_oneway[[1]]$pos)

  # With both_ways: site 1 should be kept IF control passes
  # Control has G=2,1 → 2 reps > 0, total=3, prop>0
  result_bothways <- restrict_data(
    dl, design_vector = design,
    min_samp_treat = 2, min_count = 2, min_prop = 0.01, both_ways = TRUE
  )
  expect_true(1L %in% result_bothways[[1]]$pos)
})

test_that("output structure is preserved after filtering", {
  dl <- make_test_data()
  design <- c("treat", "treat", "treat", "control", "control")

  result <- restrict_data(dl, design_vector = design)

  expect_type(result, "list")
  expect_equal(length(result), 5)
  expect_named(result, c("T1", "T2", "T3", "C1", "C2"))
  expect_true(all(c("chr", "pos", "ref", "A", "T", "C", "G", "site_id") %in% names(result[[1]])))
})


# tests/testthat/test-get_edits_vs_reference.R

library(testthat)

# -- Helper: build minimal inputs -------------------------------------------
make_fake_inputs <- function(stranded = FALSE) {
  # 4 sites, 3 samples
  if (stranded) {
    site_names <- c("chr1_100,+", "chr1_200,+", "chr1_300,-", "chr1_400,+")
  } else {
    site_names <- c("chr1_100", "chr1_200", "chr1_300", "chr1_400")
  }

  # Each sample has counts for A, T, C, G
  make_sample <- function(a, t, c, g) {
    df <- data.frame(A = a, T = t, C = c, G = g)
    row.names(df) <- site_names
    df
  }

  # Site 1: ref=A, good A->G edit (high A, high G, low T/C)
  # Site 2: ref=A, no target (high A, no G)
  # Site 3: ref=T, wrong ref for A->G
  # Site 4: ref=A, too much noise (high C)
  data_list <- list(
    s1 = make_sample(a = c(20, 25, 5, 18),
                     t = c(1,  0,  30, 0),
                     c = c(0,  0,  0,  15),
                     g = c(10, 0,  1,  8)),
    s2 = make_sample(a = c(15, 20, 3, 22),
                     t = c(0,  1,  25, 1),
                     c = c(1,  0,  1,  12),
                     g = c(8,  0,  0,  5)),
    s3 = make_sample(a = c(18, 22, 4, 20),
                     t = c(0,  0,  28, 0),
                     c = c(0,  1,  0,  10),
                     g = c(12, 0,  2,  7))
  )

  # Reference GRanges
  ref_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(100, 200, 300, 400), width = 1),
    ref      = c("A", "A", "T", "A")
  )

  design_vector <- c(s1 = "treat", s2 = "treat", s3 = "treat")

  edits_of_interest <- matrix(c("A", "G"), ncol = 2, byrow = TRUE)

  list(
    data_list         = data_list,
    ref_gr            = ref_gr,
    design_vector     = design_vector,
    edits_of_interest = edits_of_interest,
    site_names        = site_names
  )
}

# -- parse_site_positions ----------------------------------------------------
test_that("parse_site_positions extracts chr/pos unstranded", {
  dnames <- c("chr1_100", "chr2_999")
  res <- parse_site_positions(dnames, stranded = FALSE)

  expect_equal(res$chr, c("chr1", "chr2"))
  expect_equal(res$pos, c(100, 999))
  expect_equal(row.names(res), dnames)
  expect_false("strand" %in% colnames(res))
})

test_that("parse_site_positions extracts chr/pos/strand stranded", {
  dnames <- c("chr1_100,+", "chr2_999,-")
  res <- parse_site_positions(dnames, stranded = TRUE)

  expect_equal(res$chr, c("chr1", "chr2"))
  expect_equal(res$pos, c(100, 999))
  expect_equal(res$strand, c("+", "-"))
})

# -- lookup_reference_base ---------------------------------------------------
test_that("lookup_reference_base returns correct bases unstranded", {
  inp <- make_fake_inputs(stranded = FALSE)
  pos_dat <- parse_site_positions(inp$site_names, stranded = FALSE)

  ref <- lookup_reference_base(pos_dat, inp$ref_gr, stranded = FALSE)

  expect_equal(ref, c("A", "A", "T", "A"))
})

test_that("lookup_reference_base complements minus strand", {
  ref_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 100, width = 1),
    ref      = "A"
  )
  pos_dat <- data.frame(chr = "chr1", pos = 100, strand = "-",
                        row.names = "chr1_100,-",
                        stringsAsFactors = FALSE)

  ref <- lookup_reference_base(pos_dat, ref_gr, stranded = TRUE)
  expect_equal(ref, "T")  # complement of A
})

test_that("lookup_reference_base uppercases lowercase refs", {
  ref_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 100, width = 1),
    ref      = "a"  # lowercase
  )
  pos_dat <- data.frame(chr = "chr1", pos = 100,
                        row.names = "chr1_100",
                        stringsAsFactors = FALSE)

  ref <- lookup_reference_base(pos_dat, ref_gr, stranded = FALSE)
  expect_equal(ref, "A")
})

# -- get_edits_vs_reference --------------------------------------------------
test_that("get_edits_vs_reference finds expected A->G site", {
  inp <- make_fake_inputs()

  res <- get_edits_vs_reference(
    data_list         = inp$data_list,
    ref_gr            = inp$ref_gr,
    stranded          = FALSE,
    edits_of_interest = inp$edits_of_interest,
    sample_type       = "treat",
    design_vector     = inp$design_vector,
    min_ref  = 3,
    min_targ = 3,
    max_other = 2
  )

  expect_true("A:G" %in% names(res))
  # Site 1 should pass: ref=A, all 3 samples have A>0 and G>0, T/C low
  expect_true("chr1_100" %in% res[["A:G"]])
})

test_that("get_edits_vs_reference excludes site with no target", {
  inp <- make_fake_inputs()

  res <- get_edits_vs_reference(
    data_list         = inp$data_list,
    ref_gr            = inp$ref_gr,
    stranded          = FALSE,
    edits_of_interest = inp$edits_of_interest,
    sample_type       = "treat",
    design_vector     = inp$design_vector,
    min_ref  = 3,
    min_targ = 3,
    max_other = 2
  )

  # Site 2: ref=A but G=0 in all samples
  expect_false("chr1_200" %in% res[["A:G"]])
})

test_that("get_edits_vs_reference excludes noisy site", {
  inp <- make_fake_inputs()

  res <- get_edits_vs_reference(
    data_list         = inp$data_list,
    ref_gr            = inp$ref_gr,
    stranded          = FALSE,
    edits_of_interest = inp$edits_of_interest,
    sample_type       = "treat",
    design_vector     = inp$design_vector,
    min_ref  = 3,
    min_targ = 3,
    max_other = 2
  )

  # Site 4: ref=A, has G, but C is high in all 3 samples -> too much noise
  expect_false("chr1_400" %in% res[["A:G"]])
})

test_that("get_edits_vs_reference skips wrong-ref sites", {
  inp <- make_fake_inputs()

  res <- get_edits_vs_reference(
    data_list         = inp$data_list,
    ref_gr            = inp$ref_gr,
    stranded          = FALSE,
    edits_of_interest = inp$edits_of_interest,
    sample_type       = "treat",
    design_vector     = inp$design_vector
  )

  # Site 3: ref=T, never considered for A->G
  expect_false("chr1_300" %in% res[["A:G"]])
})

test_that("get_edits_vs_reference respects min_targ threshold", {
  inp <- make_fake_inputs()

  # Only 1 sample has G>0 at site 2? Actually 0.
  # Better: raise min_targ so even site 1 fails
  res <- get_edits_vs_reference(
    data_list         = inp$data_list,
    ref_gr            = inp$ref_gr,
    stranded          = FALSE,
    edits_of_interest = inp$edits_of_interest,
    sample_type       = "treat",
    design_vector     = inp$design_vector,
    min_ref  = 3,
    min_targ = 100,  # impossible threshold
    max_other = 2
  )

  expect_length(res[["A:G"]], 0)
})


# tests/testthat/test-editing_proportions_vs_reference.R

library(testthat)

# -- Helper ------------------------------------------------------------------
make_prop_inputs <- function(stranded = FALSE) {
  if (stranded) {
    site_names <- c("chr1_100,+", "chr1_200,-")
  } else {
    site_names <- c("chr1_100", "chr1_200")
  }

  make_sample <- function(a, t, c, g) {
    df <- data.frame(A = a, T = t, C = c, G = g)
    row.names(df) <- site_names
    df
  }

  # Site 1: 80A + 20G = 100 total  → edit.A=0.8, edit.G=0.2
  # Site 2: 50T + 50C = 100 total  → edit.T=0.5, edit.C=0.5
  data_list <- list(
    s1 = make_sample(a = c(40, 0), t = c(0, 25), c = c(0, 25), g = c(10, 0)),
    s2 = make_sample(a = c(40, 0), t = c(0, 25), c = c(0, 25), g = c(10, 0))
  )

  ref_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(100, 200), width = 1),
    ref      = c("A", "T")
  )

  list(data_list = data_list, ref_gr = ref_gr, site_names = site_names)
}

# -- Tests -------------------------------------------------------------------
test_that("editing_proportions returns correct proportions", {
  inp <- make_prop_inputs()

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = FALSE
  )

  expect_equal(res$total, c(100, 100))
  expect_equal(res$edit.A, c(0.8, 0.0))
  expect_equal(res$edit.G, c(0.2, 0.0))
  expect_equal(res$edit.T, c(0.0, 0.5))
  expect_equal(res$edit.C, c(0.0, 0.5))
})

test_that("editing_proportions assigns correct ref base", {
  inp <- make_prop_inputs()

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = FALSE
  )

  expect_equal(res$ref, c("A", "T"))
})

test_that("editing_proportions handles zero-coverage sites", {
  inp <- make_prop_inputs()
  # Zero out all counts for site 2
  inp$data_list <- lapply(inp$data_list, function(x) {
    x[2, ] <- 0L
    x
  })

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = FALSE
  )

  expect_equal(res$total[2], 0)
  expect_true(is.nan(res$edit.A[2]))  # 0/0
})

test_that("editing_proportions rows sum to 1 (non-zero sites)", {
  inp <- make_prop_inputs()

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = FALSE
  )

  row_sums <- res$edit.A + res$edit.T + res$edit.C + res$edit.G
  expect_equal(row_sums, c(1, 1))
})

test_that("editing_proportions works stranded", {
  inp <- make_prop_inputs(stranded = TRUE)

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = TRUE
  )

  expect_true("strand" %in% colnames(res))
  expect_equal(nrow(res), 2)
  # Proportions still valid
  expect_equal(res$total, c(100, 100))
})

test_that("editing_proportions returns expected columns", {
  inp <- make_prop_inputs()

  res <- editing_proportions_vs_reference(
    data_list = inp$data_list,
    ref_gr    = inp$ref_gr,
    stranded  = FALSE
  )

  expected_cols <- c("chr", "pos", "ref", "total",
                     "edit.A", "edit.T", "edit.C", "edit.G")
  expect_true(all(expected_cols %in% colnames(res)))
})


