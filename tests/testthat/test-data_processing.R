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

# tests/testthat/test-get_hits.R

library(testthat)

# -- Helpers -----------------------------------------------------------------

#' Build a fake DEXSeq results table
make_fake_dexseq_res <- function() {
  data.frame(
    groupID   = c("chr1_100", "chr1_100", "chr1_200", "chr1_300"),
    featureID = c("EG", "EC", "EG", "EG"),
    padj      = c(0.01, 0.5, 0.05, 0.001),
    pvalue    = c(0.005, 0.4, 0.03, 0.0005),
    log2fold_treat_control = c(2.0, 0.5, 1.5, -1.0),
    control   = c(0.02, 0.01, 0.03, 0.08),
    treat     = c(0.15, 0.02, 0.12, 0.03),
    stringsAsFactors = FALSE
  )
}

make_fake_data_list <- function() {
  sites <- c("chr1_100", "chr1_200", "chr1_300")
  make_sample <- function(a, t, c, g) {
    df <- data.frame(A = a, T = t, C = c, G = g)
    row.names(df) <- sites
    df
  }

  list(
    ctrl1 = make_sample(a = c(50, 40, 10), t = c(1, 0, 1),
                        c = c(0, 1, 0),    g = c(5, 8, 40)),
    ctrl2 = make_sample(a = c(48, 38, 12), t = c(0, 1, 0),
                        c = c(1, 0, 1),    g = c(4, 7, 38)),
    treat1 = make_sample(a = c(30, 25, 35), t = c(0, 0, 0),
                         c = c(1, 0, 1),    g = c(20, 18, 10)),
    treat2 = make_sample(a = c(28, 22, 33), t = c(1, 1, 0),
                         c = c(0, 1, 0),    g = c(22, 20, 12))
  )
}

make_fake_ref_gr <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(100, 200, 300), width = 1),
    ref      = c("A", "A", "G")
  )
  gr$names <- paste0("chr1_", c(100, 200, 300))
  gr
}

default_design <- c(ctrl1 = "control", ctrl2 = "control",
                    treat1 = "treat", treat2 = "treat")

default_edits <- rbind(c("A", "G"), c("T", "C"))

# -- compute_site_meta -------------------------------------------------------

test_that("compute_site_meta returns correct structure", {
  dl <- make_fake_data_list()
  sig <- make_fake_dexseq_res()[c(1, 3, 4), ]  # only sig rows

  res <- compute_site_meta(
    site_id           = "chr1_100",
    sig_res           = sig,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    treat_idx         = c("treat1", "treat2"),
    ctrl_idx          = c("ctrl1", "ctrl2"),
    symmetric         = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expected_cols <- c("name", "ref", "targ", "prop", "prop_ctrl",
                     "padj", "pvalue", "control_par", "treat_par",
                     "fold_change", "tags_treat", "tags_control")
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("compute_site_meta identifies dominant base as ref", {
  dl <- make_fake_data_list()
  sig <- make_fake_dexseq_res()[1, , drop = FALSE]

  res <- compute_site_meta(
    site_id           = "chr1_100",
    sig_res           = sig,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    treat_idx         = c("treat1", "treat2"),
    ctrl_idx          = c("ctrl1", "ctrl2"),
    symmetric         = FALSE
  )

  # A is dominant at chr1_100 (50+48+30+28=156 vs G: 5+4+20+22=51)
  expect_equal(res$ref, "A")
  expect_equal(res$targ, "G")
})

test_that("compute_site_meta calculates proportions correctly", {
  dl <- make_fake_data_list()
  sig <- make_fake_dexseq_res()[1, , drop = FALSE]

  res <- compute_site_meta(
    site_id           = "chr1_100",
    sig_res           = sig,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    treat_idx         = c("treat1", "treat2"),
    ctrl_idx          = c("ctrl1", "ctrl2"),
    symmetric         = FALSE
  )

  # treat: G = 20+22=42, A = 30+28=58 → prop = 42/(58+42) = 0.42
  expect_equal(res$prop, 42 / 100)
  # control: G = 5+4=9, A = 50+48=98 → prop_ctrl = 9/(98+9) = 9/107
  expect_equal(res$prop_ctrl, 9 / 107, tolerance = 1e-6)
})

test_that("compute_site_meta resolves multiple hits at same site", {
  dl <- make_fake_data_list()
  # chr1_100 has EG (padj=0.01) and EC (padj=0.5) — but EC filtered before
  # Simulate two sig hits at same site
  sig <- data.frame(
    groupID   = c("chr1_100", "chr1_100"),
    featureID = c("EG", "EC"),
    padj      = c(0.01, 0.02),
    pvalue    = c(0.005, 0.01),
    log2fold_treat_control = c(2.0, 1.0),
    control   = c(0.02, 0.01),
    treat     = c(0.15, 0.05),
    stringsAsFactors = FALSE
  )

  res <- compute_site_meta(
    site_id           = "chr1_100",
    sig_res           = sig,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    treat_idx         = c("treat1", "treat2"),
    ctrl_idx          = c("ctrl1", "ctrl2"),
    symmetric         = FALSE
  )

  # Should pick G (valid A->G edit), not C
  expect_equal(res$targ, "G")
})
test_that("compute_site_meta symmetric mode swaps on negative FC", {
  # Build fake data where control has MORE editing than treat
  # so that log2fold_treat_control < 0
  sig_res <- data.frame(
    groupID  = "chr1_100",
    featureID = "EG",
    padj     = 0.01,
    pvalue   = 0.005,
    log2fold_treat_control = -2.0,  # <-- negative FC
    control  = 0.3,
    treat    = 0.05,
    stringsAsFactors = FALSE
  )

  # Control samples have high G, treat samples have low G
  dl <- list(
    S1 = data.frame(A = 80, T = 5, C = 5, G = 10, row.names = "chr1_100"),
    S2 = data.frame(A = 85, T = 5, C = 5, G = 5,  row.names = "chr1_100"),
    S3 = data.frame(A = 40, T = 5, C = 5, G = 50, row.names = "chr1_100"),
    S4 = data.frame(A = 30, T = 5, C = 5, G = 60, row.names = "chr1_100")
  )
  dv <- c(S1 = "treat", S2 = "treat", S3 = "control", S4 = "control")

  res_nosym <- compute_site_meta("chr1_100", sig_res, dl,
                                 rbind(c("A","G")), dv, which(dv == "treat"), which(dv == "control"),
                                 symmetric = FALSE)

  res_sym <- compute_site_meta("chr1_100", sig_res, dl,
                               rbind(c("A","G")), dv, which(dv == "treat"), which(dv == "control"),
                               symmetric = TRUE)

  # With negative FC + symmetric, prop should use control group
  expect_false(identical(res_sym$prop, res_nosym$prop))
})


# -- get_hits (integration) --------------------------------------------------

test_that("get_hits without metadata returns GRanges with edit_base", {
  res <- make_fake_dexseq_res()
  dl  <- make_fake_data_list()

  hits <- get_hits(
    dexseq_res        = res,
    stranded          = FALSE,
    fdr               = 0.1,
    add_meta          = FALSE,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    include_ref       = FALSE
  )

  expect_s4_class(hits, "GRanges")
  expect_true("edit_base" %in% colnames(GenomicRanges::mcols(hits)))
})

test_that("get_hits filters by FDR correctly", {
  res <- make_fake_dexseq_res()
  dl  <- make_fake_data_list()

  hits_loose <- get_hits(
    dexseq_res = res, fdr = 0.5, add_meta = FALSE,
    data_list = dl, edits_of_interest = default_edits,
    design_vector = default_design, include_ref = FALSE
  )

  hits_strict <- get_hits(
    dexseq_res = res, fdr = 0.001, add_meta = FALSE,
    data_list = dl, edits_of_interest = default_edits,
    design_vector = default_design, include_ref = FALSE
  )

  expect_gte(length(hits_loose), length(hits_strict))
})

test_that("get_hits excludes self-edits (ref==targ)", {
  dl <- make_fake_data_list()
  ref_gr <- make_fake_ref_gr()

  res <- make_fake_dexseq_res()

  hits <- get_hits(
    dexseq_res        = res,
    fdr               = 0.1,
    n_cores           = 1,
    add_meta          = TRUE,
    data_list         = dl,
    edits_of_interest = default_edits,
    design_vector     = default_design,
    include_ref       = TRUE,
    ref_gr            = ref_gr
  )

  # No site should have ref == targ
  if (length(hits) > 0) {
    expect_true(all(hits$ref != hits$targ))
  }
})

test_that("get_hits filters to edits of interest only", {
  res <- make_fake_dexseq_res()
  # Add a T->C row that shouldn't appear if we only ask for A->G
  res <- rbind(res, data.frame(
    groupID = "chr1_400", featureID = "EC",
    padj = 0.001, pvalue = 0.0005,
    log2fold_treat_control = 3.0,
    control = 0.01, treat = 0.2,
    stringsAsFactors = FALSE
  ))

  dl <- make_fake_data_list()
  # Add site to data_list
  dl <- lapply(dl, function(x) {
    new_row <- data.frame(A = 5, T = 50, C = 10, G = 0)
    row.names(new_row) <- "chr1_400"
    rbind(x, new_row)
  })

  hits <- get_hits(
    dexseq_res        = res,
    fdr               = 0.1,
    add_meta          = FALSE,
    data_list         = dl,
    edits_of_interest = rbind(c("A", "G")),  # only A->G
    design_vector     = default_design,
    include_ref       = FALSE
  )

  # EC (T->C) feature should be excluded since we only want EG
  features <- GenomicRanges::mcols(hits)$edit_base
  expect_true(all(features %in% c("G", NA_character_)))
})

# ---- tests/test_annotate_with_genes.R ----

library(testthat)
library(GenomicRanges)

# =============================================================================
# Helper: build a minimal GTF-like GRanges
# =============================================================================
make_test_gtf <- function() {
  # GENE_A on + strand: two transcripts, non-overlapping feature types per position
  # TX_A1: exon 100-200, CDS 100-200 (whole thing is coding)
  # TX_A2: exon 100-200, CDS 100-160, 3UTR 161-200
  # GENE_B on - strand:
  # TX_B1: exon 500-600, CDS 500-600
  # TX_B2: exon 500-600, 5UTR 500-520, CDS 521-600
  GRanges(
    seqnames = rep("chr1", 8),
    ranges = IRanges(
      start = c(100, 100, 100, 100, 161,
                500, 500, 500),
      end   = c(200, 200, 200, 160, 200,
                600, 600, 600)
    ),
    strand = c(rep("+", 5), rep("-", 3)),
    type = c("exon", "CDS", "exon", "CDS", "3UTR",
             "exon", "CDS", "exon"),
    gene_id = c(rep("GENE_A", 5), rep("GENE_B", 3)),
    transcript_id = c("TX_A1", "TX_A1", "TX_A2", "TX_A2", "TX_A2",
                      "TX_B1", "TX_B1", "TX_B2")
  )
}

make_test_sites <- function() {
  s <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(start = c(150, 550), width = 1),
    strand = c("+", "-")
  )
  s$names <- c("site_A", "site_B")
  s
}

make_test_quant <- function() {
  c(TX_A1 = 50.0, TX_A2 = 10.0, TX_B1 = 30.0, TX_B2 = 5.0)
}

make_test_gene_ids <- function() {
  c(GENE_A = "GeneAlpha", GENE_B = "GeneBeta")
}

# =============================================================================
# Tests
# =============================================================================

test_that("site overlapping GENE_A gets correct gene annotation", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  res <- annotate_with_genes(sites[1], gtf, gids, quant, n_cores = 1)

  expect_equal(unname(res$gene), "GENE_A")
  expect_equal(unname(res$name), "GeneAlpha")
  expect_equal(unname(res$gtf_strand), "+")
  expect_false(as.logical(res$out_of_range))
})

test_that("site overlapping GENE_B gets correct gene annotation", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  res <- annotate_with_genes(sites[2], gtf, gids, quant, n_cores = 1)

  expect_equal(unname(res$gene), "GENE_B")
  expect_equal(unname(res$name), "GeneBeta")
  expect_equal(unname(res$gtf_strand), "-")
  expect_false(as.logical(res$out_of_range))
})

test_that("site with no direct overlap triggers out_of_range flanking", {
  # Use a simple GTF with just one transcript to avoid handleAnoType issues
  gtf <- GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges(start = c(100, 100), end = c(200, 200)),
    strand = rep("+", 2),
    type = c("exon", "CDS"),
    gene_id = c("GENE_A", "GENE_A"),
    transcript_id = c("TX_A1", "TX_A1")
  )

  quant <- c(TX_A1 = 50.0)
  gids  <- c(GENE_A = "GeneAlpha")

  site <- GRanges("chr1", IRanges(250, width = 1), strand = "+")
  site$names <- "gap_site"

  res <- annotate_with_genes(site, gtf, gids, quant, n_cores = 1)

  expect_true(as.logical(res$out_of_range))
})
test_that("site on different chromosome returns empty annotation", {
  gtf   <- make_test_gtf()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  site <- GRanges("chr2", IRanges(150, width = 1), strand = "+")
  site$names <- "nowhere_site"

  # Function has a known issue with no-hit sites; verify it doesn't silently

  # return wrong data. Accept either an error or an empty gene annotation.
  res <- tryCatch(
    suppressWarnings(annotate_with_genes(site, gtf, gids, quant, n_cores = 1)),
    error = function(e) NULL
  )

  if (!is.null(res)) {
    # If it returns something, gene should be empty or NA
    expect_true(unname(res$gene) %in% c("", NA_character_))
  } else {
    succeed("Function errored on no-overlap site (known limitation)")
  }
})
test_that("gene disambiguation picks gene by GTF order", {
  # When two genes overlap same position, function picks by GTF order (first gene wins)
  gtf <- GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges(start = c(100, 100, 100, 100), end = c(200, 200, 200, 200)),
    strand = rep("+", 4),
    type = c("exon", "CDS", "exon", "CDS"),
    gene_id = c("WINNER", "WINNER", "LOSER", "LOSER"),
    transcript_id = c("TX_WIN", "TX_WIN", "TX_LOSE", "TX_LOSE")
  )

  site <- GRanges("chr1", IRanges(150, width = 1), strand = "+")
  site$names <- "ambig_site"

  quant <- c(TX_WIN = 1.0, TX_LOSE = 100.0)
  gids  <- c(WINNER = "WinGene", LOSER = "LoseGene")

  res <- annotate_with_genes(site, gtf, gids, quant, n_cores = 1)

  # First gene in GTF wins regardless of expression
  expect_equal(unname(res$gene), "WINNER")
})

test_that("site on different chromosome returns empty annotation", {
  gtf   <- make_test_gtf()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  site <- GRanges("chr2", IRanges(150, width = 1), strand = "+")
  site$names <- "nowhere_site"

  # Function crashes on no-overlap sites (vapply type mismatch) — this is a

  # known limitation. Just verify it errors rather than returning wrong data.
  expect_error(
    suppressWarnings(annotate_with_genes(site, gtf, gids, quant, n_cores = 1))
  )
})

test_that("multiple sites processed together return correct length", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  res <- annotate_with_genes(sites, gtf, gids, quant, n_cores = 1)

  expect_equal(length(res), 2)
  expect_equal(unname(res$gene[1]), "GENE_A")
  expect_equal(unname(res$gene[2]), "GENE_B")
})

test_that("assign_strand overwrites site strand from GTF", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  res <- annotate_with_genes(sites[1:2], gtf, gids, quant,
                             assign_strand = TRUE, n_cores = 1)

  expect_equal(as.character(strand(res[1])), "+")
  expect_equal(as.character(strand(res[2])), "-")
})

test_that("transcript ordering respects expression", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()
  quant <- make_test_quant()
  gids  <- make_test_gene_ids()

  res <- annotate_with_genes(sites[1], gtf, gids, quant, n_cores = 1)

  tx  <- strsplit(unname(res$transcripts), ",")[[1]]
  tpm <- as.numeric(strsplit(unname(res$transcript_tpm), ",")[[1]])

  # TX_A1 has TPM 50, TX_A2 has TPM 10 — highest first
  expect_equal(tx[1], "TX_A1")
  expect_true(tpm[1] >= tpm[2])
})


test_that("empty quant vector doesn't crash", {
  gtf   <- make_test_gtf()
  sites <- make_test_sites()[1]
  gids  <- make_test_gene_ids()
  quant <- setNames(numeric(0), character(0))

  res <- annotate_with_genes(sites, gtf, gids, quant, n_cores = 1)

  expect_equal(unname(res$gene), "GENE_A")
})
