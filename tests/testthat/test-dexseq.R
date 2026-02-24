
library(testthat)

# -- Helper: minimal test data matching extract_count_data output ------------
make_test_data_gcf <- function() {
  sites <- data.frame(
    chr     = c("Chr1", "Chr1", "Chr2"),
    pos     = c(100L, 200L, 300L),
    ref     = c("A", "A", "C"),
    A       = c(50L, 30L, 0L),
    T       = c(0L, 1L, 10L),
    C       = c(2L, 0L, 80L),
    G       = c(5L, 3L, 1L),
    site_id = c("Chr1_100", "Chr1_200", "Chr2_300"),
    stringsAsFactors = FALSE
  )

  list(
    T1 = sites,
    T2 = within(sites, { A <- c(48L, 28L, 0L); G <- c(7L, 5L, 2L) }),
    C1 = within(sites, { A <- c(55L, 35L, 0L); G <- c(2L, 1L, 0L) }),
    C2 = within(sites, { A <- c(60L, 40L, 0L); G <- c(1L, 0L, 0L) })
  )
}

# -- Test 1: correct files are written --------------------------------------
test_that("generate_count_files creates expected output files", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  expected_files <- c(
    "countsTreat1.txt",
    "countsTreat2.txt",
    "countsControl1.txt",
    "countsControl2.txt",
    "annotation_file_made_up_for_DEXSeq.gff"
  )

  created <- list.files(tmpdir)
  for (f in expected_files) {
    expect_true(f %in% created, info = paste("Missing file:", f))
  }
})

# -- Test 2: count file structure -------------------------------------------
test_that("count files have correct number of rows (sites × 4 bases)", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  treat1 <- read.table(file.path(tmpdir, "countsTreat1.txt"),
                       header = FALSE, sep = "\t")

  # 3 sites × 4 bases = 12 rows
  expect_equal(nrow(treat1), 12L)
  # Two columns: row_id and count
  expect_equal(ncol(treat1), 2L)
})

# -- Test 3: row names format is site_id:base -------------------------------
test_that("count file row names follow site_id:base format", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  treat1 <- read.table(file.path(tmpdir, "countsTreat1.txt"),
                       header = FALSE, sep = "\t",
                       row.names = 1)

  rn <- row.names(treat1)

  # All row names should contain a colon
  expect_true(all(grepl(":", rn)))

  # Check specific entries
  expect_true("Chr1_100:A" %in% rn)
  expect_true("Chr1_200:G" %in% rn)
  expect_true("Chr2_300:C" %in% rn)
})

# -- Test 4: spot-check count values ----------------------------------------
test_that("count values match input data", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  treat1 <- read.table(file.path(tmpdir, "countsTreat1.txt"),
                       header = FALSE, sep = "\t",
                       row.names = 1)

  expect_equal(treat1["Chr1_100:A", 1], 50L)
  expect_equal(treat1["Chr1_100:G", 1], 5L)
  expect_equal(treat1["Chr2_300:C", 1], 80L)

  ctrl1 <- read.table(file.path(tmpdir, "countsControl1.txt"),
                      header = FALSE, sep = "\t",
                      row.names = 1)

  expect_equal(ctrl1["Chr1_100:A", 1], 55L)
  expect_equal(ctrl1["Chr1_100:G", 1], 2L)
})

# -- Test 5: GFF annotation structure ---------------------------------------
test_that("GFF annotation has correct structure", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  gff <- read.table(
    file.path(tmpdir, "annotation_file_made_up_for_DEXSeq.gff"),
    header = FALSE, sep = "\t", stringsAsFactors = FALSE
  )

  # 3 sites × 4 bases = 12 rows, 9 columns
  expect_equal(nrow(gff), 12L)
  expect_equal(ncol(gff), 9L)

  # col3 is always "exonic_part"
  expect_true(all(gff$V3 == "exonic_part"))

  # col7 should be "*" for unstranded
  expect_true(all(gff$V7 == "*"))
})

# -- Test 6: GFF col9 contains required DEXSeq fields -----------------------
test_that("GFF column 9 contains transcripts, exonic_part_number, gene_id", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  gff <- read.table(
    file.path(tmpdir, "annotation_file_made_up_for_DEXSeq.gff"),
    header = FALSE, sep = "\t", stringsAsFactors = FALSE
  )

  expect_true(all(grepl("transcripts ", gff$V9)))
  expect_true(all(grepl("exonic_part_number ", gff$V9)))
  expect_true(all(grepl("gene_id ", gff$V9)))

  # Spot-check one entry
  row1 <- gff$V9[1]
  expect_true(grepl("gene_id Chr1_100", row1))
})

# -- Test 7: gene_id includes strand when stranded = TRUE -------------------
test_that("stranded mode adds strand to gene_id and col7", {
  dl <- make_test_data_gcf()
  # Add strand column
  dl <- lapply(dl, function(df) {
    df$strand <- "+"
    df
  })

  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir,
                       stranded = TRUE)

  gff <- read.table(
    file.path(tmpdir, "annotation_file_made_up_for_DEXSeq.gff"),
    header = FALSE, sep = "\t", stringsAsFactors = FALSE
  )

  # col7 should be "+" not "*"
  expect_true(all(gff$V7 == "+"))

  # gene_id should contain strand
  expect_true(all(grepl("gene_id Chr[0-9]+_[0-9]+,\\+", gff$V9)))
})

# -- Test 8: treat/control mapping is correct --------------------------------
test_that("samples are mapped to correct treat/control files", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  # T2 has G=7 at site Chr1_100
  treat2 <- read.table(file.path(tmpdir, "countsTreat2.txt"),
                       header = FALSE, sep = "\t", row.names = 1)
  expect_equal(treat2["Chr1_100:G", 1], 7L)

  # C2 has A=60 at site Chr1_100
  ctrl2 <- read.table(file.path(tmpdir, "countsControl2.txt"),
                      header = FALSE, sep = "\t", row.names = 1)
  expect_equal(ctrl2["Chr1_100:A", 1], 60L)
})

# -- Test 9: function returns annotation invisibly ---------------------------
test_that("function returns annotation invisibly", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- withr::local_tempdir()

  result <- generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 12L)
  expect_equal(ncol(result), 9L)
})

# -- Test 10: output dir is created if it doesn't exist ----------------------
test_that("output directory is created if it does not exist", {
  dl     <- make_test_data_gcf()
  design <- c("treat", "treat", "control", "control")
  tmpdir <- file.path(withr::local_tempdir(), "subdir", "deep")

  expect_false(dir.exists(tmpdir))

  generate_count_files(dl, design_vector = design, out_dir = tmpdir)

  expect_true(dir.exists(tmpdir))
  expect_true(file.exists(file.path(tmpdir, "countsTreat1.txt")))
})

# tests/testthat/test-dexseq.R

library(testthat)

# -- Helper: create minimal count files + GFF --------------------------------
setup_fake_dexseq_dir <- function() {
  tmpdir <- withr::local_tempdir(.local_envir = parent.frame())

  # Minimal count data
  counts <- data.frame(site = c("gene1:001", "gene1:002"), count = c(10L, 20L))

  writeLines("countsControl1", file.path(tmpdir, "countsControl1.txt"))
  writeLines("countsTreat1",   file.path(tmpdir, "countsTreat1.txt"))
  write.table(counts, file.path(tmpdir, "countsControl1.txt"), sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(counts, file.path(tmpdir, "countsTreat1.txt"), sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  writeLines("Chr1\tDARESS\texonic_part\t100\t100\t.\t+\t.\tgene_id gene1; exonic_part_number 001",
             file.path(tmpdir, "annotation.gff"))

  tmpdir
}

# -- Test 1: file discovery --------------------------------------------------
test_that("make_test finds count files and GFF", {
  tmpdir <- setup_fake_dexseq_dir()

  count_files <- sort(list.files(tmpdir, pattern = "^counts.*\\.txt$", full.names = TRUE))
  gff_files   <- list.files(tmpdir, pattern = "\\.gff$", full.names = TRUE)

  expect_length(count_files, 2)
  expect_length(gff_files, 1)
})

# -- Test 2: design mismatch is caught --------------------------------------
test_that("make_test errors when design doesn't match count files", {
  skip_if_not_installed("DEXSeq")
  tmpdir <- setup_fake_dexseq_dir()

  # 3 samples but only 2 count files
  expect_error(
    make_test(tmpdir, design_vector = c("control", "treat", "treat")),
    "Count files don't match design"
  )
})

# -- Test 3: missing GFF is caught ------------------------------------------
test_that("make_test errors when no GFF found", {
  skip_if_not_installed("DEXSeq")
  tmpdir <- setup_fake_dexseq_dir()
  file.remove(list.files(tmpdir, pattern = "\\.gff$", full.names = TRUE))

  expect_error(
    make_test(tmpdir, design_vector = c("control", "treat")),
    "No GFF file found"
  )
})

# # -- Test 4: full integration (slow, skip on CI) ----------------------------
# test_that("make_test runs end-to-end on real data", {
#   skip_on_ci()
#   skip_if_not_installed("DEXSeq")
#   skip_if_not_installed("BiocParallel")
#   skip("Slow integration test — run manually")
#
#   # Would need real generate_count_files() output here
# })


