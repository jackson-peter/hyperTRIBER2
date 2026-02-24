# tests/testthat/test-extract_count_data.R

library(testthat)

# Fix path: when running via test_file(), working dir is project root
# but when running via test(), it's tests/testthat/
# Use here::here() or testthat::test_path()
mpileup_files <- list.files(test_path("testdata"), pattern = "txt$", full.names = TRUE)
cat("Files found:", length(mpileup_files), "\n")

mpileup_df <- do.call(rbind, lapply(mpileup_files, read.table, header = FALSE))
print(mpileup_df)
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
