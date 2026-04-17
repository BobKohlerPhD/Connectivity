# tests/test_io_utils.R
library(testthat)

# Use absolute path for robustness or correctly resolve relative to test file
# The standard in R packages is to use test_check(), but for this script-based setup
# we assume source() relative to the file.
source("../R/io_utils.R")

test_that("load_matrices works for .rds", {
  tmp_file <- tempfile(fileext = ".rds")
  test_data <- list(sub1 = matrix(1, 10, 10), sub2 = matrix(2, 10, 10))
  # Ensure symmetry for later validation tests
  test_data <- lapply(test_data, function(m) { m[lower.tri(m)] <- t(m)[lower.tri(m)]; m })
  saveRDS(test_data, tmp_file)
  
  loaded <- load_matrices(tmp_file)
  expect_equal(length(loaded), 2)
  expect_equal(names(loaded), c("sub1", "sub2"))
  expect_true(is.matrix(loaded[[1]]))
  unlink(tmp_file)
})

test_that("validate_matrices works for valid matrices", {
  mats <- list(
    sub1 = matrix(0, 5, 5),
    sub2 = matrix(0, 5, 5)
  )
  expect_true(validate_matrices(mats))
})

test_that("validate_matrices fails for non-symmetric matrices", {
  mats <- list(
    sub1 = matrix(runif(25), 5, 5)
  )
  expect_error(validate_matrices(mats), "not symmetric")
})

test_that("validate_matrices fails for inconsistent dimensions", {
  mats <- list(
    sub1 = matrix(0, 5, 5),
    sub2 = matrix(0, 6, 6)
  )
  expect_error(validate_matrices(mats), "Inconsistent node counts")
})
