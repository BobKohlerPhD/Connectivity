# tests/test_preprocessing.R
source("../R/preprocessing.R")
library(testthat)

test_that("align_data correctly intersects subjects", {
  mats <- list(sub1 = matrix(0, 5, 5), sub2 = matrix(0, 5, 5))
  meta <- data.frame(subject_id = c("sub2", "sub3"), outcome = c(10, 20))
  
  aligned <- align_data(mats, meta, "subject_id", "outcome")
  expect_equal(names(aligned$matrices), "sub2")
  expect_equal(nrow(aligned$metadata), 1)
})

test_that("clean_matrices removes matrices with NA or Inf", {
  mats <- list(
    valid = matrix(1, 2, 2),
    has_na = matrix(c(1, NA, 3, 4), 2, 2),
    has_inf = matrix(c(1, Inf, 3, 4), 2, 2)
  )
  
  cleaned <- clean_matrices(mats)
  expect_equal(names(cleaned), "valid")
})
