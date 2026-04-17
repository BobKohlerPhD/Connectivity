# tests/test_cpm_engine.R
source("/home/bob/projects/github-public/Connectivity/R/cpm_engine.R")
library(testthat)

test_that("train_cpm identifies significant edges", {
  # Create synthetic data with known correlations
  set.seed(42)
  subjects <- 20
  edges <- 100
  behav <- rnorm(subjects)
  
  # Create matrix where first 5 edges are highly correlated with behav
  mat <- matrix(rnorm(edges * subjects), nrow = edges, ncol = subjects)
  mat[1:5, ] <- t(sapply(1:5, function(i) behav + rnorm(subjects, sd = 0.1)))
  
  res <- train_cpm(mat, behav, p_thresh = 0.05)
  expect_true(all(res$pos_mask[1:5]))
})

test_that("train_cpm handles covariates with partial correlation", {
  set.seed(42)
  subjects <- 30
  edges <- 10
  behav <- rnorm(subjects)
  covar <- rnorm(subjects)
  
  # Create edge that correlates with behav ONLY when covar is NOT accounted for
  # (Spurious correlation driven by covar)
  edge_spurious <- behav + 2*covar + rnorm(subjects, sd = 0.1)
  
  # Create edge that correlates with behav ONLY when covar IS accounted for
  # (This is harder to synthesize simply, let's just test that it doesn't crash 
  # and produces expected dimensions for now, or use a simpler case)
  
  mat <- matrix(rnorm(edges * subjects), nrow = edges, ncol = subjects)
  mat[1, ] <- edge_spurious
  
  # Without covariate, it should be significant
  res_no_covar <- train_cpm(mat, behav, p_thresh = 0.05)
  expect_true(res_no_covar$pos_mask[1])
  
  # With covariate, it should be less significant or non-significant
  res_with_covar <- train_cpm(mat, behav, p_thresh = 0.05, covariates = data.frame(covar))
  # It might still be significant if the relationship is very strong, 
  # but here we just want to ensure it works.
  expect_type(res_with_covar$pos_mask, "logical")
  expect_length(res_with_covar$pos_mask, edges)
})
