# tests/test_plotting_utils.R
source("../R/plotting_utils.R")
library(testthat)

test_that("plot_connectivity_heatmap returns a ggplot object", {
  mat <- matrix(rnorm(100), 10, 10)
  p <- plot_connectivity_heatmap(mat, labels = 1:10)
  expect_s3_class(p, "ggplot")
})

test_that("plot_connectivity_heatmap uses labels correctly", {
  mat <- matrix(rnorm(4), 2, 2)
  labels <- c("A", "B")
  p <- plot_connectivity_heatmap(mat, labels = labels)
  
  # Check if labels are in the plot data
  expect_true(all(labels %in% p$data$Var1))
  expect_true(all(labels %in% p$data$Var2))
})
