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

test_that("plot_connectivity_heatmap rearranges nodes with node_order", {
  mat <- matrix(0, 2, 2)
  mat[1, 2] <- 5 # Original: [1,1]=0, [1,2]=5, [2,1]=0, [2,2]=0
  
  # Reorder [2, 1]
  # Rearranged matrix:
  #            [,1] (orig 2)  [,2] (orig 1)
  # [1,] (orig 2)    0              0
  # [2,] (orig 1)    5              0
  p <- plot_connectivity_heatmap(mat, node_order = c(2, 1))
  
  df <- p$data
  # In as.data.frame.table for a matrix:
  # Var1 corresponds to rows, Var2 to columns
  
  # Get names of first and second levels
  level1 <- levels(df$Var1)[1]
  level2 <- levels(df$Var1)[2]
  
  val_at_row2_col1 <- df$value[df$Var1 == level2 & df$Var2 == level1] 
  
  expect_equal(as.numeric(val_at_row2_col1), 5)
})
