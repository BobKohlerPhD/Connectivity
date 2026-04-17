# tests/test_config_loading.R
library(testthat)
library(yaml)

root_dir <- "/home/bob/projects/github-public/Connectivity"

test_that("default configuration can be loaded", {
  config_path <- file.path(root_dir, "config/default_config.yaml")
  
  # Ensure the config file exists
  expect_true(file.exists(config_path))
  
  # Load the config
  config <- read_yaml(config_path)
  
  # Check some key fields
  expect_equal(config$project_name, "Connectivity_Study_2026")
  expect_equal(config$atlas$name, "Shen")
  expect_equal(config$atlas$mapping_file, "atlases/Shen/shen_mapping.yaml")
})

test_that("Shen atlas mapping can be loaded", {
  mapping_path <- file.path(root_dir, "atlases/Shen/shen_mapping.yaml")
  
  # Ensure the mapping file exists
  expect_true(file.exists(mapping_path))
  
  # Load the mapping
  mapping <- read_yaml(mapping_path)
  
  # Check some key fields
  expect_equal(mapping$name, "Shen")
  expect_equal(mapping$num_nodes, 268)
})
