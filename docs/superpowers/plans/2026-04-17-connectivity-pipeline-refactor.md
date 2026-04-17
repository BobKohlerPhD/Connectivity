# Connectivity Pipeline Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor fragmented R scripts into a modular, configuration-driven, and atlas-agnostic CPM pipeline.

**Architecture:** Centralized `main.R` orchestrates modular R components (`io_utils`, `preprocessing`, `cpm_engine`, `posthoc_analysis`, `plotting_utils`) using a `config.yaml` for all parameters.

**Tech Stack:** R, `tidyverse`, `R.matlab`, `yaml`, `future`, `caret`, `ggplot2`.

---

### Task 1: Project Structure & Configuration

**Files:**
- Create: `config/default_config.yaml`
- Create: `atlases/Shen/shen_mapping.yaml`
- Create: `main.R`

- [ ] **Step 1: Create the default configuration file**

```yaml
# config/default_config.yaml
project_name: "Connectivity_Study_2026"
output_dir: "./results"

data:
  matrices: ""
  outcomes: ""
  subject_id_var: "subject_id"
  outcome_var: "outcome"
  covariates: []

cpm_params:
  k_folds: 10
  iterations: 100
  p_threshold: 0.05
  mode: "linear"
  z_score: true

atlas:
  name: "Shen"
  mapping_file: "atlases/Shen/shen_mapping.yaml"
```

- [ ] **Step 2: Create the Shen atlas mapping file**

```yaml
# atlases/Shen/shen_mapping.yaml
name: "Shen"
num_nodes: 268
network_labels:
  - "Medial Frontal"
  - "Frontoparietal"
  - "Default"
  - "Motor/Sensory"
  - "Visual"
  - "Visual B"
  - "Visual Association"
  - "Salience"
  - "Subcortical"
  - "Brainstem & Cerebellum"
node_list: "Connectomics/AtlasConfigs/ShenAtlasConfig/NodeList.txt"
canonical_nodes: "Connectomics/AtlasConfigs/ShenAtlasConfig/canonical_nodes.mat"
network_config: "Connectomics/AtlasConfigs/ShenAtlasConfig/node_to_network_config.mat"
```

- [ ] **Step 3: Create the main entry point skeleton**

```r
# main.R
library(yaml)
library(optparse)

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/default_config.yaml",
              help = "Path to configuration file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

config <- read_yaml(opt$config)
message(sprintf("Starting pipeline for project: %s", config$project_name))
```

- [ ] **Step 4: Commit**

```bash
git add config/default_config.yaml atlases/Shen/shen_mapping.yaml main.R
git commit -m "chore: initialize project structure and config"
```

---

### Task 2: I/O Utilities Module

**Files:**
- Create: `R/io_utils.R`
- Create: `tests/test_io_utils.R`

- [ ] **Step 1: Write failing tests for I/O utilities**

```r
# tests/test_io_utils.R
source("R/io_utils.R")
library(testthat)

test_that("load_matrices works for .rds", {
  tmp_file <- tempfile(fileext = ".rds")
  test_data <- list(sub1 = matrix(1, 10, 10), sub2 = matrix(2, 10, 10))
  saveRDS(test_data, tmp_file)
  
  loaded <- load_matrices(tmp_file)
  expect_equal(length(loaded), 2)
  expect_true(is.matrix(loaded[[1]]))
})
```

- [ ] **Step 2: Implement I/O utilities**

```r
# R/io_utils.R
library(R.matlab)
library(readr)

load_matrices <- function(path) {
  ext <- tools::file_ext(path)
  if (ext == "rds") {
    return(readRDS(path))
  } else if (ext == "mat") {
    data <- readMat(path)
    # Handle both list of matrices and 3D array
    if ("connectivity" %in% names(data)) {
      mats <- data$connectivity
      ids <- as.vector(data$subject_id)
      res <- lapply(seq_len(dim(mats)[3]), function(i) mats[,,i])
      names(res) <- ids
      return(res)
    }
  }
  stop("Unsupported matrix format")
}

validate_matrices <- function(mats) {
  if (length(mats) == 0) stop("No matrices found")
  node_count <- nrow(mats[[1]])
  for (i in seq_along(mats)) {
    if (!isSymmetric(mats[[i]], tol = 1e-8)) stop(sprintf("Matrix %d is not symmetric", i))
    if (nrow(mats[[i]]) != node_count) stop("Inconsistent node counts")
  }
  return(TRUE)
}
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e 'testthat::test_file("tests/test_io_utils.R")'`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add R/io_utils.R tests/test_io_utils.R
git commit -m "feat: implement I/O utilities with validation"
```

---

### Task 3: Preprocessing Module

**Files:**
- Create: `R/preprocessing.R`
- Create: `tests/test_preprocessing.R`

- [ ] **Step 1: Write failing tests for preprocessing**

```r
# tests/test_preprocessing.R
source("R/preprocessing.R")
library(testthat)

test_that("align_data correctly intersects subjects", {
  mats <- list(sub1 = matrix(0, 5, 5), sub2 = matrix(0, 5, 5))
  meta <- data.frame(subject_id = c("sub2", "sub3"), outcome = c(10, 20))
  
  aligned <- align_data(mats, meta, "subject_id", "outcome")
  expect_equal(names(aligned$matrices), "sub2")
  expect_equal(nrow(aligned$metadata), 1)
})
```

- [ ] **Step 2: Implement preprocessing**

```r
# R/preprocessing.R
library(dplyr)

align_data <- function(mats, meta, id_var, outcome_var, covariates = NULL) {
  meta_ids <- meta[[id_var]]
  mat_ids <- names(mats)
  
  common_ids <- intersect(meta_ids, mat_ids)
  
  filtered_mats <- mats[common_ids]
  filtered_meta <- meta %>%
    filter(!!sym(id_var) %in% common_ids) %>%
    filter(!is.na(!!sym(outcome_var)))
  
  # Final alignment after NA removal
  final_ids <- filtered_meta[[id_var]]
  final_mats <- filtered_mats[final_ids]
  
  list(matrices = final_mats, metadata = filtered_meta)
}

clean_matrices <- function(mats) {
  keep <- vapply(mats, function(m) !any(is.na(m)) && !any(is.infinite(m)), logical(1))
  mats[keep]
}
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e 'testthat::test_file("tests/test_preprocessing.R")'`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add R/preprocessing.R tests/test_preprocessing.R
git commit -m "feat: implement data alignment and cleaning"
```

---

### Task 4: CPM Engine Module

**Files:**
- Create: `R/cpm_engine.R`

- [ ] **Step 1: Implement core CPM functions**

```r
# R/cpm_engine.R
library(ppcor)
library(caret)

train_cpm <- function(train_mat, train_behav, p_thresh = 0.05, mode = "linear", covariates = NULL) {
  # train_mat is [edges x subjects]
  num_edges <- nrow(train_mat)
  
  if (!is.null(covariates)) {
    corrs <- apply(train_mat, 1, function(edge) pcor.test(edge, train_behav, covariates)$p.value)
    coeffs <- apply(train_mat, 1, function(edge) pcor.test(edge, train_behav, covariates)$estimate)
  } else {
    corrs <- apply(train_mat, 1, function(edge) cor.test(edge, train_behav)$p.value)
    coeffs <- apply(train_mat, 1, function(edge) cor.test(edge, train_behav)$estimate)
  }
  
  pos_mask <- (coeffs > 0) & (corrs < p_thresh)
  neg_mask <- (coeffs < 0) & (corrs < p_thresh)
  
  # Return masks and trained models (simplified for now)
  list(pos_mask = pos_mask, neg_mask = neg_mask)
}

# Add kfold_cpm with future support later
```

- [ ] **Step 2: Commit**

```bash
git add R/cpm_engine.R
git commit -m "feat: implement core CPM training logic"
```

---

### Task 5: Plotting Utilities (Atlas-Agnostic)

**Files:**
- Create: `R/plotting_utils.R`

- [ ] **Step 1: Implement heatmap plotting**

```r
# R/plotting_utils.R
library(ggplot2)
library(reshape2)

plot_connectivity_heatmap <- function(matrix, labels, title = "Consensus Network") {
  df <- melt(matrix)
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
```

- [ ] **Step 2: Commit**

```bash
git add R/plotting_utils.R
git commit -m "feat: implement atlas-agnostic plotting"
```
