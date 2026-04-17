# CPM Pipeline Integration & Completion Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bridge the gaps in the refactored connectivity pipeline to ensure full functionality, including robust cross-validation, canonical network plotting, and end-to-end orchestration.

**Architecture:** Extend `R/cpm_engine.R` with k-fold CV, `R/plotting_utils.R` with atlas-rearrangement, and `main.R` with full orchestration logic.

**Tech Stack:** R, `tidyverse`, `future`, `ggplot2`.

---

### Task 1: Complete CPM Engine (K-Fold & Prediction)

**Files:**
- Modify: `R/cpm_engine.R`
- Test: `tests/test_cpm_engine.R`

- [ ] **Step 1: Add k-fold cross-validation logic**

```r
# R/cpm_engine.R
# (Append to existing file)

kfold_cpm <- function(x_mats, y_vec, k, p_thresh = 0.05, mode = "linear", covariates = NULL) {
  num_subs <- ncol(x_mats)
  rand_inds <- sample(seq_len(num_subs))
  sample_size <- floor(num_subs / k)
  
  y_pred <- numeric(num_subs)
  
  for (fold in seq_len(k)) {
    test_inds <- if (fold != k) {
      rand_inds[((fold - 1) * sample_size + 1):(fold * sample_size)]
    } else {
      rand_inds[((fold - 1) * sample_size + 1):num_subs]
    }
    train_inds <- setdiff(rand_inds, test_inds)
    
    # Train
    res <- train_cpm(x_mats[, train_inds], y_vec[train_inds], p_thresh, mode, covariates[, train_inds, drop=F])
    
    # Predict (Simple linear sum model for now)
    pos_sum_test <- colSums(x_mats[res$pos_mask, test_inds, drop=F])
    neg_sum_test <- colSums(x_mats[res$neg_mask, test_inds, drop=F])
    combined_test <- pos_sum_test - neg_sum_test
    
    # Fit on training data
    pos_sum_train <- colSums(x_mats[res$pos_mask, train_inds, drop=F])
    neg_sum_train <- colSums(x_mats[res$neg_mask, train_inds, drop=F])
    combined_train <- pos_sum_train - neg_sum_train
    
    fit <- lm(y_vec[train_inds] ~ combined_train)
    y_pred[test_inds] <- predict(fit, newdata = data.frame(combined_train = combined_test))
  }
  
  return(y_pred)
}
```

- [ ] **Step 2: Update tests for k-fold**

```r
# tests/test_cpm_engine.R
test_that("kfold_cpm returns predictions for all subjects", {
  set.seed(42)
  subjects <- 30
  edges <- 50
  behav <- rnorm(subjects)
  mat <- matrix(rnorm(edges * subjects), nrow = edges, ncol = subjects)
  
  y_pred <- kfold_cpm(mat, behav, k = 5)
  expect_equal(length(y_pred), subjects)
  expect_false(any(is.na(y_pred)))
})
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e 'testthat::test_file("tests/test_cpm_engine.R")'`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add R/cpm_engine.R tests/test_cpm_engine.R
git commit -m "feat: complete k-fold cross-validation in cpm_engine"
```

---

### Task 2: Atlas-Aware Plotting (Canonical Rearrangement)

**Files:**
- Modify: `R/plotting_utils.R`

- [ ] **Step 1: Update plotting function to support node ordering**

```r
# R/plotting_utils.R

plot_connectivity_heatmap <- function(matrix, node_order = NULL, labels = NULL, title = "Consensus Network") {
  if (!is.null(node_order)) {
    matrix <- matrix[node_order, node_order]
  }
  
  df <- as.data.frame.table(matrix)
  colnames(df) <- c("Var1", "Var2", "value")
  
  # Map numeric indices to factor levels for plotting
  df$Var1 <- factor(df$Var1, levels = unique(df$Var1))
  df$Var2 <- factor(df$Var2, levels = rev(unique(df$Var2)))
  
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}
```

- [ ] **Step 2: Commit**

```bash
git add R/plotting_utils.R
git commit -m "feat: add support for canonical node rearrangement in plotting"
```

---

### Task 3: Pipeline Orchestration (main.R Completion)

**Files:**
- Modify: `main.R`

- [ ] **Step 1: Implement full pipeline orchestration**

```r
# main.R
# (Update existing file)

library(yaml)
library(optparse)
source("R/io_utils.R")
source("R/preprocessing.R")
source("R/cpm_engine.R")
source("R/plotting_utils.R")

# ... existing option parsing ...

# 1. Load Data
mats <- load_matrices(config$data$matrices)
meta <- readr::read_csv(config$data$outcomes)

# 2. Align and Clean
aligned <- align_data(mats, meta, config$data$subject_id_var, config$data$outcome_var)
mats_clean <- clean_matrices(aligned$matrices)
# Re-align meta after matrix cleaning
meta_clean <- aligned$metadata %>% filter(subject_id %in% names(mats_clean))

# 3. Flatten for CPM
num_nodes <- nrow(mats_clean[[1]])
upper_tri <- upper.tri(matrix(0, num_nodes, num_nodes))
x_flat <- sapply(mats_clean, function(m) m[upper_tri])
y_vec <- meta_clean[[config$data$outcome_var]]

# 4. Run CPM (Single iteration for main script example)
message("Running CPM...")
y_pred <- kfold_cpm(x_flat, y_vec, k = config$cpm_params$k_folds)

# 5. Output Results
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
readr::write_csv(data.frame(actual = y_vec, predicted = y_pred), 
                 file.path(config$output_dir, "predictions.csv"))

# 6. Plot Consensus (simplified example)
res_full <- train_cpm(x_flat, y_vec, p_thresh = config$cpm_params$p_threshold)
consensus_mat <- matrix(0, num_nodes, num_nodes)
consensus_mat[upper_tri] <- as.numeric(res_full$pos_mask)
consensus_mat <- consensus_mat + t(consensus_mat)

# Load atlas mapping for rearrangement
atlas_map <- read_yaml(config$atlas$mapping_file)
node_order <- if (!is.null(atlas_map$node_order)) atlas_map$node_order else seq_len(num_nodes)

p <- plot_connectivity_heatmap(consensus_mat, node_order = node_order, title = "Consensus Network (Positive)")
ggsave(file.path(config$output_dir, "consensus_heatmap.png"), p)

message("Pipeline completed successfully.")
```

- [ ] **Step 2: Commit**

```bash
git add main.R
git commit -m "feat: complete end-to-end pipeline orchestration in main.R"
```
