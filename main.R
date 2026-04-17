# main.R
#' @title Connectivity Pipeline Orchestration
#' @description This script orchestrates the Connectome-based Predictive Modeling (CPM) pipeline.
#' It handles data loading, preprocessing, model training via cross-validation, 
#' results persistence, and visualization.

library(yaml)
library(optparse)
library(dplyr)
library(readr)
library(ggplot2)

# Source modular components
source("R/io_utils.R")
source("R/preprocessing.R")
source("R/cpm_engine.R")
source("R/plotting_utils.R")

# 1. Option Parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/default_config.yaml",
              help = "Path to configuration file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. Load Configuration
if (!file.exists(opt$config)) {
  stop(sprintf("Configuration file not found: %s", opt$config))
}
config <- read_yaml(opt$config)
message(sprintf("Starting pipeline for project: %s", config$project_name))

# 3. Load Data
message("Loading connectivity matrices and metadata...")
if (!file.exists(config$data$matrices)) {
    stop(sprintf("Matrix file not found: %s", config$data$matrices))
}
if (!file.exists(config$data$outcomes)) {
    stop(sprintf("Outcomes file not found: %s", config$data$outcomes))
}

mats <- load_matrices(config$data$matrices)
meta <- readr::read_csv(config$data$outcomes)

# 4. Align and Clean
message("Aligning data and cleaning matrices...")
aligned <- align_data(mats, meta, config$data$subject_id_var, config$data$outcome_var)
mats_clean <- clean_matrices(aligned$matrices)

# Ensure meta is also cleaned and aligned with final set of matrices
final_subject_ids <- names(mats_clean)
meta_clean <- aligned$metadata %>% 
  filter(!!sym(config$data$subject_id_var) %in% final_subject_ids) %>%
  arrange(match(!!sym(config$data$subject_id_var), final_subject_ids))

if (nrow(meta_clean) == 0) {
  stop("No valid subjects remaining after alignment and cleaning.")
}

# 5. Flatten for CPM
message("Flattening matrices for CPM...")
num_nodes <- nrow(mats_clean[[1]])
upper_tri <- upper.tri(matrix(0, num_nodes, num_nodes))
x_flat <- sapply(mats_clean, function(m) m[upper_tri])
y_vec <- meta_clean[[config$data$outcome_var]]

# 6. Run CPM (K-fold Cross-Validation)
message(sprintf("Running CPM with %d-fold cross-validation...", config$cpm_params$k_folds))
y_pred <- kfold_cpm(
  x_mats = x_flat, 
  y_vec = y_vec, 
  k = config$cpm_params$k_folds,
  p_thresh = config$cpm_params$p_threshold,
  mode = config$cpm_params$mode
)

# 7. Output Results
message(sprintf("Saving results to %s...", config$output_dir))
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

results_df <- data.frame(
  subject_id = names(mats_clean),
  actual = y_vec, 
  predicted = y_pred
)
readr::write_csv(results_df, file.path(config$output_dir, "predictions.csv"))

# 8. Plot Consensus Network
message("Generating consensus network plot...")
# Train on full dataset to identify consensus edges
res_full <- train_cpm(x_flat, y_vec, p_thresh = config$cpm_params$p_threshold)

# Create a square matrix for plotting
consensus_mat <- matrix(0, num_nodes, num_nodes)
consensus_mat[upper_tri] <- as.numeric(res_full$pos_mask)
# Symmetrize
consensus_mat <- consensus_mat + t(consensus_mat)

# Load atlas mapping for rearrangement
atlas_map <- read_yaml(config$atlas$mapping_file)
node_order <- if (!is.null(atlas_map$node_order)) atlas_map$node_order else seq_len(num_nodes)

p <- plot_connectivity_heatmap(
  matrix = consensus_mat, 
  node_order = node_order, 
  title = "Consensus Network (Positive Edges)"
)
ggsave(file.path(config$output_dir, "consensus_heatmap.png"), p)

message("Pipeline completed successfully.")
