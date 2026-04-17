#~~~Requires thresholded masks~~~#
library(tidyverse)
library(R.matlab)
library(stringr)
library(reshape2)
library(scales)
library(ggplot2)


# Paths to CPM output folder and shen atlas config 
base_path   <- '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/mid'
output_path <- file.path(base_path, '2024-10-07_10fold_p_thresh_0.05_repeat100_iter1000_timepoint_1_z0_mode_linear_mat_mid_t1_subject.txt_3815_t1_MID_matrices')
shen_path   <- file.path('/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/shen268')

# thresholded files 
pos_thresh_file <- file.path(output_path, "net_pos_thresh_1.0.txt")
neg_thresh_file <- file.path(output_path, "net_neg_thresh_1.0.txt")

# Shen Atlas config files
lst_nodes_txt <- readLines(file.path(shen_path, "lst_nodes_orig.txt"))
num_nodes     <- length(lst_nodes_txt)

# node -> network mapping
nets  <- readMat(file.path(shen_path, "nets.mat")) %>% unlist()
nodes <- readMat(file.path(shen_path, "nodes.mat")) %>% unlist() %>% as.numeric()

#===================== HELPERS =====================#
read_binary_mat <- function(path) {
  # Handles comma OR whitespace; trims blank lines and trailing empties
  lines <- readr::read_lines(path)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("Empty file: ", path)
  delim <- if (grepl(",", lines[1])) "," else "\\s+"
  split_rows <- strsplit(lines, delim, perl = TRUE)
  split_rows <- lapply(split_rows, function(v) v[nzchar(trimws(v))])
  row_len <- vapply(split_rows, length, integer(1))
  modal_len <- as.integer(names(sort(table(row_len), decreasing = TRUE))[1])
  split_rows <- split_rows[row_len == modal_len]
  m <- suppressWarnings(matrix(as.numeric(unlist(split_rows)),
                               nrow = length(split_rows), byrow = TRUE))
  # Drop trailing all-NA row/col if present
  if (ncol(m) != nrow(m)) {
    if (ncol(m) > nrow(m) && all(is.na(m[, ncol(m)]))) m <- m[, -ncol(m), drop = FALSE]
    if (nrow(m) > ncol(m) && all(is.na(m[nrow(m), ]))) m <- m[-nrow(m), , drop = FALSE]
  }
  if (ncol(m) != nrow(m)) stop(sprintf("Non-square matrix %dx%d in %s", nrow(m), ncol(m), basename(path)))
  # Assume already 0/1; if not, coerce non-zero to 1 without changing 0/1 files
  uniqv <- unique(na.omit(as.vector(m)))
  if (!all(uniqv %in% c(0,1))) m <- (m != 0) * 1
  m
}

#===================== READ THRESHOLDED MATRICES =====================#
data_thresh_pos <- read_binary_mat(pos_thresh_file)
data_thresh_neg <- read_binary_mat(neg_thresh_file)

stopifnot(nrow(data_thresh_pos) == num_nodes,
          ncol(data_thresh_pos) == num_nodes,
          nrow(data_thresh_neg) == num_nodes,
          ncol(data_thresh_neg) == num_nodes)

#===================== WHOLE-MATRIX HEATMAPS (NO NETWORK ORDER) =====================#
pos_threshold_df <- as.data.frame(as.table(data_thresh_pos))
colnames(pos_threshold_df) <- c("Row", "Column", "Value")
pos_threshold_df <- pos_threshold_df %>%
  mutate(Row = as.numeric(Row), Column = as.numeric(Column))

positive_edge_nonetwork <- ggplot(pos_threshold_df, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile(lwd = 1.5) +
  scale_x_continuous(breaks = seq(0, num_nodes, by = 50)) +
  scale_y_continuous(breaks = seq(0, num_nodes, by = 50)) +
  coord_cartesian() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.75),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_blank())
ggsave(file.path(output_path, "positive_edge_nonetwork.jpeg"),
       positive_edge_nonetwork, dpi = 300, width = 8, height = 6, units = "in")

neg_threshold_df <- as.data.frame(as.table(data_thresh_neg))
colnames(neg_threshold_df) <- c("Row", "Column", "Value")
neg_threshold_df <- neg_threshold_df %>%
  mutate(Row = as.numeric(Row), Column = as.numeric(Column))

negative_edge_nonetwork <- ggplot(neg_threshold_df, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile(lwd = 1.5) +
  scale_x_continuous(breaks = seq(0, num_nodes, by = 50)) +
  scale_y_continuous(breaks = seq(0, num_nodes, by = 50)) +
  coord_cartesian() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.75),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_blank())
ggsave(file.path(output_path, "negative_edge_nonetwork.jpeg"),
       negative_edge_nonetwork, dpi = 300, width = 8, height = 6, units = "in")

#===================== CANONICAL REARRANGEMENT =====================#
results_dir <- file.path(output_path, "canonical_networks")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Rearrange the *thresholded* matrices by canonical node order
correlation_rearr_pos <- data_thresh_pos[nodes, nodes]
correlation_rearr_neg <- data_thresh_neg[nodes, nodes]

write.table(correlation_rearr_pos,
            file.path(results_dir, "correlation_canonical_rearr_pos_thresh.txt"),
            row.names = FALSE, col.names = FALSE)
write.table(correlation_rearr_neg,
            file.path(results_dir, "correlation_canonical_rearr_neg_thresh.txt"),
            row.names = FALSE, col.names = FALSE)

#===================== MODULE COMPOSITION / SIZES =====================#
num_module   <- max(nets)
module_sizes <- sapply(1:num_module, function(i) sum(nets == i))  # numeric
modules      <- lapply(1:num_module, function(i) lst_nodes_txt[nets == i])

# (Optional) save module compositions as in your original
for (i in 1:num_module) {
  out_file <- file.path(results_dir, sprintf("module%d_comp_realvalued.txt", i))
  write.table(as.character(modules[[i]]), file = out_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

module_files <- list.files(path = results_dir, pattern = "module.*_comp_realvalued\\.txt", full.names = TRUE)
module_names <- module_files[!str_detect(module_files, "modules_comp")]
module_names <- module_names[order(as.numeric(str_extract(basename(module_names), "\\d+")))]
cat("Number of modules:", length(module_names), "\n")
module_size <- sapply(module_names, function(fn) length(readr::read_lines(fn)))
cat("Module sizes:", module_size, "\n")

# Network heatmap setup 
network_labels <- c("Medial Frontal", "Frontoparietal", "Default", "Motor/Sensory", "Visual",
                    "Visual B", "Visual Association", "Salience", "Subcortical", "Brainstem & Cerebellum")
network_labels_reversed <- rev(network_labels)

# Use rearranged, thresholded matrices 
data_thresh_pos_canon <- as.data.frame(correlation_rearr_pos)
data_thresh_neg_canon <- as.data.frame(correlation_rearr_neg)

# Ensure numeric
data_thresh_pos_canon[] <- lapply(data_thresh_pos_canon, as.numeric)
data_thresh_neg_canon[] <- lapply(data_thresh_neg_canon, as.numeric)

# Create count tables
df_pos       <- matrix(NA_real_, nrow = length(network_labels), ncol = length(network_labels),
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()
df_ratio_pos <- matrix(NA_real_, nrow = length(network_labels), ncol = length(network_labels),
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()

df_neg       <- matrix(NA_real_, nrow = length(network_labels), ncol = length(network_labels),
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()
df_ratio_neg <- matrix(NA_real_, nrow = length(network_labels), ncol = length(network_labels),
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()

# Network counts 
module_size_num <- as.numeric(module_size)

for (i in seq_len(length(network_labels))) {
  s1 <- module_size_num[i]
  l1 <- if (i > 1) sum(module_size_num[1:(i-1)]) else 0
  h1 <- l1 + s1
  net1 <- network_labels[i]
  
  for (j in seq_len(length(network_labels))) {
    s2 <- module_size_num[j]
    l2 <- if (j > 1) sum(module_size_num[1:(j-1)]) else 0
    h2 <- l2 + s2
    net2 <- network_labels[j]
    
    # Positive network counts
    block_pos <- data_thresh_pos_canon[(l1+1):h1, (l2+1):h2]
    num_edges_pos <- sum(block_pos, na.rm = TRUE)
    if (i == j) {
      num_edges_pos <- num_edges_pos / 2
      network_size <- s1 * (s1 - 1) / 2
    } else {
      network_size <- s1 * s2
    }
    df_pos[net1, net2]       <- num_edges_pos
    df_ratio_pos[net1, net2] <- ifelse(network_size > 0, num_edges_pos / network_size, NA_real_)
    
    # Negative network counts
    block_neg <- data_thresh_neg_canon[(l1+1):h1, (l2+1):h2]
    num_edges_neg <- sum(block_neg, na.rm = TRUE)
    if (i == j) {
      num_edges_neg <- num_edges_neg / 2
      network_size_neg <- s1 * (s1 - 1) / 2
    } else {
      network_size_neg <- s1 * s2
    }
    df_neg[net1, net2]       <- num_edges_neg
    df_ratio_neg[net1, net2] <- ifelse(network_size_neg > 0, num_edges_neg / network_size_neg, NA_real_)
  }
}

# Mask LOWER triangle to DISPLAY UPPER (same as your original)
df_pos  <- df_pos  %>% mutate(across(everything(), as.numeric))
df_neg  <- df_neg  %>% mutate(across(everything(), as.numeric))
mask_pos <- lower.tri(as.matrix(df_pos), diag = FALSE)
mask_neg <- lower.tri(as.matrix(df_neg), diag = FALSE)
df_pos[mask_pos] <- NA
df_neg[mask_neg] <- NA


x_labels <- network_labels
y_labels <- rev(network_labels) 

# Heat maps 
create_heatmap_plot_positive <- function(data, x, y, fill, x_labels, y_labels) {
  plot_data <- data %>% filter(!is.na(.data[[fill]]))
  min_fill  <- min(plot_data[[fill]], na.rm = TRUE)
  max_fill  <- max(plot_data[[fill]], na.rm = TRUE)
  threshold <- (min_fill + max_fill) / 2
  
  ggplot(plot_data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(color = "black", size = 0.5) +
    geom_text(aes_string(label = fill,
                         colour = sprintf("%s > %s", fill, threshold)),
              size = 6.5) +
    scale_fill_gradientn(
      colours  = c("ivory", "#FFE4E1", "lightcoral", "salmon",
                   "indianred2", "brown3", "firebrick4"),
      limits   = c(min_fill, max_fill),
      na.value = "white"
    ) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_x_discrete(labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(labels = y_labels, expand = c(0,0)) +
    coord_fixed(expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", size = 0), 
      axis.line.x       = element_line(color = "black", size = 0.5),
      axis.line.y       = element_line(color = "black", size = 0.5),
      axis.line.x.top   = element_line(color = "black", size = 0.5),
      axis.line.y.right = element_line(color = "black", size = 0.5),
      axis.text.x       = element_text(color = "black", angle = 45, hjust = 1,size = 16),
      axis.text.x.top   = element_blank(),
      axis.text.y       = element_text(color = "black", size = 16),
      axis.text.y.right = element_blank(),
      legend.position   = "none"
    )
}

create_heatmap_plot_negative <- function(data, x, y, fill, x_labels, y_labels) {
  plot_data <- data %>% filter(!is.na(.data[[fill]]))
  min_fill  <- min(plot_data[[fill]], na.rm = TRUE)
  max_fill  <- max(plot_data[[fill]], na.rm = TRUE)
  threshold <- (min_fill + max_fill) / 2

  ggplot(plot_data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(color = "black", size = 0.5) +
    geom_text(aes_string(label = fill,
                         colour = sprintf("%s > %s", fill, threshold)),
              size = 6.5) +
    scale_fill_gradientn(
      colors   = c("ivory", "lightblue", "skyblue2",
                   "dodgerblue", "dodgerblue2", "dodgerblue3", "royalblue3"),
      limits   = c(min_fill, max_fill),
      na.value = "white"
    ) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_x_discrete(labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(labels = y_labels, expand = c(0,0)) +
    coord_fixed(expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", size = 0),
      axis.line.x       = element_line(color = "black", size = 0.5),
      axis.line.y       = element_line(color = "black", size = 0.5),
      axis.line.x.top   = element_line(color = "black", size = 0.5),
      axis.line.y.right = element_line(color = "black", size = 0.5),
      axis.text.x       = element_text(color = "black", angle = 45, hjust = 1,  size = 16),
      axis.text.x.top   = element_blank(),
      axis.text.y       = element_text(color = "black", size = 16),
      axis.text.y.right = element_blank(),
      legend.position   = "none")
}
# Melt, reverse y labels, plot, and save
df_pos_melt <- melt(as.matrix(df_pos), varnames = c("x", "y"), value.name = "value")
df_pos_melt$y <- factor(df_pos_melt$y, levels = rev(colnames(df_pos)))
positive_edge_plot <- create_heatmap_plot_positive(df_pos_melt, "x", "y", "value", x_labels, y_labels)
ggsave(file.path(output_path, "positive_edge_plot.jpeg"),
       positive_edge_plot, dpi = 300, width = 8, height = 6, units = "in")

df_neg_melt <- melt(as.matrix(df_neg), varnames = c("x", "y"), value.name = "value")
df_neg_melt$y <- factor(df_neg_melt$y, levels = rev(colnames(df_neg)))
negative_edge_plot <- create_heatmap_plot_negative(df_neg_melt, "x", "y", "value", x_labels, y_labels)
ggsave(file.path(output_path, "negative_edge_plot.jpeg"),
       negative_edge_plot, dpi = 300, width = 8, height = 6, units = "in")
