
#~~~~This script performs the following~~~~~~~~~~~#
# 1. Extracts positive and negative networks
# 2. Applies thresholds to extract significant edge 
# 3. Rearranges networks based on Shen atlas
# 4. Generates heatmaps of  significant edges
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(tidyverse)
library(R.matlab)
library(stringr)
library(reshape2)
library(scales)

base_path   <- ""
output_path <- file.path(base_path, "")
shen_path   <- file.path(base_path, "ShenAtlasConfig")



#~~~~Define CPM Predictive Model Parameters~~~~#
n_repeats <- 100   # set based on the number of CPM repeats (usually 100)
n_folds   <- 10    # number of CV folds
lst_of_folds <- 1:n_folds


#~~~~Read in atlas mapping files (Shen atlas is used here)~~~~#
lst_nodes_txt  <- readLines(file.path(shen_path, "NodeList.txt"))
num_nodes <- length(lst_nodes_txt)
cat("Number of nodes:", num_nodes, "\n")


nets <- readMat(file.path(shen_path, "node_to_network_config.mat")) %>% unlist()
nodes <- readMat(file.path(shen_path, "canonical_nodes.mat")) %>% unlist() %>% 
  as.numeric()


#~~~~Define edge threshold~~~~#
thresh    <- 1.00


#~~~~Extract Positive Network~~~~#
network_p <- list()
for (i in 1:n_repeats) {
  for (fold in lst_of_folds) {
    message(sprintf("Reading POSITIVE iter %d, fold %d", i, fold))
    file_path <- file.path(output_path, sprintf("positive_network_from_training_iter%d_fold_%d.txt", i, fold))
    if (file.exists(file_path)) {
      network_p <- append(network_p, list(as.matrix(read.table(file_path))))
    } else {
      warning(sprintf("File not found: %s", file_path))
    }
  }
}
network_p_array <- array(unlist(network_p), dim = c(nrow(network_p[[1]]), ncol(network_p[[1]]), length(network_p)))
network_p_consensus <- apply(network_p_array, c(1, 2), function(x) mean(x > 0))
write.table(network_p_consensus, file.path(output_path, "network_pos_consensus.txt"), 
            row.names = FALSE, col.names = FALSE)

#~~~~Extract Negative Network~~~~#
network_n <- list()
for (i in 1:n_repeats) {
  for (fold in lst_of_folds) {
    message(sprintf("Reading NEGATIVE iter %d, fold %d", i, fold))
    file_path <- file.path(output_path, sprintf("negative_network_from_training_iter%d_fold_%d.txt", i, fold))
    if (file.exists(file_path)) {
      network_n <- append(network_n, list(as.matrix(read.table(file_path))))
    } else {
      warning(sprintf("File not found: %s", file_path))
    }
  }
}
network_n_array <- array(unlist(network_n), dim = c(nrow(network_n[[1]]), ncol(network_n[[1]]), length(network_n)))
network_n_consensus <- apply(network_n_array, c(1, 2), function(x) mean(x > 0))
write.table(network_n_consensus, file.path(output_path, "network_neg_consensus.txt"), 
            row.names = FALSE, col.names = FALSE)


#~~~~Apply threshold to positive network and extract edges~~~~#
data_pos <- as.matrix(read.table(file.path(output_path, "network_pos_consensus.txt")))
data_thresh_pos <- (data_pos >= thresh) * 1.0
write.table(data_thresh_pos, 
            file.path(output_path, sprintf("net_pos_thresh_%.1f.txt", thresh)), 
            row.names = FALSE, col.names = FALSE)

iu_pos <- which(upper.tri(data_thresh_pos), arr.ind = TRUE)
lst_sig_edges_pos <- list()
for (idx in seq_len(nrow(iu_pos))) {
  x_pos <- iu_pos[idx, 1]
  y_pos <- iu_pos[idx, 2]
  if (data_thresh_pos[x_pos, y_pos] == 1) {
    lst_sig_edges_pos <- append(lst_sig_edges_pos, list(c(lst_rois[x_pos], lst_rois[y_pos])))
  }
}
cat("Positive significant edges count:", length(lst_sig_edges_pos), "\n")
write_conn_pos <- file(file.path(output_path, sprintf("sig_edges_pos_thresh%.1f.txt", thresh)), "w")
for (edge in lst_sig_edges_pos) {
  writeLines(paste(edge, collapse = "\t"), write_conn_pos)
}
close(write_conn_pos)

#~~~~Apply threshold to negative network and extract edges~~~~#
data_neg <- as.matrix(read.table(file.path(output_path, "network_neg_consensus.txt")))
data_thresh_neg <- (data_neg >= thresh) * 1
write.table(data_thresh_neg, 
            file.path(output_path, sprintf("net_neg_thresh_%.1f.txt", thresh)), 
            row.names = FALSE, col.names = FALSE)

iu_neg <- which(upper.tri(data_thresh_neg), arr.ind = TRUE)
lst_sig_edges_neg <- list()
for (idx in seq_len(nrow(iu_neg))) {
  x_neg <- iu_neg[idx, 1]
  y_neg <- iu_neg[idx, 2]
  if (data_thresh_neg[x_neg, y_neg] == 1) {
    lst_sig_edges_neg <- append(lst_sig_edges_neg, list(c(lst_rois[x_neg], lst_rois[y_neg])))
  }
}
cat("Negative significant edges count:", length(lst_sig_edges_neg), "\n")
write_conn_neg <- file(file.path(output_path, sprintf("sig_edges_neg_thresh%.1f.txt", thresh)), "w")
for (edge in lst_sig_edges_neg) {
  writeLines(paste(edge, collapse = "\t"), write_conn_neg)
}
close(write_conn_neg)


#~~~~Positive edge heatmap~~~~#
pos_threshold_df <- as.data.frame(as.table(data_thresh_pos))
colnames(pos_threshold_df) <- c("Row", "Column", "Value")
pos_threshold_df <- pos_threshold_df %>% mutate(Sum = as.numeric(Row) + as.numeric(Column))
pos_threshold_df$Row    <- as.numeric(pos_threshold_df$Row)
pos_threshold_df$Column <- as.numeric(pos_threshold_df$Column)
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
output_file <- file.path(output_path, "positive_edge_nonetwork.jpeg")
ggsave(filename = output_file, plot = positive_edge_nonetwork, dpi = 300, width = 8, height = 6, units = "in")

#~~~~Negative edge heatmap~~~~#
neg_threshold_df <- as.data.frame(as.table(data_thresh_neg))
colnames(neg_threshold_df) <- c("Row", "Column", "Value")
neg_threshold_df <- neg_threshold_df %>% mutate(Sum = as.numeric(Row) + as.numeric(Column))
neg_threshold_df$Row    <- as.numeric(neg_threshold_df$Row)
neg_threshold_df$Column <- as.numeric(neg_threshold_df$Column)
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
output_file <- file.path(output_path, "negative_edge_nonetwork.jpeg")
ggsave(filename = output_file, plot = negative_edge_nonetwork, dpi = 300, width = 8, height = 6, units = "in")



#~~~~Canonical Networks Rearrangement~~~~#
results_dir <- file.path(output_path, "canonical_networks")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


#~~~~Rearrange networks~~~~#
correlation_rearr_pos <- network_p_consensus[nodes, nodes]
fn_corr_rearr_mat_pos <- file.path(results_dir, "correlation_canonical_rearr_pos.txt")
write.table(correlation_rearr_pos, fn_corr_rearr_mat_pos, row.names = FALSE, col.names = FALSE)

correlation_rearr_neg <- network_n_consensus[nodes, nodes]
fn_corr_rearr_mat_neg <- file.path(results_dir, "correlation_canonical_rearr_neg.txt")
write.table(correlation_rearr_neg, fn_corr_rearr_mat_neg, row.names = FALSE, col.names = FALSE)

#~~~~Output module compositions per network based on the atlas~~~~#
num_module   <- max(nets)
module_sizes <- sapply(1:num_module, function(i) sum(nets == i)) %>% as.character()
modules      <- list()
for (i in 1:num_module) {
  tmp_module <- lst_nodes_txt[nets == i]
  modules[[i]] <- tmp_module
  out_file <- file.path(results_dir, sprintf("module%d_comp_realvalued.txt", i))
  write.table(as.character(tmp_module), file = out_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#~~~~Module File List~~~#
module_files <- list.files(path = results_dir, pattern = "module.*_comp_realvalued\\.txt", full.names = TRUE)
module_names <- module_files[!str_detect(module_files, "modules_comp")]
module_names <- module_names[order(as.numeric(str_extract(basename(module_names), "\\d+")))]
cat("Number of modules:", length(module_names), "\n")
module_size <- sapply(module_names, function(fn) length(read_lines(fn)))
cat("Module sizes:", module_size, "\n")

#~~~~Define network names (adjust as needed)~~~~#
network_labels <- c("Medial Frontal", "Frontoparietal", "Default", "Motor/Sensory", "Visual",
                    "Visual B", "Visual Association", "Salience", "Subcortical", "Brainstem & Cerebellum")

#~~~~Process Positive Canonical Network~~~~#
data_pos_canon <- read.delim(file.path(results_dir, "correlation_canonical_rearr_pos.txt"),
                             header = FALSE, stringsAsFactors = FALSE, sep = " ")
data_matrix_pos <- data.frame(lapply(data_pos_canon, as.numeric))
data_thresh_pos_canon <- ifelse(data_matrix_pos >= thresh, 1, 0)

df_pos       <- matrix(NA, nrow = length(network_labels), ncol = length(network_labels), 
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()
df_ratio_pos <- matrix(NA, nrow = length(network_labels), ncol = length(network_labels), 
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()

for (i in seq_len(length(network_labels))) {
  s1 <- module_size[i]
  l1 <- if(i > 1) sum(module_size[1:(i-1)]) else 0
  h1 <- l1 + s1
  net1 <- network_labels[i]
  for (j in seq_len(length(network_labels))) {
    s2 <- module_size[j]
    l2 <- if(j > 1) sum(module_size[1:(j-1)]) else 0
    h2 <- l2 + s2
    net2 <- network_labels[j]
    num_edges_pos <- sum(data_thresh_pos_canon[(l1+1):h1, (l2+1):h2])
    if (i == j) {
      num_edges_pos <- num_edges_pos / 2
      network_size <- s1 * (s1 - 1) / 2
    } else {
      network_size <- s1 * s2
    }
    df_pos[net1, net2]       <- num_edges_pos
    df_ratio_pos[net1, net2] <- num_edges_pos / network_size
  }
}
df_pos <- df_pos %>% mutate(across(everything(), as.numeric))
mask_pos <- lower.tri(df_pos, diag = FALSE)
df_pos[mask_pos] <- NA

#~~~~Negative Edge Networks~~~~#
data_neg_canon <- read.delim(file.path(results_dir, "correlation_canonical_rearr_neg.txt"),
                             header = FALSE, stringsAsFactors = FALSE, sep = " ")
data_matrix_neg <- data.frame(lapply(data_neg_canon, as.numeric))

data_thresh_neg_canon <- ifelse(data_matrix_neg == thresh, 1, 0)

df_neg       <- matrix(NA, nrow = length(network_labels), ncol = length(network_labels), 
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()
df_ratio_neg <- matrix(NA, nrow = length(network_labels), ncol = length(network_labels), 
                       dimnames = list(network_labels, network_labels)) %>% as.data.frame()

for (i in seq_len(length(network_labels))) {
  s1_neg <- module_size[i]
  l1_neg <- if(i > 1) sum(module_size[1:(i-1)]) else 0
  h1_neg <- l1_neg + s1_neg
  net1_neg <- network_labels[i]
  for (j in seq_len(length(network_labels))) {
    s2_neg <- module_size[j]
    l2_neg <- if(j > 1) sum(module_size[1:(j-1)]) else 0
    h2_neg <- l2_neg + s2_neg
    net2_neg <- network_labels[j]
    num_edges_neg <- sum(data_thresh_neg_canon[(l1_neg+1):h1_neg, (l2_neg+1):h2_neg])
    if (i == j) {
      num_edges_neg <- num_edges_neg / 2
      network_size_neg <- s1_neg * (s1_neg - 1) / 2
    } else {
      network_size_neg <- s1_neg * s2_neg
    }
    df_neg[net1_neg, net2_neg]       <- num_edges_neg
    df_ratio_neg[net1_neg, net2_neg] <- num_edges_neg / network_size_neg
  }
}
df_neg <- df_neg %>% mutate(across(everything(), as.numeric))
mask_neg <- lower.tri(df_neg, diag = FALSE)
df_neg[mask_neg] <- NA

#~~~~Network Heatmap Functions~~~~#
network_labels_reversed <- rev(network_labels) # Reverse for Y axis 

create_heatmap_plot <- function(data, x, y, fill, network_labels, network_labels_reversed) {
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) +
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s > 100, 'white', 'black')", fill)),
              size = 5.5) +
    scale_fill_gradientn(
      colors   = c("ivory", "#FFE4E1", "lightcoral", "salmon", "indianred2", "brown3", "firebrick4"),
      limits   = c(min_fill, max_fill),
      na.value = "white"
    ) +
    scale_color_identity() +
    scale_x_discrete(labels = network_labels) +
    scale_y_discrete(labels = network_labels_reversed) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.text.y = element_text(face = "bold", size = 12),
          legend.position = "none") +
    coord_fixed()
}

create_heatmap_plot_negative <- function(data, x, y, fill, network_labels, network_labels_reversed) {
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) +
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s < 0, 'white', 'black')", fill)),
              size = 5.5) +
    scale_fill_gradientn(
      colors   = c("ivory", "lightblue", "skyblue2", "dodgerblue", "dodgerblue2", "dodgerblue3", "royalblue3"),
      limits   = c(min_fill, max_fill),
      na.value = "white"
    ) +
    scale_color_identity() +
    scale_x_discrete(labels = network_labels) +
    scale_y_discrete(labels = network_labels_reversed) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.text.y = element_text(face = "bold", size = 12),
          legend.position = "none") +
    coord_fixed()
}

#~~~~Positive Heatmap~~~~#
df_pos_melt <- melt(as.matrix(df_pos), varnames = c("x", "y"), value.name = "value")
df_pos_melt$y <- factor(df_pos_melt$y, levels = rev(colnames(df_pos)))
positive_edge_plot <- create_heatmap_plot(df_pos_melt, "x", "y", "value", network_labels, network_labels_reversed)
output_file <- file.path(output_path, "positive_edge_plot.jpeg")
ggsave(filename = output_file, plot = positive_edge_plot, dpi = 300, width = 8, height = 6, units = "in")

#~~~~Negative Heatmap~~~~#
df_neg_melt <- melt(as.matrix(df_neg), varnames = c("x", "y"), value.name = "value")
df_neg_melt$y <- factor(df_neg_melt$y, levels = rev(colnames(df_neg)))
negative_edge_plot <- create_heatmap_plot_negative(df_neg_melt, "x", "y", "value", network_labels, network_labels_reversed)
output_file <- file.path(output_path, "negative_edge_plot.jpeg")
ggsave(filename = output_file, plot = negative_edge_plot, dpi = 300, width = 8, height = 6, units = "in")

