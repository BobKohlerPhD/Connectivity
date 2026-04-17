# R/plotting_utils.R
library(ggplot2)

#' Plot Connectivity Heatmap
#' 
#' This function takes a connectivity matrix and plots it as a heatmap.
#' 
#' @param matrix A square connectivity matrix.
#' @param node_order Optional vector of indices to rearrange nodes in a canonical order.
#' @param labels Optional labels for the nodes.
#' @param title The title of the plot.
#' @return A ggplot object.
#' @export
plot_connectivity_heatmap <- function(matrix, node_order = NULL, labels = NULL, title = "Consensus Network") {
  if (!is.null(node_order)) {
    matrix <- matrix[node_order, node_order]
  }
  
  if (!is.null(labels)) {
    dimnames(matrix) <- list(labels, labels)
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
