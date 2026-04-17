# R/plotting_utils.R
library(ggplot2)

#' Plot Connectivity Heatmap
#' 
#' This function takes a connectivity matrix and plots it as a heatmap.
#' 
#' @param matrix A square connectivity matrix.
#' @param labels Optional labels for the nodes.
#' @param title The title of the plot.
#' @return A ggplot object.
#' @export
plot_connectivity_heatmap <- function(matrix, labels = NULL, title = "Consensus Network") {
  if (!is.null(labels)) {
    dimnames(matrix) <- list(labels, labels)
  }
  
  df <- as.data.frame.table(matrix)
  colnames(df) <- c("Var1", "Var2", "value")
  
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
