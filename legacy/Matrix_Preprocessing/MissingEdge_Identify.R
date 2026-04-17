library(tidyverse)

# Identifies any instance of '0' in connectivity matrix and counts per subject and per edge
rds_matrices <- ""

connectivity_values <- readRDS(rds_matrices)
ids <- names(connectivity_values)

missing_edges_count <- sapply(connectivity_values, function(mat) {
  sum(mat == 0, na.rm = TRUE)
})

df_missing <- tibble::tibble(
  subj_id       = ids,
  missing_edges = missing_edges_count)

cat(sum(missing_edges_count > 0),
    "subjects have at least one missing edge.\n")


x_array <- simplify2array(connectivity_values)
missing_by_edge <- apply(x_array == 0, c(1,2), sum)
edge_df <- as.data.frame(missing_by_edge)

