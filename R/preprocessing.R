# R/preprocessing.R
#' @title Align matrices and metadata
#' @description Intersects subject IDs between a list of matrices and a metadata data frame, 
#' removing subjects with missing outcomes.
#' @param mats A named list of connectivity matrices.
#' @param meta A data frame containing metadata.
#' @param id_var Character string for the subject ID column name in metadata.
#' @param outcome_var Character string for the outcome column name in metadata.
#' @param covariates Optional character vector of covariate column names.
#' @return A list with 'matrices' and 'metadata' elements.
#' @export
library(dplyr)

align_data <- function(mats, meta, id_var, outcome_var, covariates = NULL) {
  meta_ids <- as.character(meta[[id_var]])
  mat_ids <- names(mats)
  
  common_ids <- intersect(meta_ids, mat_ids)
  
  filtered_mats <- mats[common_ids]
  filtered_meta <- meta %>%
    filter(!!sym(id_var) %in% common_ids) %>%
    filter(!is.na(!!sym(outcome_var)))
  
  # Final alignment after NA removal
  final_ids <- as.character(filtered_meta[[id_var]])
  final_mats <- filtered_mats[final_ids]
  
  list(matrices = final_mats, metadata = filtered_meta)
}

#' @title Clean matrices
#' @description Removes matrices containing NA or Infinite values.
#' @param mats A list of matrices.
#' @return A list of cleaned matrices.
#' @export
clean_matrices <- function(mats) {
  keep <- vapply(mats, function(m) !any(is.na(m)) && !any(is.infinite(m)), logical(1))
  mats[keep]
}
