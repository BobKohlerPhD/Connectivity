# R/cpm_engine.R
library(ppcor)
# library(caret)

#' Train Connectome-based Predictive Model (CPM)
#'
#' @param train_mat [edges x subjects] matrix of connectivity values
#' @param train_behav [subjects] vector of behavioral measures
#' @param p_thresh p-value threshold for edge selection (default: 0.05)
#' @param mode prediction mode (currently only "linear" is implemented)
#' @param covariates [subjects x n] matrix or data.frame of covariates for partial correlation
#'
#' @return A list containing:
#'   - pos_mask: Logical mask for positive edges
#'   - neg_mask: Logical mask for negative edges
#'
#' @export
train_cpm <- function(train_mat, train_behav, p_thresh = 0.05, mode = "linear", covariates = NULL) {
  # train_mat is [edges x subjects]
  num_edges <- nrow(train_mat)
  
  if (!is.null(covariates)) {
    # Using partial correlation if covariates are provided
    results <- apply(train_mat, 1, function(edge) {
      pcor_res <- pcor.test(edge, train_behav, covariates)
      c(p.value = pcor_res$p.value, estimate = pcor_res$estimate)
    })
    corrs <- results["p.value", ]
    coeffs <- results["estimate", ]
  } else {
    # Using standard Pearson correlation
    results <- apply(train_mat, 1, function(edge) {
      cor_res <- cor.test(edge, train_behav)
      # cor.test returns estimate named 'cor' for pearson
      # We extract it by index or use names
      c(p.value = cor_res$p.value, estimate = as.numeric(cor_res$estimate))
    })
    corrs <- results["p.value", ]
    coeffs <- results["estimate", ]
  }
  
  pos_mask <- (coeffs > 0) & (corrs < p_thresh)
  neg_mask <- (coeffs < 0) & (corrs < p_thresh)
  
  # Return masks and trained models (simplified for now as per Task 4 goal)
  list(pos_mask = pos_mask, neg_mask = neg_mask)
}

# Add kfold_cpm with future support later
