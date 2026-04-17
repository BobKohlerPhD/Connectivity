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

#' Run K-fold Cross-Validation for CPM
#'
#' @param x_mats [edges x subjects] matrix of connectivity values
#' @param y_vec [subjects] vector of behavioral measures
#' @param k number of folds
#' @param p_thresh p-value threshold for edge selection (default: 0.05)
#' @param mode prediction mode (default: "linear")
#' @param covariates [subjects x n] matrix or data.frame of covariates
#'
#' @return A vector of predicted values for all subjects
#' @export
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
