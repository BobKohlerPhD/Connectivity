library(caret)
library(dplyr)
library(purrr)
library(readr)
library(R.matlab)
library(jsonlite)
library(parallel)
library(fs)
library(stringr)
library(lubridate)
library(ppcor) 

#~~Model~~#
train_cpm <- function(train_mat, 
                      train_behav, 
                      num_nodes,
                      p_thresh = 0.05,
                      mode = "linear",
                      covariates = NULL) {
  # Check for covariates 
   if (!is.null(covariates)) { 
    if (is.data.frame(covariates)) {
      cov_mat <- model.matrix(~ . - 1, data = covariates) # Dummy code covariates if necessary
    } else {
      cov_mat <- covariates
    }
  } else {
    cov_mat <- NULL
  }
  
  # Correlations (partial correlations are computed when covariates are provided)
  if (!is.null(cov_mat)) {
    corr_train <- map(train_mat, ~ ppcor::pcor.test(.x, train_behav, cov_mat)) # Partial correlations
    r_lst <- sapply(corr_train, function(res) res$estimate)
    p_lst <- sapply(corr_train, function(res) res$p.value)
  } else {
    corr_train <- map(train_mat, ~ cor.test(.x, train_behav, method = "pearson")) # Pearson correlations
    r_lst <- sapply(corr_train, `[[`, "estimate")
    p_lst <- sapply(corr_train, `[[`, "p.value")
  }
  
  # Create empty matrices
  r_mat <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  p_mat <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  upper_tri <- upper.tri(r_mat)
  diag(r_mat) <- NA
  diag(p_mat) <- NA
  
  r_mat[upper_tri] <- r_lst
  p_mat[upper_tri] <- p_lst
  r_mat[lower.tri(r_mat)] <- t(r_mat)[lower.tri(r_mat)]
  p_mat[lower.tri(p_mat)] <- t(p_mat)[lower.tri(p_mat)]
  
  # Symmetry Check of r and p mats 
  if (!isSymmetric(r_mat) || !isSymmetric(p_mat)) {
    stop("ERROR: r_mat or p_mat is not symmetric. Please check your data.")
  }
  
  # Identify Significant Edges
  pos_edges <- (r_mat > 0) & (p_mat < p_thresh)
  neg_edges <- (r_mat < 0) & (p_mat < p_thresh)
  # Sum positive and negative edges 
  pos_sum <- colSums(train_mat[pos_edges[upper_tri], , drop = FALSE])
  neg_sum <- colSums(train_mat[neg_edges[upper_tri], , drop = FALSE])
  both <- pos_sum - neg_sum
  
  #  Model Train Function
  train_model <- function(y, x, mode) {
    if (mode == "ridge") {
      train(x, y, method = "ridge", tuneLength = 10)
    } else if (mode == "linear") {
      lm(y ~ x)
    } else if (mode == "logistic") {
      glm(y ~ x, family = binomial())
    } else {
      stop("Mode not implemented!")
    }
  }
  
  pos_estimator <- if (sum(pos_edges) > 0) train_model(train_behav, pos_sum, mode) else NA
  neg_estimator <- if (sum(neg_edges) > 0) train_model(train_behav, neg_sum, mode) else NA
  both_estimator <- if (sum(pos_edges | neg_edges) > 0) train_model(train_behav, both, mode) else NA
  
  
  list(pos_estimator = pos_estimator, 
       neg_estimator = neg_estimator, 
       both_estimator = both_estimator, 
       pos_edges = pos_edges, 
       neg_edges = neg_edges)
  }



kfold_cpm <- function(x, y, k,
                      p_thresh = 0.05,
                      zscore = FALSE, 
                      mode = "linear", 
                      covariates = NULL) {
  num_subs <- dim(x)[3]
  num_nodes <- dim(x)[1]
  upper_tri <- upper.tri(matrix(0, num_nodes, num_nodes))
  
  # Flatten Edge Matrix
  all_edges <- apply(x, 3, function(mat) mat[upper_tri])
  
  rand_inds <- sample(seq_len(num_subs))
  sample_size <- floor(num_subs / k)
  
  results <- list(
    y_pred_pos = numeric(num_subs),
    y_pred_neg = numeric(num_subs),
    y_pred_both = numeric(num_subs),
    fit_p = vector("list", k),
    fit_n = vector("list", k),
    fit_b = vector("list", k),
    edges_p = vector("list", k),
    edges_n = vector("list", k)
  )
  
  for (fold in seq_len(k)) {
    test_inds <- if (fold != k) {
      rand_inds[((fold - 1) * sample_size + 1):(fold * sample_size)]
    } else {
      rand_inds[((fold - 1) * sample_size + 1):num_subs]
    }
    train_inds <- setdiff(rand_inds, test_inds)
    
    train_mats <- all_edges[, train_inds]
    train_behav <- y[train_inds]
    test_mats <- all_edges[, test_inds]
    
    if (zscore) {
      train_mats <- scale(train_mats)
      test_mats <- scale(test_mats)
    }
    
    cpm <- train_cpm(train_mats, 
                     train_behav,
                     num_nodes,
                     p_thresh,
                     mode, 
                     covariates)
    pos_sum <- colSums(test_mats[cpm$pos_edges[upper_tri], , drop = FALSE])
    neg_sum <- colSums(test_mats[cpm$neg_edges[upper_tri], , drop = FALSE])
    both <- pos_sum - neg_sum
    
    results$y_pred_pos[test_inds] <- predict(cpm$pos_estimator, pos_sum)
    results$y_pred_neg[test_inds] <- predict(cpm$neg_estimator, neg_sum)
    results$y_pred_both[test_inds] <- predict(cpm$both_estimator, both)
    results$fit_p[[fold]] <- cpm$pos_estimator
    results$fit_n[[fold]] <- cpm$neg_estimator
    results$fit_b[[fold]] <- cpm$both_estimator
    results$edges_p[[fold]] <- cpm$pos_edges
    results$edges_n[[fold]] <- cpm$neg_edges
  }
  
  results
}

run_cpm_thread <- function(y, iter, x, k, out_path, p_thresh = 0.05, zscore = FALSE, mode = "linear", covariates = NULL) {
  message(sprintf("Running iteration #%d", iter))
  results <- kfold_cpm(x, y, k, p_thresh, zscore, mode, covariates)
  save_run_outputs(out_path, iter, results, y, mode)
}

#----------------------------------------#
#----------Utility Functions-------------#
#----------------------------------------#

save_run_outputs_subsample_nogrid <- function(out_path, iter, outputs, y) {
  for (fold in seq_along(outputs$edges_p)) {
    write.table(outputs$edges_p[[fold]], 
                file = file.path(out_path, sprintf("positive_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
    write.table(outputs$edges_n[[fold]], 
                file = file.path(out_path, sprintf("negative_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
  }
  
  df_y_predict <- tibble(y_pred_both = outputs$y_pred_both, y_actual = y)
  write_csv(df_y_predict, file.path(out_path, sprintf("y_prediction_iter%d.csv", iter)))
  
  df_fit <- tibble(
    both_m = map_dbl(outputs$fit_b, ~ .x$coefficients[2]),
    both_b = map_dbl(outputs$fit_b, ~ .x$coefficients[1])
  )
  write_csv(df_fit, file.path(out_path, sprintf("fit_parameters_iter%d.csv", iter)))
}

save_run_outputs_subsample <- function(out_path, iter, outputs, y) {
  save_run_outputs_subsample_nogrid(out_path, iter, outputs, y)  # Save basic outputs
  
  for (fold in seq_along(outputs$best_params)) {
    df_params <- as_tibble(outputs$best_params[[fold]])
    write_csv(df_params, file.path(out_path, sprintf("best_params_iter%d_fold_%d.csv", iter, fold)))
  }
}

y_transform <- function(y, y_norm = "id") {
  if (y_norm == "yj") {
    transformer <- preProcess(y, method = c("YeoJohnson", "center", "scale"))
  } else if (y_norm == "id") {
    transformer <- preProcess(y, method = "none")
  } else if (y_norm == "norm") {
    transformer <- preProcess(y, method = c("center", "scale"))
  } else {
    warning(sprintf("Undefined y_norm %s. Using identity function instead.", y_norm))
    transformer <- preProcess(y, method = "none")
  }
  yn <- predict(transformer, as.data.frame(y))
  list(yn = yn[, 1], transformer = transformer)
}

save_matlab_mat <- function(path, matname, x, y, lst_subj) {
  mdict <- list(x = x, y = y, subjectkey = lst_subj)
  writeMat(file.path(path, matname), mdict)
}

read_matlab_mat <- function(path, matname) {
  mdict <- readMat(file.path(path, matname))
  list(x = mdict$x, y = mdict$y, lst_subjectkey = mdict$subjectkey)
}

check_symmetric <- function(a, rtol = 1e-5, atol = 1e-8) {
  all(abs(a - t(a)) <= (atol + rtol * abs(t(a))), na.rm = TRUE)
}

generate_file_list <- function(path, lst_subj, num_roi, num_contrasts, t) {
  sprintf("%s/%s_%dROI_%dcontrasts_corr_matrix_%s.txt", path, lst_subj, num_roi, num_contrasts, t)
}

read_mats <- function(fn_list) {
  fns <- map(fn_list, ~ read.table(.x, header = FALSE))
  
  if (any(map_lgl(fns, ~ any(is.na(.x))))) {
    stop("ERROR: there are NaNs in the correlation matrices! Please check your data.")
  }
  
  array(unlist(fns), dim = c(nrow(fns[[1]]), ncol(fns[[1]]), length(fns)))
}

return_estimator_coef <- function(est, mode) {
  if (mode %in% c("linear", "logistic")) {
    if (is.null(est)) {
      c(NA, NA)
    } else {
      c(coef(est)[2], coef(est)[1])
    }
  } else if (mode == "ridge") {
    if (is.null(est)) {
      c(NA, NA, NA)
    } else {
      c(coef(est$bestTune)[2], coef(est$bestTune)[1], est$bestTune$alpha)
    }
  } else {
    stop(sprintf("ERROR: mode %s not implemented!", mode))
  }
}

save_run_outputs <- function(out_path, iter, outputs, y_run, mode = "linear") {
  for (fold in seq_along(outputs$edges_p)) {
    write.table(outputs$edges_p[[fold]], 
                file = file.path(out_path, sprintf("positive_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
    write.table(outputs$edges_n[[fold]], 
                file = file.path(out_path, sprintf("negative_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
  }
  
  df_y_predict <- tibble(
    y_pred_pos = outputs$y_pred_pos,
    y_pred_neg = outputs$y_pred_neg,
    y_pred_both = outputs$y_pred_both,
    y_actual = y_run
  )
  write_csv(df_y_predict, file.path(out_path, sprintf("y_prediction_iter%d.csv", iter)))
  
  df_fit <- tibble(
    pos_m = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[1]),
    pos_b = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[2]),
    neg_m = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[1]),
    neg_b = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[2]),
    both_m = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[1]),
    both_b = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[2])
  )
  
  if (mode == "ridge") {
    df_fit <- df_fit %>% 
      mutate(
        pos_alpha = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[3]),
        neg_alpha = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[3]),
        both_alpha = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[3])
      )
  }
  
  write_csv(df_fit, file.path(out_path, sprintf("fit_parameters_iter%d.csv", iter)))
}

#------------------------------------------------#
#---------------MAIN CPM FUNCTION----------------#
#------------------------------------------------#

run_cpm_pipeline <- function(json_file, mat_file_path, mat_file_name, output_dir, covariates = NULL) {
  
  # Load JSON configuration
  config <- fromJSON(json_file)
  k <- config$k
  p_thresh <- config$p_thresh
  zscore <- config$zscore
  mode <- config$mode
  num_iterations <- config$num_iter
  
  # Set up output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load .mat file with utility function above 
  data <- read_matlab_mat(mat_file_path, mat_file_name)
  x <- data$x
  y <- data$y
  lst_subj <- data$lst_subjectkey
  
  # Save subject ids
  writeLines(lst_subj, file.path(output_dir, "subject_keys.txt"))
  
  # Parallelize 
  cl <- makeCluster(detectCores() - 1)  # Use all but 1 core
  clusterExport(cl, c("run_cpm_thread", "x", "y", "k", "p_thresh", "zscore", "mode", "output_dir",
                      "kfold_cpm", "train_cpm", "save_run_outputs", "return_estimator_coef"))
  parLapply(cl, seq_len(num_iterations), function(iter) {
    run_cpm_thread(y, iter, x, k, output_dir, p_thresh, zscore, mode, covariates)
  })
  stopCluster(cl)
  
  # Combine output files from each iteration
  message("Combining results...")
  prediction_files <- list.files(output_dir, pattern = "y_prediction_iter.*\\.csv", full.names = TRUE)
  combined_results <- map_dfr(prediction_files, read_csv)
  write_csv(combined_results, file.path(output_dir, "combined_predictions.csv"))
  
  message("CPM pipeline completed. Results saved in ", output_dir)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~Run model by setting your file paths and running functions~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Example usage:
# json_file <- "path/to/config.json"
# mat_file_path <- "path/to/mat_directory"
# mat_file_name <- "your_data.mat"
# output_dir <- "path/to/output_directory"
# covariates <- your_covariate_dataframe_or_matrix   # covariates can be continuous, categorical or discrete
# run_cpm_pipeline(json_file, mat_file_path, mat_file_name, output_dir, covariates)
