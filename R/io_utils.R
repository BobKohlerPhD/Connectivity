#' Load Connectivity Matrices
#'
#' This function loads connectivity matrices from .rds or .mat files.
#' For .mat files, it expects a 'connectivity' field and a 'subject_id' field.
#'
#' @param path String. Path to the file to load.
#' @return A list of connectivity matrices, named by subject ID.
#' @export
load_matrices <- function(path) {
  ext <- tools::file_ext(path)
  if (ext == "rds") {
    return(readRDS(path))
  } else if (ext == "mat") {
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
      stop("Package 'R.matlab' is required for .mat files. Please install it.")
    }
    data <- R.matlab::readMat(path)
    # Handle both list of matrices and 3D array
    if ("connectivity" %in% names(data)) {
      mats <- data$connectivity
      ids <- as.vector(data$subject_id)
      
      # Handle 3D array [row, col, subject]
      if (length(dim(mats)) == 3) {
        res <- lapply(seq_len(dim(mats)[3]), function(i) mats[,,i])
      } else {
        # Assuming list or other format
        res <- mats
      }
      
      names(res) <- ids
      return(res)
    }
  }
  stop("Unsupported matrix format or missing expected fields in .mat file")
}

#' Validate Connectivity Matrices
#'
#' Checks if matrices are symmetric and have consistent dimensions.
#'
#' @param mats A list of matrices.
#' @return Logical TRUE if valid, otherwise throws an error.
#' @export
validate_matrices <- function(mats) {
  if (length(mats) == 0) stop("No matrices found")
  
  node_count <- nrow(mats[[1]])
  
  for (i in seq_along(mats)) {
    if (!is.matrix(mats[[i]])) stop(sprintf("Element %d is not a matrix", i))
    if (!isSymmetric(mats[[i]], tol = 1e-8)) stop(sprintf("Matrix %d is not symmetric", i))
    if (nrow(mats[[i]]) != node_count) stop("Inconsistent node counts")
    if (ncol(mats[[i]]) != node_count) stop("Matrices must be square")
  }
  
  return(TRUE)
}
