library(tidyverse)


# RDS with square matrices + subject ids so that names(mats) can extract IDs
mat_path <- ""  # e.g., "/path/to/matrices.rds"

# Load mats
mats <- readRDS(mat_path)

# Calculate using upper triangle (no diagonal)
n <- nrow(mats[[1]])
UT <- upper.tri(matrix(FALSE, n, n), diag = FALSE)

strength_df <- purrr::map_dfr(names(mats), function(id) {
  mat <- mats[[id]]
  vec <- mat[UT]
  pos_strength  <- sum(vec[vec > 0], na.rm = TRUE)
  neg_strength  <- sum(vec[vec < 0], na.rm = TRUE)    
  both_strength <- pos_strength - neg_strength        # adds magnitude of negatives but can change to absolute value if preferred 

   tibble(
    subject_id    = id,
    pos_strength  = pos_strength,
    neg_strength  = neg_strength,
    both_strength = both_strength
  )
})

glimpse(strength_df)
# readr::write_csv(strength_df, "network_strengths.csv")
