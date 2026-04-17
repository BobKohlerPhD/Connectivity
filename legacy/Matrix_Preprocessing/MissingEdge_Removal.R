matrices <- readRDS("")
str(matrices)


# Define a 'missing' edge
is_zero <- function(x, tol = 0) x == 0                     # exactly 0
# is_zero <- function(x, tol = 1e-12) abs(x) <= tol        # near-zero option

# Keep only matrices with NO NAs and NO zeros
keep <- vapply(
  matrices,
  function(M) !anyNA(M) && !any(is_zero(M)),
  logical(1)
)

# Filter
mats_clean <- matrices[keep]
length(mats_clean)

readr::write_rds(mats_clean, "")


dropped_ids <- names(matrices)[!keep]


