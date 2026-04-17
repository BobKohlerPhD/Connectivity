library(tidyverse)

# This script takes a csv containing an outcome of interest with corresponding subject IDs and
# pulls connectivity matrices from an .RDS file that have subject IDs matching those in the .csv outcome file

# Paths to subject lists and matrix file 
outcome_csv_path    <- '' # Assumes outcome is in .csv format 
matrices_path       <- '' # Assumes matrices are in RDS format 


# Subject and outcome variable strings
outcome_var         <- '' # will unifying naming to "subject_id" 
subject_var      <- ''


# Output path and file names for filtered data
output_path         <- ''

filtered_outcome        <- '' # saves as .rds 
filtered_matrices   <- '' # saves as .csv



filter_by_outcome <- function(matrices_path, outcome_csv_path, subject_var, outcome_var) {
  mats <- read_rds(matrices_path)
  
  out <- read_csv(outcome_csv) %>%
    select(all_of(subject_var), all_of(outcome_var)) %>%
    filter(!is.na(.data[[outcome_var]]))
  
  ids <- tibble(!!subject_var := names(mats))
  ids_to_keep <- inner_join(ids, out, by = subject_var)
  
  list(matrices = mats[ids_to_keep[[subject_var]]],
       outcome_filtered = ids_to_keep)
}

# Run filtering function
filtered_data_by_outcome <- filter_by_outcome(matrices_path,
                                              outcome_csv_path,
                                              subject_var,
                                              outcome_var)
                         

# Save filtered outcome and matrices 
write_csv(filtered_data_by_outcome$outcome_filtered, file.path(output_path, filtered_outcome))
write_rds(filtered_data_by_outcome$matrices, file.path(output_path, filtered_matrices))

