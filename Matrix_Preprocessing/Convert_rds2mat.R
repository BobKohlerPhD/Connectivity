rds_matrices <- readRDS() # Path to .rds file 


# Extract relevant components from .rds file
subject_id <- names(rds_matrices) # Extract subject list 
connectivity_matrices <- unname(rds_matrices) # Extract connectivity 

number_of_rows      <- nrow(connectivity_matrices[[1]])
number_of_columns   <- ncol(connectivity_matrices[[1]])
number_of_subjects  <- length(connectivity_matrices)

# Flatten matrices into vector and then convert to 3-D array 
connectivity_unlist <- array(
  unlist(connectivity_matrices),
  dim = c(number_of_rows, number_of_columns, number_of_subjects))

R.matlab::writeMat(
  "", # Path for writing the .mat file with connectivity and subject id 
  connectivity  = connectivity_unlist,
  subject_id = subject_id)

write.csv(subject_id, "") # Path for writing the subject list .csv file 
