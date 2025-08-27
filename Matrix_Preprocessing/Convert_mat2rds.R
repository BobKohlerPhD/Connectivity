library(R.matlab)
library(hdf5r)



mat_file_path <- "/Users/bobkohler/Desktop/tp1_MID_FD20_avg_matrix.mat"   # path to your .mat file
output_rds_path <- "/Users/bobkohler/Desktop/"  # path to save .rds

mat_data <- readMat(mat_file_path)


# Assumes that mat_data$connectivity is a 3D array [rows, cols, subjects]
connectivity_array <- mat_data$connectivity 
subject_id <- unlist(mat_data$subject_id)  # unlist to vector 


number_of_subjects <- dim(connectivity_array)[3]
connectivity_list <- vector("list", length = number_of_subjects)

for (i in seq_len(number_of_subjects)) {
  connectivity_list[[i]] <- connectivity_array[ , , i]
}

# Add subject id to each matrix 
names(connectivity_list) <- subject_id

# Save RDS
saveRDS(connectivity_list, output_rds_path)


