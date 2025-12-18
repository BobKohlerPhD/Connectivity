
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gt)
library(colorspace)
library(scales)
 sibrary(readr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~EDIT THESE PATHS~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# "Combined notation is from old code when this script was more study specific. can be considered third group but will be given the top circle in the venn diagram plot
paths <- list(positive_combined = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/mid_net_pos_thresh_1.0.txt',
              negative_combined = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/mid_net_neg_thresh_1.0.txt',
              positive_group1 = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/nback_net_pos_thresh_1.0.txt',
              negative_group1 = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/nback_net_neg_thresh_1.0.txt',
              positive_group2  = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/sst_net_pos_thresh_1.0.txt',
              negative_group2  = '/Users/bobkohler/Desktop/Current Work/CumulativeAdversityCPM/output/ca_method_1/net_thresh_meth1/sst_net_neg_thresh_1.0.txt')

# ---------- Helper functions ---------- #
read_mask_txt <- function(path) {
  stopifnot(file.exists(path))
  as.matrix(read.table(path, header = FALSE))
}

# Convert pos/neg matrices to a single logical union mask (>0 = present)
# Also zero the diagonal and keep only upper triangle when extracting edges
make_union_mask <- function(pos_group1at, neg_group1at) {
  stopifnot(is.matrix(pos_group1at), is.matrix(neg_group1at))
  stopifnot(nrow(pos_group1at) == ncol(pos_group1at), nrow(neg_group1at) == ncol(neg_group1at))
  stopifnot(all(dim(pos_group1at) == dim(neg_group1at)))
  mask <- (pos_group1at > 0) | (neg_group1at > 0)
  diag(mask) <- FALSE
  mask
}

# Build edge names ("i_j") for TRUE entries of upper triangle
upper_edge_names <- function(mask) {
  UT <- upper.tri(mask)
  idx <- which(UT & mask, arr.ind = TRUE)
  sprintf("%d_%d", idx[,1], idx[,2])
}

# Count positive and negative edges separately (>0) in upper triangle
count_pos_neg <- function(pos_group1at, neg_group1at) {
  UT <- upper.tri(pos_group1at)
  c(total_pos = sum(pos_group1at[UT] > 0, na.rm = TRUE),
    total_neg = sum(neg_group1at[UT] > 0, na.rm = TRUE))
}

# Convert back from "i_j" to (i, j) data frame for any set of edges
edges_to_ij <- function(edge_names) {
  if (length(edge_names) == 0) return(tibble(i = integer(), j = integer()))
  ij <- do.call(rbind, strsplit(edge_names, "_", fixed = TRUE))
  tibble(i = as.integer(ij[,1]), j = as.integer(ij[,2]))
}

# Create 3-circle Venn geometry using sf
make_venn_circles <- function(labels = c("Combined","Group 1","Group 2"), r = 1.25, s = 1.6) {
  stopifnot(length(labels) == 3)
  h <- sqrt(3)/2 * s
  centers <- tibble(
    name = labels,
    x    = c(0, -s/2, s/2),
    y    = c(h/2, 0, 0))
  geoms <- mapply(function(x, y) st_buffer(st_sfc(st_point(c(x, y))), dist = r),
                  centers$x, centers$y, SIMPLIFY = FALSE)
  sfc_circles <- do.call(c, geoms)
  st_crs(sfc_circles) <- NA
  list(circles = st_sf(name = centers$name, geometry = sfc_circles), centers = centers)
}

# Build region polygons + counts for a,b,c sets
build_regions <- function(circles, set_combined, set_group1, set_group2,
                          labels = c("Combined","Group 1","Group 2")) {
  L1 <- labels[1]; L2 <- labels[2]; L3 <- labels[3]
  
  only_1  <- st_difference(circles$geometry[1], st_union(circles$geometry[2:3]))[[1]]
  only_2  <- st_difference(circles$geometry[2], st_union(circles$geometry[c(1,3)]))[[1]]
  only_3  <- st_difference(circles$geometry[3], st_union(circles$geometry[1:2]))[[1]]
  
  int_12  <- st_difference(st_intersection(circles$geometry[1], circles$geometry[2]), circles$geometry[3])[[1]]
  int_13  <- st_difference(st_intersection(circles$geometry[1], circles$geometry[3]), circles$geometry[2])[[1]]
  int_23  <- st_difference(st_intersection(circles$geometry[2], circles$geometry[3]), circles$geometry[1])[[1]]
  int_123 <- st_intersection(st_intersection(circles$geometry[1], circles$geometry[2]), circles$geometry[3])[[1]]
  
  regions <- st_sf(
    region   = c(paste("Only", L1), paste("Only", L2), paste("Only", L3),
                 paste(L1, "&", L2), paste(L1, "&", L3), paste(L2, "&", L3), "All three"),
    geometry = st_sfc(only_1, only_2, only_3, int_12, int_13, int_23, int_123))
  
  regions$count <- c(
    length(setdiff(set_combined, union(set_group1, set_group2))),
    length(setdiff(set_group1,  union(set_combined, set_group2))),
    length(setdiff(set_group2,  union(set_combined, set_group1))),
    length(setdiff(intersect(set_combined, set_group1), set_group2)),
    length(setdiff(intersect(set_combined, set_group2), set_group1)),
    length(setdiff(intersect(set_group1, set_group2), set_combined)),
    length(Reduce(intersect, list(set_combined, set_group1, set_group2))))
  
  regions
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~Create Venn Diagram~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
make_venn_plot <- function(regions, circles, centers, totals,
                           labels = c("Combined","Group 1","Group 2"),
                           green_base = "#2B7A0B", purple_base = "#8E44AD",
                           orange_base = "#FF6A00", white_base = "grey") {
  
  total_unique_edges <- length(unique(unlist(totals$sets)))
  
  regions <- regions %>%
    mutate(
      fill_color = case_when(
        region == paste(labels[1], "&", labels[2]) ~ lighten(green_base, 0.1),
        region == paste(labels[1], "&", labels[3]) ~ lighten(purple_base, 0.1),
        region == "All three"                      ~ lighten(orange_base, 0.1),
        region == paste(labels[2], "&", labels[3]) ~ lighten(white_base,  0.1)))
  
  cent <- st_point_on_surface(st_geometry(regions))
  coords <- st_coordinates(cent)
  regions <- bind_cols(regions, as_tibble(coords))
  
  # label positions by index (avoid name matching)
  centers$idx <- seq_len(nrow(centers))
  label_pos <- centers %>%
    left_join(totals %>% select(name, total), by = c("name" = "name")) %>%
    mutate(
      label = paste(name, paste0("Total = ", total), sep = "\n"),
      label_x = case_when(idx == 1 ~ x,
                          idx == 2 ~ x - (1.25 - 0.1),
                          idx == 3 ~ x + (1.25 - 0.1)),
      label_y = case_when(idx == 1 ~ y + (1.25 + 0.4),
                          idx == 2 ~ y - (1.25 + 0.1),
                          idx == 3 ~ y - (1.25 + 0.1)))
  ggplot() +
    geom_sf(data = regions, aes(fill = fill_color), color = NA) +
    geom_sf(data = circles, fill = NA, color = "black", linewidth = 3) +
    geom_text(data = regions, aes(X, Y, label = count), size = 12, family = "Arial") +
    geom_text(data = label_pos, aes(label_x, label_y, label = label),
              size = 14, family = "Arial", lineheight = 0.9) +
    scale_fill_identity() +
    coord_sf(datum = NA, clip = "off", expand = TRUE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin     = margin(1, 1, 1, 1, "cm"),
      plot.title      = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 12))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~MAIN CODE~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
run_cpm_overlap <- function(paths,
                            labels = c("Combined", "Group 1", "Group 2"),
                            save_plot = NULL,
                            save_table = NULL) {
  stopifnot(length(labels) == 3) 
  # Read masks
  pos_combined <- read_mask_txt(paths$positive_combined); neg_combined <- read_mask_txt(paths$negative_combined)
  pos_group1 <- read_mask_txt(paths$positive_group1); neg_group1 <- read_mask_txt(paths$negative_group1)
  pos_group2 <- read_mask_txt(paths$positive_group2);  neg_group2 <- read_mask_txt(paths$negative_group2)
  
  # Dimensions check
  dims <- rbind(dim(pos_combined), dim(neg_combined), dim(pos_group1), dim(neg_group1), dim(pos_group2), dim(neg_group2))
  if (length(unique(apply(dims, 1, paste, collapse = "x"))) != 1)
    stop("All input matrices must have the same dimensions.")
  
  # Union masks
  mask_combined <- make_union_mask(pos_combined, neg_combined)
  mask_group1 <- make_union_mask(pos_group1, neg_group1)
  mask_group2 <- make_union_mask(pos_group2, neg_group2)
  
  # Edge sets (upper triangle only)
  set_combined <- upper_edge_names(mask_combined)
  set_group1 <- upper_edge_names(mask_group1)
  set_group2 <- upper_edge_names(mask_group2)
  
  # Totals per circle for labels
  counts_c <- count_pos_neg(pos_combined, neg_combined)
  counts_m <- count_pos_neg(pos_group1, neg_group1)
  counts_f <- count_pos_neg(pos_group2, neg_group2)
  
  totals_tbl <- tibble(
    name  = labels,
    total = c(length(set_combined), length(set_group1), length(set_group2)),
    sets  = list(set_combined, set_group1, set_group2))
  
  # Summary tables for each group
  tbl_summary_combined <- tibble(
    Metric = c("Total Positive", "Total Negative", "Overlap with Group 1", "Overlap with Group 2"),
    Count  = c(counts_c["total_pos"], counts_c["total_neg"],
               length(intersect(set_combined, set_group1)), length(intersect(set_combined, set_group2))))
  
  tbl_summary_group1 <- tibble(
    Metric = c("Total Positive", "Total Negative", "Overlap with Combined", "Overlap with Group 2"),
    Count  = c(counts_m["total_pos"], counts_m["total_neg"],
               length(intersect(set_group1, set_combined)), length(intersect(set_group1, set_group2))))
  
  tbl_summary_group2 <- tibble(
    Metric = c("Total Positive", "Total Negative", "Overlap with Combined", "Overlap with Group 1"),
    Count  = c(counts_f["total_pos"], counts_f["total_neg"],
               length(intersect(set_group2, set_combined)), length(intersect(set_group2, set_group1))))
  
  # Build Venn circle + regions
  vg      <- make_venn_circles(labels = labels, r = 1.25, s = 1.6)
  regions <- build_regions(vg$circles, set_combined, set_group1, set_group2, labels = labels)
  
  # Plot
  venn_plot <- make_venn_plot(regions, vg$circles, vg$centers, totals_tbl, labels = labels)
  
  # save 
  if (!is.null(save_plot)) {
    ggsave(filename = save_plot, plot = venn_plot, width = 8, height = 6, dpi = 300)
  }
  
  if (!is.null(save_table)) {
    write_csv(regions %>% st_drop_geometry(), save_table)
  }
  
  # Intersections for (i,j) pairs
  inter_combined_group1          <- intersect(set_combined, set_group1)
  inter_combined_group2          <- intersect(set_combined, set_group2)
  inter_group1_group2            <- intersect(set_group1,  set_group2)
  inter_combined_group1_group2   <- Reduce(intersect, list(set_combined, set_group1, set_group2))
  
  list(
    sets = setNames(
      list(set_combined, set_group1, set_group2),
      labels),
    intersections = setNames(
      list(
        edges_to_ij(inter_combined_group1),
        edges_to_ij(inter_combined_group2),
        edges_to_ij(inter_group1_group2),
        edges_to_ij(inter_combined_group1_group2)),
      c(paste(labels[1],"&",labels[2]),
        paste(labels[1],"&",labels[3]),
        paste(labels[2],"&",labels[3]),
        "All three")),
    summaries = setNames(
      list(
        tibble(Metric = c("Total Positive","Total Negative",
                          paste("Overlap with", labels[2]),
                          paste("Overlap with", labels[3])),
               Count  = c(counts_c["total_pos"], counts_c["total_neg"],
                          length(intersect(set_combined, set_group1)),
                          length(intersect(set_combined, set_group2)))),
        tibble(Metric = c("Total Positive","Total Negative",
                          paste("Overlap with", labels[1]),
                          paste("Overlap with", labels[3])),
               Count  = c(counts_m["total_pos"], counts_m["total_neg"],
                          length(intersect(set_group1, set_combined)),
                          length(intersect(set_group1, set_group2)))),
        tibble(Metric = c("Total Positive","Total Negative",
                          paste("Overlap with", labels[1]),
                          paste("Overlap with", labels[2])),
               Count  = c(counts_f["total_pos"], counts_f["total_neg"],
                          length(intersect(set_group2, set_combined)),
                          length(intersect(set_group2, set_group1))))),
      labels),
    regions_table = regions %>% st_drop_geometry(),
    venn_plot = venn_plot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~Run CPM Overlap Function~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
result <- run_cpm_overlap(
  paths,
  labels = c("MID", "NBACK", "SST"),   # update as needed 
  save_plot  = NULL,
  save_table = NULL)
result$venn_plot
