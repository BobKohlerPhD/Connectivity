
#~~~~This script performs the following~~~~#
# 1. Reads in CPM model output files 
# 2. Compute performance metrics
# 3. Creates visualizations of performance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(gt)
options(scipen = 999)

#~~~~Path to CPM output folder~~~~#
read_path <-''


#~~~~CPM Parameters that were used~~~~#
regression_type <- "linear" # Type of Regression ("linear" or "logistic" currently implemented)


n_repeats    <- 100 # Number of true repeats
n_iterations <- 1 # Null (permutation) iterations
lst_of_i     <- 1:(n_repeats + n_iterations)

#~~~~Vectors to store output~~~~#
r_pos  <- numeric()
r_neg  <- numeric()
r_both <- numeric()
data_list <- vector("list", length(lst_of_i))

#~~~~Loop through output folder for y_prediction csv files~~~~#
for (i in lst_of_i) {
  file_path <- file.path(read_path, sprintf("y_prediction_iter%d.csv", i))
  if (file.exists(file_path)) {
    df_tmp <- read.csv(file_path)
    data_list[[i]] <- list(
      y_actual = as.numeric(df_tmp$y_actual),
      y_pos    = as.numeric(df_tmp$y_pred_pos),
      y_neg    = as.numeric(df_tmp$y_pred_neg),
      y_both   = as.numeric(df_tmp$y_pred_both)
    )
  } else {
    warning(sprintf("File not found: %s", file_path))
  }
}

#~~~~Calculate Model Performance~~~~#  
for (data in data_list) {
  if (is.null(data)) next 
  if (regression_type == "linear") {
    r_pos  <- c(r_pos,  cor(data$y_pos,  data$y_actual, method = "spearman", use = "complete.obs"))
    r_neg  <- c(r_neg,  cor(data$y_neg,  data$y_actual, method = "spearman", use = "complete.obs"))
    r_both <- c(r_both, cor(data$y_both, data$y_actual, method = "spearman", use = "complete.obs"))
  } else if (regression_type == "logistic") {
    accuracy_pos  <- mean(data$y_pos == data$y_actual, na.rm = TRUE)
    accuracy_neg  <- mean(data$y_neg == data$y_actual, na.rm = TRUE)
    accuracy_both <- mean(data$y_both == data$y_actual, na.rm = TRUE)
    r_pos  <- c(r_pos, accuracy_pos)
    r_neg  <- c(r_neg, accuracy_neg)
    r_both <- c(r_both, accuracy_both)
  }
}

#~~~~Separate true (first n_repeats) and null (remaining iterations) performance~~~~#
r_pos_true  <- r_pos[1:n_repeats]
r_neg_true  <- r_neg[1:n_repeats]
r_both_true <- r_both[1:n_repeats]
r_pos_null  <- r_pos[(n_repeats + 1):length(r_pos)]
r_neg_null  <- r_neg[(n_repeats + 1):length(r_neg)]
r_both_null <- r_both[(n_repeats + 1):length(r_both)]


#~~~~Compute p-values~~~~#
p_value_one_tail <- function(null, true_mean) {
  if (true_mean >= 0) {
    mean(null >= true_mean, na.rm = TRUE)
  } else {
    warning("Metric is not positive.")
    mean(null <= true_mean, na.rm = TRUE)
  }
}

p_value_two_tail <- function(null, true_mean) {
  mean(abs(null) >= abs(true_mean), na.rm = TRUE)
}

#~~~~One tail~~~~#
p_pos_onetail  <- p_value_one_tail(r_pos_null,  mean(r_pos_true))
p_neg_onetail  <- p_value_one_tail(r_neg_null,  mean(r_neg_true))
p_both_onetail <- p_value_one_tail(r_both_null, mean(r_both_true))

#~~~~Two tail~~~~#
p_pos_twotail  <- p_value_two_tail(r_pos_null,  mean(r_pos_true))
p_neg_twotail  <- p_value_two_tail(r_neg_null,  mean(r_neg_true))
p_both_twotail <- p_value_two_tail(r_both_null, mean(r_both_true))


#~~~~Summary Tables of CPM Performance~~~#
(df_boxplot <- tibble(
  `Null Model Size` = c(length(r_pos_null), length(r_neg_null), length(r_both_null)),
  `True Model Size` = c(length(r_pos_true), length(r_neg_true), length(r_both_true)),
  `Metric (Mean)` = c(mean(r_pos_true, na.rm = TRUE), mean(r_neg_true, na.rm = TRUE), mean(r_both_true, na.rm = TRUE)),
  `p-value (One-Tail)` = c(p_pos_onetail, p_neg_onetail, p_both_onetail),
  `Metric (Median)` = c(median(r_pos_true, na.rm = TRUE), median(r_neg_true, na.rm = TRUE), median(r_both_true, na.rm = TRUE))))

df_p_table <- tibble(
  Edge = c("Positive Edges", "Negative Edges", "Both Edges"),
  #`Null Model Iterations` = c(length(r_pos_null), length(r_neg_null), length(r_both_null)),
  `True Metric (Mean)` = c(mean(r_pos_true, na.rm = TRUE),
                           mean(r_neg_true, na.rm = TRUE),
                           mean(r_both_true, na.rm = TRUE)),
  `True Metric (SD)` = c(sd(r_pos_true, na.rm = TRUE),
                         sd(r_neg_true, na.rm = TRUE),
                         sd(r_both_true, na.rm = TRUE)),
  `p-value (One-Tail)` = c(p_pos_onetail, p_neg_onetail, p_both_onetail),
  `p-value (Two-Tail)` = c(p_pos_twotail, p_neg_twotail, p_both_twotail))


#~~~~Accuracy across iterations~~~~#
df_p_table_fmt <- df_p_table %>%
  mutate(
    # One‐tailed
    `One-Tailed p-value` = if_else(
      `p-value (One-Tail)` < 0.001,
      "*p* < 0.001",
      sprintf("*p* = %.3f", `p-value (One-Tail)`)
    ),
    # Two‐tailed
    `Two-Tailed p-value` = if_else(
      `p-value (Two-Tail)` < 0.001,
      "*p* < 0.001",
      sprintf("*p* = %.3f", `p-value (Two-Tail)`)
    )
  ) %>%
  # drop the original numeric p‐value cols if you like
  select(-`p-value (One-Tail)`, -`p-value (Two-Tail)`)

df_gt <- df_p_table_fmt %>%
  gt(rowname_col = "Edge") %>%
  tab_header(
    title    = "Sex-Agnostic ERT Performance Summary",
  #  subtitle = "Summary for Positive, Negative, and Both Edges"
  ) %>%
  fmt_markdown(
    columns = c(`One-Tailed p-value`, `Two-Tailed p-value`)
  ) %>%
  fmt_number(
    columns  = c(`True Metric (Mean)`, `True Metric (SD)`),
    decimals = 4
  ) %>%
  cols_label(
    `True Metric (Mean)`    = "Mean",
    `True Metric (SD)`      = "SD",
    `One-Tailed p-value`    = "One‐Tailed p",
    `Two-Tailed p-value`    = "Two‐Tailed p"
  ) %>%
  cols_align(
    align   = "center",
    columns = everything()
  ) %>%
  tab_style(
    style = cell_borders(
      sides  = "all",
      color  = "black",
      weight = px(2)
    ),
    locations = list(
      cells_body(),
      cells_column_labels(),
      cells_stub()
    )
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgrey"),
      cell_borders(
        sides  = c("left","bottom","top"), 
        color  = "black", 
        weight = px(2)
      )
    ),
    locations = cells_stubhead()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_title(groups = "title")
  ) %>%
  tab_options(
    table.border.top.style      = "solid",
    table.border.top.width      = px(2),
    table.border.top.color      = "black",
    table.border.bottom.style   = "solid",
    table.border.bottom.width   = px(2),
    table.border.bottom.color   = "black",
    heading.border.bottom.style = "solid",
    heading.border.bottom.width = px(2),
    heading.border.bottom.color = "black",
    table_body.hlines.color     = "black",
    table_body.vlines.color     = "black",
    column_labels.border.bottom.color = "black")
print(df_gt)


# df_gt <- df_p_table %>%
#   dplyr::select(-contains("Null Model")) %>% 
#   gt(rowname_col = "Edge") %>%
#   tab_header(
#     title = "Performance Summary Table",
#     subtitle = "Summary for Positive, Negative, and Both Edges") %>%
#   fmt_number(
#     columns = c(`True Metric (Mean)`, `True Metric (SD)`,
#                 `p-value (One-Tail)`, `p-value (Two-Tail)`),
#     decimals = 4) %>%
#   fmt_number(
#     columns = c(`p-value (One-Tail)`, `p-value (Two-Tail)`),
#     decimals = 0) %>%
#   cols_label(
#    # `Null Model Iterations` = "Null Iterations",
#     `True Metric (Mean)` = "Mean",
#     `True Metric (SD)` = "SD",
#     `p-value (One-Tail)` = "One-Tailed p-value",
#     `p-value (Two-Tail)` = "Two-Tailed p-value") %>%
#   cols_align(
#     align = "center",
#     columns = everything()) %>%
#   tab_style(
#     style = cell_borders(
#       sides = "all",
#       color = "black",
#       weight = px(2)),
#     locations = list(
#       cells_body(),
#       cells_column_labels(),
#       cells_stub())) %>%
#   tab_options(
#     table.border.top.style = "solid",
#     table.border.top.width = px(2),
#     table.border.top.color = "black",
#     table.border.bottom.style = "solid",
#     table.border.bottom.width = px(2),
#     table.border.bottom.color = "black",
#     heading.border.bottom.style = "solid",
#     heading.border.bottom.width = px(2),
#     heading.border.bottom.color = "black",
#     table_body.hlines.color = "black",
#     table_body.vlines.color = "black",
#     column_labels.border.bottom.color = "black")
# print(df_gt)  




#~~~~Boxplot of Performance~~~~#
level_order <- c('Positive', 'Negative', 'Both') 
create_plot <- function(data, x_var, y_var, edge_var, level_order, fill_values, color_values) {
  ggplot(data, aes(x = factor(!!rlang::sym(x_var), levels = level_order), 
                   y = !!rlang::sym(y_var), color = !!rlang::sym(edge_var), fill = !!rlang::sym(edge_var))) +
    geom_boxplot(aes(fill = !!rlang::sym(edge_var),
                     fill = after_scale(colorspace::lighten(fill, .2))),
                 color = "black",
                 width = .55,
                 fatten = .75,
                 linewidth = 1.25,
                 outlier.shape = NA)+
    geom_point(position = position_jitter(width = .2, seed = 0), size = 4, alpha = .4) +
    geom_point(position = position_jitter(width = .2, seed = 0), size = 4, stroke = .4, shape = 1) +
    scale_fill_manual(values = fill_values) +
    scale_color_manual(values = color_values) +
    theme_minimal() +
    labs(x = "", y = "Accuracy" , title = "Average Cross-Validated Performance by Edge Type")+
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.line = element_line(linewidth = 1.25),
          axis.ticks.length = unit(.25, "cm"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.text.x = element_text(size = 22, color = "black", face = "bold"),
          legend.text = element_text(size = 0, color = "black", face = "bold"),
          legend.title = element_text(size = 0, color = "black", face = "bold"),
          legend.position = "none")
}

r_true_df <- data.frame(Positive = r_pos_true, 
                        Negative = r_neg_true, 
                        Both = r_both_true) %>% 
  pivot_longer(cols = everything(),
               names_to = "edge",
               values_to = "r_value",
               values_drop_na = FALSE)
                  
(r_true_boxplot <- 
    create_plot(r_true_df, 
                x_var = "edge", 
                y_var = "r_value", 
                edge_var = "edge", 
                level_order = level_order, 
                fill_values = c('Positive' = 'red', 
                                'Negative' = 'blue', 
                                'Both' = 'purple'), 
                color_values = c('Positive' = 'red', 
                                 'Negative' = 'blue',
                                 'Both' = 'purple')))

output_file <- file.path(read_path, "boxplot_performance.jpeg")
#ggsave(filename = output_file, plot = r_true_boxplot, dpi = 300, width = 8, height = 6, units = "in")

                              