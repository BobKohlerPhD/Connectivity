library(tidyverse)
library(MASS)
library(lme4)
library(lmerTest)
library(gt)
library(pROC)


# File paths
matrices_path    <- ""   # .rds of named matrices list: names(mats) = subject IDs
outcome_csv_path <- ""   # .csv with subject ID, outcome, and any covariates
pos_raw_path     <- ""   # positive mask (CSV or TXT of 0/1)
neg_raw_path     <- ""   # negative mask (CSV or TXT of 0/1)
output_dir       <- ""   # if you want to write files
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Variables (as they appear in outcome CSV)
subject_id_var   <- ""   
outcome_var      <- ""   
covariate_vars   <- c() 


# Model choice: "lm", "glm", "lmer", or "glmer"
model_type <- "glm"

# Family for glm/glmer:
#   "auto"      → binomial() if outcome is binary; otherwise gaussian()
#   "binomial", "poisson", or "gaussian"
glm_family <- "auto"

# random-effect grouping variable for lmer/glmer model selection (leave "" if there is none)
random_effect_var <- ""  



# Functions for inputs
is_binary_vector <- function(x) {
  ux <- unique(na.omit(x))
  length(ux) == 2
}

pick_family <- function(df_outcome, glm_family) {
  if (identical(glm_family, "auto")) {
    if (is_binary_vector(df_outcome)) return(binomial())
    return(gaussian())
  }
  switch(glm_family,
         "binomial" = binomial(),
         "poisson"  = poisson(link = "log"),
         "gaussian" = gaussian(),
         stop("glm_family must be one of 'auto','binomial','poisson','gaussian'"))
}

coef_extract <- function(model_obj, term_name) {
  s <- summary(model_obj)$coefficients
  if (is.null(dim(s)) || !(term_name %in% rownames(s))) {
    return(tibble(
      term      = term_name,
      estimate  = NA_real_,
      std_error = NA_real_,
      t_value   = NA_real_,
      z_value   = NA_real_,
      p_value   = NA_real_))
  }
  row   <- s[term_name, , drop = FALSE]
  p_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)"
  else if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)"
  else NA_character_
  
  tibble(
    term      = term_name,
    estimate  = unname(row[1, "Estimate",   drop = TRUE]),
    std_error = unname(row[1, "Std. Error", drop = TRUE]),
    t_value   = unname(if ("t value" %in% colnames(s)) row[1, "t value"] else NA_real_),
    z_value   = unname(if ("z value" %in% colnames(s)) row[1, "z value"] else NA_real_),
    p_value   = unname(if (!is.na(p_col)) row[1, p_col, drop = TRUE] else NA_real_))
}


# Load masks 
read_mask <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) == "csv") {
    as.matrix(read.csv(path, header = TRUE, check.names = FALSE))
  } else {
    as.matrix(read.table(path, header = TRUE, check.names = FALSE))
  }
}
pos_edges <- read_mask(pos_raw_path) == 1
neg_edges <- read_mask(neg_raw_path) == 1


# Join matrices with outcomes
filter_by_outcome <- function(matrices_path, outcome_csv,
                              subject_id_var, outcome_var,
                              covariate_vars = NULL,
                              random_effect_var = "") {
  mats <- read_rds(matrices_path)
  
  sel_vars <- c(subject_id_var, outcome_var, covariate_vars)
  out <- read_csv(outcome_csv, show_col_types = FALSE) %>%
    dplyr::select(all_of(sel_vars)) %>%
    rename(
      subject_id = !!rlang::sym(subject_id_var),
      outcome    = !!rlang::sym(outcome_var)
    ) %>%
    filter(!is.na(outcome))
  
  # Keep subjects with matrix AND outcome
  ids  <- tibble(subject_id = names(mats))
  keep <- inner_join(ids, out, by = "subject_id")
  
  list(matrices = mats[keep$subject_id],
       model_df = keep)
}

filtered <- filter_by_outcome(
  matrices_path,
  outcome_csv_path,
  subject_id_var,
  outcome_var,
  covariate_vars,
  random_effect_var)


# Masking, sum strengths and model

apply_external_mask <- function(mats_list, model_df,
                                pos_edges, neg_edges,
                                covariate_vars = NULL,
                                model_type = c("lm","glm","lmer","glmer"),
                                glm_family = "auto",
                                random_effect_var = "") {
  model_type <- match.arg(model_type)
  
  # Check if lmer/glmer is needed 
  if (random_effect_var == "" && model_type %in% c("lmer","glmer")) {
    message("No random_effect_var given → using fixed-effects ",
            ifelse(model_type == "lmer", "lm", "glm"), " instead.")
    model_type <- ifelse(model_type == "lmer", "lm", "glm")
  }
  
  ut     <- upper.tri(pos_edges)
  mask_p <- pos_edges[ut]
  mask_n <- neg_edges[ut]
  
  strength_df <- purrr::map_dfr(names(mats_list), function(id) {
    vec <- mats_list[[id]][ut]
    tibble(
      subject_id    = id,
      pos_strength  = sum(vec[mask_p], na.rm = TRUE),
      neg_strength  = sum(vec[mask_n], na.rm = TRUE),
      both_strength = sum(vec[mask_p], na.rm = TRUE) - sum(vec[mask_n], na.rm = TRUE)
    )
  })
  
  df <- strength_df %>%
    left_join(model_df, by = "subject_id")
  
  # Scale net strengths and covariates
  cols_to_scale <- intersect(
    c("pos_strength","neg_strength","both_strength", covariate_vars),
    names(df)
  )
  if (length(cols_to_scale) > 0) {
    df <- df %>% mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.x))))
  }
  
  # Drop rows with missing covariates
  if (length(covariate_vars) > 0) {
    df <- df %>% filter(if_all(all_of(covariate_vars), ~ !is.na(.x)))
  }
  
  # Outcome handling
  is_bin <- is_binary_vector(df$outcome)
  fam <- if (model_type %in% c("glm","glmer")) pick_family(df$outcome, glm_family) else NULL
  
  # Cast outcome for model type
  if (model_type %in% c("lm","lmer")) {
    df$outcome <- as.numeric(df$outcome)
  } else if (model_type %in% c("glm","glmer")) {
    if (fam$family == "binomial") {
      # factor with two levels & keep "0","1" ordering the same 
      lvls <- sort(unique(as.character(na.omit(df$outcome))))
      if (all(lvls %in% c("0","1"))) {
        df$outcome <- factor(as.character(df$outcome), levels = c("0","1"))
      } else {
        df$outcome <- factor(df$outcome)
      }
    } else if (fam$family == "gaussian") {
      df$outcome <- as.numeric(df$outcome)
    }
  }
  
  strength_cols <- c("pos_strength","neg_strength","both_strength")
  cov_part <- if (length(covariate_vars) > 0) paste(covariate_vars, collapse = " + ") else ""
  
  make_formula <- function(sc) {
    rhs <- paste(sc, if (cov_part != "") paste("+", cov_part) else "")
    if (random_effect_var != "" && model_type %in% c("lmer","glmer")) {
      as.formula(paste("outcome ~", rhs, "+ (1|", random_effect_var, ")", sep=""))
    } else {
      as.formula(paste("outcome ~", rhs))
    }
  }
  
  # fit models 
  coef_tbl <- map_dfr(strength_cols, function(sc) {
    form <- make_formula(sc)
    
    fit <- switch(model_type,
                  "lm"   = lm(form,  data = df),
                  "glm"  = glm(form, data = df, family = fam),
                  "lmer" = lmer(form, data = df, control = lmerControl(optimizer = "bobyqa")),
                  "glmer"= glmer(form, data = df, family = fam, control = glmerControl(optimizer = "bobyqa"))
    )
    
    row <- coef_extract(fit, sc)
    
    # dispersion check if Poisson
    if (model_type %in% c("glm","glmer") && fam$family == "poisson") {
      pearson_resid <- residuals(fit, type = "pearson")
      row$dispersion <- sum(pearson_resid^2) / df.residual(fit)
    }
    
    # AUC for binomial models (glm/glmer)
    if (model_type %in% c("glm","glmer") && fam$family == "binomial") {
      # If prediction fails (RE or separation), catch and return NA
      prob <- tryCatch(predict(fit, type = "response"), error = function(e) rep(NA_real_, nrow(df)))
      # Only compute if 2 classes w/ probabilities 
      if (all(!is.na(prob)) && nlevels(df$outcome) == 2) {
        roc_obj <- tryCatch(pROC::roc(df$outcome, prob), error = function(e) NULL)
        row$AUC <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
      } else {
        row$AUC <- NA_real_
      }
    }
    
    row
  })
  
  list(strength = df, coefficients = coef_tbl)
}


# Run above code
results <- apply_external_mask(
  mats_list         = filtered$matrices,
  model_df          = filtered$model_df,
  pos_edges         = pos_edges,
  neg_edges         = neg_edges,
  covariate_vars    = covariate_vars,
  model_type        = model_type,
  glm_family        = glm_family,
  random_effect_var = random_effect_var)
print(results$coefficients)

# Table of results
label_map <- c(
  pos_strength  = "Positive Edge Strength",
  neg_strength  = "Negative Edge Strength",
  both_strength = "Combined Strength"
)

results$coefficients %>%
  mutate(term = recode(term, !!!label_map)) %>%
  gt() %>%
  tab_header(title = md("**Model Summary**")) %>%
  fmt_number(columns = c(estimate, std_error, p_value, AUC, dispersion), decimals = 4) %>%
  cols_label(
    term      = "",
    estimate  = "Estimate",
    std_error = "Std. Error",
    t_value   = "t",
    z_value   = "z",
    p_value   = "p",
    AUC       = "AUC (binom.)",
    dispersion= "Dispersion (Poisson)"
  ) %>%
  tab_style(
    style = cell_borders(sides = "right", weight = px(2), color = "gray"),
    locations = cells_body(columns = "term")
  ) %>%
  cols_align(align = "left",   columns = "term") %>%
  cols_align(align = "center", columns = everything()) %>%
  opt_table_lines(extent = "none")


