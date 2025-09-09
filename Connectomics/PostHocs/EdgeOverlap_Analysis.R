
library(gt)
#~~~~~~~~~~~~~~~~~~~~~#
#~~~Two-Way Overlap~~~#
#~~~~~~~~~~~~~~~~~~~~~#

## Inputs
A <- 1975     # observed overlap (present in BOTH sets)
C <- 2850     # size of set 1 
D <- 5494    # size of set 2 
B <- 35778   # total edges 

## Expected 
expected <- C * D / B

## Strength 
strength <- A / expected

## Similarity
cosine_similarity <- A / sqrt(C * D) # sameness of overall networks
jaccard_similarity <-  A / (C + D - A) # similarity metric more sensitive to size of sets 

overlap_coefficient <- A / min(C, D) # how much the smaller set is explained by the overlap with larger set

## p-values using phyper for hypergeometric 
p_right_tail <- phyper(q = A - 1, m = D, n = B - D, k = C, lower.tail = FALSE)  # enrichment- prob of observing more than expected
p_left_tail   <- phyper(q = A, m = D, n = B - D, k = C, lower.tail = TRUE)   # depletion- prob of observing fewer than expected
p_two_sided <- min(1, 2 * min(p_right_tail, p_left_tail))

results_twoway <- tibble::tibble(
        Metric = c("Expected Overlap", "Overlap Strength", "Overlap Coefficient",
                   "Cosine", "Jaccard", 
                   "p (right tail)", "p (left tail)", "p (two-sided)"),
        Value = c(expected, strength, overlap_coefficient, 
                  cosine_similarity,jaccard_similarity, 
                  p_right_tail, p_left_tail, p_two_sided))   %>%
        gt() %>%
        fmt_number(columns = "Value", decimals = 4) %>% 
        tab_header(
                title = "Female ∩ Sex-Agnostic") #update as needed
print(results_twoway)

     
#~~~~~~~~~~~~~~~~~~~~~~~#
#~~~Three-Way Overlap~~~#
#~~~~~~~~~~~~~~~~~~~~~~~#

# Inputs
All_overlap <- 211  # observed overlap (ALL THREE)
B <- 35778     # total edges
C <- 2850       # set 1
D <- 549       # set 2
E <- 5494      # set 3  

# Expected 
expected_triple <- (C * D * E) / (B^2)

# Strength  
strength_triple <- All_overlap / expected_triple

## Similarity

# Jaccard needs pairwise overlaps for union. May not be know so can leave blank 
A_CD <- 211        # |1 ∩ 2|
A_CE <- 1975      # |1 ∩ 3|
A_DE <- 362       # |2 ∩ 3|


jaccard_triple <- if (!any(is.na(c(A_CD, A_CE, A_DE)))) {
        union3 <- C + D + E - (A_CD + A_CE + A_DE) + All_overlap
        if (union3 > 0) All_overlap / union3 else NA_real_
} else NA_real_


cosine_CD <- if (!is.na(A_CD)) A_CD / sqrt(C * D) else NA_real_
cosine_CE <- if (!is.na(A_CE)) A_CE / sqrt(C * E) else NA_real_
cosine_DE <- if (!is.na(A_DE)) A_DE / sqrt(D * E) else NA_real_

cosine_triple <- All_overlap / sqrt(C * D * E)
overlap_coefficient_triple <- All_overlap / min(C, D, E)


# p-value (binomial) = (C/B)*(D/B)*(E/B)
p_binomial <- (C/B) * (D/B) * (E/B)
p_right_binomial <- pbinom(q = All_overlap - 1, size = B, prob = p_binomial, lower.tail = FALSE)
p_left_binomial <- pbinom(q = All_overlap - 1, size = B, prob = p_binomial, lower.tail = TRUE)
p_twoside_binomial   <- min(1, 2 * min(p_left_binomial, p_right_binomial))


list('Expected Overlap' = expected_triple,
     'Overlap Strength' = strength_triple,
     'Overlap Coefficient (3-Way)' = overlap_coefficient_triple,
     'Cosine C ∩ D' = cosine_CD,
     'Cosine C ∩ E' = cosine_CE,
     'Cosine D ∩ E' = cosine_DE,
     'Cosine (3-Way)' = cosine_triple,
     'Jaccard (3-Way)' = jaccard_triple, 
     'p (right tail)' = p_right_binomial,
     'p (left tail)' = p_left_binomial,
     'p (two-sided)' = p_twoside_binomial)

results_threeway <- tibble::tibble(
        Metric = c("Expected Overlap", 
                   "Overlap Strength", 
                   "Overlap Coefficient (3-Way)", 
                   "Cosine C ∩ D",
                   "Cosine C ∩ E",
                   "Cosine D ∩ E", 
                   'Cosine (3-Way)',
                   "Jaccard (3-Way)",
                   "p (right tail)", 
                   "p (left tail)", 
                   "p (two-sided)"),
        Value = c(expected_triple,
                  strength_triple, 
                  overlap_coefficient_triple,
                  cosine_CD,
                  cosine_CE,
                  cosine_DE, 
                  cosine_triple, 
                  jaccard_triple,
                  p_right_binomial,
                  p_left_binomial, 
                  p_twoside_binomial))   %>%
        gt() %>%
        fmt_number(columns = "Value", decimals = 4) %>%
        tab_header(
                title = "Female ∩ Male ∩ Sex-Agnostic") #update as needed
print(results_threeway)
