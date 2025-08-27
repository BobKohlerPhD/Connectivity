#########
# Given #
#########

#~~~~Can find these values from output of EdgeOverlap_Plotting.R~~~~#
B <- 35778   # all possible edges
C <- 549     # number of edges in set 1
D <- 2850    # number of edges in  2
A <- 14      # observed overlap in edges between set 1 and set 2

############
# Expected #
############
expected <- C * D / B
enrichment <- A / expected
log2_enrichment <- log2(enrichment) # for clarity 

#########################################
### One-sided hypergeometric p-values ###
#########################################
# Enrichment P(X >= A) -- more overlap than expected 
p_enrich <- phyper(q = A - 1, m = D, n = B - D, k = C, lower.tail = FALSE)

# Depletion P(X <= A) -- less overlap than expected
p_deplete <- phyper(q = A,     m = D, n = B - D, k = C, lower.tail = TRUE)

######################################
### two-sided w/ Fisher exact test ###
######################################
tab <- matrix(c(A, C - A, D - A, B - C - D + A), nrow = 2, byrow = TRUE)
p_two_sided <- fisher.test(tab, alternative = "two.sided")$p.value

###############
### Results ###
###############
list(
  expected = expected,
  enrichment = enrichment, # for looking at more than expected overlap 
  log2_enrichment = log2_enrichment, # easier to describe after logtransform
  p_enrich_right_tail = p_enrich, # if significant, then there are more overlaps than expected
  p_deplete_left_tail = p_deplete, # if significant, then there are fewer overlaps than expected
  p_two_sided = p_two_sided # if significant, then there are significantly different number of overlaps than expected 
)
