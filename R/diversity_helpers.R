# Shared helpers for the extraction scripts (R/extract_forward_scenarios.R,
# R/extract_information_null.R). The simulated loci are biallelic with
# genotypes coded "aa:bb" (alleles "01"/"02"), and every deme holds the same
# number of individuals, so the metapopulation frequency is the unweighted
# mean of the deme frequencies.

# Frequency of allele "01" for each deme (rows) and locus (columns).
allele_freqs <- function(genotypes, loci) {
  pops <- genotypes$Population
  sapply(loci, function(l) {
    s   <- as.character(genotypes[[l]])
    n01 <- (substr(s, 1, 2) == "01") + (substr(s, 4, 5) == "01")
    tapply(n01, pops, sum) / (2 * tapply(n01, pops, length))
  })
}

# Diversity summaries over the given loci:
#   n_poly_global  loci still polymorphic in the metapopulation
#   n_poly_pop     mean number of loci polymorphic within a deme
#   He             mean within-deme expected heterozygosity
#   sumAe_total    sum over loci of the metapopulation effective number of
#                  alleles, A_e = 1 / sum(p_i^2); 1 per fixed locus
#   sumAe_within   the same sum computed within each deme, averaged over demes
# sumAe_total - sumAe_within is the among-deme component of diversity.
diversity_summary <- function(genotypes, loci) {
  p    <- allele_freqs(genotypes, loci)
  pbar <- colMeans(p)
  ae   <- function(q) 1 / (q^2 + (1 - q)^2)
  c(n_poly_global = sum(pbar > 0 & pbar < 1),
    n_poly_pop    = mean(rowSums(p > 0 & p < 1)),
    He            = mean(2 * p * (1 - p)),
    sumAe_total   = sum(ae(pbar)),
    sumAe_within  = mean(rowSums(ae(p))))
}
