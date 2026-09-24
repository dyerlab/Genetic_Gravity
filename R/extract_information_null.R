# Builds data/information_null.csv: the isotropic locus-subsample null used by
# R/information_loss.R.
#
# Without mutation, drift and migration fix loci, and under strong directional
# flow the metapopulation approaches fixation quickly, so late-phase graphs
# are estimated from few informative loci. Fewer informative loci make edge
# weights noisier, which inflates |Delta| regardless of direction. This null
# measures that inflation directly: isotropic snapshots (no direction) at
# gens 2504 and 2999 are rebuilt from k randomly chosen metapopulation-
# polymorphic loci (k = 2, 3, 4, 5, 6, 8, 10, 12; 3 random subsets per k and
# snapshot), exactly as the original graphs were built
# (popgraph(to_mv(genotypes), Population) with default settings), and the graph-mean
# signed and absolute Delta_ij recorded alongside the diversity of the loci
# used (R/diversity_helpers.R).
#
# This script needs the raw simulation output and so only runs where the
# private research repo is available; everything downstream reads the CSV.
#
#   Rscript R/extract_information_null.R

suppressPackageStartupMessages({
  library(gstudio)
  library(igraph)
  library(parallel)
})

source("R/diversity_helpers.R")

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
KS      <- c(2, 3, 4, 5, 6, 8, 10, 12)
NSUB    <- 3
GENS    <- c(2504, 2999)

jobs <- expand.grid(replicate = 1:50, generation = GENS)

null_for_snapshot <- function(i) {
  r <- jobs$replicate[i]; g <- jobs$generation[i]
  e <- new.env()
  load(file.path(SIM_DIR, sprintf("replicate%d/rep%d-genotypes-scenario1-%d.rda", r, r, g)),
       envir = e)
  x    <- e$genotypes
  loci <- setdiff(names(x), c("ID", "Population"))
  pbar <- colMeans(allele_freqs(x, loci))
  poly <- loci[pbar > 0 & pbar < 1]

  set.seed(r * 10000 + g)
  out <- list()
  for (k in KS[KS <= length(poly)]) for (s in seq_len(NSUB)) {
    L   <- sample(poly, k)
    res <- tryCatch({
      gph <- popgraph(to_mv(x[, L, drop = FALSE]), groups = as.factor(x$Population))
      d   <- E(graph_asymmetries(gph))$delta
      c(ok = 1, n_edges = ecount(gph), mean_delta = mean(d), mean_abs_delta = mean(abs(d)))
    }, error = function(err) c(ok = 0, n_edges = NA, mean_delta = NA, mean_abs_delta = NA))
    out[[length(out) + 1]] <- data.frame(replicate = r, generation = g, k = k, subset = s,
                                         t(diversity_summary(x, L)), t(res))
  }
  do.call(rbind, out)
}

nl <- do.call(rbind, mclapply(seq_len(nrow(jobs)), null_for_snapshot,
                              mc.cores = max(1L, detectCores() - 2L)))

write.csv(nl, "data/information_null.csv", row.names = FALSE)
cat(sprintf("Wrote %d null graphs (%d failed) to data/information_null.csv\n",
            nrow(nl), sum(nl$ok == 0)))
