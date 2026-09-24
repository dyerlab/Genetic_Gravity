# Builds data/forward_scenarios.csv from the stored simulation output.
#
# For every N=100 forward-phase snapshot of the three base-study scenarios
# (1 = Isotropic, 2 = Redistributed, 3 = Obstructed; 50 replicates x 200
# snapshots, gens 2004-2999) the stored Population Graph is loaded, directional
# weights are recomputed with gstudio::graph_asymmetries(), and the matching
# genotype snapshot is summarized. Columns:
#   graph      n_edges, CD (C_D-bar = 2|E| / (K(K-1))), diameter (weighted by
#              conditional genetic distance), mean_delta and mean_abs_delta
#              (graph means of signed and absolute Delta_ij), sd_delta (SD of
#              Delta_ij across edges), n_pendant (degree-one nodes)
#   diversity  n_poly_global, n_poly_pop, He, sumAe_total, sumAe_within
#              (see R/diversity_helpers.R); sumAe_total is the lineage clock
#              used in R/lineage_clock.R
#
# The graph summaries are recomputed from the graph files rather than read from
# the private repo's data/graph_summary.rda because that summary silently omits
# 31 late Redistributed snapshots (reps 2, 4, 18, 21, 22, 27, 34, 44, 45;
# gens 2729-2999) -- the most asymmetric graphs in the study. Values for every
# snapshot it does contain agree with the recomputation to machine precision.
#
# This script needs the raw simulation output and so only runs where the
# private research repo is available; everything downstream reads the CSV.
#
#   Rscript R/extract_forward_scenarios.R

suppressPackageStartupMessages({
  library(gstudio)
  library(igraph)
  library(parallel)
})

source("R/diversity_helpers.R")

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
K       <- 25
REPS    <- 1:50
SCN     <- c("1" = "Isotropic", "2" = "Redistributed", "3" = "Obstructed")

files <- unlist(lapply(REPS, function(r)
  list.files(file.path(SIM_DIR, sprintf("replicate%d", r)),
             pattern = sprintf("^rep%d-graph-scenario[123]-[0-9]+[.]rda$", r),
             full.names = TRUE)))

summarise_snapshot <- function(path) {
  e <- new.env()
  load(path, envir = e)                                   # 'graph'
  load(sub("-graph-", "-genotypes-", path), envir = e)    # 'genotypes'
  g <- graph_asymmetries(e$graph)
  d <- E(g)$delta
  b <- basename(path)
  cbind(
    data.frame(
      Scenario       = SCN[sub(".*-scenario([0-9]+)-.*", "\\1", b)],
      replicate      = as.integer(sub("rep([0-9]+)-.*", "\\1", b)),
      generation     = as.integer(sub(".*-([0-9]+)[.]rda$", "\\1", b)),
      n_edges        = ecount(g),
      CD             = 2 * ecount(g) / (K * (K - 1)),
      diameter       = diameter(e$graph),
      mean_delta     = mean(d),
      mean_abs_delta = mean(abs(d)),
      sd_delta       = sd(d),
      n_pendant      = sum(degree(g) == 1),
      stringsAsFactors = FALSE),
    t(diversity_summary(e$genotypes, setdiff(names(e$genotypes), c("ID", "Population")))))
}

out <- do.call(rbind, mclapply(files, summarise_snapshot,
                               mc.cores = max(1L, detectCores() - 2L)))
out <- out[order(out$Scenario, out$replicate, out$generation), ]
rownames(out) <- NULL

write.csv(out, "data/forward_scenarios.csv", row.names = FALSE)
cat(sprintf("Wrote %d snapshots to data/forward_scenarios.csv\n", nrow(out)))
