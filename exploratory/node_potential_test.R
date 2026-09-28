# exploratory/node_potential_test.R
#
# Node-potential covariate test with a graph-spectral null
# (local/node_potential_test_instructions.md, Steps 1-5; Step 0 is
# exploratory/node_potential_step0.R).
#
# Question: is the node-level asymmetry field (phi, div, or the degree-removed
# phi_res) associated with a node covariate (here chain position x) more
# strongly than a field with the same spatial smoothness on the Population
# Graph would be by chance?
#
# Nulls (B = 999 per census; two-sided p on |r|):
#   N0       field permuted across nodes;
#   N1       Moran spectral randomization (adespatial::msr, method "singleton") on
#            the census's own Population Graph, binary adjacency (spdep listw);
#   N1w      same, adjacency weighted exp(-e^2 / 2 bbar^2), bbar = mean edge weight;
#   N1true   same, on the true 24-link chain (simulation-only diagnostic).
# Each test is also recomputed on interior nodes only (k > 1 and 3 <= x <= 23),
# with the same null draws restricted to those nodes.
#
# Samples: burn-in 1904-1999 every 5 gens (1,000 censuses); forward 2004-2954
# every 50 gens for Isotropic / Redistributed / Obstructed and the symmetric
# gradient sym-mid / sym-low / sym-verylow (scenarios 4-6), 50 replicates.
#
# Run from the repository root:  Rscript exploratory/node_potential_test.R
# Outputs: data/derived/node_potential_nodes.rds (node table, Step 1 cache);
#          data/derived/node_potential_test.rda (census table, tests).

source("exploratory/node_potential_helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(parallel) })

B_DRAWS <- 999L
SEED <- 20260928L
NODES_RDS <- "data/derived/node_potential_nodes.rds"
OUT <- "data/derived/node_potential_test.rda"
scen <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L, `sym-mid` = 4L, `sym-low` = 5L,
          `sym-verylow` = 6L)

jobs <- rbind(
  expand.grid(Replicate = 1:50, scenario = "Burn-in", generation = seq(1904L, 1999L, 5L),
              stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, scenario = names(scen), generation = seq(2004L, 2954L, 50L),
              stringsAsFactors = FALSE))

chain_adj <- function(nodes) {
  x <- deme_index(nodes); A <- 1 * (abs(outer(x, x, "-")) == 1); dimnames(A) <- list(nodes, nodes); A
}

census_one <- function(i) {
  j <- jobs[i, ]
  set.seed(SEED + i)
  g <- load_graph(j$Replicate, j$generation, if (j$scenario == "Burn-in") NULL else scen[[j$scenario]])
  if (is.null(g)) return(NULL)
  nf <- tryCatch(node_fields(g), error = function(e) NULL)
  if (is.null(nf)) return(NULL)
  nd <- nf$nodes; nodes <- nd$node
  gl <- tryCatch(node_fields(g, bandwidth = mean(E(g)$weight))$census[["grad_frac"]],
                 error = function(e) NA_real_)

  Abin <- as_adjacency_matrix(nf$ga, sparse = FALSE)[nodes, nodes]
  bbar <- mean(E(nf$ga)$weight)
  Aw <- as_adjacency_matrix(nf$ga, attr = "weight", sparse = FALSE)[nodes, nodes]
  Aw[Aw > 0] <- exp(-Aw[Aw > 0]^2 / (2 * bbar^2))
  listws <- list(N1 = as_listw(Abin), N1w = as_listw(Aw), N1true = as_listw(chain_adj(nodes)))

  ## interior-only tests restrict the statistic (observed and null) to interior nodes
  interior <- nd$k > 1 & nd$x >= 3 & nd$x <= 23
  tests <- list()
  for (fld in c("phi", "div", "phi_res")) {
    f <- nd[[fld]]
    draws <- null_draws(f, listws, B_DRAWS)
    tests[[fld]] <- rbind(cbind(field = fld, nodes = "all", field_tests(f, nd$x, draws)),
                          cbind(field = fld, nodes = "interior", field_tests(f, nd$x, draws, interior)))
  }
  tests <- do.call(rbind, tests)
  ext <- c(argmax_leaf = nd$leaf[which.max(nd$phi)], argmin_leaf = nd$leaf[which.min(nd$phi)])
  list(nodes = cbind(j, nd),
       census = cbind(j, t(nf$census), grad_frac_global = gl, has_leaf = any(nd$leaf), t(ext)),
       tests = cbind(j, tests, row.names = NULL))
}

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), census_one, mc.cores = max(1L, detectCores() - 1L),
                mc.preschedule = FALSE)
ok <- !vapply(res, is.null, logical(1))
nodes_tbl <- do.call(rbind, lapply(res[ok], `[[`, "nodes"))
census_tbl <- do.call(rbind, lapply(res[ok], `[[`, "census"))
tests_tbl <- do.call(rbind, lapply(res[ok], `[[`, "tests"))
wall <- difftime(Sys.time(), t0, units = "mins")
dir.create("data/derived", showWarnings = FALSE)
saveRDS(nodes_tbl, NODES_RDS)
meta <- list(gstudio = as.character(packageVersion("gstudio")), seed = SEED, B = B_DRAWS,
             wall_min = as.numeric(wall), censuses = sum(ok), missing = sum(!ok))
save(census_tbl, tests_tbl, meta, file = OUT)
cat(sprintf("%d censuses (%d missing), %.1f min\n", sum(ok), sum(!ok), as.numeric(wall)))
