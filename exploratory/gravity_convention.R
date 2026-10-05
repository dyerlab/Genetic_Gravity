# exploratory/gravity_convention.R
#
# Project I/O and the manuscript's sign convention. Every gravity quantity is
# computed with stock gstudio (>= 1.15.0); this file only loads the simulation
# output and converts archived old-sign summaries. Higher = more of a source,
# for every quantity.
#
# Manuscript notation and the gstudio function that returns it, for an edge
# between populations i and j (gamma = bandwidth multiplier, b_i = gamma s_i / k_i):
#   w_{j|i}        neighbourhood_weights(g, gamma)[i, j]
#   Delta_{i->j}   gravity_edges(g, gamma, x)$Delta, = w_{i|j} - w_{j|i}, > 0 when
#                  i is the source; with x = deme index each edge is written with
#                  i upstream (lower index), so mean(Delta) = Delta-bar
#   Delta0_{i->j}  gravity_edges()$Delta0, = 1/k_j - 1/k_i (degree-only term)
#   g_i, g_i - 1   source_sink_scores()$gravity, $balance
#   S_i            source_sink_scores()$S: least-squares Delta_{i->j} ~ S_i - S_j,
#                  centred within each connected component
#   S0_i           source_sink_scores()$S0: -1/k_i, centred within each component
#   n_i            source_sink_scores(g, gamma, x)$n: edges with i upstream minus
#                  edges with i downstream
#   r(S, x) test   source_sink_test(g, x, gamma)
#   pGD, Delta R^2 pgd(), ibgd(g, x, mode = "pgd", gamma)
#
# Stored old-sign values: data/forward_scenarios.csv (mean_delta),
# data/information_null.csv (mean_delta), data/isotropic_baseline.csv,
# data/boundary_summary_edges.csv, data/graph_summary.csv (mean_delta) and
# data/divmigrate_timecourse_timecourse.csv (gravity_signed) hold the graph mean
# of gstudio::graph_asymmetries()'s `delta` over edges stored lo -> hi, which is
# -Delta-bar (graph_asymmetries() predates the convention: its delta is
# w_{v|u} - w_{u|v} for an edge stored (u, v)). Read them only through
# gc_read_*(), which return `dbar` and drop the old column.
#
# Check that gstudio reproduces the archived edge and score values:
#   Rscript exploratory/gravity_convention.R

suppressPackageStartupMessages({ library(igraph); library(gstudio) })

GC_SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"

#' Deme index from population names (Pop01 -> 1)
gc_x <- function(nodes) as.integer(sub("\\D+", "", nodes))

#' Deme index for every population in a graph, named by population (the `x`
#' argument of gravity_edges(), source_sink_scores(), source_sink_test(), ibgd())
gc_xg <- function(g) { nodes <- V(g)$name; setNames(gc_x(nodes), nodes) }

#' Load the Population Graph saved for one census (NULL if not saved).
#' Burn-in files are zero-padded to three digits (rep1-graph-004.rda), as written by
#' the simulation driver; %03d leaves generations >= 100 unchanged.
gc_load_graph <- function(rep, gen, scenario = NULL) {
  f <- if (is.null(scenario)) sprintf("%s/replicate%d/rep%d-graph-%03d.rda", GC_SIM_DIR, rep, rep, gen)
       else sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", GC_SIM_DIR, rep, rep, scenario, gen)
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); e$graph
}

## ---- stored (old-sign) summaries ------------------------------------------------

#' Convert a stored graph-mean of graph_asymmetries() `delta` (lo -> hi) to Dbar
gc_dbar_from_stored <- function(x) -x

#' Any stored CSV whose `col` is a graph mean of graph_asymmetries() delta (old sign):
#' returns the table with `dbar` (new sign) and without `col`. Used for
#' data/isotropic_baseline.csv, boundary_summary_edges.csv, graph_summary.csv.
gc_read_stored <- function(path, col = "mean_delta") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$dbar <- gc_dbar_from_stored(d[[col]]); d[[col]] <- NULL; d
}

gc_read_forward <- function(path = "data/forward_scenarios.csv") gc_read_stored(path)
gc_read_information_null <- function(path = "data/information_null.csv") gc_read_stored(path)

#' divMigrate time course: gravity_signed -> gravity_dbar (new sign). The
#' divMigrate column divm_signed (M[lo,hi] - M[hi,lo], row = source) is
#' already source-positive (local/divMigrate_audit.md) and is returned as is.
gc_read_divmigrate_timecourse <- function(path = "data/divmigrate_timecourse_timecourse.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$gravity_dbar <- gc_dbar_from_stored(d$gravity_signed); d$gravity_signed <- NULL; d
}

## ---- check: gstudio against the archived values ----------------------------------
if (sys.nframe() == 0L) {
  set.seed(20260930)
  cat("gstudio", as.character(packageVersion("gstudio")), "against the archived values\n\n")
  fu <- new.env(); load("data/derived/node_potential_followup.rda", envir = fu)
  nd <- readRDS("data/derived/node_potential_nodes.rds")
  # 5 Redistributed censuses at the peak block (2504-2749), gamma = 1 (the archive's bandwidth)
  cand <- unique(fu$edges[fu$edges$scenario == "Redistributed" & fu$edges$block == "2504-2749",
                          c("Replicate", "generation")])
  pick <- cand[sample(nrow(cand), 5), ]
  out <- list()
  for (r in seq_len(nrow(pick))) {
    rep <- pick$Replicate[r]; gen <- pick$generation[r]
    g <- gc_load_graph(rep, gen, 2L); x <- gc_xg(g)
    ed <- gravity_edges(g, gamma = 1, x = x); sc <- source_sink_scores(g, gamma = 1, x = x)
    old <- fu$edges[fu$edges$scenario == "Redistributed" & fu$edges$Replicate == rep &
                      fu$edges$generation == gen, ]
    ok <- match(paste(old$lo, old$hi), paste(ed$from, ed$to))
    ph <- nd[nd$scenario == "Redistributed" & nd$Replicate == rep & nd$generation == gen, ]
    stored <- gc_read_forward(); s1 <- stored[stored$Scenario == "Redistributed" & stored$replicate == rep &
                                                stored$generation == gen, ]
    out[[r]] <- data.frame(
      rep = rep, gen = gen,
      max_err_Delta_vs_archive = max(abs(ed$Delta[ok] + old$delta_lohi)),        # archived delta = -Delta
      max_err_S_vs_minus_phi = if (nrow(ph)) max(abs(sc$S[match(ph$node, sc$Stratum)] + ph$phi)) else NA,
      err_dbar_vs_csv = abs(mean(ed$Delta) - s1$dbar),
      sum_balance = sum(sc$balance),
      telescoping_err = abs(sum(ed$Delta) - sum(sc$S * sc$n) - attr(sc, "dbar_decomposition")[["residual"]]),
      frac_pos = mean(ed$Delta > 0),
      r_S_x = cor(signif(sc$S, 12), x, method = "spearman"))
  }
  chk <- do.call(rbind, out); print(chk, digits = 3, row.names = FALSE)
  stopifnot(all(chk$max_err_Delta_vs_archive < 1e-12), all(chk$err_dbar_vs_csv < 1e-12),
            all(abs(chk$sum_balance) < 1e-10), all(chk$telescoping_err < 1e-10),
            all(is.na(chk$max_err_S_vs_minus_phi) | chk$max_err_S_vs_minus_phi < 1e-10))
  cat("\nAll checks passed: gstudio's Delta, S and Delta-bar match the archived values.\n")
}
