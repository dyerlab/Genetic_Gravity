# exploratory/gravity_convention.R
#
# The one place the manuscript's sign convention is applied. Every script that
# computes or plots a signed edge asymmetry, graph mean, chain-oriented sum or
# population score goes through these functions (local/node_potential_reanalysis.md,
# Task 1). Higher = more of a source, for every quantity here.
#
# Manuscript notation, for an edge between populations i and j:
#   w_{j|i}        j's weight in i's neighbourhood (Gaussian kernel over i's
#                  retained neighbours, bandwidth b_i); sum_j w_{j|i} = 1.
#   Delta_{i->j}   w_{i|j} - w_{j|i}; > 0 when i is the source and j the sink.
#   g_i            gravity, sum_{j in N(i)} w_{i|j}; g_i - 1 is the source-sink
#                  balance and sums to 0 over a graph without isolated nodes.
#   S_i            source-sink score: least-squares solution of
#                  Delta_{i->j} ~ S_i - S_j, S = L^+ (g - 1) with L the
#                  unweighted Laplacian; centred within each component.
#   S0_i           topological score, -1/k_i (Delta0_{i->j} = 1/k_j - 1/k_i).
#   n_i            orientation imbalance: edges on which i is upstream minus
#                  edges on which it is downstream.
#   Dbar           |E|^-1 sum Delta_{i->j}, each edge written with i upstream
#                  (lower deme index) of j.
#
# gstudio mapping [Certain, from graph_asymmetries() source]: for an edge
# stored as (u, v), w_away = w_{v|u}, w_to = w_{u|v}, and
# delta = w_away - w_to = -Delta_{u->v}. gstudio::asymmetric_weights() has
# pij = w_{j|i}, so its pij - pji is also -Delta_{i->j}.
#
# Stored old-sign values: data/forward_scenarios.csv (mean_delta),
# data/information_null.csv (mean_delta) and
# data/divmigrate_timecourse_timecourse.csv (gravity_signed) hold the mean of
# gstudio's delta over edges stored lo -> hi, i.e. -Dbar (checked on 8 random
# censuses to machine precision, 2026-09-28). Read them only through
# gc_read_*(), which returns `dbar` and drops the old column.
#
# Run the unit checks with:  Rscript exploratory/gravity_convention.R

suppressPackageStartupMessages({ library(igraph); library(gstudio) })

GC_SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"

#' Deme index from population names (Pop01 -> 1)
gc_x <- function(nodes) as.integer(sub("\\D+", "", nodes))

#' Load the Population Graph saved for one census (NULL if not saved)
gc_load_graph <- function(rep, gen, scenario = NULL) {
  f <- if (is.null(scenario)) sprintf("%s/replicate%d/rep%d-graph-%d.rda", GC_SIM_DIR, rep, rep, gen)
       else sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", GC_SIM_DIR, rep, rep, scenario, gen)
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); e$graph
}

#' Neighbourhood-weight matrix, W[i, j] = w_{j|i} (0 off-graph, NA diagonal)
#'
#' @param g Undirected Population Graph (edge attribute `weight` = e_ij).
#' @param bandwidth,scale Passed to gstudio::graph_asymmetries() (NULL = local
#'   bandwidth s_i / k_i); ignored when `perplexity` is set.
#' @param perplexity NULL, or a target perplexity for
#'   gstudio::asymmetric_weights()'s per-node sigma search.
#' @param nodes Row/column order (default: the graph's nodes, sorted).
gc_weights <- function(g, bandwidth = NULL, scale = 1, perplexity = NULL,
                       nodes = sort(V(g)$name)) {
  W <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  if (is.null(perplexity)) {
    ga <- graph_asymmetries(g, bandwidth = bandwidth, scale = scale)
    el <- as_edgelist(ga, names = TRUE)
    W[el] <- E(ga)$w_away                    # w_{v|u} at [u, v]
    W[el[, 2:1, drop = FALSE]] <- E(ga)$w_to # w_{u|v} at [v, u]
  } else {
    aw <- asymmetric_weights(g, perplexity = perplexity)
    W[cbind(aw$i, aw$j)] <- aw$pij           # w_{j|i}
    W[cbind(aw$j, aw$i)] <- aw$pji
  }
  diag(W) <- NA
  W
}

#' Edge table in the manuscript convention, one row per retained edge,
#' oriented so that i is upstream (lower deme index) of j
#'
#' @param W From gc_weights().
#' @param g The graph (for e_ij).
#' @return data.frame: i, j, x_i, x_j, e, w_ij (= w_{i|j}), w_ji (= w_{j|i}),
#'   delta (= Delta_{i->j}).
gc_edges <- function(W, g) {
  el <- as_edgelist(g, names = TRUE)
  up <- gc_x(el[, 1]) < gc_x(el[, 2])
  i <- ifelse(up, el[, 1], el[, 2]); j <- ifelse(up, el[, 2], el[, 1])
  w_ij <- W[cbind(j, i)]                     # i's weight in j's neighbourhood
  w_ji <- W[cbind(i, j)]
  data.frame(i = i, j = j, x_i = gc_x(i), x_j = gc_x(j), e = E(g)$weight,
             w_ij = w_ij, w_ji = w_ji, delta = w_ij - w_ji, stringsAsFactors = FALSE)
}

#' Population-level quantities in the manuscript convention
#'
#' @param ed From gc_edges().
#' @param nodes Node order (all populations in the graph).
#' @return list: `nodes` data.frame (node, x, k, g, balance = g - 1, S, S0, n,
#'   component), `census` named vector (edges, dbar, dbar_grad, grad_frac,
#'   resid_frac, n_comp, r_S_x), and the incidence matrix B.
gc_fields <- function(ed, nodes) {
  B <- matrix(0, nrow(ed), length(nodes), dimnames = list(NULL, nodes))
  B[cbind(seq_len(nrow(ed)), match(ed$i, nodes))] <- 1   # upstream endpoint
  B[cbind(seq_len(nrow(ed)), match(ed$j, nodes))] <- -1
  gi <- as.numeric(tapply(c(ed$w_ij, ed$w_ji), factor(c(ed$i, ed$j), nodes), sum))
  gi[is.na(gi)] <- 0
  bal <- as.numeric(crossprod(B, ed$delta))              # sum_j Delta_{i->j} = g_i - 1
  L <- crossprod(B)
  comp <- components(graph_from_edgelist(cbind(ed$i, ed$j), directed = FALSE) +
                       vertices(setdiff(nodes, c(ed$i, ed$j))))$membership[nodes]
  S <- numeric(length(nodes))
  for (cc in unique(comp)) {                             # least squares within each component
    m <- comp == cc
    if (sum(m) > 1) { s <- as.numeric(MASS::ginv(L[m, m, drop = FALSE]) %*% bal[m]); S[m] <- s - mean(s) }
  }
  k <- as.numeric(diag(L)); n <- as.numeric(colSums(B))
  fit <- as.numeric(B %*% S)
  list(nodes = data.frame(node = nodes, x = gc_x(nodes), k = k, g = gi, balance = bal,
                          S = S, S0 = ifelse(k > 0, -1 / k, NA), n = n,
                          component = as.integer(comp), stringsAsFactors = FALSE),
       census = c(edges = nrow(ed), dbar = mean(ed$delta), dbar_grad = sum(S * n) / nrow(ed),
                  grad_frac = sum(fit^2) / sum(ed$delta^2), resid_frac = 1 - sum(fit^2) / sum(ed$delta^2),
                  n_comp = length(unique(comp)),
                  r_S_x = suppressWarnings(cor(S, gc_x(nodes), method = "spearman"))),
       B = B)
}

#' Everything for one graph: weights, edges, fields
gc_census <- function(g, ...) {
  W <- gc_weights(g, ...)
  ed <- gc_edges(W, g)
  c(list(W = W, edges = ed), gc_fields(ed, rownames(W)))
}

## ---- stored (old-sign) summaries ------------------------------------------------

#' Convert a stored graph-mean of gstudio `delta` (lo -> hi) to Dbar
gc_dbar_from_stored <- function(x) -x

#' Any stored CSV whose `col` is a graph mean of gstudio's delta (old sign):
#' returns the table with `dbar` (new sign) and without `col`. Used for
#' data/isotropic_baseline.csv, boundary_summary_edges.csv, graph_summary.csv.
gc_read_stored <- function(path, col = "mean_delta") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$dbar <- gc_dbar_from_stored(d[[col]]); d[[col]] <- NULL; d
}

gc_read_forward <- function(path = "data/forward_scenarios.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$dbar <- gc_dbar_from_stored(d$mean_delta); d$mean_delta <- NULL; d
}
gc_read_information_null <- function(path = "data/information_null.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$dbar <- gc_dbar_from_stored(d$mean_delta); d$mean_delta <- NULL; d
}
#' divMigrate time course: gravity_signed -> gravity_dbar (new sign). The
#' divMigrate column divm_signed (M[lo,hi] - M[hi,lo], row = source) is
#' already source-positive (local/divMigrate_audit.md) and is returned as is.
gc_read_divmigrate_timecourse <- function(path = "data/divmigrate_timecourse_timecourse.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$gravity_dbar <- gc_dbar_from_stored(d$gravity_signed); d$gravity_signed <- NULL; d
}

## ---- unit checks ----------------------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(20260930)
  cat("gravity_convention.R unit checks\n\n")
  fu <- new.env(); load("data/derived/node_potential_followup.rda", envir = fu)
  nd <- readRDS("data/derived/node_potential_nodes.rds")
  # 5 Redistributed censuses at the peak block (2504-2749), present in round 2
  cand <- unique(fu$edges[fu$edges$scenario == "Redistributed" & fu$edges$block == "2504-2749",
                          c("Replicate", "generation")])
  pick <- cand[sample(nrow(cand), 5), ]
  out <- list()
  for (r in seq_len(nrow(pick))) {
    rep <- pick$Replicate[r]; gen <- pick$generation[r]
    cs <- gc_census(gc_load_graph(rep, gen, 2L))
    old <- fu$edges[fu$edges$scenario == "Redistributed" & fu$edges$Replicate == rep &
                      fu$edges$generation == gen, ]
    key <- paste(cs$edges$i, cs$edges$j); ok <- match(paste(old$lo, old$hi), key)
    ph <- nd[nd$scenario == "Redistributed" & nd$Replicate == rep & nd$generation == gen, ]
    out[[r]] <- data.frame(
      rep = rep, gen = gen,
      max_err_vs_round2 = max(abs(cs$edges$delta[ok] + old$delta_lohi)),     # Delta = -old
      frac_pos_adj_int = mean(cs$edges$delta[ok][old$adjacent & old$interior] > 0),
      frac_pos_all = mean(cs$edges$delta > 0),
      max_err_S_vs_minus_phi = if (nrow(ph)) max(abs(cs$nodes$S[match(ph$node, cs$nodes$node)] + ph$phi)) else NA,
      max_err_balance = max(abs(cs$nodes$balance - (cs$nodes$g - 1))),
      sum_balance = sum(cs$nodes$g - 1),
      telescoping_err = abs(sum(cs$B %*% cs$nodes$S) - sum(cs$nodes$S * cs$nodes$n)),
      r_S_x = cs$census[["r_S_x"]], dbar = cs$census[["dbar"]])
  }
  chk <- do.call(rbind, out); print(chk, digits = 3, row.names = FALSE)
  stopifnot(all(chk$max_err_vs_round2 < 1e-12),
            all(chk$max_err_balance < 1e-12), all(abs(chk$sum_balance) < 1e-10),
            all(chk$telescoping_err < 1e-10),
            all(is.na(chk$max_err_S_vs_minus_phi) | chk$max_err_S_vs_minus_phi < 1e-10))
  cat(sprintf("\nDelta_{lo->hi} > 0 on the majority of all edges: %d/5 censuses; adjacent interior: %d/5\n",
              sum(chk$frac_pos_all > 0.5), sum(chk$frac_pos_adj_int > 0.5)))
  cat(sprintf("r(S, x) < 0: %d/5 censuses\n", sum(chk$r_S_x < 0)))
  # stored CSV convention
  fw <- gc_read_forward(); s1 <- fw[fw$Scenario == "Redistributed" & fw$replicate == pick$Replicate[1] &
                                      fw$generation == pick$generation[1], ]
  cat(sprintf("CSV dbar vs recomputed (rep %d, gen %d): %.3g\n", pick$Replicate[1], pick$generation[1],
              abs(s1$dbar - chk$dbar[1])))
  cat("All identity checks passed.\n")
}
