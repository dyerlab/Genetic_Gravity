# exploratory/node_potential_helpers.R
#
# Helpers for exploratory/node_potential_test.R (see
# local/node_potential_test_instructions.md for the definitions used here).
#
# Orientation: an edge stored by graph_asymmetries() as (u, v) has
# delta = w_{u->v} - w_{v->u}. The incidence matrix B has +1 at u and -1 at v,
# so the gradient model is Delta = B phi, i.e. Delta_uv ~ phi_u - phi_v
# (high phi = source), div = t(B) Delta, and phi = L^+ div with L = t(B) B.

suppressPackageStartupMessages({ library(igraph); library(gstudio) })

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"

#' Load the Population Graph saved for one census
load_graph <- function(rep, gen, scenario = NULL) {
  f <- if (is.null(scenario)) sprintf("%s/replicate%d/rep%d-graph-%d.rda", SIM_DIR, rep, rep, gen)
       else sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", SIM_DIR, rep, rep, scenario, gen)
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if (!exists("graph", envir = e)) stop("no object named 'graph' in ", f)
  e$graph
}

deme_index <- function(nodes) as.integer(sub("\\D+", "", nodes))

#' Edge-node incidence matrix (+1 at the stored first endpoint, -1 at the second)
incidence <- function(el, nodes) {
  B <- matrix(0, nrow(el), length(nodes), dimnames = list(NULL, nodes))
  B[cbind(seq_len(nrow(el)), match(el[, 1], nodes))] <- 1
  B[cbind(seq_len(nrow(el)), match(el[, 2], nodes))] <- -1
  B
}

#' Node-level fields for one census graph
#'
#' @param g Undirected Population Graph (edge attribute `weight` = cGD).
#' @param bandwidth Passed to graph_asymmetries() (NULL = local bandwidth).
#' @param gamma Bandwidth multiplier (graph_asymmetries(scale = gamma)); 1 = local mean.
#' @return List: `nodes` (data.frame, one row per node) and `census` (named
#'   numeric vector), plus the edge list, delta, B and L for reuse.
node_fields <- function(g, bandwidth = NULL, gamma = 1) {
  ga <- graph_asymmetries(g, bandwidth = bandwidth, scale = gamma)
  nodes <- V(ga)$name
  el <- as_edgelist(ga, names = TRUE)
  delta <- E(ga)$delta
  B <- incidence(el, nodes); L <- crossprod(B); Lp <- MASS::ginv(L)
  k <- degree(ga)[nodes]
  u <- match(el[, 1], nodes); v <- match(el[, 2], nodes)
  # sum_j w_{i->j} (edge (u,v): u sends w_away, v sends w_to) and in-weight c_i
  rowsum_w <- as.numeric(tapply(c(E(ga)$w_away, E(ga)$w_to), factor(c(u, v), seq_along(nodes)), sum))
  c_in     <- as.numeric(tapply(c(E(ga)$w_to, E(ga)$w_away), factor(c(u, v), seq_along(nodes)), sum))
  div <- as.numeric(crossprod(B, delta))
  phi <- as.numeric(Lp %*% div)
  delta0 <- 1 / k[u] - 1 / k[v]
  phi0 <- as.numeric(Lp %*% crossprod(B, delta0))
  phi_res <- if (sd(phi0) > 0) as.numeric(resid(lm(phi ~ phi0))) else phi - mean(phi)
  grad <- sum((B %*% phi)^2) / sum(delta^2)
  comp <- components(ga)$membership[nodes]
  list(nodes = data.frame(node = nodes, x = deme_index(nodes), k = as.integer(k),
                          strength = as.numeric(strength(ga)[nodes]), b = V(ga)$bandwidth,
                          div = div, c = c_in, rowsum = rowsum_w, phi = phi, phi0 = phi0,
                          phi_res = phi_res, leaf = k == 1, component = as.integer(comp)),
       census = c(n_nodes = length(nodes), edges = nrow(el), n_comp = max(comp),
                  grad_frac = grad, cycle_frac = 1 - grad,
                  cor_phi_phi0 = if (sd(phi0) > 0) cor(phi, phi0) else NA),
       el = el, delta = delta, B = B, L = L, ga = ga)
}

#' Spatial-weights list from a symmetric weights matrix (spdep; style "B"
#' keeps the supplied weights unstandardized)
as_listw <- function(A) spdep::mat2listw(A, style = "B", zero.policy = TRUE)

#' Null draws of a node field
#'
#' N0 permutes the field across nodes; every other entry is Moran spectral
#' randomization (adespatial::msr, default method "pair", which preserves the
#' field's Moran's I on the supplied weights).
#'
#' @param f Numeric field, one value per node.
#' @param listws Named list of listw objects (same node order as f).
#' @param B Number of null draws.
#' @return Named list of n x B matrices.
null_draws <- function(f, listws, B = 999L) {
  out <- list(N0 = replicate(B, sample(f)))
  for (nm in names(listws))
    out[[nm]] <- as.matrix(adespatial::msr(f, listws[[nm]], nrepet = B, method = "singleton", simplify = TRUE))
  out
}

#' Two-sided p-values for r(field, x), optionally on a subset of nodes
#'
#' @param f Numeric field; x covariate; draws from null_draws().
#' @param keep Logical; nodes entering the correlation.
#' @return data.frame, one row per null x correlation type: r_obs, p.
field_tests <- function(f, x, draws, keep = rep(TRUE, length(f))) {
  B <- ncol(draws[[1]]); out <- list()
  f <- signif(f, 12)     # tie rounding (gravity_convention.R gc_rank_S); pass the same rounded f to null_draws()
  for (m in c("spearman", "pearson")) {
    r_obs <- cor(f[keep], x[keep], method = m)
    for (nm in names(draws)) {
      r_null <- as.numeric(cor(draws[[nm]][keep, , drop = FALSE], x[keep], method = m))
      out[[length(out) + 1]] <- data.frame(null = nm, cor = m, r_obs = r_obs,
                                           p = (1 + sum(abs(r_null) >= abs(r_obs) - 1e-12)) / (B + 1))
    }
  }
  do.call(rbind, out)
}
