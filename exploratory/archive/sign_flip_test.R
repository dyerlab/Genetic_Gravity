# exploratory/sign_flip_test.R
#
# Prototype of the graph-level sign-flip test for existence of asymmetry
# (local/sign_test_spec.md), before implementing it in gstudio.
#
# Null: sign(Delta_ij) is a fair coin, independently per edge; |Delta_ij| and
# the topology are held fixed. Three statistics are compared:
#   T_spec  sum(Delta^2)                   -- the spec as written. Invariant
#                                             under sign flips, so p = 1 always.
#   S_axis  |sum(o_ij * Delta_ij)|         -- Delta oriented along the chain
#                                             (o = +1 when the edge runs toward
#                                             the higher-numbered deme). Needs a
#                                             known axis.
#   G_grad  ||P_grad Delta||^2             -- energy of the node-potential
#                                             (gradient) part of Delta, the
#                                             least-squares fit Delta_ij ~
#                                             phi_j - phi_i. Axis-free.
#   G_res   G_grad on Delta - (1/k_i - 1/k_j), the uniform-weight degree
#           baseline removed.
#
# (a) False-positive rate: the Location-test calibration sample, burn-in
#     generations 1904-1999, 50 replicates (Delta taken from the saved window
#     run; Delta does not depend on the permutation).
# (b) Chain-end split of (a).
# (c) Power: forward-phase censuses of all three treatments, every 50
#     generations from 2004, 50 replicates, Delta from graph_asymmetries().
#
# Run from the repository root:  Rscript exploratory/sign_flip_test.R

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
OUT_DIR <- "data/derived"
B <- 9999L
set.seed(20260927)

#' Sign-flip p-values for one census
#'
#' @param from,to Character; edge endpoints, in the orientation Delta refers to.
#' @param delta Numeric; Delta_ij = w_{from->to} - w_{to->from}.
#' @return Named numeric vector of p-values and observed statistics.
sign_flip_one <- function(from, to, delta) {
  m <- length(delta)
  nodes <- sort(unique(c(from, to)))
  k <- table(factor(c(from, to), nodes))
  inc <- matrix(0, m, length(nodes), dimnames = list(NULL, nodes))
  inc[cbind(seq_len(m), match(from, nodes))] <- -1
  inc[cbind(seq_len(m), match(to, nodes))]   <-  1
  P <- inc %*% MASS::ginv(crossprod(inc)) %*% t(inc)      # projection onto gradients
  o <- sign(as.integer(sub("\\D+", "", to)) - as.integer(sub("\\D+", "", from)))
  base <- as.numeric(1 / k[from] - 1 / k[to])
  res  <- delta - base

  S <- matrix(sample(c(-1, 1), m * B, replace = TRUE), m, B)
  pv <- function(obs, null) (1 + sum(null >= obs - 1e-12)) / (B + 1)
  grad <- function(x) colSums((P %*% x)^2)

  t_obs <- sum(delta^2);            t_null <- colSums((S * abs(delta))^2)
  a_obs <- abs(sum(o * delta));     a_null <- abs(colSums(o * S * abs(delta)))
  g_obs <- grad(matrix(delta));     g_null <- grad(S * abs(delta))
  r_obs <- grad(matrix(res));       r_null <- grad(S * abs(res))
  c(edges = m, p_spec = pv(t_obs, t_null), p_axis = pv(a_obs, a_null),
    p_grad = pv(g_obs, g_null), p_res = pv(r_obs, r_null),
    axis_signed = sum(o * delta), grad_frac = g_obs / t_obs, res_grad_frac = r_obs / sum(res^2))
}

run_censuses <- function(df, keys) {
  groups <- split(df, df[keys], drop = TRUE)
  out <- mclapply(groups, function(d) {
    cbind(d[1, keys, drop = FALSE], t(sign_flip_one(d$from, d$to, d$delta)))
  }, mc.cores = max(1L, detectCores() - 1L), mc.set.seed = TRUE)
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}

## ---- (a, b) calibration sample ----------------------------------------------
e <- new.env()
load(file.path(SIM_DIR, "divMigrate/firstsig/gravity.location_window.burnin_1900_2000.rda"), envir = e)
win <- e$window_edges |> filter(is.finite(delta))
cal <- run_censuses(win, c("Replicate", "generation"))

## chain-end split: rerun each census on chain-end and interior edges separately
pos <- function(x) as.integer(sub("\\D+", "", x))
win$chain_end <- pmin(pos(win$from), pos(win$to)) <= 2 | pmax(pos(win$from), pos(win$to)) >= 24
cal_end <- run_censuses(win, c("Replicate", "generation", "chain_end"))

## ---- (c) forward phase --------------------------------------------------------
scen <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
jobs <- expand.grid(Replicate = 1:50, scenario = names(scen), generation = seq(2004L, 2954L, 50L),
                    stringsAsFactors = FALSE)
fwd_edges <- mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]
  f <- sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", SIM_DIR, j$Replicate, j$Replicate,
               scen[[j$scenario]], j$generation)
  if (!file.exists(f)) return(NULL)
  ee <- new.env(); load(f, envir = ee)
  ga <- graph_asymmetries(ee$graph)
  el <- as_edgelist(ga, names = TRUE)
  data.frame(j, from = el[, 1], to = el[, 2], delta = E(ga)$delta)
}, mc.cores = max(1L, detectCores() - 1L))
fwd_edges <- do.call(rbind, fwd_edges)
fwd <- run_censuses(fwd_edges, c("Replicate", "scenario", "generation"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save(cal, cal_end, fwd, file = file.path(OUT_DIR, "sign_flip_test.rda"))

## ---- report --------------------------------------------------------------------
alphas <- c(0.01, 0.05, 0.10)
fpr <- function(d) sapply(c("p_spec", "p_axis", "p_grad", "p_res"),
                          function(v) sapply(alphas, function(a) mean(d[[v]] <= a)))
cat(sprintf("(a) FPR, %d burn-in censuses (rows = alpha 0.01, 0.05, 0.10)\n", nrow(cal)))
print(round(`rownames<-`(fpr(cal), alphas), 3))
cat("\np-value histograms (uniform = 0.10):\n")
br <- seq(0, 1, 0.1)
print(round(t(sapply(c("p_axis", "p_grad", "p_res"), function(v)
  table(cut(cal[[v]], br, include.lowest = TRUE)) / nrow(cal))), 3))
cat("\nGradient share of sum(Delta^2), burn-in: median", round(median(cal$grad_frac), 3),
    "| degree-removed:", round(median(cal$res_grad_frac), 3),
    "| random-sign expectation (n-1)/m:", round(median(24 / cal$edges), 3), "\n")

cat("\n(b) FPR at 0.05 by chain end:\n")
for (ce in c(FALSE, TRUE)) {
  d <- cal_end[cal_end$chain_end == ce, ]
  cat(sprintf("  chain_end = %-5s censuses %4d  median edges %3.0f  axis %.3f  grad %.3f  res %.3f\n",
              ce, nrow(d), median(d$edges), mean(d$p_axis <= .05), mean(d$p_grad <= .05), mean(d$p_res <= .05)))
}

cat("\n(c) Rejection rate at 0.05, forward phase, by treatment and generation block:\n")
fwd$block <- cut(fwd$generation, c(2000, 2250, 2500, 2750, 3000),
                 labels = c("2004-2249", "2254-2499", "2504-2749", "2754-2999"))
print(fwd |> group_by(scenario, block) |>
        summarise(n = n(), axis = mean(p_axis <= .05), grad = mean(p_grad <= .05), res = mean(p_res <= .05),
                  axis_signed = median(axis_signed), .groups = "drop") |> as.data.frame(), digits = 3)
