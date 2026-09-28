# exploratory/ibd_cgd_vs_pgd.R
#
# Isolation by distance on the Population Graph topology, with symmetric (cGD)
# versus directional (pGD) edge weights. Reframes the significance question as
# the downstream use it serves: does the directional partition of the edge
# weights explain spatial structure that the symmetric graph distance does not?
#
# Landscape: 1-D stepping stone, unit spacing, x_i = deme index (Pop01 = 1).
#
# Per census:
#   cGD   undirected shortest-path distance, arc cost e_ij (symmetric).
#   pGD   directed shortest-path distance, arc i->j cost
#         p_ij = e_ij * w_{i->j} / (w_{i->j} + w_{j->i}), so p_ij + p_ji = e_ij
#         (as in R/divmigrate_truth.R).
#   X     |x_i - x_j|.
# Statistics:
#   r_cgd     Pearson r(cGD, X) over unordered pairs          (classic IBD)
#   r_pgd     Pearson r(pGD, X) over ordered pairs i != j      (directional IBD)
#   r_dir     Pearson r(pGD_ij - pGD_ji, x_j - x_i) over i < j. Zero expected
#             under symmetric migration; negative if moving "forward" (toward
#             higher x) is cheaper than moving back.
#   b_fwd, b_rev   slopes of pGD on X for forward (x_j > x_i) and reverse
#             ordered pairs, fitted jointly: pGD ~ X + X:reverse.
#   rel_fwd   b_fwd / b_rev (1 under symmetry).
#   r2_cgd        R^2 of cGD ~ X (standard IBGD).
#   r2_pgd_blind  adjusted R^2 of pGD ~ X (direction ignored).
#   r2_pgd_dir    adjusted R^2 of pGD ~ X + X:reverse (direction-aware).
#
# Samples (50 replicates each): symmetric burn-in 1904-1999 (every 5 gens, the
# false-positive reference used for the Location test) and the forward phase of
# all three treatments, every 50 generations 2004-2954.
#
# Run from the repository root:  Rscript exploratory/ibd_cgd_vs_pgd.R
# Output: data/derived/ibd_cgd_vs_pgd.rda (object `ibd`, one row per census).

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
OUT     <- "data/derived/ibd_cgd_vs_pgd.rda"

#' IBD statistics for one Population Graph
#'
#' @param graph Undirected igraph Population Graph with a `weight` attribute
#'   (cGD) and node names PopNN.
#' @return Named numeric vector of statistics (see header).
ibd_one <- function(graph) {
  ga <- graph_asymmetries(graph)
  el <- as_edgelist(ga, names = TRUE)
  e  <- E(ga)$weight; wa <- E(ga)$w_away; wt <- E(ga)$w_to
  p_fwd <- e * wa / (wa + wt); p_rev <- e * wt / (wa + wt)

  nodes <- sort(V(ga)$name)
  x <- setNames(as.integer(sub("\\D+", "", nodes)), nodes)
  C <- distances(ga, v = nodes, to = nodes, weights = e)
  dg <- graph_from_data_frame(data.frame(from = c(el[, 1], el[, 2]), to = c(el[, 2], el[, 1]),
                                         weight = c(p_fwd, p_rev)), directed = TRUE,
                              vertices = data.frame(name = nodes))
  P <- distances(dg, v = nodes, to = nodes, mode = "out")
  X <- abs(outer(x, x, "-")); S <- outer(x, x, function(a, b) b - a)   # S_ij = x_j - x_i

  up  <- upper.tri(C); off <- row(P) != col(P)
  ok_u <- up & is.finite(C) & is.finite(P) & is.finite(t(P)); ok_o <- off & is.finite(P)
  d_o <- data.frame(p = P[ok_o], X = X[ok_o], rev = as.numeric(S[ok_o] < 0))
  fit <- lm(p ~ X + X:rev, d_o)
  cf <- coef(fit)
  c(nodes = length(nodes), edges = ecount(ga), unreachable = sum(!is.finite(C[up])),
    r_cgd = cor(C[ok_u], X[ok_u]), r_pgd = cor(P[ok_o], X[ok_o]),
    r_dir = cor((P - t(P))[ok_u], S[ok_u]),
    b_fwd = unname(cf["X"]), b_rev = unname(cf["X"] + cf["X:rev"]),
    rel_fwd = unname(cf["X"] / (cf["X"] + cf["X:rev"])),
    r2_cgd = cor(C[ok_u], X[ok_u])^2,                        # standard IBGD
    r2_pgd_blind = summary(lm(p ~ X, d_o))$adj.r.squared,     # pGD, one slope
    r2_pgd_dir = summary(fit)$adj.r.squared)                  # pGD, forward/reverse slopes
}

if (sys.nframe() == 0L) {
scen <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
jobs <- rbind(
  expand.grid(Replicate = 1:50, Treatment = "Burn-in", generation = seq(1904L, 1999L, 5L),
              stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, Treatment = names(scen), generation = seq(2004L, 2954L, 50L),
              stringsAsFactors = FALSE))
jobs$file <- with(jobs, ifelse(Treatment == "Burn-in",
  sprintf("%s/replicate%d/rep%d-graph-%d.rda", SIM_DIR, Replicate, Replicate, generation),
  sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", SIM_DIR, Replicate, Replicate,
          scen[Treatment], generation)))

res <- mclapply(seq_len(nrow(jobs)), function(i) {
  if (!file.exists(jobs$file[i])) return(NULL)
  ee <- new.env(); load(jobs$file[i], envir = ee)
  out <- tryCatch(ibd_one(ee$graph), error = function(err) NULL)
  if (is.null(out)) NULL else cbind(jobs[i, 1:3], t(out))
}, mc.cores = max(1L, detectCores() - 1L))
ibd <- do.call(rbind, res)
ibd$Treatment <- factor(ibd$Treatment, c("Burn-in", names(scen)))
ibd$block <- cut(ibd$generation, c(1900, 2000, 2250, 2500, 2750, 3000),
                 labels = c("1904-1999", "2004-2249", "2254-2499", "2504-2749", "2754-2954"))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
save(ibd, file = OUT)

cat(sprintf("%d of %d censuses (missing/failed: %d)\n\n", nrow(ibd), nrow(jobs), nrow(jobs) - nrow(ibd)))
cat("Medians by treatment and generation block:\n")
print(ibd |> group_by(Treatment, block) |>
        summarise(n = n(), edges = median(edges), r_cgd = median(r_cgd), r_pgd = median(r_pgd),
                  r_dir = median(r_dir), rel_fwd = median(rel_fwd), .groups = "drop") |> as.data.frame(),
      digits = 3)

## Symmetric reference: burn-in distribution of r_dir, and how often each
## forward-phase census falls outside its central 95%.
ref <- quantile(ibd$r_dir[ibd$Treatment == "Burn-in"], c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("\nBurn-in r_dir 95%% range: %.3f to %.3f\n", ref[1], ref[2]))
cat("Share of censuses outside the burn-in 95% range of r_dir:\n")
print(ibd |> filter(Treatment != "Burn-in") |> group_by(Treatment, block) |>
        summarise(outside = mean(r_dir < ref[1] | r_dir > ref[2]), below = mean(r_dir < ref[1]),
                  .groups = "drop") |> as.data.frame(), digits = 3)
cat("\nPaired within census: share with r_pgd > r_cgd\n")
print(ibd |> group_by(Treatment) |> summarise(pgd_better = mean(r_pgd > r_cgd)) |> as.data.frame(), digits = 3)
}
