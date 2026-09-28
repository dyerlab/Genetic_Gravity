# exploratory/node_potential_followup.R
#
# Round 2 of the node-potential analysis (local/node_potential_followup_instructions.md).
#   Task 1  direction semantics: sign of Delta on edges oriented along the imposed
#           migration (lo -> hi), node in-weight c vs position, and the
#           bandwidth mechanism described in the supplement.
#   Task 2  decomposition of the chain-oriented sum A = sum_edges o_ij Delta_ij
#           into gradient (chain ends / interior) and cycle parts through the
#           Redistributed sign reversal, with Obstructed for comparison.
# Simulated data only.
#
# Orientation: for an edge with endpoints lo < hi (deme index), migration in the
# asymmetric scenarios runs lo -> hi (m[i, i+1] = m_fwd, rows = source in
# gstudio::migrate). Delta_lohi = w_{lo->hi} - w_{hi->lo}.
#
# Run from the repository root:  Rscript exploratory/node_potential_followup.R
# Output: data/derived/node_potential_followup.rda

source("exploratory/node_potential_helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(parallel); library(mgcv) })
SEED <- 20260929L
set.seed(SEED)
OUT <- "data/derived/node_potential_followup.rda"
scen <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
cores <- max(1L, detectCores() - 1L)
blocks <- function(g) cut(g, c(2000, 2250, 2500, 2750, 3000),
                          labels = c("2004-2249", "2254-2499", "2504-2749", "2754-2999"))
t0 <- Sys.time()

## ---- Task 1: edges oriented lo -> hi ---------------------------------------------
edge_one <- function(rep, s, gen) {
  g <- load_graph(rep, gen, scen[[s]]); if (is.null(g)) return(NULL)
  ga <- graph_asymmetries(g)
  el <- as_edgelist(ga, names = TRUE); x <- deme_index(V(ga)$name); names(x) <- V(ga)$name
  k <- degree(ga); b <- setNames(V(ga)$bandwidth, V(ga)$name)
  xu <- x[el[, 1]]; xv <- x[el[, 2]]; flip <- xu > xv
  lo <- ifelse(flip, el[, 2], el[, 1]); hi <- ifelse(flip, el[, 1], el[, 2])
  w_lohi <- ifelse(flip, E(ga)$w_to, E(ga)$w_away); w_hilo <- ifelse(flip, E(ga)$w_away, E(ga)$w_to)
  data.frame(Replicate = rep, scenario = s, generation = gen, lo = lo, hi = hi,
             x_lo = x[lo], x_hi = x[hi], k_lo = k[lo], k_hi = k[hi], b_lo = b[lo], b_hi = b[hi],
             delta_lohi = w_lohi - w_hilo, row.names = NULL)
}
jobs1 <- expand.grid(rep = 1:50, s = names(scen), gen = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
edges <- do.call(rbind, mclapply(seq_len(nrow(jobs1)), function(i)
  edge_one(jobs1$rep[i], jobs1$s[i], jobs1$gen[i]), mc.cores = cores))
edges <- edges |> mutate(adjacent = x_hi - x_lo == 1,
                         interior = k_lo > 1 & k_hi > 1 & x_lo >= 3 & x_hi <= 23,
                         block = blocks(generation), scenario = factor(scenario, names(scen)))

## ---- Task 2: decomposition of A --------------------------------------------------
decomp_one <- function(rep, s, gen) {
  g <- load_graph(rep, gen, scen[[s]]); if (is.null(g)) return(NULL)
  nf <- node_fields(g); nd <- nf$nodes
  xu <- nd$x[match(nf$el[, 1], nd$node)]; xv <- nd$x[match(nf$el[, 2], nd$node)]
  o <- ifelse(xu < xv, 1, -1)
  n_i <- as.numeric(crossprod(nf$B, o))
  A <- sum(o * nf$delta); A_grad <- sum(nd$phi * n_i)
  ends <- nd$x %in% c(1, 25)
  int_e <- pmin(xu, xv) >= 3 & pmax(xu, xv) <= 23
  data.frame(Replicate = rep, scenario = s, generation = gen, edges = length(o), A = A, A_grad = A_grad,
             A_cycle = A - A_grad, A_ends = sum((nd$phi * n_i)[ends]), A_interior = sum((nd$phi * n_i)[!ends]),
             A_int_edges = sum((o * nf$delta)[int_e]))
}
jobs2 <- expand.grid(rep = 1:50, s = c("Redistributed", "Obstructed"), gen = seq(2254L, 2954L, 5L),
                     stringsAsFactors = FALSE)
decomp <- do.call(rbind, mclapply(seq_len(nrow(jobs2)), function(i)
  decomp_one(jobs2$rep[i], jobs2$s[i], jobs2$gen[i]), mc.cores = cores))

wall <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
meta <- list(seed = SEED, gstudio = as.character(packageVersion("gstudio")),
             adespatial = as.character(packageVersion("adespatial")), wall_min = wall,
             n_edge_censuses = nrow(jobs1), n_decomp_censuses = nrow(decomp))
save(edges, decomp, meta, file = OUT)

## ---- Reports -----------------------------------------------------------------------
fmt <- function(d) print(as.data.frame(d), digits = 3, row.names = FALSE)
cat(sprintf("Wall time %.1f min\n", wall))
summ <- function(d) d |> group_by(scenario, block) |>
  summarise(edges = n(), frac_neg = mean(delta_lohi < 0), med_delta = median(delta_lohi),
            rep_frac_neg = {r <- tapply(delta_lohi < 0, Replicate, mean)
              sprintf("%.2f (%.2f-%.2f)", median(r), quantile(r, .25), quantile(r, .75))}, .groups = "drop")
cat("\n== 1a: chain-adjacent interior edges, Delta_{lo,hi} (per-replicate frac_neg: median (IQR)) ==\n")
fmt(summ(filter(edges, adjacent, interior)))
cat("\n== 1b: all edges, oriented lo -> hi ==\n")
fmt(summ(edges))

cat("\n== 1c: median Spearman r(c_i, x) by scenario x block (cached node table) ==\n")
nodes <- readRDS("data/derived/node_potential_nodes.rds")
fmt(nodes |> filter(scenario %in% names(scen)) |> group_by(scenario, Replicate, generation) |>
      summarise(r = cor(c, x, method = "spearman"), .groups = "drop") |>
      mutate(block = blocks(generation), scenario = factor(scenario, names(scen))) |>
      group_by(scenario, block) |> summarise(med_r_c_x = median(r), frac_neg = mean(r < 0), .groups = "drop"))

cat("\n== 1d: bandwidth mechanism on chain-adjacent interior edges ==\n")
bd <- filter(edges, adjacent, interior)
fmt(bd |> group_by(scenario, block) |>
      summarise(frac_b_lo_gt_b_hi = mean(b_lo > b_hi),
                `b_lo>b_hi & D>0` = mean(b_lo > b_hi & delta_lohi > 0),
                `b_lo>b_hi & D<0` = mean(b_lo > b_hi & delta_lohi < 0),
                `b_lo<b_hi & D>0` = mean(b_lo < b_hi & delta_lohi > 0),
                `b_lo<b_hi & D<0` = mean(b_lo < b_hi & delta_lohi < 0),
                agree_sign = mean(sign(b_lo - b_hi) == sign(delta_lohi)), .groups = "drop"))

cat("\n== Task 2: median trajectories (every 50 generations shown) ==\n")
comps <- c("A", "A_grad", "A_cycle", "A_ends", "A_interior", "A_int_edges")
fmt(decomp |> filter(generation %% 50 == 4) |> group_by(scenario, generation) |>
      summarise(across(all_of(comps), median), .groups = "drop") |> arrange(scenario, generation))

## zero crossings of each lineage's smoothed component after the Delta-bar minimum
grid <- seq(2254, 2954, 5)
cross <- function(g, y) {
  f <- as.numeric(predict(gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML"),
                          data.frame(g = grid)))
  s <- sign(f); i <- which(s[-1] != s[-length(s)] & grid[-1] > 2400)
  if (length(i)) grid[i[1] + 1] else NA_real_
}
zc <- decomp |> group_by(scenario, Replicate) |>
  summarise(across(all_of(comps), ~ cross(generation, .x)), .groups = "drop")
cat("\n== Task 2: first zero crossing after 2400 (smoothed per lineage): n crossing, median (IQR) ==\n")
fmt(zc |> pivot_longer(all_of(comps), names_to = "component", values_to = "gen") |>
      group_by(scenario, component) |>
      summarise(n_cross = sum(!is.na(gen)),
                median_IQR = if (any(!is.na(gen))) sprintf("%.0f (%.0f-%.0f)", median(gen, na.rm = TRUE),
                  quantile(gen, .25, na.rm = TRUE), quantile(gen, .75, na.rm = TRUE)) else "-", .groups = "drop"))
cat("\n== Task 2: among lineages whose A crosses, share where the crossing coincides (+/-50 gens) ==\n")
fmt(zc |> filter(!is.na(A)) |> group_by(scenario) |>
      summarise(n = n(), with_ends = mean(!is.na(A_ends) & abs(A_ends - A) <= 50),
                with_interior = mean(!is.na(A_interior) & abs(A_interior - A) <= 50),
                with_int_edges = mean(!is.na(A_int_edges) & abs(A_int_edges - A) <= 50),
                with_cycle = mean(!is.na(A_cycle) & abs(A_cycle - A) <= 50),
                interior_no_cross = mean(is.na(A_interior)), int_edges_no_cross = mean(is.na(A_int_edges)),
                .groups = "drop"))
save(edges, decomp, zc, meta, file = OUT)
