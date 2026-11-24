# R/simulation/half_bandwidth_compute.R
#
# Census-level inputs for the main Results at a given bandwidth multiplier
# (local/half_bandwidth_rerun_instructions.md). Run once at GAMMA = 1 and once at
# GAMMA = 0.5 (the manuscript's bandwidth), with identical code, graphs, phases,
# null adjacency and seeds, so the two differ only by bandwidth. Every gravity
# quantity comes from stock gstudio (>= 1.15.0) with b_i = GAMMA * s_i / k_i; see
# R/common/gravity_convention.R for the notation.
#
# Parts (each cached; resumable):
#   grid   every saved census on the 5-generation grid: burn-in 0-1999 and the six
#          forward scenarios 2004-2999 (80,050 censuses). Per census: edge count,
#          Delta-bar, mean |Delta|, SD; cor(S, S0) (the topological-floor statistic);
#          leaf statistics; IBGD fits (ibgd(mode = "dcGD"); burn-in from 1904 only);
#          truth comparison of w for Symmetric / Redistributed / Obstructed; the
#          decomposition of |E| Delta-bar; and node S and degree on the
#          Fig_SourceSinkScore grid (burn-in 1904-1999 every 5, forward every 10).
#   tests  the source-sink test, source_sink_test(), at 7,000 censuses (burn-in
#          1904-1999 every 5, forward every 50): the manuscript's null (Moran
#          spectral randomization on the graph, null = "graph"), the same on interior
#          nodes, and three comparison nulls (plain permutation of S; spectral
#          randomization on the unweighted graph; on the true stepping-stone chain).
#          One seed per census (20260928 + job index), set before the first call.
#   bw     (GAMMA = 0.5 only) Tab_BandwidthSensitivity and
#          Tab_BandwidthSensitivitySourceSink relative to GAMMA: global bandwidth,
#          local x 2 (= gamma 1), local x 1/2 (= gamma 0.25), perplexity 4 (seed
#          20260930 + i, reference first, every row on fresh draws).
#   mirror (GAMMA = 0.5 only) the mirroring check: the 50 generation-1999 burn-in
#          graphs rebuilt from genotypes with deme labels reversed; r(S, x) compared.
#
# Run from the repository root:
#   GAMMA=0.5 PART=grid Rscript R/simulation/half_bandwidth_compute.R   (PART = grid | tests | bw | mirror)
# Output: data/derived/half_bandwidth_<part>_g<tag>.rda (tag 1 or 05)

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })
source("R/common/gravity_convention.R")
GAMMA <- as.numeric(Sys.getenv("GAMMA", "1")); TAG <- gsub("\\.", "", format(GAMMA))
PART  <- Sys.getenv("PART", "grid")
CORES <- max(1L, detectCores() - 1L)
SCN   <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-mid` = 4L, `sym-low` = 5L, `sym-verylow` = 6L)
RATES <- list(Symmetric = c(0.025, 0.025), Redistributed = c(0.040, 0.010), Obstructed = c(0.025, 0.010))
POPS  <- sprintf("Pop%02d", 1:25); OFF <- row(diag(25)) != col(diag(25))
out_file <- function(part) sprintf("data/derived/half_bandwidth_%s_g%s.rda", part, TAG)
load_g <- function(rep, smp, gen) gc_load_graph(rep, gen, if (smp == "Burn-in") NULL else SCN[[smp]])

## ---- grid ---------------------------------------------------------------------------------
#' Neighbourhood weights against the true migration matrix (Spearman rho and RMSE, both
#' rescaled to a maximum of 1), and how often w favours the true forward direction on
#' retained stepping-stone links
truth_w <- function(W, smp) {
  Wf <- matrix(0, 25, 25, dimnames = list(POPS, POPS)); Wf[rownames(W), colnames(W)] <- W; diag(Wf) <- NA
  r <- RATES[[smp]]; T0 <- matrix(0, 25, 25, dimnames = list(POPS, POPS)); i <- 1:24
  T0[cbind(i, i + 1L)] <- r[1]; T0[cbind(i + 1L, i)] <- r[2]; diag(T0) <- NA
  Tr <- T0 / max(T0[OFF]); TrT <- t(Tr); Ws <- Wf / max(Wf[OFF])
  fw <- Wf[cbind(i + 1, i)]; rv <- Wf[cbind(i, i + 1)]; ret <- fw > 0 & rv > 0
  c(w_rho_run = suppressWarnings(cor(Ws[OFF], Tr[OFF], method = "spearman")), w_rmse_run = sqrt(mean((Ws[OFF] - Tr[OFF])^2)),
    w_rho_tr = suppressWarnings(cor(Ws[OFF], TrT[OFF], method = "spearman")), w_rmse_tr = sqrt(mean((Ws[OFF] - TrT[OFF])^2)),
    links = sum(ret), w_fwd_gt = mean(fw[ret] > rv[ret]))
}
#' Adjusted R^2 of pGD on separation with one slope (direction ignored), beside
#' ibgd(mode = "dcGD")'s forward/reverse fit
pgd_blind_r2 <- function(g, x) {
  P <- dcgd(g, GAMMA, output = "matrix"); X <- abs(outer(x[rownames(P)], x[colnames(P)], "-"))
  ok <- row(P) != col(P) & is.finite(P)
  summary(lm(P[ok] ~ X[ok]))$adj.r.squared
}
grid_census <- function(rep, smp, gen) {
  g <- load_g(rep, smp, gen); if (is.null(g)) return(NULL)
  x <- gc_xg(g)
  ed <- gravity_edges(g, GAMMA, x)                        # each edge written upstream (from) -> downstream (to)
  sc <- source_sink_scores(g, GAMMA, x); sc <- sc[order(gc_x(sc$Stratum)), ]
  xs <- gc_x(sc$Stratum); xf <- x[ed$from]; xt <- x[ed$to]
  ends <- xs %in% c(1, 25); int_e <- pmin(xf, xt) >= 3 & pmax(xf, xt) <= 23
  Sn <- sc$S * sc$n; leaf <- sc$degree == 1
  v <- c(n_edges = nrow(ed), dbar = mean(ed$Delta), mean_abs_delta = mean(abs(ed$Delta)), sd_delta = sd(ed$Delta),
         n_pendant = sum(leaf), n_comp = max(sc$component), grad_frac = attr(sc, "gradient_share"),
         cor_S_S0f = attr(sc, "cor_S_S0"), has_leaf = any(leaf),
         src_leaf = leaf[which.max(sc$S)], sink_leaf = leaf[which.min(sc$S)],
         A = sum(ed$Delta), A_grad = sum(Sn), A_cycle = sum(ed$Delta) - sum(Sn), A_ends = sum(Sn[ends]),
         A_interior = sum(Sn[!ends]), A_int_edges = sum(ed$Delta[int_e]),
         r2_cgd = NA, r2_pgd_blind = NA, r2_pgd_dir = NA, rel_fwd = NA,
         w_rho_run = NA, w_rmse_run = NA, w_rho_tr = NA, w_rmse_tr = NA, links = NA, w_fwd_gt = NA)
  if (gen >= 1904) { ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = GAMMA), error = function(e) NULL)
    if (!is.null(ib)) v[c("r2_cgd", "r2_pgd_dir", "rel_fwd", "r2_pgd_blind")] <-
      c(ib$estimate[c("R2 cGD", "adj R2 dcGD", "slope ratio")], pgd_blind_r2(g, x)) }
  if (smp %in% names(RATES)) { tw <- truth_w(neighbourhood_weights(g, GAMMA), smp); v[names(tw)] <- tw }
  heat <- (smp == "Burn-in" && gen >= 1904) || (smp %in% names(RATES) && (gen - 2004) %% 10 == 0)
  list(row = data.frame(Replicate = rep, sample = smp, generation = gen, t(v)),
       nodes = if (heat) data.frame(Replicate = rep, scenario = smp, generation = gen, deme = xs, S = sc$S, k = sc$degree) else NULL)
}
if (PART == "grid" && !file.exists(out_file("grid"))) {
  cache <- sprintf("data/derived/half_bandwidth_cache_g%s/grid", TAG); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  units <- rbind(expand.grid(Replicate = 1:50, sample = "Burn-in", stringsAsFactors = FALSE),
                 expand.grid(Replicate = 1:50, sample = names(SCN), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  invisible(mclapply(seq_len(nrow(units)), function(u) {
    f <- file.path(cache, sprintf("%s_%02d.rds", gsub("[^A-Za-z]", "", units$sample[u]), units$Replicate[u])); if (file.exists(f)) return(NULL)
    gens <- if (units$sample[u] == "Burn-in") c(0L, seq(4L, 1999L, 5L)) else seq(2004L, 2999L, 5L)
    res <- lapply(gens, function(gg) tryCatch(grid_census(units$Replicate[u], units$sample[u], gg), error = function(e)
      list(row = data.frame(Replicate = units$Replicate[u], sample = units$sample[u], generation = gg, err = conditionMessage(e)), nodes = NULL)))
    saveRDS(list(rows = bind_rows(lapply(res, `[[`, "row")), nodes = bind_rows(lapply(res, `[[`, "nodes"))), f); NULL
  }, mc.cores = CORES, mc.preschedule = FALSE))
  fs <- list.files(cache, full.names = TRUE)
  if (length(fs) == nrow(units)) {
    parts <- lapply(fs, readRDS); grid <- bind_rows(lapply(parts, `[[`, "rows")); heat <- bind_rows(lapply(parts, `[[`, "nodes"))
    meta <- list(gamma = GAMMA, date = as.character(Sys.Date()), gstudio = as.character(packageVersion("gstudio")),
                 wall_min_last = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    save(grid, heat, meta, file = out_file("grid")); cat(sprintf("Wrote %s: %d censuses\n", out_file("grid"), nrow(grid)))
  } else cat(sprintf("%d of %d units cached; rerun to continue\n", length(fs), nrow(units)))
}

## ---- tests ---------------------------------------------------------------------------------
r1_jobs <- rbind(
  expand.grid(Replicate = 1:50, scenario = "Burn-in", generation = seq(1904L, 1999L, 5L), stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, scenario = c("Isotropic", names(SCN)[-1]), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE))
#' Adjacency of the true stepping-stone chain (neighbouring demes)
chain_adj <- function(nodes) { x <- gc_x(nodes); A <- 1 * (abs(outer(x, x, "-")) == 1); dimnames(A) <- list(nodes, nodes); A }
test_census <- function(i) {
  j <- r1_jobs[i, ]; set.seed(20260928L + i)
  g <- gc_load_graph(j$Replicate, j$generation, if (j$scenario == "Burn-in") NULL else c(Isotropic = 1L, SCN[-1])[[j$scenario]])
  if (is.null(g)) return(NULL)
  x <- gc_xg(g)
  t_w    <- source_sink_test(g, x, GAMMA)                                           # the manuscript's test
  t_int  <- source_sink_test(g, x, GAMMA, interior = TRUE)                          # demes 3-23, degree > 1
  t_perm <- source_sink_test(g, x, GAMMA, null = "permute")
  t_bin  <- source_sink_test(g, x, GAMMA, null = "adjacency", adjacency = as_adjacency_matrix(g, sparse = FALSE))
  t_true <- source_sink_test(g, x, GAMMA, null = "adjacency", adjacency = chain_adj(names(x)))
  sc <- t_w$scores; leaf <- sc$degree == 1
  data.frame(Replicate = j$Replicate, scenario = j$scenario, generation = j$generation, r_S = unname(t_w$statistic),
             p_N0 = t_perm$p.value, p_N1 = t_bin$p.value, p = t_w$p.value, p_N1true = t_true$p.value,
             r_S_int = unname(t_int$statistic), p_int = t_int$p.value,
             n_comp = max(sc$component), cor_phi_phi0 = attr(sc, "cor_S_S0"), has_leaf = any(leaf),
             src_leaf = leaf[which.max(sc$S)], sink_leaf = leaf[which.min(sc$S)])
}
if (PART == "tests" && !file.exists(out_file("tests"))) {
  t0 <- Sys.time()
  tests <- bind_rows(mclapply(seq_len(nrow(r1_jobs)), function(i) tryCatch(test_census(i), error = function(e)
    data.frame(r1_jobs[i, ], err = conditionMessage(e))), mc.cores = CORES, mc.preschedule = FALSE))
  meta <- list(gamma = GAMMA, date = as.character(Sys.Date()), gstudio = as.character(packageVersion("gstudio")),
               wall_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  save(tests, meta, file = out_file("tests")); cat(sprintf("Wrote %s: %d censuses\n", out_file("tests"), nrow(tests)))
}

## ---- bandwidth tables relative to GAMMA ---------------------------------------------------------
# Global bandwidth and perplexity are not gamma kernels, so their weights come from
# graph_asymmetries(bandwidth =) and asymmetric_weights(), and their scores from the same
# least squares as source_sink_scores() (checked against it at the reference bandwidth).
#' W[i, j] = w_{j|i} for one alternative bandwidth
W_alt <- function(g, a) {
  if (!is.null(a$gamma)) return(neighbourhood_weights(g, a$gamma))
  nodes <- V(g)$name; W <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  if (isTRUE(a$global)) {
    ga <- graph_asymmetries(g, bandwidth = mean(E(g)$weight)); el <- as_edgelist(ga, names = TRUE)
    W[el] <- E(ga)$w_away; W[el[, 2:1, drop = FALSE]] <- E(ga)$w_to          # w_{v|u} at [u, v]; w_{u|v} at [v, u]
  } else {
    aw <- asymmetric_weights(g, perplexity = a$perplexity)
    W[cbind(aw$from, aw$to)] <- aw$pij; W[cbind(aw$to, aw$from)] <- aw$pji   # pij = w_{j|i}
  }
  diag(W) <- NA; W
}
#' Delta_{from->to} on gravity_edges()'s oriented edges, from a weight matrix
delta_W <- function(W, ed) W[cbind(ed$to, ed$from)] - W[cbind(ed$from, ed$to)]
#' Source-sink score S from Delta on oriented edges (least squares, centred within components)
scores_from_delta <- function(g, ed, delta) {
  nodes <- V(g)$name; B <- matrix(0, nrow(ed), length(nodes), dimnames = list(NULL, nodes))
  B[cbind(seq_len(nrow(ed)), match(ed$from, nodes))] <- 1; B[cbind(seq_len(nrow(ed)), match(ed$to, nodes))] <- -1
  L <- crossprod(B); bal <- as.numeric(crossprod(B, delta)); comp <- components(g)$membership[nodes]; S <- numeric(length(nodes))
  for (cc in unique(comp)) { m <- comp == cc
    if (sum(m) > 1) { s <- as.numeric(MASS::ginv(L[m, m, drop = FALSE]) %*% bal[m]); S[m] <- s - mean(s) } }
  setNames(S, nodes)
}
if (PART == "bw" && !file.exists(out_file("bw"))) {
  alts <- list(global = list(global = TRUE), `local x2` = list(gamma = 2 * GAMMA), `local x1/4`  = list(gamma = GAMMA / 2),
               perplexity = list(perplexity = 4))
  ej <- expand.grid(Replicate = 1:50, Scenario = c("Redistributed", "Obstructed"), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
  edge_bw <- bind_rows(mclapply(seq_len(nrow(ej)), function(i) {
    g <- load_g(ej$Replicate[i], ej$Scenario[i], ej$generation[i]); if (is.null(g)) return(NULL)
    ref <- gravity_edges(g, GAMMA, gc_xg(g)); d0 <- ref$Delta; out <- list()
    for (nm in names(alts)) { W <- tryCatch(W_alt(g, alts[[nm]]), error = function(e) NULL); if (is.null(W)) next
      d1 <- delta_W(W, ref); nz <- d0 != 0 | d1 != 0
      out[[nm]] <- data.frame(ej[i, ], bandwidth = nm, concord = mean(sign(d1[nz]) == sign(d0[nz])),
                              rho_signed = suppressWarnings(cor(d1, d0, method = "spearman")),
                              rho_abs = suppressWarnings(cor(abs(d1), abs(d0), method = "spearman")), dbar_kept = sign(mean(d1)) == sign(mean(d0))) }
    bind_rows(out) }, mc.cores = CORES, mc.preschedule = FALSE))
  sj <- expand.grid(Replicate = 1:50, scenario = c("Isotropic", "Redistributed", "Obstructed"), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
  ss_bw <- bind_rows(mclapply(seq_len(nrow(sj)), function(i) {
    set.seed(20260930L + i)
    g <- gc_load_graph(sj$Replicate[i], sj$generation[i], c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)[[sj$scenario[i]]]); if (is.null(g)) return(NULL)
    x <- gc_xg(g); ed <- gravity_edges(g, GAMMA, x); out <- list()
    Aw <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE); bbar <- mean(E(g)$weight)
    Aw[Aw > 0] <- exp(-Aw[Aw > 0]^2 / (2 * bbar^2))                           # source_sink_test()'s null = "graph" kernel
    for (nm in c("reference", names(alts))) {
      a <- if (nm == "reference") list(gamma = GAMMA) else alts[[nm]]
      row <- tryCatch(if (!is.null(a$gamma)) { t <- source_sink_test(g, x, a$gamma); c(r_S = unname(t$statistic), p = t$p.value) } else {
        S <- signif(scores_from_delta(g, ed, delta_W(W_alt(g, a), ed)), 12)
        r_obs <- cor(S, x, method = "spearman"); r_null <- as.numeric(cor(msr_null(S, Aw, 999L), x, method = "spearman"))
        c(r_S = r_obs, p = (1 + sum(abs(r_null) >= abs(r_obs) - 1e-12)) / 1000) }, error = function(e) NULL)
      out[[nm]] <- if (is.null(row)) data.frame(sj[i, ], bandwidth = nm, r_S = NA, p = NA, failed = TRUE) else
        data.frame(sj[i, ], bandwidth = nm, r_S = row[["r_S"]], p = row[["p"]], failed = FALSE) }
    bind_rows(out) }, mc.cores = CORES, mc.preschedule = FALSE))
  save(edge_bw, ss_bw, file = out_file("bw")); cat("Wrote", out_file("bw"), "\n")
}

## ---- mirroring check -----------------------------------------------------------------------------
if (PART == "mirror" && !file.exists(out_file("mirror"))) {
  mirror <- bind_rows(mclapply(1:50, function(r) {
    e <- new.env(); load(sprintf("%s/replicate%d/rep%d-genotypes-1999.rda", GC_SIM_DIR, r, r), envir = e); gt <- e$genotypes
    rS <- function(pop) { gph <- popgraph(to_mv(gt[, setdiff(names(gt), c("ID", "Population"))]), groups = as.factor(pop))
      x <- gc_xg(gph); sc <- source_sink_scores(gph, GAMMA)
      c(r = cor(signif(sc$S, 12), x[sc$Stratum], method = "spearman"), dbar = mean(gravity_edges(gph, GAMMA, x)$Delta)) }
    flip <- sprintf("Pop%02d", 26L - as.integer(sub("\\D+", "", gt$Population)))
    o <- rS(gt$Population); m <- rS(flip)
    data.frame(Replicate = r, r = o[["r"]], r_mirror = m[["r"]], dbar = o[["dbar"]], dbar_mirror = m[["dbar"]])
  }, mc.cores = CORES))
  save(mirror, file = out_file("mirror")); cat("Wrote", out_file("mirror"), "\n")
}
