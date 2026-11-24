# R/case_study/araptus_null_ext_helpers.R  (2026-11-23)
# Local re-implementation of gstudio::popgraph_null() (1.15.8) that
#  (a) accepts a per-site Ne vector (the internal .pgn_wf() already recycles 2*Ne by
#      stratum, but popgraph_null() documents Ne as a scalar, and its tuning/print code
#      assumes one), and
#  (b) records extra per-null-graph statistics: cor(S, S0), degree variance, max/mean
#      degree, cor(k, site Hs), and per-site simulated Hs.
# Uses gstudio internals read-only; gstudio itself is not modified.  The simulation steps
# (pooled starting frequencies, Sinkhorn-balanced exp(-d/s) kernel, WF drift, HW genotypes
# with the observed missingness, popgraph with default alpha/tol, Gst of Nei & Chesser)
# are those of popgraph_null().
suppressPackageStartupMessages({ library(gstudio); library(igraph) })
pg <- function(f) getFromNamespace(f, "gstudio")
.pgn_freqs <- pg(".pgn_freqs"); .pgn_gst <- pg(".pgn_gst"); .pgn_wf <- pg(".pgn_wf")
.pgn_mix <- pg(".pgn_mix"); .pgn_balance <- pg(".pgn_balance"); .pgn_genotypes <- pg(".pgn_genotypes")
.pgn_sample_freqs <- pg(".pgn_sample_freqs"); .pgn_stats <- pg(".pgn_stats")

# Design objects exactly as in popgraph_null()
null_design <- function(x, axis, stratum = "Population") {
  loci <- column_class(x, "locus")
  grp <- factor(as.character(x[[stratum]])); strata <- levels(grp)
  mv <- lapply(loci, function(l) { m <- to_mv(x[[l]]); m[rowSums(m) == 0, ] <- NA; m }); names(mv) <- loci
  mv <- mv[vapply(mv, function(m) ncol(m) > 1, TRUE)]
  n_i <- c(table(grp))
  p0 <- lapply(mv, function(m) { p <- colMeans(m, na.rm = TRUE); p / sum(p) })
  fr_obs <- lapply(mv, function(m) .pgn_freqs(m, grp))
  ord <- order(grp)
  miss <- lapply(mv, function(m) !stats::complete.cases(m)[ord])
  n_l <- lapply(mv, function(m) c(table(factor(grp[stats::complete.cases(m)], levels = strata))))
  axis <- axis[strata]
  D <- abs(outer(axis, axis, "-")); dimnames(D) <- list(strata, strata)
  K <- length(strata); nn <- apply(D + diag(Inf, K), 1, min); s <- max(nn[is.finite(nn) & nn > 0])
  B <- .pgn_balance(exp(-D / s))
  list(strata = strata, n_i = n_i, p0 = p0, miss = miss, n_l = n_l, axis = axis, B = B, scale = s,
       obs_gst = .pgn_gst(fr_obs, n_i), hs_site_obs = site_hs(fr_obs))
}

# Per-site unbiased He, averaged over loci (as popgraph_diagnostics / Table B1a)
site_hs <- function(fr) {
  he <- vapply(fr, function(f) (2 * f$n / pmax(2 * f$n - 1, 1)) * (1 - rowSums(f$p^2)), numeric(nrow(fr[[1]]$p)))
  setNames(rowMeans(he, na.rm = TRUE), rownames(fr[[1]]$p))
}

# One simulated sample: per-site sampled frequencies (for Hs / Gst) from a WF run
sim_sample <- function(des, M, Ne, ngen) {
  P <- .pgn_wf(des$p0, M, Ne[des$strata], ngen)
  fr <- .pgn_sample_freqs(P, des$n_l)
  fr <- lapply(fr, function(f) { rownames(f$p) <- des$strata; names(f$n) <- des$strata; f })
  list(gst = .pgn_gst(fr, des$n_i)[["gst"]], hs = site_hs(fr))
}

# Bisection on m (12 steps, pilot replicates), as popgraph_null()
tune_m <- function(des, Ne, ngen, target, pilot = 8L) {
  # NA replicates (Gst undefined when every locus is fixed in the sample) are dropped;
  # if all are NA, differentiation is treated as above target (more migration needed).
  sim_gst <- function(m) mean(replicate(pilot, sim_sample(des, .pgn_mix(des$B, m), Ne, ngen)$gst), na.rm = TRUE)
  lo <- log(1e-4); hi <- log(0.95)
  for (it in seq_len(12)) { mid <- (lo + hi) / 2; g <- sim_gst(exp(mid))
    if (is.na(g) || g > target) lo <- mid else hi <- mid }
  exp((lo + hi) / 2)
}

# Null replicates with extended statistics
run_null_ext <- function(x, des, graph, distance, forward, M, Ne, ngen, nsim, gamma = 0.5, cores = 1L) {
  stat_ext <- function(g, hs_site) {
    base <- .pgn_stats(g, des$axis, distance, forward, gamma)
    k <- degree(g)[des$strata]
    sc <- tryCatch(source_sink_scores(g, gamma), error = function(e) NULL)
    c(base,
      cor_S_S0 = if (is.null(sc)) NA else attr(sc, "cor_S_S0"),
      cor_S_k  = if (is.null(sc)) NA else attr(sc, "cor_S_degree"),
      var_k = var(k), maxmean_k = max(k) / mean(k), mean_k = mean(k),
      cor_k_hs = suppressWarnings(cor(k, hs_site[des$strata], method = "spearman")))
  }
  hs_obs <- des$hs_site_obs
  observed <- stat_ext(graph, hs_obs)
  one <- function(b) {
    P <- .pgn_wf(des$p0, M, Ne[des$strata], ngen)
    X <- .pgn_genotypes(P, des$n_i, des$miss)
    fr <- lapply(X$by_locus, function(m) .pgn_freqs(m, X$groups))
    gs <- .pgn_gst(fr, des$n_i); hs <- site_hs(fr)
    g <- tryCatch(suppressWarnings(popgraph(X$mv, groups = X$groups)), error = function(e) NULL)
    st <- if (is.null(g)) rep(NA_real_, length(observed)) else stat_ext(g, hs)
    list(stats = c(st, gst = gs[["gst"]], hs = gs[["hs"]]), hs_site = hs)
  }
  res <- if (cores > 1L) parallel::mclapply(seq_len(nsim), one, mc.cores = cores, mc.set.seed = TRUE) else lapply(seq_len(nsim), one)
  null <- as.data.frame(do.call(rbind, lapply(res, `[[`, "stats")))
  names(null)[seq_along(observed)] <- names(observed)
  hs_site <- do.call(rbind, lapply(res, `[[`, "hs_site"))
  list(observed = observed, null = null, hs_site = hs_site, Ne = Ne, M = M, ngen = ngen, nsim = nsim)
}

# Table A1-style summary of a run_null_ext() result
summ_null <- function(z, alpha = 0.05) {
  nl <- z$null; o <- z$observed; conn <- nl$n_comp == 1
  q <- function(v) max(0, unname(quantile(v, 1 - alpha, na.rm = TRUE)))
  pct <- function(v, x) mean(v <= x, na.rm = TRUE)
  c(n_valid = sum(!is.na(nl$dR2)), frag = mean(nl$n_comp > 1, na.rm = TRUE),
    cutoff = q(nl$dR2), cutoff_conn = q(nl$dR2[conn]), fpr_rule = mean(nl$dR2 > 0, na.rm = TRUE),
    dR2_obs = o[["dR2"]], p_dR2 = (1 + sum(nl$dR2 >= o[["dR2"]], na.rm = TRUE)) / (1 + sum(!is.na(nl$dR2))),
    rS_obs = o[["r_S"]], p_rS = (1 + sum(abs(nl$r_S) >= abs(o[["r_S"]]) - 1e-12, na.rm = TRUE)) / (1 + sum(!is.na(nl$r_S))),
    frac_rS_neg = mean(nl$r_S < 0, na.rm = TRUE), frac_ratio_lt1 = mean(nl$slope_ratio < 1, na.rm = TRUE),
    ds_obs = o[["degree_share"]], ds_med = median(nl$degree_share, na.rm = TRUE),
    ds_q05 = unname(quantile(nl$degree_share, 0.05, na.rm = TRUE)), ds_q95 = unname(quantile(nl$degree_share, 0.95, na.rm = TRUE)),
    ds_pct = pct(nl$degree_share, o[["degree_share"]]),
    cSS0_obs = o[["cor_S_S0"]], cSS0_med = median(nl$cor_S_S0, na.rm = TRUE),
    cSS0_q05 = unname(quantile(nl$cor_S_S0, 0.05, na.rm = TRUE)), cSS0_q95 = unname(quantile(nl$cor_S_S0, 0.95, na.rm = TRUE)),
    cSS0_pct = pct(nl$cor_S_S0, o[["cor_S_S0"]]),
    vark_obs = o[["var_k"]], vark_med = median(nl$var_k, na.rm = TRUE), vark_pct = pct(nl$var_k, o[["var_k"]]),
    mmk_obs = o[["maxmean_k"]], mmk_med = median(nl$maxmean_k, na.rm = TRUE), mmk_pct = pct(nl$maxmean_k, o[["maxmean_k"]]),
    ckhs_obs = o[["cor_k_hs"]], ckhs_med = median(nl$cor_k_hs, na.rm = TRUE), ckhs_pct = pct(nl$cor_k_hs, o[["cor_k_hs"]]),
    gst_null = mean(nl$gst, na.rm = TRUE), hs_null = mean(nl$hs, na.rm = TRUE))
}
