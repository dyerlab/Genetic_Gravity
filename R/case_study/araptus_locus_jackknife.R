# R/case_study/araptus_locus_jackknife.R  (2026-11-23)
# Task B of local/araptus_checks_instructions.md (S11 item 3; S10.5 consistency
# across loci subsets): information content per site, leave-one- and leave-two-
# loci-out rebuilds of the peninsular Population Graph, and design-matched
# recalibration of Delta R^2 for each leave-one-out panel (central component).
# Construction, separation, forward rule and gamma = 1/2 follow
# R/case_study/case_study_arapat.R and R/case_study/case_study_arapat_calibration.R.
suppressPackageStartupMessages({ library(gstudio); library(igraph) })
stopifnot(packageVersion("gstudio") >= "1.15.8")
load("data/case_study_arapat.rda")
loci  <- column_class(peninsula, "locus")
cores <- max(1L, parallel::detectCores() - 1L)
comp_sites <- lapply(results, function(r) V(r$graph)$name)   # South = central, North = northern
t0 <- Sys.time()

## 1. Information content per site -------------------------------------------------------------
info_table <- function(sites) {
  x <- droplevels(peninsula[peninsula$Population %in% sites, ])
  grp <- factor(as.character(x$Population))
  per <- lapply(loci, function(l) {
    m <- to_mv(x[[l]]); m[rowSums(m) == 0, ] <- NA
    f <- gstudio:::.pgn_freqs(m, grp)
    he <- (2 * f$n / pmax(2 * f$n - 1, 1)) * (1 - rowSums(f$p^2))
    data.frame(site = levels(grp), locus = l, n = f$n, A = rowSums(f$p > 0, na.rm = TRUE), Hs = he, nHs = f$n * he)
  })
  per <- do.call(rbind, per)
  out <- do.call(rbind, lapply(split(per, per$site), function(z) data.frame(
    site = z$site[1], n = max(z$n), poly = sum(z$A > 1), A_mean = mean(z$A), Hs_mean = mean(z$Hs, na.rm = TRUE),
    nHs_mean = mean(z$nHs, na.rm = TRUE), nHs_min = min(z$nHs, na.rm = TRUE))))
  out$d_km <- round(d[out$site]); out$lat <- V(graph)$Latitude[match(out$site, V(graph)$name)]
  out <- out[order(out$d_km), ]; rownames(out) <- NULL
  list(sites = out, per_locus = per)
}
info <- lapply(comp_sites, info_table)

## 2-3. Rebuilds ---------------------------------------------------------------------------------
component_stats <- function(g, sites, with_p) {
  memb <- components(g)$membership
  m_sites <- unique(memb[sites])
  same <- length(m_sites) == 1 && sum(memb == m_sites) == length(sites)
  sub <- induced_subgraph(g, sites)
  dk  <- d[V(sub)$name]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
  fp  <- tryCatch(ibgd(sub, mode = "dcGD", distance = sep, forward = fwd, gamma = 0.5), error = function(e) NULL)
  sc  <- tryCatch(source_sink_scores(sub, 0.5), error = function(e) NULL)
  ga  <- tryCatch(graph_asymmetries(sub, gamma = 0.5), error = function(e) NULL)
  rS  <- if (is.null(sc)) NA else suppressWarnings(cor(gstudio:::.merge_ties(sc$S), dk[sc$Stratum], method = "spearman"))
  pS  <- if (with_p) tryCatch(source_sink_test(sub, x = dk, gamma = 0.5, nperm = 999)$p.value, error = function(e) NA) else NA
  o   <- if (is.null(sc)) NULL else sc[order(-sc$S), ]
  data.frame(same_component = same, sub_edges = ecount(sub), sub_components = components(sub)$no,
             dR2 = if (is.null(fp)) NA else unname(fp$statistic),
             ratio = if (is.null(fp)) NA else unname(fp$estimate[["forward slope"]] / fp$estimate[["reverse slope"]]),
             r_S = rS, p_S = pS,
             degree_share = if (is.null(sc)) NA else attr(sc, "degree_share"),
             mean_abs_delta = if (is.null(ga)) NA else mean(abs(E(ga)$delta)),
             sources = if (is.null(o)) NA else paste(head(o$Stratum, 3), collapse = ", "),
             sinks   = if (is.null(o)) NA else paste(rev(tail(o$Stratum, 3)), collapse = ", "))
}
rebuild <- function(drop, with_p) {
  x <- peninsula[, setdiff(names(peninsula), drop)]
  g <- population_graph(x, stratum = "Population", decorate = TRUE)
  cm <- components(g)
  rows <- lapply(names(comp_sites), function(h) cbind(component = h, component_stats(g, comp_sites[[h]], with_p)))
  cbind(dropped = paste(drop, collapse = "+"), edges = ecount(g), components = cm$no,
        sizes = paste(sort(cm$csize, decreasing = TRUE), collapse = ","), do.call(rbind, rows))
}
set.seed(2013)
full <- rebuild(character(0), TRUE)
loo  <- do.call(rbind, lapply(loci, rebuild, with_p = TRUE))
pairs <- combn(loci, 2, simplify = FALSE)
lto  <- do.call(rbind, lapply(pairs, rebuild, with_p = FALSE))
t1 <- Sys.time()
cat(sprintf("rebuilds: %.1f s\n", as.numeric(difftime(t1, t0, units = "secs"))))

lto_summary <- do.call(rbind, lapply(split(lto, lto$component), function(z) data.frame(
  component = z$component[1], n = nrow(z), same_component = mean(z$same_component),
  frac_rS_neg = mean(z$r_S < 0, na.rm = TRUE), frac_ratio_lt1 = mean(z$ratio < 1, na.rm = TRUE),
  frac_dR2_pos = mean(z$dR2 > 0, na.rm = TRUE), dR2_range = paste(sprintf("%.3f", range(z$dR2, na.rm = TRUE)), collapse = " to "),
  rS_range = paste(sprintf("%.2f", range(z$r_S, na.rm = TRUE)), collapse = " to "))))

## 4. Recalibration of each leave-one-out panel (central component) ------------------------------
recal <- function(drop) {
  x <- peninsula[, setdiff(names(peninsula), drop)]
  g <- population_graph(x, stratum = "Population", decorate = TRUE)
  sites <- comp_sites$South
  sub <- induced_subgraph(g, sites); nm <- V(sub)$name
  dd  <- droplevels(x[x$Population %in% nm, ])
  dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
  set.seed(20261129L)
  z <- popgraph_null(x = dd, axis = dk, distance = sep, forward = fwd, graph = sub, nsim = 999L, cores = cores)
  data.frame(dropped = paste(drop, collapse = "+"), dR2 = z$observed[["dR2"]], cutoff = z$cutoff[["all"]],
             cutoff_conn = z$cutoff[["connected"]], p_dR2 = z$p.value[["dR2"]], fpr_rule = z$fpr_rule,
             p_rS_design = z$p.value[["r_S"]], frag = mean(z$null$n_comp > 1, na.rm = TRUE),
             gst_target = z$tuning[["gst_target"]], gst_null = z$tuning[["gst_null"]], m = z$tuning[["m"]])
}
t2 <- Sys.time()
calib <- do.call(rbind, lapply(loci, recal))
t3 <- Sys.time()
cat(sprintf("recalibration (8 x 999): %.1f s\n", as.numeric(difftime(t3, t2, units = "secs"))))

## 5. Inflation check ----------------------------------------------------------------------------
infl <- merge(loo[, c("dropped", "component", "mean_abs_delta")],
              full[, c("component", "mean_abs_delta")], by = "component", suffixes = c("", "_full"))
infl$ratio <- infl$mean_abs_delta / infl$mean_abs_delta_full

options(width = 200)
for (h in names(info)) { cat("\n== info:", h, "==\n"); print(info[[h]]$sites, digits = 3) }
cat("\n== full ==\n"); print(full, digits = 3)
cat("\n== leave-one-out ==\n"); print(loo, digits = 3)
cat("\n== leave-two-out summary ==\n"); print(lto_summary, digits = 3)
cat("\n== recalibration ==\n"); print(calib, digits = 3)
cat("\n== inflation ==\n"); print(infl[order(infl$component, -infl$ratio), ], digits = 3)
run_times <- c(rebuilds = as.numeric(difftime(t1, t0, units = "secs")), recalibration = as.numeric(difftime(t3, t2, units = "secs")))
save(info, full, loo, lto, lto_summary, calib, infl, run_times, file = "data/derived/araptus_locus_jackknife.rda")
