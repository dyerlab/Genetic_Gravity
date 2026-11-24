# R/case_study/araptus_degree_structure.R  (2026-11-23)
# Task C1 of local/araptus_degree_check_instructions.md: degree structure of the
# Araptus components, and the central component's observed degree statistics against
# the default (equal-Ne) design-matched null.  The saved default null
# (data/case_study_arapat_calibration.rda) does not store per-graph cor(S, S0) or degree,
# so the null is regenerated with the same settings (Ne = 100, 400 generations, the saved
# tuned m = 0.0209, 999 data sets, seed 20261129) via R/case_study/araptus_null_ext_helpers.R.
source("R/case_study/araptus_null_ext_helpers.R")
load("data/case_study_arapat.rda")
load("data/derived/araptus_locus_jackknife.rda")      # info: per-site Hs (Table B1a/b)
cores <- max(1L, parallel::detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG")

## 1-3. Per-site table and correlations ----------------------------------------------------------
site_table <- function(h) {
  sub <- results[[h]]$graph; sc <- source_sink_scores(sub, 0.5)
  tb <- data.frame(site = sc$Stratum, d_km = round(d[sc$Stratum]), k = sc$degree,
                   s = strength(sub)[sc$Stratum], Hs = info[[h]]$sites$Hs_mean[match(sc$Stratum, info[[h]]$sites$site)],
                   n = info[[h]]$sites$n[match(sc$Stratum, info[[h]]$sites$site)],
                   S = sc$S, S0 = sc$S0, inv_k = -1 / sc$degree, balance = sc$balance)
  tb <- tb[order(tb$d_km), ]; rownames(tb) <- NULL; tb
}
sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(rho = unname(ct$estimate), p = ct$p.value) }
tabs <- lapply(c(South = "South", North = "North"), site_table)
cors <- do.call(rbind, lapply(names(tabs), function(h) { t <- tabs[[h]]
  rbind(data.frame(component = h, pair = "k ~ d", t(sp(t$k, t$d_km))),
        data.frame(component = h, pair = "k ~ Hs", t(sp(t$k, t$Hs))),
        data.frame(component = h, pair = "S ~ S0", t(sp(t$S, t$S0))),
        data.frame(component = h, pair = "S ~ Hs", t(sp(t$S, t$Hs))),
        data.frame(component = h, pair = "s ~ Hs", t(sp(t$s, t$Hs))),
        data.frame(component = h, pair = "Hs ~ d", t(sp(t$Hs, t$d_km)))) }))
rank_info <- lapply(tabs, function(t) { o <- t[order(-t$S), ]
  list(sources = head(o[, c("site", "k", "S")], 3), sinks = tail(o[, c("site", "k", "S")], 3),
       k_rank_desc = setNames(rank(-t$k, ties.method = "min"), t$site)) })

## 4. Central null comparison (regenerated equal-Ne default null) --------------------------------
load("data/case_study_arapat_calibration.rda")    # cal$default: saved tuning
sub <- results$South$graph; nm <- V(sub)$name
dd  <- droplevels(peninsula[peninsula$Population %in% nm, ])
dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
des <- null_design(dd, dk)
m0  <- cal$default$tuning[["m"]]
M0  <- .pgn_mix(des$B, m0); Ne0 <- setNames(rep(100, length(des$strata)), des$strata)
set.seed(20261129L)
t_pilot <- system.time(run_null_ext(dd, des, sub, sep, fwd, M0, Ne0, 400L, 20L, cores = 1L))[["elapsed"]]
cat(sprintf("pilot: 20 data sets, 1 core: %.1f s (999 on %d cores ~ %.0f s)\n", t_pilot, cores, t_pilot / 20 * 999 / cores))
set.seed(20261129L)
t_full <- system.time(null_eq <- run_null_ext(dd, des, sub, sep, fwd, M0, Ne0, 400L, 999L, cores = cores))[["elapsed"]]
cat(sprintf("full: %.1f s\n", t_full))
eq_summary <- summ_null(null_eq)

options(width = 200)
for (h in names(tabs)) { cat("\n==", h, "==\n"); print(tabs[[h]], digits = 3) }
cat("\n"); print(cors, digits = 3)
for (h in names(rank_info)) { cat("\n--", h, "sources / sinks (k rank among sites, 1 = highest) --\n")
  r <- rank_info[[h]]; print(cbind(r$sources, k_rank = r$k_rank_desc[r$sources$site]), digits = 3)
  print(cbind(r$sinks, k_rank = r$k_rank_desc[r$sinks$site]), digits = 3) }
cat("\n-- central, equal-Ne null (regenerated) --\n"); print(round(eq_summary, 4))
cat("\nsaved default run for comparison: cutoff", cal$default$cutoff[["all"]], " p", cal$default$p.value[["dR2"]],
    " ds median", median(cal$default$null$degree_share, na.rm = TRUE), "\n")
save(tabs, cors, rank_info, null_eq, eq_summary, t_pilot, t_full, file = "data/derived/araptus_degree_structure.rda")
