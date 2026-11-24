# R/case_study/araptus_northern_calibration.R  (2026-11-23)
# Task A of local/araptus_checks_instructions.md: design-matched calibration of the
# northern Araptus component (S11 item 7), following the default run of
# R/case_study/case_study_arapat_calibration.R (999 data sets, Ne = 100, 400 generations,
# seed 20261129, gamma = 1/2, migration tuned to the component's own Gst).
# Also the null degree-share summary for the central component from its saved run.
suppressPackageStartupMessages({ library(gstudio); library(igraph) })
stopifnot(packageVersion("gstudio") >= "1.15.8")
load("data/case_study_arapat.rda")
sub <- results$North$graph; nm <- V(sub)$name
dd  <- droplevels(peninsula[peninsula$Population %in% nm, ])
dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
cores <- max(1L, parallel::detectCores() - 1L)

# pilot timing: 20 data sets
set.seed(20261129L)
t_pilot <- system.time(popgraph_null(x = dd, axis = dk, distance = sep, forward = fwd, graph = sub,
                                     nsim = 20L, cores = 1L))[["elapsed"]]
cat(sprintf("pilot (tuning + 20 data sets, 1 core): %.1f s\n", t_pilot))

diag_north <- popgraph_diagnostics(dd, graph = sub)
set.seed(20261129L)
t_full <- system.time(
  cal_north <- popgraph_null(x = dd, axis = dk, distance = sep, forward = fwd, graph = sub,
                             nsim = 999L, cores = cores))[["elapsed"]]
cat(sprintf("full run (999 data sets, %d cores): %.1f s\n", cores, t_full))
print(cal_north); print(diag_north)

null_summary <- function(z) {
  nl <- z$null; ok <- !is.na(nl$dR2)
  ds <- nl$degree_share
  c(n_valid = sum(ok), frag = mean(nl$n_comp > 1, na.rm = TRUE),
    cutoff = z$cutoff[["all"]], cutoff_conn = z$cutoff[["connected"]],
    fpr_rule = z$fpr_rule, dR2_obs = z$observed[["dR2"]], p_dR2 = z$p.value[["dR2"]],
    rS_obs = z$observed[["r_S"]], p_rS = z$p.value[["r_S"]],
    frac_rS_neg = mean(nl$r_S < 0, na.rm = TRUE), frac_ratio_lt1 = mean(nl$slope_ratio < 1, na.rm = TRUE),
    ds_obs = z$observed[["degree_share"]], ds_med = median(ds, na.rm = TRUE),
    ds_q05 = unname(quantile(ds, 0.05, na.rm = TRUE)), ds_q95 = unname(quantile(ds, 0.95, na.rm = TRUE)),
    ds_pctile = mean(ds <= z$observed[["degree_share"]], na.rm = TRUE),
    ds_med_conn = median(ds[nl$n_comp == 1], na.rm = TRUE),
    m = z$tuning[["m"]], gst_target = z$tuning[["gst_target"]], gst_null = z$tuning[["gst_null"]],
    hs_obs = z$tuning[["hs_observed"]], hs_null = z$tuning[["hs_null"]])
}
load("data/case_study_arapat_calibration.rda")   # cal (central), summ, diag_central
tab <- rbind(North = null_summary(cal_north), Central = null_summary(cal$default))
print(t(round(tab, 4)))
save(cal_north, diag_north, tab, t_pilot, t_full, file = "data/derived/araptus_northern_calibration.rda")
