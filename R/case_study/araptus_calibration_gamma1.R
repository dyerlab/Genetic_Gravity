# R/case_study/araptus_calibration_gamma1.R  (2026-11-23)
# Design-matched calibration of Delta R^2 for the central Araptus component at
# gamma = 1, mirroring the default run of R/case_study/case_study_arapat_calibration.R
# (same seed and settings; the migration tuning lands slightly differently, m 0.0215 vs 0.0209).
suppressPackageStartupMessages({ library(gstudio); library(igraph) })
load("data/case_study_arapat.rda")
sub <- results$South$graph; nm <- V(sub)$name
dd  <- droplevels(peninsula[peninsula$Population %in% nm, ])
dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
cores <- max(1L, parallel::detectCores() - 1L)
set.seed(20261129L)
cal_g1 <- popgraph_null(x = dd, axis = dk, distance = sep, forward = fwd, graph = sub,
                        nsim = 999L, cores = cores, gamma = 1)
print(cal_g1)
cat(sprintf("\ngamma = 1: dR2 obs %.4f, cutoff %.4f (connected %.4f), p %.3f, rule FPR %.3f, p r(S,d) design %.3f\n",
            cal_g1$observed[["dR2"]], cal_g1$cutoff[["all"]], cal_g1$cutoff[["connected"]],
            cal_g1$p.value[["dR2"]], cal_g1$fpr_rule, cal_g1$p.value[["r_S"]]))
save(cal_g1, file = "data/derived/araptus_calibration_gamma1.rda")
