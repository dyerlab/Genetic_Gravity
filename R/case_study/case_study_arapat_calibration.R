#' Case study: design-specific calibration of the directional isolation comparison
#' for the central component of Araptus attenuatus (R/case_study/case_study_arapat.R).
#'
#' gstudio::popgraph_null() (gstudio >= 1.15.8) simulates data sets with the
#' central component's design -- its 14 sites, the number of genotyped beetles at
#' each, and its 8 loci with their alleles -- under symmetric migration through a
#' distance kernel along the expansion axis, tuned so that the simulated samples
#' match the observed differentiation (Gst). For each, the Population Graph is rebuilt
#' and Delta R^2 computed exactly as for the observed data (separation |d_i - d_j|,
#' forward when j lies farther from the refuge). The cutoff is the 95th percentile of
#' the null Delta R^2; the one-sided p is the rank of the observed value.
#' Sensitivity: Ne (50, 200), generations before sampling (200, 1000), and a
#' two-dimensional great-circle kernel in place of the axis.
#'
#' Output: data/case_study_arapat_calibration.rda and the summary printed below.
suppressPackageStartupMessages({ library(gstudio); library(igraph); library(dplyr) })
stopifnot(packageVersion("gstudio") >= "1.15.8")
load("data/case_study_arapat.rda")
sub <- results$South$graph; nm <- V(sub)$name
dd  <- droplevels(peninsula[peninsula$Population %in% nm, ])
dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
xy  <- data.frame(Stratum = nm, Longitude = V(sub)$Longitude, Latitude = V(sub)$Latitude)
geo <- strata_distance(xy, mode = "Circle"); geo <- as.matrix(geo)[nm, nm]
cores <- max(1L, parallel::detectCores() - 1L)
runs <- list(default = list(), Ne50 = list(Ne = 50L, ngen = 200L), Ne200 = list(Ne = 200L, ngen = 800L),
             ngen200 = list(ngen = 200L), ngen1000 = list(ngen = 1000L), geographic = list(kernel_distance = geo))
diag_central <- popgraph_diagnostics(dd, graph = sub)
cal <- list()
for (r in names(runs)) {
  set.seed(20261129L)
  cal[[r]] <- do.call(popgraph_null, c(list(x = dd, axis = dk, distance = sep, forward = fwd, graph = sub,
                                            nsim = 999L, cores = cores), runs[[r]]))
  cat("==", r, "==\n"); print(cal[[r]])
}
summ <- bind_rows(lapply(names(cal), function(r) { z <- cal[[r]]; tu <- z$tuning
  data.frame(run = r, Ne = tu[["Ne"]], ngen = tu[["ngen"]], m = tu[["m"]], gst_target = tu[["gst_target"]], gst_null = tu[["gst_null"]],
             hs_obs = tu[["hs_observed"]], hs_null = tu[["hs_null"]], frag = mean(z$null$n_comp > 1, na.rm = TRUE),
             fpr_rule = z$fpr_rule, cutoff = z$cutoff[["all"]], cutoff_conn = z$cutoff[["connected"]],
             dR2_obs = z$observed[["dR2"]], p_dR2 = z$p.value[["dR2"]], p_rS_design = z$p.value[["r_S"]]) }))
print(diag_central); print(summ, digits = 3, row.names = FALSE)
save(cal, summ, diag_central, file = "data/case_study_arapat_calibration.rda")
