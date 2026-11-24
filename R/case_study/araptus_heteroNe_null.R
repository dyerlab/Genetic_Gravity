# R/case_study/araptus_heteroNe_null.R  (2026-11-23)
# Task C2 of local/araptus_degree_check_instructions.md: a symmetric design-matched null
# with per-site Ne tuned so that each site's simulated Hs tracks its observed Hs.
# popgraph_null() documents Ne as a scalar (its tuning summary and print assume one),
# so this uses the local wrapper in R/case_study/araptus_null_ext_helpers.R, verified to
# reproduce popgraph_null() exactly at equal Ne (same seed: identical dR2 and degree share).
# Same kernel (Sinkhorn-balanced exp(-|d_i - d_j|/s)), 400 generations, mean Ne 100,
# migration tuned to the component's observed Gst, 999 data sets, seed 20261129.
#
# Three nulls, each with m tuned by bisection to the observed Gst and mean Ne = 100:
#  equalNe    -- the default (reference);
#  gradientNe -- primary: Ne log-linear in distance from the refuge, slope chosen so the
#                null slope of site Hs on d matches the observed slope;
#  siteNe     -- per-site Ne fitted (damped) to each site's relative Hs.
# Relative Hs is matched because the equal-Ne null already sits below the observed mean Hs
# at Ne = 100 (0.296 vs 0.337) and the mean Ne is fixed by the instructions.
source("R/case_study/araptus_null_ext_helpers.R")
load("data/case_study_arapat.rda")
cores <- max(1L, parallel::detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG")

# Mean per-site Hs and Gst of R simulated samples
eval_Ne <- function(des, Ne, m, ngen, R) {
  M <- .pgn_mix(des$B, m)
  sims <- replicate(R, sim_sample(des, M, Ne, ngen), simplify = FALSE)
  list(hs = Reduce(`+`, lapply(sims, `[[`, "hs")) / R, gst = mean(vapply(sims, `[[`, 1, "gst"), na.rm = TRUE))
}
fit_stats <- function(des, ev) {
  hs_t <- des$hs_site_obs[des$strata]; hs <- ev$hs[des$strata]; dk <- des$axis[des$strata]
  c(cor_hs = cor(hs, hs_t), rmse_rel = sqrt(mean((hs / mean(hs) - hs_t / mean(hs_t))^2)),
    slope_null = unname(coef(lm(hs ~ dk))[2]), slope_obs = unname(coef(lm(hs_t ~ dk))[2]), gst = ev$gst)
}

# (a) Per-site, damped: log Ne_i += beta * log(rel_target / rel_sim); mean 100; Ne in [lo, hi].
#     Undamped (beta = 3, bounds 5-5000) diverged: Ne collapsed to 5 at most sites and the
#     Gst bisection failed by iteration 12.  Keep the iterate with the smallest rmse_rel.
tune_site <- function(des, ngen = 400L, iters = 10L, R = 40L, beta = 1, bounds = c(10, 1000), mean_Ne = 100) {
  target <- des$obs_gst[["gst"]]; rel_t <- des$hs_site_obs / mean(des$hs_site_obs)
  Ne <- setNames(rep(mean_Ne, length(des$strata)), des$strata)
  hist <- list(); best <- NULL
  for (it in seq_len(iters)) {
    m <- tune_m(des, Ne, ngen, target); ev <- eval_Ne(des, Ne, m, ngen, R); fs <- fit_stats(des, ev)
    hist[[it]] <- data.frame(iter = it, m = m, t(fs), Ne_min = min(Ne), Ne_max = max(Ne))
    if (is.null(best) || fs[["rmse_rel"]] < best$fs[["rmse_rel"]]) best <- list(Ne = Ne, m = m, fs = fs, iter = it)
    rel <- ev$hs / mean(ev$hs)
    Ne <- Ne * exp(beta * log(rel_t[des$strata] / rel[des$strata]))
    for (k in 1:5) { Ne <- pmin(pmax(Ne, bounds[1]), bounds[2]); Ne <- Ne * mean_Ne / mean(Ne) }
    Ne <- round(pmin(pmax(Ne, bounds[1]), bounds[2]))
  }
  c(best, list(history = do.call(rbind, hist), gst_target = target))
}

# (b) Smooth gradient: Ne_i = 100 * exp(b z_i) / mean(exp(b z)), z = standardized d; b chosen
#     on a grid so that the null slope of Hs on d matches the observed slope.
tune_gradient <- function(des, ngen = 400L, R = 40L, grid = seq(0, -2, by = -0.25), mean_Ne = 100) {
  target <- des$obs_gst[["gst"]]; dk <- des$axis[des$strata]; z <- (dk - mean(dk)) / sd(dk)
  mk <- function(b) { Ne <- mean_Ne * exp(b * z) / mean(exp(b * z)); setNames(pmax(1, round(Ne)), des$strata) }
  hist <- lapply(grid, function(b) { Ne <- mk(b); m <- tune_m(des, Ne, ngen, target)
    fs <- fit_stats(des, eval_Ne(des, Ne, m, ngen, R)); data.frame(b = b, m = m, t(fs), Ne_min = min(Ne), Ne_max = max(Ne)) })
  hist <- do.call(rbind, hist)
  # first crossing: the smallest |b| whose null slope reaches the observed slope, or the
  # grid point just before it, whichever is closer (the null slope is not monotone in b:
  # at large |b| the low-Ne sites fix and the slope flattens again, giving spurious matches)
  j <- which(hist$slope_null <= hist$slope_obs)[1]
  i <- if (is.na(j)) which.min(abs(hist$slope_null - hist$slope_obs)) else
         c(j - 1L, j)[which.min(abs(hist$slope_null[c(j - 1L, j)] - hist$slope_obs[1]))]
  b <- hist$b[i]; Ne <- mk(b)
  list(Ne = Ne, m = hist$m[i], b = b, fs = unlist(hist[i, c("cor_hs", "rmse_rel", "slope_null", "slope_obs", "gst")]),
       history = hist, gst_target = target)
}

run_component <- function(h) {
  sub <- results[[h]]$graph; nm <- V(sub)$name
  dd  <- droplevels(peninsula[peninsula$Population %in% nm, ])
  dk  <- d[nm]; sep <- abs(outer(dk, dk, "-")); fwd <- outer(dk, dk, function(di, dj) dj > di)
  des <- null_design(dd, dk)
  target <- des$obs_gst[["gst"]]
  set.seed(20261129L); t_g <- system.time(tg <- tune_gradient(des))[["elapsed"]]
  set.seed(20261129L); t_s <- system.time(ts <- tune_site(des))[["elapsed"]]
  set.seed(20261129L)
  Ne0 <- setNames(rep(100, length(des$strata)), des$strata); m0 <- tune_m(des, Ne0, 400L, target)
  set.seed(20261129L)
  t_pilot <- system.time(run_null_ext(dd, des, sub, sep, fwd, .pgn_mix(des$B, tg$m), tg$Ne, 400L, 20L, cores = 1L))[["elapsed"]]
  cat(sprintf("%s: tuning gradient %.0f s, per-site %.0f s; pilot 20 data sets %.1f s (999 on %d cores ~ %.0f s)\n",
              h, t_g, t_s, t_pilot, cores, t_pilot / 20 * 999 / cores))
  runs <- list(equalNe = list(Ne = Ne0, m = m0), gradientNe = tg[c("Ne", "m")], siteNe = ts[c("Ne", "m")])
  nulls <- lapply(runs, function(r) { set.seed(20261129L)
    run_null_ext(dd, des, sub, sep, fwd, .pgn_mix(des$B, r$m), r$Ne, 400L, 999L, cores = cores) })
  site <- data.frame(site = des$strata, d_km = round(dk[des$strata]), n = des$n_i,
                     Ne_gradient = tg$Ne[des$strata], Ne_site = ts$Ne[des$strata], Hs_obs = des$hs_site_obs[des$strata],
                     Hs_null_equal = colMeans(nulls$equalNe$hs_site, na.rm = TRUE)[des$strata],
                     Hs_null_gradient = colMeans(nulls$gradientNe$hs_site, na.rm = TRUE)[des$strata],
                     Hs_null_site = colMeans(nulls$siteNe$hs_site, na.rm = TRUE)[des$strata])
  site <- site[order(site$d_km), ]; rownames(site) <- NULL
  hs_fit <- sapply(c(equalNe = "Hs_null_equal", gradientNe = "Hs_null_gradient", siteNe = "Hs_null_site"), function(v)
    c(pearson = cor(site$Hs_obs, site[[v]]), spearman = cor(site$Hs_obs, site[[v]], method = "spearman"),
      slope = unname(coef(lm(site[[v]] ~ site$d_km))[2])))
  list(gradient = tg, site_tune = ts, m_eq = m0, site = site, hs_fit = hs_fit, nulls = nulls,
       slope_obs = unname(coef(lm(site$Hs_obs ~ site$d_km))[2]),
       summary = do.call(rbind, lapply(nulls, summ_null)), gst_target = target,
       times = c(tune_gradient = t_g, tune_site = t_s, pilot = t_pilot))
}

out <- lapply(c(South = "South", North = "North"), run_component)
options(width = 220)
for (h in names(out)) {
  o <- out[[h]]
  cat("\n=====", h, "=====  Gst target", round(o$gst_target, 3), "\n")
  cat("-- gradient tuning --\n"); print(o$gradient$history, digits = 3); cat("chosen b =", o$gradient$b, "\n")
  cat("-- per-site tuning (best iterate", o$site_tune$iter, ") --\n"); print(o$site_tune$history, digits = 3)
  cat(sprintf("m: equal %.4g, gradient %.4g, site %.4g\n", o$m_eq, o$gradient$m, o$site_tune$m))
  print(o$site, digits = 3); cat("observed Hs~d slope", signif(o$slope_obs, 3), "\n"); print(round(o$hs_fit, 5))
  print(t(round(o$summary, 4)))
}
save(out, file = "data/derived/araptus_heteroNe_null.rda")
