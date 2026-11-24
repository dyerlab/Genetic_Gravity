# R/analysis/half_bandwidth_results.R
#
# The main Results at gamma = 1/2 beside gamma = 1
# (local/half_bandwidth_rerun_instructions.md). Both columns come from
# R/simulation/half_bandwidth_compute.R run with identical code, graphs, phases, null
# adjacency and seeds (S rounded to 12 significant digits before ranking in both), so
# old and new differ only by bandwidth. divMigrate's Nm and every graph-only quantity
# (edge counts, degree, C_D, diameter, phase boundaries) do not depend on the bandwidth.
#
# Inputs: data/derived/half_bandwidth_{grid,tests}_g{1,05}.rda, half_bandwidth_{bw,mirror}_g05.rda,
#   data/information_null.csv and data/information_null_g05.csv, node_potential_reanalysis.rda
#   (Nm), node_potential_test.rda (unrounded gamma = 1 tests), phase_partitioned_results.rda and
#   ss_bandwidth_by_phase.rda (current bandwidth tables), data/forward_scenarios.csv (diversity).
# Outputs: local/half_bandwidth_rerun_results.md (table; the summary is written by hand above it),
#   data/derived/half_bandwidth_results.rda, and the gamma = 1/2 data files the figure scripts read
#   with GG_TAG=g05: data/forward_scenarios_g05.csv, data/isotropic_baseline_g05.csv,
#   data/derived/ibgd_trajectory_fit_g05.rda, data/derived/node_gravity_heatmap_g05.rds,
#   data/derived/sstest_g05.rds, data/derived/truth_g05.rds.
#
# Run from the repository root:  Rscript R/analysis/half_bandwidth_results.R

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv) })
source("R/common/gravity_convention.R")
source("R/common/lineage_phases.R")
ld <- function(f) { e <- new.env(); load(f, envir = e); e }
G1 <- ld("data/derived/half_bandwidth_grid_g1.rda"); G5 <- ld("data/derived/half_bandwidth_grid_g05.rda")
T1 <- ld("data/derived/half_bandwidth_tests_g1.rda")$tests; T5 <- ld("data/derived/half_bandwidth_tests_g05.rda")$tests
BW <- ld("data/derived/half_bandwidth_bw_g05.rda"); MIR <- ld("data/derived/half_bandwidth_mirror_g05.rda")$mirror
rean <- ld("data/derived/node_potential_reanalysis.rda")$rean
r1u <- ld("data/derived/node_potential_test.rda")$tests_tbl |> filter(field == "phi", nodes == "all", cor == "spearman", null == "N1w")
ep <- ld("data/derived/phase_partitioned_results.rda")$ppr; es <- ld("data/derived/ss_bandwidth_by_phase.rda")
fwd <- gc_read_forward() |> rename(Replicate = replicate)
AS <- c("Redistributed", "Obstructed"); SYMS <- c("Burn-in", "Symmetric", "sym-mid", "sym-low", "sym-verylow")
PH <- c(Reorganize = "ascent", Plateau = "plateau", `Fall apart` = "decay")
for (x in c("G1", "G5")) { g <- get(x); stopifnot(!"err" %in% names(g$grid) || all(is.na(g$grid$err))) }
stopifnot(!"err" %in% names(T1) || all(is.na(T1$err)), !"err" %in% names(T5) || all(is.na(T5$err)))
fixT <- function(t) t |> mutate(sample = recode(scenario, Isotropic = "Symmetric"))
T1 <- fixT(T1); T5 <- fixT(T5)
ph <- function(d, s = "sample", gen = "generation", rcol = "phase_policy_R_gen")
  d |> mutate(phase = as.character(lp_phase(.data[[s]], Replicate, .data[[gen]], rcol = rcol)))

## ---- Step 1: gamma = 1 test (gstudio source_sink_test(), S rounded, new null draws) against the round-1 test (unrounded) ---
tie <- T1 |> inner_join(r1u |> transmute(Replicate, sample = recode(as.character(scenario), Isotropic = "Symmetric"), generation, p_unr = p, r_unr = -r_obs),
                        by = c("Replicate", "sample", "generation"))
tie_changed <- tie |> filter((p < .05) != (p_unr < .05)) |> select(sample, Replicate, generation, r_unr, r_S, p_unr, p)
tie_sum <- c(censuses = nrow(tie), r_changed = sum(abs(tie$r_S - tie$r_unr) > 1e-12), p_changed = sum(abs(tie$p - tie$p_unr) > 1e-12), rej_changed = nrow(tie_changed))

## ---- reproduction of the cached manuscript inputs by the gamma = 1 grid ---------------------------
repro <- c(
  dbar = with(inner_join(fwd |> mutate(sample = recode(Scenario, Isotropic = "Symmetric")), G1$grid, by = c("sample", "Replicate", "generation")),
              max(abs(dbar.x - dbar.y), abs(mean_abs_delta.x - mean_abs_delta.y))),
  ibgd = { ef <- ld("data/derived/ibgd_trajectory_fit.rda")
    o <- bind_rows(ef$traj |> mutate(Treatment = as.character(Treatment)), rean$ibgd_sym$traj) |> filter(status == "ok") |>
      transmute(sample = recode(Treatment, Isotropic = "Symmetric"), Replicate, generation, r2c = r2_cgd, r2d = r2_pgd_dir, rf = rel_fwd)
    with(inner_join(o, G1$grid, by = c("sample", "Replicate", "generation")), max(abs(r2c - r2_cgd), abs(r2d - r2_pgd_dir), abs(rf - rel_fwd), na.rm = TRUE)) },
  truth = with(inner_join(rean$truth$d |> transmute(sample = recode(Treatment, Isotropic = "Symmetric"), Replicate, generation = Generation, wr = w_rho_tr),
                          G1$grid, by = c("sample", "Replicate", "generation")), max(abs(wr - w_rho_tr))),
  floor = with(inner_join(ld("data/derived/node_potential_test.rda")$census_tbl |> transmute(sample = recode(scenario, Isotropic = "Symmetric"), Replicate, generation, c0 = cor_phi_phi0),
                          T1, by = c("sample", "Replicate", "generation")), max(abs(c0 - cor_phi_phi0), na.rm = TRUE)))

## ---- information horizon -----------------------------------------------------------------------------
BREAKS <- c(0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 6, 8, 10)
fit_exp <- function(t, y) tryCatch(nls(y ~ yinf + (y0 - yinf) * exp(-t / tau), start = list(yinf = mean(y[t >= 1500]), y0 = y[1], tau = 300),
                                       algorithm = "port", lower = c(-Inf, -Inf, 1), upper = c(Inf, Inf, 1e5), control = nls.control(maxiter = 500)), error = function(e) NULL)
iso_fits <- function(grid) bind_rows(lapply(1:50, function(r) { d <- grid |> filter(sample == "Burn-in", Replicate == r) |> arrange(generation)
  m <- fit_exp(d$generation, d$mean_abs_delta); if (is.null(m)) return(NULL); cf <- coef(m)
  data.frame(Replicate = r, y0 = cf[["y0"]], yinf = cf[["yinf"]], tau = cf[["tau"]],
             R2 = 1 - sum(resid(m)^2) / sum((d$mean_abs_delta - mean(d$mean_abs_delta))^2)) }))
F1 <- iso_fits(G1$grid); F5 <- iso_fits(G5$grid)
horizon <- function(nullcsv, fits) {
  nl <- gc_read_information_null(nullcsv) |> filter(ok == 1); plat <- median(fits$yinf)
  nb <- nl |> mutate(bin = cut(n_poly_pop, BREAKS)) |> filter(!is.na(bin)) |> group_by(bin) |> summarise(abs_med = median(mean_abs_delta), .groups = "drop")
  infl <- nb |> filter(abs_med > 1.1 * plat)
  thr <- BREAKS[max(as.integer(infl$bin)) + 1]
  red <- fwd |> filter(Scenario == "Redistributed") |> group_by(Replicate) |> summarise(gen = min(c(generation[n_poly_pop < thr], NA), na.rm = TRUE), .groups = "drop")
  list(thr = thr, plateau = plat, bins = nb |> mutate(ratio = abs_med / plat), horizon_gen = median(red$gen[is.finite(red$gen)]))
}
H1 <- horizon("data/information_null.csv", F1); H5 <- horizon("data/information_null_g05.csv", F5)

## ---- assemble rows -----------------------------------------------------------------------------------
rows <- list()
add <- function(item, loc, old, new, type = "pct", txt_old = NA_character_, txt_new = NA_character_, group = "")
  rows[[length(rows) + 1]] <<- data.frame(group = group, item = item, location = loc, old = as.numeric(old)[1], new = as.numeric(new)[1], type = type,
                                          old_txt = txt_old, new_txt = txt_new, stringsAsFactors = FALSE)
both <- function(f) list(old = f(G1$grid, T1), new = f(G5$grid, T5))

## horizon
grp <- "Informative horizon"
add("Informative horizon (polymorphic loci per deme; null |Δ| > 1.1 × Symmetric plateau below it)", "index.md Phases; supplement Information loss", H1$thr, H5$thr, "stat2", group = grp)
add("Symmetric plateau of mean |Δ| (burn-in exponential fit, median y∞)", "supplement Information loss; Fig_SymmetricBaseline", H1$plateau, H5$plateau, "stat4", group = grp)
for (b in levels(H1$bins$bin)[levels(H1$bins$bin) %in% c("(3.5,4.5]", "(4.5,6]", "(2.75,3.5]")])
  add(sprintf("Null median |Δ| / plateau, polymorphism bin %s", b), "supplement Information loss", H1$bins$ratio[H1$bins$bin == b], H5$bins$ratio[H5$bins$bin == b], "stat2", group = grp)

## Symmetric baseline
grp <- "Symmetric baseline"
ci <- function(grid) { m <- grid |> filter(sample == "Symmetric") |> group_by(Replicate) |> summarise(m = mean(dbar)) |> pull(m); tt <- t.test(m); c(mean(m), tt$conf.int) }
b1 <- ci(G1$grid); b5 <- ci(G5$grid)
add("Symmetric forward-phase Δ̄ (mean of lineage means, 95% CI)", "index.md Symmetric baseline", b1[1], b5[1], "stat4",
    sprintf("%.4f (%.4f to %.4f)", b1[1], b1[2], b1[3]), sprintf("%.4f (%.4f to %.4f)", b5[1], b5[2], b5[3]), group = grp)
qq <- function(v, d) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), median(v), quantile(v, .025), quantile(v, .975))
for (v in c("y0", "yinf", "tau")) add(sprintf("Mean |Δ| burn-in relaxation: %s, median (2.5–97.5%%) over lineages", v), "Fig_SymmetricBaseline; Tab_IsotropicBaseline",
                                      median(F1[[v]]), median(F5[[v]]), if (v == "tau") "count" else "stat4",
                                      qq(F1[[v]], if (v == "tau") 0 else 4), qq(F5[[v]], if (v == "tau") 0 else 4), group = grp)
fd <- function(grid, fits) { m <- grid |> filter(sample == "Symmetric") |> group_by(Replicate) |> summarise(m = mean(mean_abs_delta)) |>
  inner_join(fits, by = "Replicate"); 100 * (m$m - m$yinf) / m$yinf }
add("Mean |Δ| forward deviation from own plateau (%), median", "Tab_IsotropicBaseline", median(fd(G1$grid, F1)), median(fd(G5$grid, F5)), "stat2", group = grp)
symabs <- function(grid) grid |> filter(sample == "Symmetric") |> group_by(Replicate) |> summarise(m = mean(mean_abs_delta)) |> pull(m) |> median()
add("Symmetric mean |Δ| (median of lineage means, forward)", "index.md Consequences", symabs(G1$grid), symabs(G5$grid), "stat3", group = grp)
lean_r <- function(t, s) { m <- t |> filter(sample == s) |> group_by(Replicate) |> summarise(m = mean(r_S)) |> pull(m)
  c(k = sum(m > 0), n = length(m), p = binom.test(sum(m > 0), sum(m != 0))$p.value) }
lean_d <- function(grid, s) { m <- grid |> filter(sample == s, if (s == "Burn-in") generation >= 1904 else TRUE) |> group_by(Replicate) |> summarise(m = mean(dbar)) |> pull(m)
  c(k = sum(m < 0), n = length(m), p = binom.test(sum(m < 0), sum(m != 0))$p.value) }
lt <- function(v) sprintf("%d of %d (p = %.3f)", v[["k"]], v[["n"]], v[["p"]])
for (s in SYMS) {
  o <- lean_r(T1, s); n <- lean_r(T5, s)
  add(sprintf("Lineage lean in r(S,x), %s: lineages with mean r > 0 (sign test)", s), "index.md Symmetric baseline; supplement Differentiation without Direction", o[["k"]], n[["k"]], "count", lt(o), lt(n), group = grp)
  o <- lean_d(G1$grid, s); n <- lean_d(G5$grid, s)
  add(sprintf("Lineage lean in Δ̄, %s: lineages with mean Δ̄ < 0 (sign test)", s), "index.md Symmetric baseline", o[["k"]], n[["k"]], "count", lt(o), lt(n), group = grp)
}
add("Mirroring check (γ = ½; 50 generation-1999 graphs rebuilt with deme labels reversed): max |r + r_mirror|", "index.md Symmetric baseline",
    NA, max(abs(MIR$r + MIR$r_mirror)), "sci", "reverses exactly at γ = 1 (current text)",
    sprintf("%.1e; Δ̄ max |Δ̄ + Δ̄_mirror| %.1e", max(abs(MIR$r + MIR$r_mirror)), max(abs(MIR$dbar + MIR$dbar_mirror))), group = grp)

## topological floor
grp <- "Topological floor"
fl <- function(t) t |> group_by(sample) |> summarise(c = median(cor_phi_phi0, na.rm = TRUE), leaf = mean(has_leaf), src = mean(src_leaf), sink = mean(sink_leaf),
                                                      rint = cor(r_S, r_S_int, use = "complete.obs"), dint = median(abs(r_S - r_S_int), na.rm = TRUE), .groups = "drop")
FL <- list(old = fl(T1), new = fl(T5))
fv <- function(d, s, v) d[[v]][d$sample == s]
for (s in c("Burn-in", "Symmetric", "Redistributed", "Obstructed", "sym-mid", "sym-low", "sym-verylow"))
  add(sprintf("cor(S, S⁰-potential), %s: median over censuses", s), "index.md topological floor (≈0.8 / 0.49 / 0.61)", fv(FL$old, s, "c"), fv(FL$new, s, "c"), "stat", group = grp)
for (s in c("Symmetric", "sym-verylow", "Redistributed", "Obstructed")) {
  add(sprintf("Leaves: %s censuses with a leaf", s), "index.md topological floor", fv(FL$old, s, "leaf"), fv(FL$new, s, "leaf"), group = grp)
  add(sprintf("Leaves: %s censuses where a leaf is the strongest source", s), "index.md topological floor", fv(FL$old, s, "src"), fv(FL$new, s, "src"), group = grp)
  add(sprintf("Leaves: %s censuses where a leaf is the strongest sink", s), "index.md topological floor", fv(FL$old, s, "sink"), fv(FL$new, s, "sink"), group = grp)
}
for (s in c("Symmetric", "Redistributed", "Obstructed"))
  add(sprintf("Interior-only r(S,x), %s: correlation with full r / median |difference|", s), "index.md topological floor", fv(FL$old, s, "rint"), fv(FL$new, s, "rint"), "stat",
      sprintf("%.3f / %.3f", fv(FL$old, s, "rint"), fv(FL$old, s, "dint")), sprintf("%.3f / %.3f", fv(FL$new, s, "rint"), fv(FL$new, s, "dint")), group = grp)

## source-sink map
grp <- "Source–sink map"
hm <- function(heat) heat |> group_by(scenario, generation, deme) |> summarise(S = median(S), k = median(k), .groups = "drop")
H <- list(old = hm(G1$heat), new = hm(G5$heat))
kterm <- G1$heat |> filter(scenario == "Symmetric") |> group_by(term = deme %in% c(1, 2, 24, 25)) |> summarise(k = median(k))
add("Terminal-deme degree (demes 1–2, 24–25) vs interior, Symmetric (bandwidth-free)", "index.md source–sink map", NA, NA, "txt",
    sprintf("%g vs %g", kterm$k[kterm$term], kterm$k[!kterm$term]), sprintf("%g vs %g", kterm$k[kterm$term], kterm$k[!kterm$term]), group = grp)
band <- function(h, s, lo, hi) { x <- h |> filter(scenario == s, generation >= lo, generation < hi) |> group_by(deme) |> summarise(S = median(S)); x$deme[which.max(x$S)] }
for (s in AS) add(sprintf("%s: deme with the highest median S, generations 2300–2600", s), "index.md source–sink map (demes 4–5)", band(H$old, s, 2300, 2600), band(H$new, s, 2300, 2600), "count", group = grp)
for (w in list(c(2500, 2600), c(2600, 2700), c(2700, 2800), c(2800, 2900), c(2900, 3000)))
  add(sprintf("Redistributed: highest-S deme, generations %d–%d", w[1], w[2]), "index.md source–sink map (7 / 9 / 11)", band(H$old, "Redistributed", w[1], w[2]),
      band(H$new, "Redistributed", w[1], w[2]), "count", group = grp)
for (s in AS) add(sprintf("%s: deme with the lowest median S, generations 2300–2600", s), "index.md source–sink map",
                  { x <- H$old |> filter(scenario == s, generation >= 2300, generation < 2600) |> group_by(deme) |> summarise(S = median(S)); x$deme[which.min(x$S)] },
                  { x <- H$new |> filter(scenario == s, generation >= 2300, generation < 2600) |> group_by(deme) |> summarise(S = median(S)); x$deme[which.min(x$S)] }, "count", group = grp)
up <- function(t, g) mean(t$r_S[t$sample == "Redistributed" & t$generation == g] < 0)
for (g in c(2704, 2804, 2904, 2954)) add(sprintf("Redistributed censuses with sources upstream at generation %d", g), "index.md source–sink map (≈100% → 50%)", up(T1, g), up(T5, g), group = grp)
add("Symmetric forward median r(S,x)", "index.md source–sink map (about 0.10)", median(T1$r_S[T1$sample == "Symmetric"]), median(T5$r_S[T5$sample == "Symmetric"]), "stat", group = grp)

## descriptive by phase
grp <- "Graph metrics by phase"
desc <- function(grid) grid |> filter(sample %in% AS) |> ph() |> group_by(sample, phase, Replicate) |>
  summarise(across(c(mean_abs_delta, dbar), mean), .groups = "drop") |> group_by(sample, phase) |>
  summarise(lin = n(), abs = median(mean_abs_delta), pos = mean(dbar > 0), dbar = median(dbar), .groups = "drop")
D <- list(old = desc(G1$grid), new = desc(G5$grid))
dv <- function(d, s, p, v) d[[v]][d$sample == s & d$phase == p]
for (s in AS) for (p in names(PH)) {
  if (s == "Obstructed" && p == "Fall apart") next
  add(sprintf("%s %s: mean |Δ| (median of lineage means)", s, PH[p]), "index.md Consequences", dv(D$old, s, p, "abs"), dv(D$new, s, p, "abs"), "stat3", group = grp)
  add(sprintf("%s %s: Δ̄ positive in (lineages)", s, PH[p]), "index.md Consequences; supplement Directional balance", dv(D$old, s, p, "pos"), dv(D$new, s, p, "pos"), group = grp)
  add(sprintf("%s %s: Δ̄ median", s, PH[p]), "index.md Consequences; supplement Directional balance", dv(D$old, s, p, "dbar"), dv(D$new, s, p, "dbar"), "stat4", group = grp)
}
add("Redistributed decay: Δ̄ negative in (lineages)", "index.md Consequences; supplement", 1 - dv(D$old, "Redistributed", "Fall apart", "pos"), 1 - dv(D$new, "Redistributed", "Fall apart", "pos"), group = grp)

## source-sink test: calibration and power
grp <- "Source–sink test"
cal <- function(t) t |> filter(sample %in% SYMS) |> group_by(sample) |>
  summarise(a01 = mean(p < .01), a05 = mean(p < .05), a10 = mean(p < .10), n0 = mean(p_N0 < .05), nt = mean(p_N1true < .05),
            dec = { b <- table(cut(p, seq(0, 1, .1), include.lowest = TRUE)) / n(); sprintf("%.0f–%.0f%%", 100 * min(b), 100 * max(b)) }, .groups = "drop")
C <- list(old = cal(T1), new = cal(T5)); cv <- function(d, s, v) d[[v]][d$sample == s]
for (s in SYMS) {
  add(sprintf("Tab_SourceSinkCalibration %s: graph spectral null α = 0.01 / 0.05 / 0.10", s), "index.md Tab_SourceSinkCalibration", cv(C$old, s, "a05"), cv(C$new, s, "a05"), "pct",
      sprintf("%.3f / %.3f / %.3f", cv(C$old, s, "a01"), cv(C$old, s, "a05"), cv(C$old, s, "a10")), sprintf("%.3f / %.3f / %.3f", cv(C$new, s, "a01"), cv(C$new, s, "a05"), cv(C$new, s, "a10")), group = grp)
  add(sprintf("Tab_SourceSinkCalibration %s: plain permutation, 0.05", s), "index.md Tab_SourceSinkCalibration", cv(C$old, s, "n0"), cv(C$new, s, "n0"), group = grp)
  add(sprintf("Tab_SourceSinkCalibration %s: stepping-stone spectral null, 0.05", s), "index.md Tab_SourceSinkCalibration", cv(C$old, s, "nt"), cv(C$new, s, "nt"), group = grp)
}
add("Burn-in p-value deciles (range of shares)", "index.md source–sink calibration (8–12%)", NA, NA, "txt", cv(C$old, "Burn-in", "dec"), cv(C$new, "Burn-in", "dec"), group = grp)
mc <- function(t) { x <- t |> filter(sample %in% SYMS); c(all = mean(x$p < .05), excl = mean(x$p[x$n_comp == 1] < .05), ncomp = sum(x$n_comp > 1)) }
add("Censuses with > 1 component (symmetric samples) / largest effect of excluding them", "index.md calibration paragraph", mc(T1)[["ncomp"]], mc(T5)[["ncomp"]], "count",
    sprintf("%d; %.2f pp", mc(T1)[["ncomp"]], 100 * abs(mc(T1)[["all"]] - mc(T1)[["excl"]])), sprintf("%d; %.2f pp", mc(T5)[["ncomp"]], 100 * abs(mc(T5)[["all"]] - mc(T5)[["excl"]])), group = grp)
sst <- function(t, rcol = "phase_policy_R_gen") { d <- t |> filter(sample %in% c("Symmetric", AS)) |> ph(rcol = rcol) |> mutate(phase = ifelse(sample == "Symmetric", "all", phase))
  d |> group_by(sample, phase) |> summarise(cens = n(), lin = n_distinct(Replicate), rej = mean(p < .05), medr = median(r_S), up = mean(r_S < 0), .groups = "drop") |>
    left_join(d |> group_by(sample, phase, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |> group_by(sample, phase) |>
                summarise(lin_up = mean(m < 0), nlin_up = sum(m < 0), .groups = "drop"), by = c("sample", "phase")) }
S <- list(old = sst(T1), new = sst(T5)); sv <- function(d, s, p, v) d[[v]][d$sample == s & d$phase == p]
SSROWS <- c(list(c("Symmetric", "all")), lapply(names(PH), function(p) c("Redistributed", p)), lapply(names(PH)[1:2], function(p) c("Obstructed", p)))
lab <- function(s, p) if (p == "all") (if (s == "Burn-in") "Burn-in 1904–1999" else paste(s, "forward phase")) else paste(s, PH[p])
for (sp in SSROWS) { s <- sp[1]; p <- sp[2]; loc <- "index.md source–sink paragraph; supplement Tab_SourceSinkPower"
  add(sprintf("Tab_SourceSinkPower %s: censuses (lineages)", lab(s, p)), loc, sv(S$old, s, p, "cens"), sv(S$new, s, p, "cens"), "count",
      sprintf("%d (%d)", sv(S$old, s, p, "cens"), sv(S$old, s, p, "lin")), sprintf("%d (%d)", sv(S$new, s, p, "cens"), sv(S$new, s, p, "lin")), group = grp)
  add(sprintf("Tab_SourceSinkPower %s: rejection rate", lab(s, p)), loc, sv(S$old, s, p, "rej"), sv(S$new, s, p, "rej"), group = grp)
  add(sprintf("Tab_SourceSinkPower %s: median r(S,x)", lab(s, p)), loc, sv(S$old, s, p, "medr"), sv(S$new, s, p, "medr"), "stat", group = grp)
  add(sprintf("Tab_SourceSinkPower %s: censuses with sources upstream", lab(s, p)), loc, sv(S$old, s, p, "up"), sv(S$new, s, p, "up"), group = grp)
  add(sprintf("Tab_SourceSinkPower %s: lineages with sources upstream", lab(s, p)), loc, sv(S$old, s, p, "lin_up"), sv(S$new, s, p, "lin_up"), "pct",
      sprintf("%.0f%% (%d of %d)", 100 * sv(S$old, s, p, "lin_up"), sv(S$old, s, p, "nlin_up"), sv(S$old, s, p, "lin")),
      sprintf("%.0f%% (%d of %d)", 100 * sv(S$new, s, p, "lin_up"), sv(S$new, s, p, "nlin_up"), sv(S$new, s, p, "lin")), group = grp)
}

## directional IBGD (5-generation grid)
grp <- "Directional isolation (IBGD)"
ibd <- function(grid, rcol = "phase_policy_R_gen") grid |> filter(!is.na(r2_pgd_dir), generation >= 1904) |> ph(rcol = rcol) |>
  mutate(phase = ifelse(sample %in% AS, phase, "all"), det = r2_pgd_dir > r2_cgd)
ibg <- function(d) d |> group_by(sample, phase) |>
  summarise(cens = n(), det_rate = mean(det), lo = binom.test(sum(det), n())$conf.int[1], hi = binom.test(sum(det), n())$conf.int[2],
            med_dR2 = median(r2_pgd_dir - r2_cgd), lt1 = mean(rel_fwd[det] < 1), gt1 = mean(rel_fwd[det] > 1), med_ratio = median(rel_fwd[det]),
            ratio_all = median(rel_fwd), ratio_q1 = quantile(rel_fwd, .25), ratio_q3 = quantile(rel_fwd, .75),
            cgd = median(r2_cgd), blind = median(r2_pgd_blind), dir = median(r2_pgd_dir), gain = median(r2_pgd_dir - r2_pgd_blind),
            gain_pos = mean(r2_pgd_dir > r2_pgd_blind), posthoc = mean(det & rel_fwd < 1),
            lin_max = max(tapply(det, Replicate, mean)), .groups = "drop")
ID <- list(old = ibd(G1$grid), new = ibd(G5$grid)); I <- lapply(ID, ibg); iv <- function(d, s, p, v) { x <- d[[v]][d$sample == s & d$phase == p]; if (length(x)) x else NA }
IBROWS <- c(lapply(SYMS, function(s) c(s, "all")), lapply(names(PH), function(p) c("Redistributed", p)), lapply(names(PH), function(p) c("Obstructed", p)))
for (sp in IBROWS) { s <- sp[1]; p <- sp[2]; nm <- lab(s, p); L1 <- "Tab_IBGDFitComparison; index.md IBGD text"; L2 <- "Tab_IBGDDetection; index.md ΔR² paragraphs"
  add(sprintf("IBGD %s: censuses", nm), "Tab_IBGDDetection", iv(I$old, s, p, "cens"), iv(I$new, s, p, "cens"), "count", group = grp)
  add(sprintf("IBGD %s: cGD R² (median)", nm), L1, iv(I$old, s, p, "cgd"), iv(I$new, s, p, "cgd"), "stat", group = grp)
  add(sprintf("IBGD %s: pGD direction-blind R²", nm), L1, iv(I$old, s, p, "blind"), iv(I$new, s, p, "blind"), "stat", group = grp)
  add(sprintf("IBGD %s: pGD direction-aware R²", nm), L1, iv(I$old, s, p, "dir"), iv(I$new, s, p, "dir"), "stat", group = grp)
  add(sprintf("IBGD %s: directional gain (aware − blind)", nm), "index.md IBGD text; Fig_DirectionalIBGD", iv(I$old, s, p, "gain"), iv(I$new, s, p, "gain"), "stat", group = grp)
  add(sprintf("IBGD %s: ΔR² > 0 (%s; 95%% CI)", nm, if (p == "all") "false-positive rate" else "detection"), L2, iv(I$old, s, p, "det_rate"), iv(I$new, s, p, "det_rate"), "pct",
      sprintf("%.1f%% (%.1f–%.1f)", 100 * iv(I$old, s, p, "det_rate"), 100 * iv(I$old, s, p, "lo"), 100 * iv(I$old, s, p, "hi")),
      sprintf("%.1f%% (%.1f–%.1f)", 100 * iv(I$new, s, p, "det_rate"), 100 * iv(I$new, s, p, "lo"), 100 * iv(I$new, s, p, "hi")), group = grp)
  add(sprintf("IBGD %s: median ΔR²", nm), L2, iv(I$old, s, p, "med_dR2"), iv(I$new, s, p, "med_dR2"), "stat3", group = grp)
  add(sprintf("IBGD %s: detections with slope ratio < 1", nm), L2, iv(I$old, s, p, "lt1"), iv(I$new, s, p, "lt1"), group = grp)
  add(sprintf("IBGD %s: median slope ratio among detections", nm), L2, iv(I$old, s, p, "med_ratio"), iv(I$new, s, p, "med_ratio"), "stat", group = grp)
  if (p == "all") {
    add(sprintf("IBGD %s: post hoc rule rate (ΔR² > 0 and slope ratio < 1)", nm), "index.md post hoc rule (1.4–5.2%)", iv(I$old, s, p, "posthoc"), iv(I$new, s, p, "posthoc"), group = grp)
    add(sprintf("IBGD %s: highest per-lineage false-positive rate", nm), "index.md (39–74%)", iv(I$old, s, p, "lin_max"), iv(I$new, s, p, "lin_max"), group = grp)
    add(sprintf("IBGD %s: detections with slope ratio > 1", nm), "index.md (59–67%)", iv(I$old, s, p, "gt1"), iv(I$new, s, p, "gt1"), group = grp)
  } else add(sprintf("IBGD %s: censuses favouring cGD", nm), "index.md one-sided paragraph", 1 - iv(I$old, s, p, "det_rate"), 1 - iv(I$new, s, p, "det_rate"), group = grp)
}
add("IBGD Symmetric: slope ratio over all censuses, median (IQR)", "index.md (1.11, IQR 0.87–1.40)", iv(I$old, "Symmetric", "all", "ratio_all"), iv(I$new, "Symmetric", "all", "ratio_all"), "stat",
    sprintf("%.2f (%.2f–%.2f)", iv(I$old, "Symmetric", "all", "ratio_all"), iv(I$old, "Symmetric", "all", "ratio_q1"), iv(I$old, "Symmetric", "all", "ratio_q3")),
    sprintf("%.2f (%.2f–%.2f)", iv(I$new, "Symmetric", "all", "ratio_all"), iv(I$new, "Symmetric", "all", "ratio_q1"), iv(I$new, "Symmetric", "all", "ratio_q3")), group = grp)
add("IBGD Symmetric: within-pGD gain positive (aware > blind)", "index.md (94%)", iv(I$old, "Symmetric", "all", "gain_pos"), iv(I$new, "Symmetric", "all", "gain_pos"), group = grp)
a99 <- function(d) { x <- d |> filter(sample %in% AS, det); mean(x$rel_fwd < 1) }
add("IBGD: detections under asymmetric migration with slope ratio < 1 (all phases)", "index.md (99%)", a99(ID$old), a99(ID$new), group = grp)

## test overlap (50-generation grid)
grp <- "Test overlap"
ovd <- function(grid, t) inner_join(t |> select(sample, Replicate, generation, r_S, p), grid |> select(sample, Replicate, generation, r2_cgd, r2_pgd_dir, rel_fwd),
                                    by = c("sample", "Replicate", "generation")) |> filter(!is.na(r2_pgd_dir)) |>
  ph() |> mutate(phase = ifelse(sample %in% AS, phase, "all"), ss = p < .05, ib = r2_pgd_dir > r2_cgd)
ovl <- function(d) d |> group_by(sample, phase) |> summarise(cens = n(), both = sum(ss & ib), sso = sum(ss & !ib), ibo = sum(!ss & ib), neither = sum(!ss & !ib),
                                                            boe = both / max(1, both + sso + ibo), same = sum((r_S < 0 & rel_fwd < 1)[ss & ib]), .groups = "drop")
OD <- list(old = ovd(G1$grid, T1), new = ovd(G5$grid, T5)); O <- lapply(OD, ovl); ov <- function(d, s, p, v) d[[v]][d$sample == s & d$phase == p]
for (sp in c(lapply(names(PH), function(p) c("Redistributed", p)), lapply(names(PH), function(p) c("Obstructed", p)), list(c("Symmetric", "all"), c("Burn-in", "all")))) {
  s <- sp[1]; p <- sp[2]; if (!length(ov(O$old, s, p, "cens"))) next
  f5 <- function(d) with(d[d$sample == s & d$phase == p, ], sprintf("%d / %d / %d / %d / %d", cens, both, sso, ibo, neither))
  add(sprintf("Tab_TestOverlap %s: censuses / both / S-only / ΔR²-only / neither", lab(s, p)), "supplement Tab_TestOverlap", ov(O$old, s, p, "both"), ov(O$new, s, p, "both"), "count", f5(O$old), f5(O$new), group = grp)
  add(sprintf("Tab_TestOverlap %s: both / either", lab(s, p)), "index.md + Tab_TestOverlap", ov(O$old, s, p, "boe"), ov(O$new, s, p, "boe"), group = grp)
  add(sprintf("Tab_TestOverlap %s: same direction among joint", lab(s, p)), "Tab_TestOverlap", ov(O$old, s, p, "same"), ov(O$new, s, p, "same"), "count",
      sprintf("%d/%d", ov(O$old, s, p, "same"), ov(O$old, s, p, "both")), sprintf("%d/%d", ov(O$new, s, p, "same"), ov(O$new, s, p, "both")), group = grp)
}
jt <- function(o) { x <- o |> filter(sample %in% AS, !(sample == "Obstructed" & phase == "Fall apart")); c(sum(x$same), sum(x$both)) }
add("Joint detections agreeing on direction (ascent, plateau, Redistributed decay)", "index.md (274 of 274)", jt(O$old)[2], jt(O$new)[2], "count",
    sprintf("%d of %d", jt(O$old)[1], jt(O$old)[2]), sprintf("%d of %d", jt(O$new)[1], jt(O$new)[2]), group = grp)

## divMigrate comparison
grp <- "divMigrate comparison"
nm <- rean$truth$d |> transmute(sample = recode(Treatment, Isotropic = "Symmetric"), Replicate, generation = Generation, Exclusion, nm_rho_run, nm_rho_tr, nm_rmse_run)
nc <- rean$corridor$d |> transmute(sample = recode(Treatment, Isotropic = "Symmetric"), Replicate, generation = Generation, nm_ok, nm_fwd_gt)
trd <- function(grid) grid |> filter(sample %in% c("Symmetric", AS)) |> inner_join(nm, by = c("sample", "Replicate", "generation")) |>
  inner_join(nc, by = c("sample", "Replicate", "generation")) |> ph() |> mutate(phase = ifelse(sample == "Symmetric", "all", phase))
TD <- list(old = trd(G1$grid), new = trd(G5$grid))
tru <- function(d) d |> group_by(sample, phase) |>
  summarise(w = median(w_rho_tr), w_p = median(w_rho_tr[Exclusion == ""]), nmr = median(nm_rho_run, na.rm = TRUE), wbeat = mean((w_rho_tr > nm_rho_run)[Exclusion == ""]),
            wf = mean(w_fwd_gt, na.rm = TRUE), wf_p = mean(w_fwd_gt[nm_ok], na.rm = TRUE), nf = mean(nm_fwd_gt, na.rm = TRUE),
            tr_shift = abs(median(w_rho_tr) - median(w_rho_run)), .groups = "drop")
TR <- lapply(TD, tru); tv <- function(d, s, p, v) d[[v]][d$sample == s & d$phase == p]
add("Truth ρ, Symmetric: w (Nm 0.45, bandwidth-free)", "index.md divMigrate", tv(TR$old, "Symmetric", "all", "w"), tv(TR$new, "Symmetric", "all", "w"), "stat", group = grp)
add("Transposing the truth: largest change in w's median ρ (Symmetric / Redistributed / Obstructed, all censuses)", "index.md divMigrate (≤ 0.001)",
    max(TR$old$tr_shift), max(TR$new$tr_shift), "stat4", group = grp)
for (s in AS) {
  add(sprintf("Truth ρ, %s plateau: w (Nm %.2f)", s, tv(TR$old, s, "Plateau", "nmr")), "index.md divMigrate", tv(TR$old, s, "Plateau", "w"), tv(TR$new, s, "Plateau", "w"), "stat", group = grp)
  add(sprintf("%s plateau: w higher than Nm (paired censuses)", s), "index.md divMigrate (99%)", tv(TR$old, s, "Plateau", "wbeat"), tv(TR$new, s, "Plateau", "wbeat"), group = grp)
  for (p in c("Reorganize", "Plateau")) add(sprintf("Corridor direction, %s %s: w, same censuses as Nm (Nm %.0f%%)", s, PH[p], 100 * tv(TR$old, s, p, "nf")),
                                            "index.md divMigrate", tv(TR$old, s, p, "wf_p"), tv(TR$new, s, p, "wf_p"), group = grp)
  lb <- function(d) { x <- d |> filter(sample == s, phase == "Plateau", nm_ok) |> group_by(Replicate) |> summarise(dd = mean(nm_fwd_gt - w_fwd_gt), .groups = "drop"); c(sum(x$dd > 0), nrow(x)) }
  add(sprintf("Corridor direction, %s plateau: lineages where Nm is better", s), "index.md divMigrate", lb(TD$old)[1], lb(TD$new)[1], "count",
      sprintf("%d of %d", lb(TD$old)[1], lb(TD$old)[2]), sprintf("%d of %d", lb(TD$new)[1], lb(TD$new)[2]), group = grp)
}
add("Corridor direction, Symmetric: w", "index.md divMigrate (50%)", tv(TR$old, "Symmetric", "all", "wf"), tv(TR$new, "Symmetric", "all", "wf"), group = grp)
add("Redistributed decay: truth ρ of w, all censuses", "index.md divMigrate (0.31)", tv(TR$old, "Redistributed", "Fall apart", "w"), tv(TR$new, "Redistributed", "Fall apart", "w"), "stat", group = grp)
add(sprintf("Redistributed decay: truth ρ of w, censuses where Nm is defined (Nm %.2f)", tv(TR$old, "Redistributed", "Fall apart", "nmr")), "index.md divMigrate (0.43 vs 0.33)",
    tv(TR$old, "Redistributed", "Fall apart", "w_p"), tv(TR$new, "Redistributed", "Fall apart", "w_p"), "stat", group = grp)
late <- function(grid) { x <- grid |> filter(sample == "Redistributed", generation >= 2879) |> group_by(generation) |> summarise(m = median(w_rho_tr)); sprintf("%.2f–%.2f", min(x$m), max(x$m)) }
add("Redistributed w ρ, all censuses, over the final 25 generations (range of medians)", "Fig_DivMigrateTruth caption (0.14–0.22)", NA, NA, "txt", late(G1$grid), late(G5$grid), group = grp)

## Delta-bar sign change after the horizon
grp <- "Δ̄ sign change (Redistributed)"
cross <- function(g, y) { gg <- seq(2254, 2954, 5)
  f <- as.numeric(predict(gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML"), data.frame(g = gg)))
  s <- sign(f); i <- which(s[-1] != s[-length(s)] & gg[-1] > 2400); if (length(i)) gg[i[1] + 1] else NA_real_ }
comps <- c("A", "A_grad", "A_cycle", "A_ends", "A_interior", "A_int_edges")
zcf <- function(grid) grid |> filter(sample %in% AS, generation >= 2254, generation <= 2954) |> group_by(sample, Replicate) |>
  summarise(across(all_of(comps), ~ cross(generation, .x)), .groups = "drop")
ZC <- list(old = zcf(G1$grid), new = zcf(G5$grid))
qs <- function(v) { v <- v[!is.na(v)]; if (!length(v)) "—" else sprintf("%d — %.0f (%.0f–%.0f)", length(v), median(v), quantile(v, .25), quantile(v, .75)) }
for (s in AS) for (cp in comps) add(sprintf("Crossing of %s, %s: lineages — median (IQR) generation", c(A = "|E|Δ̄", A_grad = "part reproduced by S", A_cycle = "non-directional residual",
                                       A_ends = "terminal demes", A_interior = "interior", A_int_edges = "interior edges only")[cp], s),
                                    "supplement Tab_DeltaBarCrossings; index.md (39 lineages, 2684)", sum(!is.na(ZC$old[[cp]][ZC$old$sample == s])),
                                    sum(!is.na(ZC$new[[cp]][ZC$new$sample == s])), "count", qs(ZC$old[[cp]][ZC$old$sample == s]), qs(ZC$new[[cp]][ZC$new$sample == s]), group = grp)
coin <- function(z) { x <- z |> filter(sample == "Redistributed", !is.na(A))
  sapply(c("A_int_edges", "A_interior", "A_cycle", "A_ends"), function(cp) mean(!is.na(x[[cp]]) & abs(x[[cp]] - x$A) <= 50)) }
co <- list(old = coin(ZC$old), new = coin(ZC$new))
for (cp in names(co$old)) add(sprintf("Redistributed crossings of |E|Δ̄ coinciding (± 50 generations) with %s", c(A_int_edges = "interior-edge sum", A_interior = "interior part of reproduced sum",
                                                                                                                A_cycle = "non-directional residual", A_ends = "terminal-deme terms")[cp]),
                              "index.md + supplement (92% / 59% / 26% / 0%)", co$old[[cp]], co$new[[cp]], group = grp)
dec_tab <- function(grid) grid |> filter(sample == "Redistributed", generation %in% c(2254, 2454, 2604, 2654, 2704, 2754, 2804, 2904, 2954)) |>
  group_by(generation) |> summarise(across(all_of(comps), median), .groups = "drop")
DT <- list(old = dec_tab(G1$grid), new = dec_tab(G5$grid))
for (gg in DT$old$generation) add(sprintf("Tab_DeltaBarDecomposition, generation %d: |E|Δ̄ / reproduced / residual / terminal / interior / interior edges", gg), "supplement Tab_DeltaBarDecomposition",
                                  DT$old$A[DT$old$generation == gg], DT$new$A[DT$new$generation == gg], "stat",
                                  do.call(sprintf, c("%.2f / %.2f / %.2f / %.2f / %.2f / %.2f", as.list(unlist(DT$old[DT$old$generation == gg, comps])))),
                                  do.call(sprintf, c("%.2f / %.2f / %.2f / %.2f / %.2f / %.2f", as.list(unlist(DT$new[DT$new$generation == gg, comps])))), group = grp)
obs_pos <- function(grid) { x <- grid |> filter(sample == "Obstructed", generation >= 2254, generation <= 2954) |> group_by(generation) |> summarise(A = median(A)); c(min(x$A), x$A[1], x$A[nrow(x)]) }
add("Obstructed median |E|Δ̄, 2254–2954: minimum (first → last)", "supplement (stays positive, 0.47 → 1.47)", obs_pos(G1$grid)[1], obs_pos(G5$grid)[1], "stat",
    do.call(sprintf, c("%.2f (%.2f → %.2f)", as.list(obs_pos(G1$grid)))), do.call(sprintf, c("%.2f (%.2f → %.2f)", as.list(obs_pos(G5$grid)))), group = grp)
trend <- function(grid, s) { d <- grid |> filter(sample == s) |> mutate(rep = factor(Replicate))
  m <- gam(dbar ~ s(generation, k = 20) + s(rep, bs = "re"), data = d, method = "REML"); gg <- 2000:2999
  f <- as.numeric(predict(m, data.frame(generation = gg, rep = d$rep[1]), exclude = "s(rep)"))
  i <- which(sign(f[-1]) != sign(f[-length(f)]) & gg[-1] > 2400); c(peak = gg[which.max(f)], cross = if (length(i)) gg[i[1] + 1] else NA, lo = min(f), hi = max(f)) }
TRD <- list(old = lapply(c(Redistributed = "Redistributed", Obstructed = "Obstructed", Symmetric = "Symmetric"), function(s) trend(G1$grid, s)),
            new = lapply(c(Redistributed = "Redistributed", Obstructed = "Obstructed", Symmetric = "Symmetric"), function(s) trend(G5$grid, s)))
add("Redistributed Δ̄ trend: peak generation", "index.md + supplement (2435)", TRD$old$Redistributed[["peak"]], TRD$new$Redistributed[["peak"]], "count", group = grp)
add("Redistributed Δ̄ trend: zero crossing", "supplement (near 2700)", TRD$old$Redistributed[["cross"]], TRD$new$Redistributed[["cross"]], "count", group = grp)
add("Obstructed Δ̄ trend: peak generation", "index.md + supplement (2820)", TRD$old$Obstructed[["peak"]], TRD$new$Obstructed[["peak"]], "count", group = grp)
add("Obstructed Δ̄ trend: zero crossing after 2400", "supplement (none)", TRD$old$Obstructed[["cross"]], TRD$new$Obstructed[["cross"]], "count", group = grp)
add("Symmetric Δ̄ trend range", "Fig_SignedAsymmetry caption (−0.0016 to −0.0010)", TRD$old$Symmetric[["lo"]], TRD$new$Symmetric[["lo"]], "stat4",
    sprintf("%.4f to %.4f", TRD$old$Symmetric[["lo"]], TRD$old$Symmetric[["hi"]]), sprintf("%.4f to %.4f", TRD$new$Symmetric[["lo"]], TRD$new$Symmetric[["hi"]]), group = grp)

## bandwidth section (old: current tables, reference gamma = 1; new: reference gamma = 1/2)
grp <- "Bandwidth section"
eb_old <- ep$bw_raw |> mutate(phase = as.character(lp_phase(Scenario, Replicate, generation)))
eb_new <- BW$edge_bw |> mutate(phase = as.character(lp_phase(Scenario, Replicate, generation)))
ebs <- function(d) d |> group_by(bandwidth, phase) |> summarise(cens = n(), concord = mean(concord), rs = mean(rho_signed, na.rm = TRUE),
                                                                 ra = mean(rho_abs, na.rm = TRUE), kept = mean(dbar_kept), .groups = "drop")
E0 <- ebs(eb_old); E5 <- ebs(eb_new)
cell <- function(d, b, v, f) { g <- function(p) { x <- d[[v]][d$bandwidth == b & d$phase == p]; if (length(x)) f(x) else "—" }
  sprintf("%s (%s, %s)", g("Plateau"), g("Reorganize"), g("Fall apart")) }
pf <- function(x) sprintf("%.0f%%", 100 * x); nf <- function(x) sprintf("%.2f", x)
MAP <- list(c("Fixed global", "global", "global"), c("Local × ½ (old) / local × ¼ = γ 0.25 (new)", "local x1/2", "local x1/4"),
            c("Local × 2 (old: γ 2) / local × 2 = local mean, γ 1 (new)", "local x2", "local x2"), c("Perplexity", "perplexity", "perplexity"))
for (m in MAP) for (v in list(c("concord", "same sign"), c("rs", "rank ρ signed"), c("ra", "rank ρ |Δ|"), c("kept", "Δ̄ sign kept"))) {
  f <- if (v[1] %in% c("concord", "kept")) pf else nf
  add(sprintf("Tab_BandwidthSensitivity %s: %s, plateau (ascent, decay)", m[1], v[2]), "index.md Tab_BandwidthSensitivity",
      E0[[v[1]]][E0$bandwidth == m[2] & E0$phase == "Plateau"], E5[[v[1]]][E5$bandwidth == m[3] & E5$phase == "Plateau"], if (v[1] %in% c("concord", "kept")) "pct" else "stat",
      cell(E0, m[2], v[1], f), cell(E5, m[3], v[1], f), group = grp) }
cc <- function(d) { x <- d |> filter(bandwidth == "global"); sprintf("%d / %d / %d", sum(x$phase == "Reorganize"), sum(x$phase == "Plateau"), sum(x$phase == "Fall apart")) }
add("Tab_BandwidthSensitivity censuses ascent / plateau / decay", "index.md caption", NA, NA, "txt", cc(eb_old), cc(eb_new), group = grp)
pfail <- function(d) { k <- d |> distinct(Scenario, Replicate, generation, phase); h <- d |> filter(bandwidth == "perplexity") |> distinct(Scenario, Replicate, generation)
  x <- anti_join(k, h, by = c("Scenario", "Replicate", "generation")); sprintf("%d plateau, %d decay", sum(x$phase == "Plateau"), sum(x$phase == "Fall apart")) }
add("Perplexity failures (edge-level table)", "index.md caption (31 plateau, 131 decay)", NA, NA, "txt", pfail(eb_old), pfail(eb_new), group = grp)
sb_old <- es$bwp |> mutate(sample = recode(scenario, Isotropic = "Symmetric")); sb_new <- BW$ss_bw |> mutate(sample = recode(scenario, Isotropic = "Symmetric"))
sbs <- function(d) { d <- d |> filter(!failed) |> ph() |> mutate(phase = ifelse(sample == "Symmetric", "all", phase))
  d |> group_by(bandwidth, sample, phase) |> summarise(rej = mean(p < .05), medr = median(r_S), .groups = "drop") |>
    left_join(d |> group_by(bandwidth, sample, phase, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |> group_by(bandwidth, sample, phase) |>
                summarise(lin_up = mean(m < 0), .groups = "drop"), by = c("bandwidth", "sample", "phase")) }
B0 <- sbs(sb_old); B5 <- sbs(sb_new)
SMAP <- list(c("Reference (old: local mean γ 1; new: γ ½)", "canonical", "reference"), c("Local × ½ (old: γ ½) / local × ¼ (new: γ 0.25)", "local x1/2", "local x1/4"),
             c("Local × 2 (old: γ 2) / local × 2 (new: γ 1, the local mean)", "local x2", "local x2"), c("Perplexity", "perplexity", "perplexity"), c("Fixed global", "global", "global"))
bv <- function(d, b, s, p, v) { x <- d[[v]][d$bandwidth == b & d$sample == s & d$phase == p]; if (length(x)) x else NA }
for (m in SMAP) {
  add(sprintf("Tab_BandwidthSensitivitySourceSink %s: Symmetric FPR", m[1]), "index.md Tab_BandwidthSensitivitySourceSink", bv(B0, m[2], "Symmetric", "all", "rej"), bv(B5, m[3], "Symmetric", "all", "rej"), group = grp)
  for (s in AS) add(sprintf("Tab_BandwidthSensitivitySourceSink %s: %s plateau rejection / median r / lineages upstream", m[1], s), "index.md Tab_BandwidthSensitivitySourceSink",
                    bv(B0, m[2], s, "Plateau", "rej"), bv(B5, m[3], s, "Plateau", "rej"), "pct",
                    sprintf("%.0f%% / %.2f / %.0f%%", 100 * bv(B0, m[2], s, "Plateau", "rej"), bv(B0, m[2], s, "Plateau", "medr"), 100 * bv(B0, m[2], s, "Plateau", "lin_up")),
                    sprintf("%.0f%% / %.2f / %.0f%%", 100 * bv(B5, m[3], s, "Plateau", "rej"), bv(B5, m[3], s, "Plateau", "medr"), 100 * bv(B5, m[3], s, "Plateau", "lin_up")), group = grp)
  add(sprintf("Tab_BandwidthSensitivitySourceSink %s: Redistributed decay lineages upstream", m[1]), "index.md bandwidth paragraph", bv(B0, m[2], "Redistributed", "Fall apart", "lin_up"),
      bv(B5, m[3], "Redistributed", "Fall apart", "lin_up"), group = grp)
}

## phase-delineation sensitivity (unflagged-lineage policy)
grp <- "Phase delineation"
pol <- function(grid, t) {
  r <- function(rcol) { a <- sst(t, rcol) |> filter(sample %in% AS) |> select(sample, phase, rej, up, lin_up)
    b <- ibg(ibd(grid, rcol) |> filter(sample %in% AS)) |> select(sample, phase, det_rate, lt1)
    full_join(a, b, by = c("sample", "phase")) }
  x <- full_join(r("phase_policy_R_gen"), r("R_gen"), by = c("sample", "phase"), suffix = c(".n", ".o"))
  d <- function(p) { y <- x |> filter(phase == p); max(abs(c(y$rej.n - y$rej.o, y$up.n - y$up.o, y$lin_up.n - y$lin_up.o, y$det_rate.n - y$det_rate.o, y$lt1.n - y$lt1.o)), na.rm = TRUE) }
  c(plateau = d("Plateau"), ascent = d("Reorganize")) }
P1 <- pol(G1$grid, T1); P5 <- pol(G5$grid, T5)
add("Largest change in a plateau rate when unflagged lineages are phased from their own 90% point", "supplement Phase delineation (2.4 points)", P1[["plateau"]], P5[["plateau"]], "pp", group = grp)
add("Largest change in an ascent rate, same comparison", "supplement Phase delineation (3.4 points)", P1[["ascent"]], P5[["ascent"]], "pp", group = grp)

## ---- table -------------------------------------------------------------------------------------------
tab <- bind_rows(rows)
fmt <- function(x, type) ifelse(is.na(x), "—", ifelse(type == "pct", sprintf("%.1f%%", 100 * x), ifelse(type == "stat", sprintf("%.2f", x),
         ifelse(type == "stat2", sprintf("%.2f", x), ifelse(type == "stat3", sprintf("%.3f", x), ifelse(type == "stat4", sprintf("%.4f", x),
         ifelse(type == "pp", sprintf("%.1f pp", 100 * x), ifelse(type == "sci", sprintf("%.1e", x), sprintf("%g", round(x, 1))))))))))
tab <- tab |> mutate(old_s = ifelse(is.na(old_txt), fmt(old, type), old_txt), new_s = ifelse(is.na(new_txt), fmt(new, type), new_txt),
  moved = case_when(type %in% c("txt") ~ ifelse(old_s != new_s, "changed", "no"),
                    is.na(old) | is.na(new) ~ ifelse(old_s != new_s, "changed", "no"),
                    type %in% c("pct", "pp") ~ ifelse(abs(new - old) > 0.01, "**yes**", "no"),
                    type %in% c("stat", "stat2", "stat3", "stat4") ~ ifelse(abs(new - old) > 0.01 | (type == "stat4" & sign(new) != sign(old) & abs(old) > 1e-4), "**yes**", "no"),
                    TRUE ~ ifelse(old_s != new_s, "changed", "no")))
save(tab, tie_changed, tie_sum, repro, H1, H5, F1, F5, D, S, I, O, TR, ZC, DT, TRD, E0, E5, B0, B5, P1, P5, C, FL, file = "data/derived/half_bandwidth_results.rda")
md <- c("<!-- TABLE BELOW IS WRITTEN BY R/analysis/half_bandwidth_results.R -->", "",
  sprintf("Reproduction: the γ = 1 run reproduces the cached manuscript inputs (max |diff|: Δ̄ and mean |Δ| %.1g; IBGD %.1g; truth ρ %.1g; cor(φ, φ⁰) %.1g).",
          repro[["dbar"]], repro[["ibgd"]], repro[["truth"]], repro[["floor"]]), "",
  sprintf("γ = 1 source_sink_test() (S rounded to 12 significant digits, gstudio %s null draws) against the round-1 test (S unrounded, earlier draws), %d censuses: r(S,x) changed in %d (tied scores), rejection at α = 0.05 changed in %d (Monte Carlo error in p).",
          as.character(packageVersion("gstudio")), tie_sum[["censuses"]], tie_sum[["r_changed"]], tie_sum[["rej_changed"]]), "",
  "- **Old** = γ = 1 (local mean), recomputed with tie rounding; **New** = γ = ½. Same code, graphs, phases (`phase_policy_R_gen`, `F_gen`), null adjacency and seeds.",
  "- **Moved** = more than 1 percentage point (rates) or more than 0.01 (ρ, R², Δ, r); Δ̄ medians (4 decimals) are also flagged on a sign change; \"changed\" for counts and text.", "")
for (gname in unique(tab$group)) { x <- tab |> filter(group == gname)
  md <- c(md, sprintf("### %s", gname), "", "| Item | Old (γ = 1) | New (γ = ½) | Moved | Location |", "| :--- | :--- | :--- | :---: | :--- |",
          sprintf("| %s | %s | %s | %s | %s |", x$item, x$old_s, x$new_s, x$moved, x$location), "") }
res <- "output/half_bandwidth_rerun_results.md"
summ <- if (file.exists(res)) { x <- readLines(res); i <- grep("^<!-- TABLE BELOW", x); if (length(i)) x[seq_len(i[1] - 1)] else x } else
  c("# Half-bandwidth rerun: results", "", "_Summary to be written._", "")
writeLines(c(summ, md), res)
cat(sprintf("Wrote %s: %d items, %d moved\n", res, nrow(tab), sum(tab$moved %in% c("**yes**", "changed"))))

## ---- gamma = 1/2 data files for the figure scripts (GG_TAG=g05) ----------------------------------------
g <- G5$grid
f5 <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE) |> select(-mean_delta, -mean_abs_delta, -sd_delta) |>
  inner_join(g |> filter(sample %in% c("Symmetric", AS)) |> transmute(Scenario = recode(sample, Symmetric = "Isotropic"), replicate = Replicate, generation,
                                                                     mean_delta = gc_dbar_from_stored(dbar), mean_abs_delta, sd_delta), by = c("Scenario", "replicate", "generation"))
write.csv(f5[, names(read.csv("data/forward_scenarios.csv", nrows = 1))], "data/forward_scenarios_g05.csv", row.names = FALSE)
ib0 <- read.csv("data/isotropic_baseline.csv", stringsAsFactors = FALSE) |> select(-mean_delta, -mean_abs_delta) |>
  inner_join(g |> filter(sample %in% c("Burn-in", "Symmetric")) |> transmute(replicate = Replicate, phase = ifelse(sample == "Burn-in", "burnin", "forward"), generation,
                                                                          mean_delta = gc_dbar_from_stored(dbar), mean_abs_delta), by = c("replicate", "phase", "generation"))
write.csv(ib0[, names(read.csv("data/isotropic_baseline.csv", nrows = 1))], "data/isotropic_baseline_g05.csv", row.names = FALSE)
traj <- g |> filter(sample %in% c("Burn-in", "Symmetric", AS), generation >= 1904, !is.na(r2_pgd_dir)) |>
  transmute(Replicate, Treatment = factor(recode(sample, Symmetric = "Isotropic"), c("Burn-in", "Isotropic", "Redistributed", "Obstructed")), generation,
            r2_cgd, r2_pgd_blind, r2_pgd_dir, rel_fwd, status = "ok")
save(traj, file = "data/derived/ibgd_trajectory_fit_g05.rda")
saveRDS(G5$heat |> mutate(scenario = recode(scenario, Symmetric = "Isotropic")), "data/derived/node_gravity_heatmap_g05.rds")
saveRDS(T5 |> transmute(Replicate, scenario, generation, null = "N1w", r_S, p), "data/derived/sstest_g05.rds")
saveRDS(rean$truth$d |> select(-w_rho_tr, -w_rmse_tr, -w_rho_run, -w_rmse_run) |>
          inner_join(g |> filter(sample %in% c("Symmetric", AS)) |> transmute(Treatment = recode(sample, Symmetric = "Isotropic"), Replicate, Generation = generation,
                                                                             w_rho_tr, w_rmse_tr, w_rho_run, w_rmse_run), by = c("Treatment", "Replicate", "Generation")),
        "data/derived/truth_g05.rds")
cat("Wrote the gamma = 1/2 figure inputs\n")
