# R/analysis/phase_core_rerun.R
#
# Phase-core summaries at gamma = 1/2, with horizon and ascent-end sensitivity
# (local/phase_core_rerun_instructions.md). No new simulation or graph computation: every
# census value comes from the gamma = 1/2 (and, for run 4, gamma = 1) outputs of
# R/simulation/half_bandwidth_compute.R; C_D and diameter (bandwidth-free) and the
# polymorphism marker from data/forward_scenarios.csv; the He-gradient marker from
# data/derived/he_gradient.rds; divMigrate's Nm from node_potential_reanalysis.rda.
#
# Phase intervals per lineage (asymmetric scenarios), with R = end of ascent and F = onset
# of decay; lineages not flagged as reorganized are held in ascent until F (phase policy):
#   ascent [2000, min(R', F)), R' = R if reorganized else F; censored if neither boundary is reached
#   plateau [R', F) (or [R', 2999] if no F, censored); empty if R' >= F
#   decay [F, 2999] (always censored)
# Core = middle 50% of the interval: [a + (b - a)/4, b - (b - a)/4], b = the next boundary
# (or 2999). Variants: F at < 4.5 (current F_gen), < 6 and < 8 polymorphic loci per deme;
# R at 80%, 90% (current R_gen) and 95% of the He-gradient rise (same smoothing as
# R/simulation/he_gradient_boundary.R; the reorganized flag is unchanged).
#
# Run from the repository root:  Rscript R/analysis/phase_core_rerun.R
# Outputs: local/phase_core_rerun_results.md (tables; the summary is written by hand above
#   them), data/derived/phase_core_rerun.rda; adds F_gen6 and F_gen8 to
#   data/derived/lineage_phases.csv (F_gen unchanged).

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv) })
source("R/common/gravity_convention.R")
source("R/common/lineage_phases.R")
ld <- function(f) { e <- new.env(); load(f, envir = e); e }
G5 <- ld("data/derived/half_bandwidth_grid_g05.rda")$grid; T5 <- ld("data/derived/half_bandwidth_tests_g05.rda")$tests
G1 <- ld("data/derived/half_bandwidth_grid_g1.rda")$grid;  T1 <- ld("data/derived/half_bandwidth_tests_g1.rda")$tests
BW <- ld("data/derived/half_bandwidth_bw_g05.rda"); rean <- ld("data/derived/node_potential_reanalysis.rda")$rean
T5 <- T5 |> mutate(sample = recode(scenario, Isotropic = "Symmetric")); T1 <- T1 |> mutate(sample = recode(scenario, Isotropic = "Symmetric"))
fwd <- gc_read_forward() |> rename(Replicate = replicate) |> mutate(sample = recode(Scenario, Isotropic = "Symmetric"))
AS <- c("Redistributed", "Obstructed"); PHN <- c(ascent = "ascent", plateau = "plateau", decay = "decay")

## ---- boundary variants ----------------------------------------------------------------------
lp <- read.csv("data/derived/lineage_phases.csv", stringsAsFactors = FALSE)
Fv <- fwd |> filter(Scenario %in% AS) |> arrange(generation) |> group_by(Scenario, Replicate) |>
  summarise(F45 = generation[n_poly_pop < 4.5][1], F6 = generation[n_poly_pop < 6][1], F8 = generation[n_poly_pop < 8][1], .groups = "drop")
hg <- readRDS("data/derived/he_gradient.rds"); grid5 <- seq(2004, 2999, 5)
sm <- function(g, y) as.numeric(predict(gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML"), data.frame(g = grid5)))
base <- hg |> filter(scenario == "Burn-in") |> group_by(Replicate) |> summarise(b = mean(slope_int), .groups = "drop")
Rv <- bind_rows(lapply(AS, function(s) bind_rows(lapply(1:50, function(r) {
  x <- hg |> filter(scenario == s, Replicate == r) |> arrange(generation); f <- sm(x$generation, x$slope_int)
  b <- base$b[base$Replicate == r]; amp <- max(f) - b
  q <- function(p) { i <- which(f >= b + p * amp)[1]; if (is.na(i)) NA else grid5[i] }
  data.frame(Scenario = s, Replicate = r, R80 = q(0.8), R90 = q(0.9), R95 = q(0.95)) }))))
B <- lp |> select(Scenario, Replicate, R_gen, F_gen, reorganized) |> left_join(Fv, by = c("Scenario", "Replicate")) |> left_join(Rv, by = c("Scenario", "Replicate"))
chk_bounds <- c(F45 = sum(xor(is.na(B$F45), is.na(B$F_gen)) | (!is.na(B$F45) & !is.na(B$F_gen) & B$F45 != B$F_gen)),
                R90 = sum(xor(is.na(B$R90), is.na(B$R_gen)) | (!is.na(B$R90) & !is.na(B$R_gen) & B$R90 != B$R_gen)))
if (!all(c("F_gen6", "F_gen8") %in% names(lp))) {
  lp2 <- lp |> left_join(Fv |> transmute(Scenario, Replicate, F_gen6 = F6, F_gen8 = F8), by = c("Scenario", "Replicate"))
  write.csv(lp2, "data/derived/lineage_phases.csv", row.names = FALSE)
}

## intervals and core windows for one variant
intervals <- function(Rcol, Fcol) B |> transmute(Scenario, Replicate, reorganized,
    R = .data[[Rcol]], F = .data[[Fcol]]) |>
  mutate(Rp = ifelse(reorganized, R, F), Fi = ifelse(is.na(F), Inf, F), Rpi = ifelse(is.na(Rp), Inf, Rp),
         asc_a = 2000, asc_b = pmin(Rpi, Fi), asc_cens = !is.finite(asc_b), asc_b = ifelse(is.finite(asc_b), asc_b, 2999),
         pl_ok = Rpi < Fi, pl_a = ifelse(pl_ok, Rpi, NA), pl_b = ifelse(pl_ok, ifelse(is.finite(Fi), Fi, 2999), NA), pl_cens = pl_ok & !is.finite(Fi),
         de_ok = is.finite(Fi), de_a = ifelse(de_ok, Fi, NA), de_b = ifelse(de_ok, 2999, NA))
assign_phase <- function(d, iv, core = TRUE, s = "sample", gen = "generation") {
  x <- d |> left_join(iv |> rename(sample = Scenario), by = c("sample", "Replicate"))
  g <- x[[gen]]; inn <- function(a, b, right_closed) { if (!core) g >= a & (g < b | (right_closed & g <= b))
    else { lo <- a + (b - a) / 4; hi <- b - (b - a) / 4; g >= lo & g <= hi } }
  ph <- ifelse(!x$sample %in% AS, ifelse(x$sample == "Burn-in" | g < 2000, "Burn-in", "Symmetric"),
         ifelse(inn(x$asc_a, x$asc_b, x$asc_cens) %in% TRUE & (core | g < x$asc_b | x$asc_cens), "ascent",
         ifelse(x$pl_ok & inn(x$pl_a, x$pl_b, x$pl_cens) %in% TRUE, "plateau",
         ifelse(x$de_ok & inn(x$de_a, x$de_b, TRUE) %in% TRUE, "decay", NA))))
  d$phase <- ph; d
}

## ---- summaries (as R/analysis/half_bandwidth_results.R) --------------------------------------
summarise_all <- function(G, T, iv, core) {
  ga <- assign_phase(G, iv, core); ta <- assign_phase(T, iv, core)
  fa <- assign_phase(fwd |> filter(Scenario %in% AS), iv, core)
  out <- list()
  out$desc <- ga |> filter(sample %in% AS, !is.na(phase)) |>
    inner_join(fa |> select(sample, Replicate, generation, CD, diameter), by = c("sample", "Replicate", "generation")) |>
    group_by(sample, phase, Replicate) |> summarise(across(c(CD, diameter, mean_abs_delta, dbar), mean), n = n(), .groups = "drop") |>
    group_by(sample, phase) |> summarise(lin = n(), cens = sum(n), CD = median(CD), diameter = median(diameter), abs = median(mean_abs_delta),
                                         pos = mean(dbar > 0), dbar = median(dbar), .groups = "drop")
  t2 <- ta |> filter(sample %in% AS, !is.na(phase))
  out$ss <- t2 |> group_by(sample, phase) |> summarise(cens = n(), lin = n_distinct(Replicate), rej = mean(p < .05), medr = median(r_S), up = mean(r_S < 0), .groups = "drop") |>
    left_join(t2 |> group_by(sample, phase, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |> group_by(sample, phase) |>
                summarise(lin_up = mean(m < 0), .groups = "drop"), by = c("sample", "phase"))
  ib <- ga |> filter(sample %in% AS, !is.na(phase), !is.na(r2_pgd_dir)) |> mutate(det = r2_pgd_dir > r2_cgd)
  out$ib <- ib |> group_by(sample, phase) |> summarise(cens = n(), lin = n_distinct(Replicate), cgd = median(r2_cgd), blind = median(r2_pgd_blind), dir = median(r2_pgd_dir),
                                                       gain = median(r2_pgd_dir - r2_pgd_blind), det_rate = mean(det), med_dR2 = median(r2_pgd_dir - r2_cgd),
                                                       lt1 = mean(rel_fwd[det] < 1), med_ratio = median(rel_fwd[det]), .groups = "drop")
  ov <- inner_join(t2 |> select(sample, Replicate, generation, phase, r_S, p), G |> select(sample, Replicate, generation, r2_cgd, r2_pgd_dir, rel_fwd),
                   by = c("sample", "Replicate", "generation")) |> filter(!is.na(r2_pgd_dir)) |> mutate(ss = p < .05, ibp = r2_pgd_dir > r2_cgd)
  out$ov <- ov |> group_by(sample, phase) |> summarise(cens = n(), both = sum(ss & ibp), either = sum(ss | ibp), boe = both / max(1, either),
                                                       same = sum((r_S < 0 & rel_fwd < 1)[ss & ibp]), .groups = "drop")
  nm <- rean$truth$d |> transmute(sample = Treatment, Replicate, generation = Generation, Exclusion, nm_rho_run)
  nc <- rean$corridor$d |> transmute(sample = Treatment, Replicate, generation = Generation, nm_ok, nm_fwd_gt)
  tr <- ga |> filter(sample %in% AS, !is.na(phase)) |> inner_join(nm, by = c("sample", "Replicate", "generation")) |> inner_join(nc, by = c("sample", "Replicate", "generation"))
  out$tr <- tr |> group_by(sample, phase) |> summarise(cens = n(), w = median(w_rho_tr), nmr = median(nm_rho_run, na.rm = TRUE), wbeat = mean((w_rho_tr > nm_rho_run)[Exclusion == ""]),
                                                       wf = mean(w_fwd_gt[nm_ok], na.rm = TRUE), nf = mean(nm_fwd_gt, na.rm = TRUE), nm_def = mean(Exclusion == ""), .groups = "drop") |>
    left_join(tr |> filter(nm_ok) |> group_by(sample, phase, Replicate) |> summarise(dd = mean(nm_fwd_gt - w_fwd_gt), .groups = "drop") |>
                group_by(sample, phase) |> summarise(nm_better = sum(dd > 0), nm_lin = n(), .groups = "drop"), by = c("sample", "phase"))
  eb <- assign_phase(BW$edge_bw |> rename(sample = Scenario), iv, core) |> filter(!is.na(phase))
  out$ebw <- eb |> group_by(bandwidth, phase) |> summarise(cens = n(), concord = mean(concord), rs = mean(rho_signed, na.rm = TRUE), kept = mean(dbar_kept), .groups = "drop")
  sb <- assign_phase(BW$ss_bw |> mutate(sample = recode(scenario, Isotropic = "Symmetric")) |> filter(!failed), iv, core) |> filter(sample %in% AS, !is.na(phase))
  out$sbw <- sb |> group_by(bandwidth, sample, phase) |> summarise(rej = mean(p < .05), medr = median(r_S), .groups = "drop") |>
    left_join(sb |> group_by(bandwidth, sample, phase, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |> group_by(bandwidth, sample, phase) |>
                summarise(lin_up = mean(m < 0), .groups = "drop"), by = c("bandwidth", "sample", "phase"))
  out$counts <- list(t50 = t2 |> group_by(sample, phase) |> summarise(lin = n_distinct(Replicate), cens = n(), .groups = "drop"),
                     g5 = ga |> filter(sample %in% AS, !is.na(phase)) |> group_by(sample, phase) |> summarise(lin = n_distinct(Replicate), cens = n(), .groups = "drop"))
  out
}

## long form: one row per item
flatten <- function(S) {
  r <- list(); a <- function(sec, s, p, item, v, type = "pct") r[[length(r) + 1]] <<- data.frame(section = sec, sample = s, phase = p, item = item, value = v, type = type)
  for (i in seq_len(nrow(S$desc))) with(S$desc[i, ], { a("Graph metrics", sample, phase, "C̄_D", CD, "cd"); a("Graph metrics", sample, phase, "diameter", diameter, "len")
    a("Graph metrics", sample, phase, "mean |Δ|", abs, "stat3"); a("Graph metrics", sample, phase, "Δ̄ positive (lineages)", pos); a("Graph metrics", sample, phase, "Δ̄ median", dbar, "stat4") })
  for (i in seq_len(nrow(S$ss))) with(S$ss[i, ], { a("Tab_SourceSinkPower", sample, phase, "rejection", rej); a("Tab_SourceSinkPower", sample, phase, "median r(S,x)", medr, "stat2")
    a("Tab_SourceSinkPower", sample, phase, "censuses with sources upstream", up); a("Tab_SourceSinkPower", sample, phase, "lineages with sources upstream", lin_up) })
  for (i in seq_len(nrow(S$ib))) with(S$ib[i, ], { for (v in c("cgd", "blind", "dir", "gain")) a("Tab_IBGDFitComparison", sample, phase,
      c(cgd = "cGD R²", blind = "pGD direction-blind R²", dir = "pGD direction-aware R²", gain = "directional gain")[v], get(v), "stat2")
    a("Tab_IBGDDetection", sample, phase, "ΔR² > 0 (detection)", det_rate); a("Tab_IBGDDetection", sample, phase, "median ΔR²", med_dR2, "stat3")
    a("Tab_IBGDDetection", sample, phase, "detections with slope ratio < 1", lt1); a("Tab_IBGDDetection", sample, phase, "median slope ratio (detections)", med_ratio, "stat2") })
  for (i in seq_len(nrow(S$ov))) with(S$ov[i, ], { a("Tab_TestOverlap", sample, phase, "both / either", boe); a("Tab_TestOverlap", sample, phase, "joint agreeing on direction", same / max(1, both)) })
  for (i in seq_len(nrow(S$tr))) with(S$tr[i, ], { a("divMigrate", sample, phase, "truth ρ, w", w, "stat2"); a("divMigrate", sample, phase, "truth ρ, Nm", nmr, "stat2")
    a("divMigrate", sample, phase, "w higher than Nm (paired)", wbeat); a("divMigrate", sample, phase, "corridor direction, w (Nm-defined censuses)", wf)
    a("divMigrate", sample, phase, "corridor direction, Nm", nf); a("divMigrate", sample, phase, "Nm defined", nm_def)
    a("divMigrate", sample, phase, "lineages where Nm orients better (share)", nm_better / nm_lin) })
  for (i in which(S$ebw$phase == "plateau")) with(S$ebw[i, ], { a("Tab_BandwidthSensitivity (plateau, pooled)", bandwidth, phase, "same sign", concord)
    a("Tab_BandwidthSensitivity (plateau, pooled)", bandwidth, phase, "rank ρ signed", rs, "stat2"); a("Tab_BandwidthSensitivity (plateau, pooled)", bandwidth, phase, "Δ̄ sign kept", kept) })
  for (i in which(S$sbw$phase == "plateau")) with(S$sbw[i, ], { a("Tab_BandwidthSensitivitySourceSink (plateau)", paste(bandwidth, sample), phase, "rejection", rej)
    a("Tab_BandwidthSensitivitySourceSink (plateau)", paste(bandwidth, sample), phase, "median r(S,x)", medr, "stat2")
    a("Tab_BandwidthSensitivitySourceSink (plateau)", paste(bandwidth, sample), phase, "lineages upstream", lin_up) })
  bind_rows(r)
}

## ---- runs ---------------------------------------------------------------------------------------
IV <- list(full = intervals("R_gen", "F_gen"), F45 = intervals("R_gen", "F45"), F6 = intervals("R_gen", "F6"), F8 = intervals("R_gen", "F8"),
           R80 = intervals("R80", "F45"), R95 = intervals("R95", "F45"))
S_full <- summarise_all(G5, T5, IV$full, core = FALSE)
S <- lapply(IV[c("F45", "F6", "F8", "R80", "R95")], function(iv) summarise_all(G5, T5, iv, core = TRUE))
S1 <- summarise_all(G1, T1, IV$F45, core = TRUE); S1_full <- summarise_all(G1, T1, IV$full, core = FALSE)
L <- lapply(c(list(full = S_full), S), flatten)
key <- c("section", "sample", "phase", "item", "type")
W <- Reduce(function(a, b) full_join(a, b, by = key), lapply(names(L), function(n) L[[n]] |> rename(!!n := value)))
# reproduction of the full-phase gamma = 1/2 values (R/analysis/half_bandwidth_results.R)
hb <- ld("data/derived/half_bandwidth_results.rda")
rep_chk <- c(ss_plateau_R = S_full$ss$rej[S_full$ss$sample == "Redistributed" & S_full$ss$phase == "plateau"] - hb$S$new$rej[hb$S$new$sample == "Redistributed" & hb$S$new$phase == "Plateau"],
             ib_plateau_O = S_full$ib$det_rate[S_full$ib$sample == "Obstructed" & S_full$ib$phase == "plateau"] - hb$I$new$det_rate[hb$I$new$sample == "Obstructed" & hb$I$new$phase == "Plateau"],
             abs_decay_R = S_full$desc$abs[S_full$desc$sample == "Redistributed" & S_full$desc$phase == "decay"] - hb$D$new$abs[hb$D$new$sample == "Redistributed" & hb$D$new$phase == "Fall apart"])

## ---- lengths, censoring and counts -----------------------------------------------------------------
cnt <- function(iv, nm) iv |> group_by(Scenario) |> summarise(variant = nm, reach_decay = sum(is.finite(Fi)),
  plateau_lineages = sum(pl_ok), plateau_censored = sum(pl_cens), plateau_len_med = median((pl_b - pl_a)[pl_ok]),
  ascent_censored = sum(asc_cens), ascent_len_med = median(asc_b - asc_a), decay_len_med = median((de_b - de_a)[de_ok]), .groups = "drop")
counts <- bind_rows(lapply(names(IV)[-1], function(n) cnt(IV[[n]], n)))
lost50 <- function(iv) { g50 <- seq(2004, 2954, 50); x <- iv |> rowwise() |> mutate(
  asc = any(g50 >= asc_a + (asc_b - asc_a) / 4 & g50 <= asc_b - (asc_b - asc_a) / 4),
  pl = if (pl_ok) any(g50 >= pl_a + (pl_b - pl_a) / 4 & g50 <= pl_b - (pl_b - pl_a) / 4) else NA,
  de = if (de_ok) any(g50 >= de_a + (de_b - de_a) / 4 & g50 <= de_b - (de_b - de_a) / 4) else NA) |> ungroup()
  x |> group_by(Scenario) |> summarise(asc_lost = sum(!asc), pl_lost = sum(!pl, na.rm = TRUE), de_lost = sum(!de, na.rm = TRUE), .groups = "drop") }
LOST <- bind_rows(lapply(names(IV)[-1], function(n) lost50(IV[[n]]) |> mutate(variant = n)))

## ---- tables ---------------------------------------------------------------------------------------
fmt <- function(x, type) ifelse(is.na(x), "—", ifelse(type == "pct", sprintf("%.1f%%", 100 * x), ifelse(type %in% c("stat2", "len"), sprintf("%.2f", x),
         ifelse(type %in% c("stat3", "cd"), sprintf("%.3f", x), sprintf("%.4f", x)))))
# rates: 1 / 2 percentage points; rho, R^2, r, Delta (incl. Delta-bar): 0.01 / 0.02; C_D and diameter (graph-only) reported, not flagged
thr <- function(type, small) ifelse(type %in% c("cd", "len"), Inf, if (small) 0.01 else 0.02)
md <- function(df) { df <- as.data.frame(df); c(paste("|", paste(names(df), collapse = " | "), "|"), paste("|", paste(c(":---", ":---", ":---", rep(":---:", ncol(df) - 3)), collapse = " | "), "|"),
  apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")), "") }
ord <- c("ascent", "plateau", "decay")
W <- W |> mutate(phase = factor(phase, ord)) |> arrange(section, sample, phase, item)
out <- c("<!-- TABLES BELOW ARE WRITTEN BY R/analysis/phase_core_rerun.R -->", "",
  sprintf("Checks: the boundary recomputation reproduces `F_gen` and `R_gen` exactly (%d and %d mismatches of 100); full-phase values reproduce the γ = ½ rerun (largest difference %.1g).",
          chk_bounds[["F45"]], chk_bounds[["R90"]], max(abs(rep_chk))), "",
  "Flags: **yes** = an *ascent* or *plateau* value moves by more than 2 percentage points (rates) or 0.02 (ρ, R², r, Δ, including Δ̄) between the 4.5 and the 6 or 8 horizon (run 2) or between 90% and 80% or 95% (run 3). *Decay* rows are not flagged. C̄_D and diameter are graph-only and are reported but not flagged. In run 1 the flag marks a change of more than 1 point / 0.01 between full phase and core.", "")
out <- c(out, "## Run 1. Full phase vs core, γ = ½, F₄.₅", "")
r1 <- W |> mutate(change = F45 - full, flag = ifelse(phase != "decay" & abs(change) %in% abs(change) & !is.na(change) & abs(change) > thr(type, TRUE), "**yes**", ""))
for (sec in unique(r1$section)) { x <- r1 |> filter(section == sec)
  out <- c(out, sprintf("### %s", sec), "", md(x |> transmute(Scenario = sample, Phase = phase, Item = item, `Full phase` = fmt(full, type), Core = fmt(F45, type),
                                                              Change = ifelse(type == "pct", sprintf("%+.1f pp", 100 * change), sprintf("%+.3f", change)), Flag = flag))) }
out <- c(out, "## Run 2. Horizon sensitivity of the cores (F₄.₅ / F₆ / F₈), γ = ½", "")
r2 <- W |> mutate(d = F8 - F45, d6 = F6 - F45, flag = ifelse(phase != "decay" & (abs(d) > thr(type, FALSE) | abs(d6) > thr(type, FALSE)) %in% TRUE, "**yes**", ""))
for (sec in unique(r2$section)) { x <- r2 |> filter(section == sec)
  out <- c(out, sprintf("### %s", sec), "", md(x |> transmute(Scenario = sample, Phase = phase, Item = item, `Core 4.5` = fmt(F45, type), `Core 6` = fmt(F6, type), `Core 8` = fmt(F8, type), Flag = flag))) }
out <- c(out, "### Phase lengths and lineages reaching decay, by horizon", "",
         md(counts |> filter(variant %in% c("F45", "F6", "F8")) |> transmute(Scenario, Horizon = recode(variant, F45 = "4.5", F6 = "6", F8 = "8"), Variant = "",
            `Reach decay` = reach_decay, `Lineages with a plateau` = plateau_lineages, `Plateau censored` = plateau_censored, `Median plateau length` = plateau_len_med,
            `Median ascent length` = ascent_len_med, `Ascent censored` = ascent_censored, `Median decay length` = decay_len_med) |> select(-Variant)))
out <- c(out, "## Run 3. Ascent-end sensitivity of the plateau core (80% / 90% / 95% of the H_e-gradient rise), γ = ½, F₄.₅", "")
r3 <- W |> filter(phase == "plateau") |> mutate(flag = ifelse((abs(R80 - F45) > thr(type, FALSE) | abs(R95 - F45) > thr(type, FALSE)) %in% TRUE, "**yes**", ""))
for (sec in unique(r3$section)) { x <- r3 |> filter(section == sec)
  out <- c(out, sprintf("### %s", sec), "", md(x |> transmute(Scenario = sample, Phase = phase, Item = item, `80%` = fmt(R80, type), `90% (current)` = fmt(F45, type), `95%` = fmt(R95, type), Flag = flag))) }
out <- c(out, "### Plateau lineages and lengths by ascent-end rule", "",
         md(counts |> filter(variant %in% c("R80", "F45", "R95")) |> transmute(Scenario, Rule = recode(variant, R80 = "80%", F45 = "90%", R95 = "95%"), x = "",
            `Lineages with a plateau` = plateau_lineages, `Plateau censored` = plateau_censored, `Median plateau length` = plateau_len_med, `Median ascent length` = ascent_len_med) |> select(-x)))
h <- function(S, s) c(rej = S$ss$rej[S$ss$sample == s & S$ss$phase == "plateau"], medr = S$ss$medr[S$ss$sample == s & S$ss$phase == "plateau"],
                      lin = S$ss$lin_up[S$ss$sample == s & S$ss$phase == "plateau"], det = S$ib$det_rate[S$ib$sample == s & S$ib$phase == "plateau"],
                      cor = S$tr$wf[S$tr$sample == s & S$tr$phase == "plateau"], cor_nm = S$tr$nf[S$tr$sample == s & S$tr$phase == "plateau"])
r4 <- bind_rows(lapply(AS, function(s) bind_rows(
  data.frame(Scenario = s, Bandwidth = "Local mean (γ = 1), full phase", t(h(S1_full, s))), data.frame(Scenario = s, Bandwidth = "Local mean (γ = 1), core", t(h(S1, s))),
  data.frame(Scenario = s, Bandwidth = "Degree-neutral (γ = ½), full phase", t(h(S_full, s))), data.frame(Scenario = s, Bandwidth = "Degree-neutral (γ = ½), core", t(h(S$F45, s))))))
out <- c(out, "## Run 4. Local mean vs degree-neutral, plateau headline values (F₄.₅)", "",
         md(r4 |> transmute(Scenario, Phase = "plateau", Bandwidth, `Source–sink rejection` = sprintf("%.1f%%", 100 * rej), `Median r(S,x)` = sprintf("%.2f", medr),
                            `Lineages upstream` = sprintf("%.0f%%", 100 * lin), `ΔR² detection` = sprintf("%.1f%%", 100 * det),
                            `Corridor direction, w (Nm)` = sprintf("%.1f%% (%.1f%%)", 100 * cor, 100 * cor_nm))))
cc <- function(S, nm) bind_rows(S$counts$t50 |> mutate(grid = "50-generation"), S$counts$g5 |> mutate(grid = "5-generation")) |> mutate(variant = nm)
CC <- bind_rows(cc(S_full, "full phase"), lapply(names(S), function(n) cc(S[[n]], n))) |> mutate(phase = factor(phase, ord)) |> arrange(variant, grid, sample, phase)
out <- c(out, "## Lineage and census counts per core", "",
         md(CC |> transmute(Variant = variant, Grid = grid, Scenario = sample, Phase = phase, Lineages = lin, Censuses = cens)),
         "**Censoring and lineages lost on the 50-generation grid** (a lineage is lost when its core window holds no census on that grid):", "",
         md(counts |> left_join(LOST, by = c("Scenario", "variant")) |> transmute(Variant = variant, Scenario, x = "", `Ascent censored` = ascent_censored,
            `Plateau censored / with plateau` = sprintf("%d / %d", plateau_censored, plateau_lineages), `Decay (all censored)` = reach_decay,
            `Lost: ascent / plateau / decay` = sprintf("%d / %d / %d", asc_lost, pl_lost, de_lost)) |> select(-x)))
res <- "output/phase_core_rerun_results.md"
summ <- if (file.exists(res)) { x <- readLines(res); i <- grep("^<!-- TABLES BELOW", x); if (length(i)) x[seq_len(i[1] - 1)] else x } else c("# Phase-core rerun: results", "", "_Summary to be written._", "")
writeLines(c(summ, out), res)
save(B, IV, S_full, S, S1, S1_full, W, counts, LOST, CC, rep_chk, chk_bounds, file = "data/derived/phase_core_rerun.rda")
cat(sprintf("Wrote %s; flags run 2: %d, run 3: %d\n", res, sum(r2$flag != ""), sum(r3$flag != "")))
