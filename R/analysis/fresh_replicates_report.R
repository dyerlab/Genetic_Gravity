# R/analysis/fresh_replicates_report.R
#
# Tables for Part 3 of local/gstudio_1145_verification_and_gamma_rule_instructions.md
# from data/derived/fresh_replicates.rda. Writes local/fresh_replicates_tables.md.
#
# Run from the repository root:  Rscript R/analysis/fresh_replicates_report.R

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("R/common/lineage_phases.R")
load("data/derived/fresh_replicates.rda")                       # fr_stats, fr_bounds, fr_tests, meta
OUTMD <- "output/fresh_replicates_tables.md"
md <- function(df) {
  df <- as.data.frame(df); hdr <- paste("|", paste(names(df), collapse = " | "), "|")
  sep <- paste("|", paste(c(":----", rep(":---:", ncol(df) - 1)), collapse = " | "), "|")
  c(hdr, sep, apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")), "")
}
pct <- function(v, dig = 1) ifelse(is.na(v), "—", sprintf(paste0("%.", dig, "f%%"), 100 * v))
num <- function(v, dig = 2) ifelse(is.na(v), "—", formatC(v, format = "f", digits = dig))
L <- character(0); add <- function(...) L <<- c(L, ...)
PH_LAB <- c(Reorganize = "ascent", Plateau = "plateau", `Fall apart` = "decay")
stopifnot(!"err" %in% names(fr_tests) || all(is.na(fr_tests$err)))

t <- fr_tests |> mutate(
  phase = ifelse(sample %in% c("Redistributed", "Obstructed"), PH_LAB[as.character(lp_phase(sample, Replicate, generation, bounds = fr_bounds))], "all"),
  col = factor(ifelse(phase == "all", sample, paste(sample, phase)),
               c("Symmetric", "sym-low", "Redistributed ascent", "Redistributed plateau", "Redistributed decay", "Obstructed ascent", "Obstructed plateau", "Obstructed decay")),
  setting = factor(recode(setting, `1` = "γ = 1", `0.5` = "γ = 0.5", star = "γ*"), c("γ = 1", "γ = 0.5", "γ*")))

add("## Part 3. Fresh replicates (50 new lineages, new seeds)", "",
    sprintf("Simulated %s with the manuscript's driver functions and parameters (lineages %d–%d; burn-in, then Symmetric, Redistributed, Obstructed, sym-low). %d census graphs on the 50-generation grid (%d expected); %d test rows. Simulation %.0f min, total %.0f min.",
            meta$date, min(meta$reps), max(meta$reps), n_distinct(paste(t$Replicate, t$sample, t$generation)), 50 * 4 * 20, nrow(t), meta$sim_min, meta$total_min), "")
b <- fr_bounds |> filter(Scenario != "Symmetric")
add(sprintf("**Phases** (current rule; reorganized threshold = 95th percentile of the fresh Symmetric rises, %.2e; the manuscript's threshold is %.2e):", meta$thr_fresh, meta$thr_manuscript), "",
    md(b |> group_by(Scenario) |> summarise(`reorganized (fresh thr)` = sprintf("%d / %d", sum(reorganized), n()),
                                            `reorganized (manuscript thr)` = sprintf("%d / %d", sum(reorganized_ms_thr), n()),
                                            `end of ascent, median (IQR), flagged` = sprintf("%.0f (%.0f–%.0f)", median(R_gen[reorganized], na.rm = TRUE), quantile(R_gen[reorganized], .25, na.rm = TRUE), quantile(R_gen[reorganized], .75, na.rm = TRUE)),
                                            `reached decay` = sum(!is.na(F_gen)), `decay onset, median` = num(median(F_gen, na.rm = TRUE), 0), .groups = "drop")),
    md(t |> filter(setting == "γ = 1") |> count(col, name = "censuses") |> rename(`Sample / phase` = col)))

add("**γ\\***: median (IQR); flags", "",
    md(t |> filter(setting == "γ*") |> group_by(col) |> summarise(n = n(), `γ*` = sprintf("%.2f (%.2f–%.2f)", median(gamma), quantile(gamma, .25), quantile(gamma, .75)),
                                                               no_root = sum(flag == "no_root"), multi_root = sum(flag == "multi_root"), .groups = "drop") |> rename(`Sample / phase` = col)))
add("**Source–sink false-positive rate (α = 0.01 / 0.05 / 0.10) and ΔR² > 0 rate, symmetric samples**", "",
    md(t |> filter(col %in% c("Symmetric", "sym-low")) |> group_by(setting, col) |>
         summarise(v = sprintf("%s / %s / %s; ΔR² %s", pct(mean(p < .01)), pct(mean(p < .05)), pct(mean(p < .10)), pct(mean(dR2 > 0, na.rm = TRUE))), .groups = "drop") |>
         pivot_wider(names_from = col, values_from = v) |> rename(Bandwidth = setting)))
ph <- t |> filter(!col %in% c("Symmetric", "sym-low"), !is.na(col))
add("**Asymmetric scenarios: source–sink rejection / median r(S,x) / lineages with mean r < 0 (lineages); ΔR² detection / detections with slope ratio < 1**", "",
    md(ph |> group_by(setting, col) |> summarise(n = n(), rej = pct(mean(p < .05)), med = num(median(r_S)), det = pct(mean(dR2 > 0, na.rm = TRUE)),
                                                 lt1 = pct(mean(rel_fwd[dR2 > 0 & !is.na(dR2)] < 1), 0), .groups = "drop") |>
         left_join(ph |> group_by(setting, col, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |> group_by(setting, col) |>
                     summarise(lin = sprintf("%s (%d)", pct(mean(m < 0), 0), n()), .groups = "drop"), by = c("setting", "col")) |>
         transmute(setting, col = sprintf("%s (n=%d)", col, n), v = sprintf("%s / %s / %s; %s / %s", rej, med, lin, det, lt1)) |>
         pivot_wider(names_from = col, values_from = v) |> rename(Bandwidth = setting)))
add("**Symmetric lineage lean (lineages with mean r(S,x) > 0; sign test)**", "",
    md(t |> filter(col %in% c("Symmetric", "sym-low")) |> group_by(setting, col, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
         group_by(setting, col) |> summarise(v = sprintf("%d / %d (p = %.2f)", sum(m > 0), n(), binom.test(sum(m > 0), sum(m != 0))$p.value), .groups = "drop") |>
         pivot_wider(names_from = col, values_from = v) |> rename(Bandwidth = setting)))
writeLines(L, OUTMD); cat("Wrote", OUTMD, "\n")
