# Tab_GravityVsDivMigrate and Fig_DivMigrateTimecourse.
#
# Reproduces the final-generation (gen 2999) comparison of the genetic-gravity
# asymmetry index against divMigrate's directional relative migration, and the
# descriptive time-course of both indices' scenario-minus-isotropic signal
# over the forward phase. Source data (data/divmigrate_compare_summary.csv,
# data/divmigrate_timecourse_timecourse.csv) are the per-replicate/
# per-scenario/per-generation summaries already computed from the individual-
# based simulation output by the private research repo; no raw simulation data
# or the diveRsity package is required to run this script. The compare summary
# is from R/divmigrate_compare.R (divMigrate stat = "d", Jost's D). The
# time-course table uses divMigrate stat = "Nm": its divm_* columns are the
# mean (signed / absolute) edge asymmetry A_ij = m_i->j - m_j->i over each
# snapshot's retained graph edges, taken from the saved per-census Nm matrices
# (R/divmigrate_matrix_census.R), at the 21 generations of
# R/divmigrate_timecourse.R (2004, every 50th to 2954, and 2999). The time-course
# detection logic below reimplements the private repo's R/divmigrate_detection.R
# (paired t-test of scenario-minus-isotropic signal per generation, plus the
# first sustained-significant generation) directly on the cached per-replicate
# time-course table.
#
# Run from the repository root:
#   Rscript R/divmigrate_comparison.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---- Tab_GravityVsDivMigrate --------------------------------------------

cmp <- read.csv("data/divmigrate_compare_summary.csv", stringsAsFactors = FALSE)

scn <- c("Isotropic", "Redistributed", "Obstructed")
cmp <- cmp |> mutate(scenario = factor(scenario, levels = scn)) |> arrange(scenario)

tbl_divmigrate <- cmp |>
  transmute(Scenario = scenario,
            `gravity |Delta|` = gravity_abs_mean,
            `divM |A|`        = divm_abs_mean,
            `median rho (Delta, A)` = median_rho)

cat("=== Tab_GravityVsDivMigrate ===\n")
print(as.data.frame(tbl_divmigrate), row.names = FALSE)

.dm_get <- function(scn_name, col) cmp[[col]][cmp$scenario == scn_name][1]
rv_dm_rho_redistributed   <- .dm_get("Redistributed", "median_rho")
rv_dm_rho_obstructed   <- .dm_get("Obstructed", "median_rho")
rv_dm_iso_grav   <- .dm_get("Isotropic", "gravity_abs_mean")
rv_dm_iso_divm   <- .dm_get("Isotropic", "divm_abs_mean")
rv_dm_spec_ratio <- round(rv_dm_iso_divm / rv_dm_iso_grav, 1)

cat(sprintf(
  "\nrv_dm_rho_redistributed   = %s\nrv_dm_rho_obstructed   = %s\nrv_dm_iso_grav   = %s\nrv_dm_iso_divm   = %s\nrv_dm_spec_ratio = %sx\n",
  rv_dm_rho_redistributed, rv_dm_rho_obstructed, rv_dm_iso_grav, rv_dm_iso_divm, rv_dm_spec_ratio))

## ---- Fig_DivMigrateTimecourse --------------------------------------------

# divMigrate breakdown: relative Nm lies in [0, 1], so a valid snapshot's mean
# |A_ij| cannot exceed 1. Snapshots whose mean is non-finite (NaN Nm on an
# edge) or exceeds 1 (a finite but negative Nm on an edge, of order -1e10 to
# -1e13) are the same differentiation-collapse failure and are set to NA for
# divMigrate only, so they drop out of means, tests, and AUC alike.
tc <- read.csv("data/divmigrate_timecourse_timecourse.csv", stringsAsFactors = FALSE) |>
  mutate(scenario  = factor(scenario, levels = scn),
         dm_broken = !is.finite(divm_abs) | divm_abs > 1,
         divm_signed = ifelse(dm_broken, NA_real_, divm_signed),
         divm_abs    = ifelse(dm_broken, NA_real_, divm_abs))

# Reimplements R/divmigrate_detection.R: pair each asymmetric-scenario
# replicate's signed signal with the same replicate's isotropic signal at the
# same generation, then a two-sided paired t-test of the differences against
# zero at each generation. Detection generation = first grid generation from
# which the test stays significant at alpha for every later generation.
long <- tc |>
  select(replicate, scenario, generation,
         gravity = gravity_signed, divMigrate = divm_signed) |>
  pivot_longer(c(gravity, divMigrate), names_to = "method", values_to = "signal")

ref <- long |>
  filter(scenario == "Isotropic") |>
  select(replicate, generation, method, ref_signal = signal)

paired <- long |>
  filter(scenario != "Isotropic") |>
  left_join(ref, by = c("replicate", "generation", "method")) |>
  mutate(diff = signal - ref_signal, scenario = droplevels(scenario))

signal <- paired |>
  group_by(method, scenario, generation) |>
  summarise(
    n         = sum(is.finite(diff)),
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff   = sd(diff, na.rm = TRUE) / sqrt(sum(is.finite(diff))),
    p_value   = tryCatch(t.test(diff)$p.value, error = function(e) NA_real_),
    .groups   = "drop") |>
  arrange(method, scenario, generation)

alpha <- 0.05
detection <- signal |>
  group_by(method, scenario) |>
  group_modify(function(df, key) {
    df  <- arrange(df, generation)
    sig <- !is.na(df$p_value) & df$p_value < alpha
    n   <- nrow(df)
    sustained <- vapply(seq_len(n), function(i) all(sig[i:n]), logical(1))
    idx <- which(sustained)[1]
    tibble::tibble(
      detect_gen    = if (length(idx) && !is.na(idx)) df$generation[idx] else NA_real_,
      first_sig_gen = if (any(sig)) df$generation[which(sig)[1]] else NA_real_)
  }) |>
  ungroup() |>
  mutate(detect_delta = detect_gen - 2000, first_sig_delta = first_sig_gen - 2000)

cat("\n=== Detection generations (scenario minus isotropic, sustained p<0.05) ===\n")
print(as.data.frame(detection), row.names = FALSE)

vlines <- detection |> filter(!is.na(detect_gen))

# Winsorize using the gravity-method 1st/99th percentile range with a 25%
# margin, applied to both methods (divMigrate is prone to a D_max near-zero
# blow-up at late asymmetric generations; gravity is immune).
grav_lo <- quantile(signal$mean_diff[signal$method == "gravity"], 0.01) * 1.25
grav_hi <- quantile(signal$mean_diff[signal$method == "gravity"], 0.99) * 1.25
sig_w <- signal |>
  mutate(
    mean_diff_w = pmax(grav_lo, pmin(grav_hi, mean_diff)),
    ymin_w      = pmax(grav_lo, mean_diff_w - se_diff),
    ymax_w      = pmin(grav_hi, mean_diff_w + se_diff)
  )

fig_dm_timecourse <- ggplot(sig_w, aes(generation, mean_diff_w, colour = method, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = ymin_w, ymax = ymax_w), colour = NA, alpha = 0.15) +
  geom_line() +
  geom_point(data = dplyr::filter(sig_w, p_value < 0.05), size = 1.6) +
  geom_vline(data = vlines, aes(xintercept = detect_gen, colour = method),
             linetype = "dotted", show.legend = FALSE) +
  facet_wrap(~ scenario) +
  labs(x = "Generation", y = "Scenario minus isotropic (mean signed asymmetry)",
       colour = NULL, fill = NULL) +
  theme_minimal() + theme(legend.position = "top")

out_png <- "data/derived/fig-divmigrate-timecourse_reproduced.png"
ggsave(out_png, fig_dm_timecourse, width = 8, height = 4, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_DivMigrateTimecourse to %s\n", out_png))
