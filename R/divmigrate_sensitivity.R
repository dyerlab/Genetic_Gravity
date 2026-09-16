# Fig_DivMigrateSensitivity and Tab_SensitivityAUC.
#
# Reproduces the single-snapshot sensitivity analysis comparing genetic
# gravity and divMigrate as detectors of asymmetry: each replicate's snapshot
# is treated as one independent trial and scored, at each forward generation,
# against the empirical null supplied by the isotropic replicates at that same
# generation. Source data (data/divmigrate_timecourse_timecourse.csv)
# is the per-replicate/scenario/generation index summary already computed
# from the individual-based simulation output by the private research repo's
# R/divmigrate_timecourse.R; no raw simulation data or the diveRsity package
# is required to run this script. The AUC/power logic below reimplements the
# private repo's R/divmigrate_sensitivity.R verbatim (the rank-based
# Mann-Whitney-style .auc() helper is shared with R/detection_analysis.R for
# consistency across the manuscript's detection analyses).
#
# Run from the repository root:
#   Rscript R/divmigrate_sensitivity.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Rank-based AUC = P(positive > null) + 0.5 P(tie). Non-finite values dropped.
.auc <- function(pos, null) {
  pos  <- pos[is.finite(pos)]
  null <- null[is.finite(null)]
  n1 <- length(pos); n0 <- length(null)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(c(pos, null))
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

scn <- c("Isotropic", "Flux-conserved", "Rate-conserved")
tc <- read.csv("data/divmigrate_timecourse_timecourse.csv", stringsAsFactors = FALSE) |>
  mutate(scenario = factor(scenario, levels = scn))

## ---- divmigrate_sensitivity(): per method x scenario x generation --------
## statistic = "abs" (mean |Delta| / mean |A|), spec = 0.95 (5% FPR threshold)

spec <- 0.95
reference <- "Isotropic"

long <- tc |>
  transmute(replicate, scenario, generation,
            gravity = gravity_abs, divMigrate = divm_abs) |>
  pivot_longer(c(gravity, divMigrate), names_to = "method", values_to = "stat")

null_df <- long |>
  filter(scenario == reference) |>
  select(method, generation, null_stat = stat)

pos_df <- long |>
  filter(scenario != reference) |>
  mutate(scenario = droplevels(scenario))

sens <- pos_df |>
  group_by(method, scenario, generation) |>
  group_modify(function(df, key) {
    nl  <- null_df$null_stat[null_df$method == key$method &
                              null_df$generation == key$generation]
    pp  <- df$stat
    thr <- stats::quantile(nl[is.finite(nl)], spec, names = FALSE, type = 7)
    tibble::tibble(
      auc       = .auc(pp, nl),
      power     = mean(pp[is.finite(pp)] > thr),
      threshold = thr,
      fpr_check = mean(nl[is.finite(nl)] > thr),
      n_pos     = sum(is.finite(pp)),
      n_null    = sum(is.finite(nl)))
  }) |>
  ungroup() |>
  # Anchored to the first forward census (gen 2004): asymmetry switches on at
  # gen 2000 but the first snapshot is gen 2004, hence the +4 offset.
  mutate(generation_delta = generation - 2004) |>
  arrange(method, scenario, generation)

## ---- Stabilized generation (grid point closest to 200 gens after onset) --

.stab_gen <- sens$generation[which.min(abs(sens$generation_delta - 200))][1]
rv_sens_stab_delta <- sens$generation_delta[which.min(abs(sens$generation_delta - 200))][1]

.sens_at <- function(meas, meth, scn_name) {
  v <- sens[[meas]][sens$method == meth & sens$scenario == scn_name &
                     sens$generation == .stab_gen]
  if (length(v) == 0) NA_real_ else v[1]
}
rv_auc_grav_flux <- round(.sens_at("auc", "gravity",    "Flux-conserved"), 3)
rv_auc_divm_flux <- round(.sens_at("auc", "divMigrate", "Flux-conserved"), 3)
rv_auc_grav_rate <- round(.sens_at("auc", "gravity",    "Rate-conserved"), 3)
rv_auc_divm_rate <- round(.sens_at("auc", "divMigrate", "Rate-conserved"), 3)
rv_pow_grav_flux <- round(.sens_at("power", "gravity",    "Flux-conserved"), 2)
rv_pow_divm_flux <- round(.sens_at("power", "divMigrate", "Flux-conserved"), 2)
rv_pow_grav_rate <- round(.sens_at("power", "gravity",    "Rate-conserved"), 2)
rv_pow_divm_rate <- round(.sens_at("power", "divMigrate", "Rate-conserved"), 2)

cat(sprintf("rv_sens_stab_delta (stabilized window, generations after onset) = %s\n\n", rv_sens_stab_delta))
cat("=== Sensitivity at stabilized window ===\n")
cat(sprintf("Flux-conserved: gravity AUC=%s power=%s | divMigrate AUC=%s power=%s\n",
            rv_auc_grav_flux, rv_pow_grav_flux, rv_auc_divm_flux, rv_pow_divm_flux))
cat(sprintf("Rate-conserved: gravity AUC=%s power=%s | divMigrate AUC=%s power=%s\n",
            rv_auc_grav_rate, rv_pow_grav_rate, rv_auc_divm_rate, rv_pow_divm_rate))

## ---- Tab_SensitivityAUC ---------------------------------------------------

targets <- c(50, 100, 150, 200, 250, 300, 350, 400)
grid_d  <- sort(unique(sens$generation_delta))
pick    <- sort(unique(vapply(targets,
                              function(t) grid_d[which.min(abs(grid_d - t))],
                              numeric(1))))
ord_cols <- paste0("+", pick)

.panel <- function(scn_name) {
  sens |>
    filter(scenario == scn_name, generation_delta %in% pick) |>
    mutate(cell = sprintf("%.2f (%.2f)", auc, power),
           col  = paste0("+", generation_delta),
           Method = recode(method, gravity = "Genetic gravity (|Delta|)",
                            divMigrate = "divMigrate (|A|)")) |>
    select(Method, col, cell) |>
    pivot_wider(names_from = col, values_from = cell) |>
    arrange(Method) |>
    select(Method, all_of(ord_cols))
}

cat("\n=== Tab_SensitivityAUC (A) Flux-conserved ===\n")
print(as.data.frame(.panel("Flux-conserved")), row.names = FALSE)
cat("\n=== Tab_SensitivityAUC (B) Rate-conserved ===\n")
print(as.data.frame(.panel("Rate-conserved")), row.names = FALSE)

## ---- Fig_DivMigrateSensitivity -------------------------------------------

ref <- tibble::tibble(measure = c("AUC", "Power (FPR = 0.05)"), y = c(0.5, 0.05))
sens_long <- sens |>
  filter(generation_delta <= 400) |>
  pivot_longer(c(auc, power), names_to = "measure", values_to = "value") |>
  mutate(measure = recode(measure, auc = "AUC", power = "Power (FPR = 0.05)"),
         method  = recode(method, gravity = "Genetic gravity", divMigrate = "divMigrate"))

fig_sensitivity <- ggplot(sens_long, aes(generation_delta, value, colour = method)) +
  geom_hline(data = ref, aes(yintercept = y), linetype = "dashed", colour = "grey60") +
  geom_line() + geom_point(size = 1.5) +
  facet_grid(measure ~ scenario, scales = "free_y") +
  labs(x = "Generations since onset of asymmetry", y = NULL, colour = NULL) +
  theme_minimal() + theme(legend.position = "top")

out_png <- "data/derived/fig-divmigrate-sensitivity_reproduced.png"
ggsave(out_png, fig_sensitivity, width = 8, height = 5.5, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_DivMigrateSensitivity to %s\n", out_png))
