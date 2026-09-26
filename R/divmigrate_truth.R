# Fig_DivMigrateTruth.
#
# Agreement of each directional estimator with the true simulated migration
# matrix. For every scenario the true matrix M_true has m_{i->j} equal to the
# scenario's forward or reverse rate on each stepping-stone link (Pop01..Pop25
# in chain order; Tab_ScenarioParameters) and zero on every non-adjacent pair.
# The two estimates, per replicate and forward census, are
#   gravity     the pGD directional weights recast as similarities,
#               s_ij = 1 / (1 + p_ij) on retained Population Graph edges and 0
#               elsewhere, with p_ij = e_ij * w_{i->j} / (w_{i->j} + w_{j->i});
#   divMigrate  relative Nm (diveRsity::divMigrate, stat = "Nm").
# Each matrix (and M_true) is rescaled by its own maximum off-diagonal value;
# diagonals are excluded throughout. A census whose Nm matrix has any
# non-finite or negative off-diagonal entry is excluded for both methods.
# Source data (data/divmigrate_truth_compare.csv) are the per-replicate/
# census metrics computed from the individual-based simulation output by the
# private research repo's R/divmigrate_truth_compare.R (score = "similarity");
# no raw simulation data or the diveRsity package is required to run this
# script. The figure shows rho_all: Spearman rho between the rescaled estimate
# and rescaled truth over all K(K-1) = 600 ordered pairs.
#
# Run from the repository root:
#   Rscript R/divmigrate_truth.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

scn <- c("Isotropic", "Redistributed", "Obstructed")
tr <- read.csv("data/divmigrate_truth_compare.csv", stringsAsFactors = FALSE) |>
  mutate(scenario = factor(scenario, levels = scn),
         method   = factor(method, levels = c("gravity", "divMigrate"),
                           labels = c("Genetic gravity (pGD similarity)", "divMigrate (Nm)")))

## ---- Exclusions ------------------------------------------------------------

excl <- tr |>
  filter(method == "divMigrate (Nm)") |>
  group_by(scenario) |>
  summarise(censuses = n(), nonfinite = sum(Exclusion == "nonfinite"),
            negative = sum(Exclusion == "negative"), excluded = sum(Excluded),
            frac = round(mean(Excluded), 3), .groups = "drop")
cat("=== Censuses excluded (non-finite or negative Nm) ===\n")
print(as.data.frame(excl), row.names = FALSE)

## ---- Per-census summaries ----------------------------------------------------

med <- tr |>
  filter(!Excluded) |>
  group_by(scenario, method, generation) |>
  summarise(n = n(), rho = median(rho_all),
            lo = quantile(rho_all, 0.25), hi = quantile(rho_all, 0.75), .groups = "drop")

paired <- tr |>
  filter(!Excluded) |>
  select(replicate, scenario, generation, method, rho_all) |>
  tidyr::pivot_wider(names_from = method, values_from = rho_all) |>
  mutate(gap = `Genetic gravity (pGD similarity)` - `divMigrate (Nm)`)
cat("\n=== Paired comparison (same replicate and census) ===\n")
print(as.data.frame(paired |> group_by(scenario) |>
  summarise(censuses = n(), gravity_higher = round(mean(gap > 0), 3),
            median_gap = round(median(gap), 3), .groups = "drop")), row.names = FALSE)

cat("\n=== Median rho, early (2004-2249) vs late (2754-2999) ===\n")
print(as.data.frame(med |>
  mutate(window = case_when(generation < 2250 ~ "early", generation >= 2754 ~ "late")) |>
  filter(!is.na(window)) |>
  group_by(scenario, method, window) |>
  summarise(rho = round(median(rho), 3), reps = round(mean(n), 1), .groups = "drop")), row.names = FALSE)

## ---- Fig_DivMigrateTruth -----------------------------------------------------

fig_truth <- ggplot(med, aes(generation, rho, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.18) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Genetic gravity (pGD similarity)" = "firebrick",
                                 "divMigrate (Nm)" = "steelblue")) +
  scale_fill_manual(values = c("Genetic gravity (pGD similarity)" = "firebrick",
                               "divMigrate (Nm)" = "steelblue")) +
  facet_wrap(~ scenario) +
  labs(x = "Generation", y = expression("Spearman " * rho * " with true migration matrix"),
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

out_png <- "data/derived/fig-divmigrate-truth_reproduced.png"
ggsave(out_png, fig_truth, width = 9, height = 3.8, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_DivMigrateTruth to %s\n", out_png))
