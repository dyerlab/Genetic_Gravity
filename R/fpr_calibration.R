# Tab_FPRCalibration and Fig_FPRCalibration.
#
# Reproduces the calibration check of the fixed-topology permutation test
# across the Stream C Arm A symmetric differentiation gradient (Isotropic,
# sym-mid, sym-low, sym-verylow; N=100): whether the empirical edge-level
# false-positive rate stays flat as neutral differentiation rises, at three
# nominal alpha. Source data (data/fpr_calibration_summary.csv) is
# the per-scenario/per-alpha summary already computed from the individual-
# based simulation output by the private research repo's
# R/specificity_analysis.R (fpr_calibration()); no raw simulation data is
# required to run this script.
#
# Run from the repository root:
#   Rscript R/fpr_calibration.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

fc_summary <- read.csv("data/fpr_calibration_summary.csv", stringsAsFactors = FALSE)

scn_levels <- c("Isotropic", "sym-mid", "sym-low", "sym-verylow")
fc_summary <- fc_summary |> mutate(scenario = factor(scenario, levels = scn_levels))

## ---- Tab_FPRCalibration ---------------------------------------------------

tbl_fpr <- fc_summary |>
  mutate(acol = sprintf("FPR (alpha=%.2f)", alpha)) |>
  select(scenario, mean_cGD, median_cGD, acol, fpr_mean) |>
  pivot_wider(names_from = acol, values_from = fpr_mean) |>
  mutate(mean_cGD = round(mean_cGD, 3), median_cGD = round(median_cGD, 3),
         across(starts_with("FPR"), ~ round(.x, 3))) |>
  rename(Scenario = scenario, `Mean cGD` = mean_cGD, `Median cGD` = median_cGD) |>
  arrange(Scenario)

cat("=== Tab_FPRCalibration ===\n")
print(as.data.frame(tbl_fpr), row.names = FALSE)

## ---- Inline text values ---------------------------------------------------

meta <- read.csv("data/_meta.csv", stringsAsFactors = FALSE)
nperm          <- meta$value[meta$dataset == "fpr_calibration" & meta$key == "nperm"]
n_per_scenario <- meta$value[meta$dataset == "fpr_calibration" & meta$key == "n_per_scenario"]

s5 <- filter(fc_summary, abs(alpha - 0.05) < 1e-9)
s1 <- filter(fc_summary, abs(alpha - 0.01) < 1e-9)

cgd_fold  <- round(max(s5$mean_cGD) / min(s5$mean_cGD), 1)
ratio5    <- round(mean(s5$fpr_mean) / 0.05, 1)
ratio1    <- round(mean(s1$fpr_mean) / 0.01, 1)
int_ratio <- round(mean(s5$fpr_int_mean) / 0.05, 1)
mm        <- round(mean(s5$mean_cGD) / mean(s5$median_cGD), 2)
bw_lo     <- min(s5$fpr_bw_mean)
bw_hi     <- max(s5$fpr_bw_mean)
bw_ratio  <- round(mean(s5$fpr_bw_mean) / 0.05, 1)

cat(sprintf("\nnperm = %s, n_per_scenario (n_snapshots) = %s\n", nperm, n_per_scenario))
cat(sprintf("cgd_fold (mean_cGD fold-increase across s5 scenarios) = %s\n", cgd_fold))
cat(sprintf("FPR at alpha=0.05 range = %.3f-%.3f, ratio to nominal = %sx\n",
            min(s5$fpr_mean), max(s5$fpr_mean), ratio5))
cat(sprintf("FPR at alpha=0.01 ratio to nominal = %sx\n", ratio1))
cat(sprintf("Interior-only FPR range at alpha=0.05 = %.3f-%.3f, ratio = %sx\n",
            min(s5$fpr_int_mean), max(s5$fpr_int_mean), int_ratio))
cat(sprintf("mean/median cGD ratio = %.2f\n", mm))
cat(sprintf("Bandwidth-alignment null FPR range at alpha=0.05 = %.3f-%.3f, ratio = %sx\n",
            bw_lo, bw_hi, bw_ratio))

## ---- Fig_FPRCalibration ---------------------------------------------------

fig_fpr <- fc_summary |>
  ggplot(aes(alpha, fpr_mean, colour = scenario)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = pmax(0, fpr_mean - fpr_se), ymax = fpr_mean + fpr_se),
                width = 0.004, alpha = 0.7) +
  geom_point(size = 2.4) +
  scale_x_continuous(breaks = unique(fc_summary$alpha)) +
  coord_equal() +
  labs(x = "Nominal alpha", y = "Empirical false-positive rate", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

out_png <- "data/derived/fig-fpr-calibration_reproduced.png"
ggsave(out_png, fig_fpr, width = 5.5, height = 4.2, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_FPRCalibration to %s\n", out_png))
