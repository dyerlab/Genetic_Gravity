# Tab_DesignDrivenTest.
#
# Reproduces the design-driven directional test applied per replicate to a
# single final-generation graph (N=100, generation 2999): a one-tailed
# Wilcoxon signed-rank test on arc weights (centered on symmetry, 0.5),
# edges oriented by ascending stepping-stone position, testing for a
# downstream bias. Source data (data/directional_sweep_summary.csv,
# data/directional_sweep_per_snapshot.csv) are the per-scenario and
# per-replicate summaries already computed from the individual-based
# simulation output by the private research repo's R/directional_sweep.R;
# no raw simulation data is required to run this script.
#
# Run from the repository root:
#   Rscript R/directional_sweep.R

suppressPackageStartupMessages({
  library(dplyr)
})

ds_summary      <- read.csv("data/directional_sweep_summary.csv", stringsAsFactors = FALSE)
ds_per_snapshot <- read.csv("data/directional_sweep_per_snapshot.csv", stringsAsFactors = FALSE)

scn <- c("Isotropic", "Flux-conserved", "Rate-conserved")
ds_summary <- ds_summary |> mutate(scenario = factor(scenario, levels = scn)) |> arrange(scenario)

## ---- Tab_DesignDrivenTest --------------------------------------------------

tbl_dirsweep <- ds_summary |>
  transmute(
    Scenario = scenario,
    Replicates = n,
    `Significant (alpha=0.05)` = sprintf("%.0f%%", 100 * frac_sig),
    `Median forward bias` = sprintf("%+.3f", median_estimate))

cat("=== Tab_DesignDrivenTest ===\n")
print(as.data.frame(tbl_dirsweep), row.names = FALSE)

## ---- Inline text values ---------------------------------------------------

alpha <- 0.05  # data/_meta.csv: directional_sweep alpha = 0.05
g <- function(cc) ds_summary[ds_summary$cond == cc, , drop = FALSE]
iso  <- g(1); flux <- g(2); rate <- g(3)
ne <- round(mean(ds_per_snapshot$n_edges, na.rm = TRUE))

cat(sprintf(
  "\nIsotropic: median forward bias = %+.3f, over-rejection = %.0f%% (nominal %.0f%%)\n",
  iso$median_estimate, 100 * iso$frac_sig, 100 * alpha))
cat(sprintf(
  "Flux-conserved: median forward bias = %+.3f, detected in %.0f%% of snapshots\n",
  flux$median_estimate, 100 * flux$frac_sig))
cat(sprintf(
  "Rate-conserved: median forward bias = %+.3f, detected in %.0f%% of snapshots\n",
  rate$median_estimate, 100 * rate$frac_sig))
cat(sprintf("Mean edges per graph (ne) = %d\n", ne))
