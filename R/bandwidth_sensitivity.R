# Tab_BandwidthSensitivity.
#
# Reproduces the manuscript's bandwidth-sensitivity table: how well the
# canonical local bandwidth estimator's edge-level asymmetry sign/rank and
# graph-mean sign survive four non-local bandwidth alternatives, pooled over
# the two directional scenarios (Flux-conserved, Rate-conserved; N=100,
# generation 2999). Source data (data/bandwidth_sensitivity_per_snapshot.csv)
# is the per-replicate/per-alternative comparison already computed from the
# individual-based simulation output by the private research repo's
# R/bandwidth_sensitivity.R; no raw simulation data is required to run this
# script. Pooling logic (filter to the two directional scenarios, mean
# sign_concord/dbar_pres, median rho_signed/rho_abs, grouped by alternative)
# matches paper1.qmd's `fit-bandwidth-sensitivity` chunk exactly.
#
# Run from the repository root:
#   Rscript R/bandwidth_sensitivity.R

suppressPackageStartupMessages({
  library(dplyr)
})

per_snapshot <- read.csv("data/bandwidth_sensitivity_per_snapshot.csv", stringsAsFactors = FALSE)

alt_levels <- c("global", "scale x0.5", "scale x2", "perplexity")

# Pool the two directional scenarios (cond 2 = Flux-conserved, cond 3 =
# Rate-conserved), where Delta carries real signal; the isotropic null (cond 1)
# is excluded because under symmetric migration the sign of Delta is noise.
bw_tab <- per_snapshot |>
  filter(cond %in% c(2L, 3L)) |>
  group_by(alternative) |>
  summarise(sign_concord = mean(sign_concord, na.rm = TRUE),
            rho_signed   = median(rho_signed, na.rm = TRUE),
            rho_abs      = median(rho_abs, na.rm = TRUE),
            dbar_pres    = mean(dbar_pres, na.rm = TRUE), .groups = "drop") |>
  mutate(alternative = factor(alternative, levels = alt_levels)) |>
  arrange(alternative)

tbl_bandwidth <- bw_tab |>
  transmute(
    `Non-local bandwidth` = recode(as.character(alternative),
        global = "Fixed global", `scale x0.5` = "Local x 1/2",
        `scale x2` = "Local x 2", perplexity = "Perplexity (sigma search)"),
    `Sign concordance` = sprintf("%.1f%%", 100 * sign_concord),
    `Rank rho (signed)` = sprintf("%.3f", rho_signed),
    `Rank rho (|Delta|)` = sprintf("%.3f", rho_abs),
    `Delta-bar sign kept` = sprintf("%.0f%%", 100 * dbar_pres))

cat("=== Tab_BandwidthSensitivity ===\n")
print(as.data.frame(tbl_bandwidth), row.names = FALSE)

## ---- Inline text values ---------------------------------------------------

sc_lo <- min(bw_tab$sign_concord); sc_hi <- max(bw_tab$sign_concord)
rs_lo <- min(bw_tab$rho_signed);   rs_hi <- max(bw_tab$rho_signed)
flux_dbar <- mean(per_snapshot$dbar_pres[per_snapshot$cond == 2L], na.rm = TRUE)
rate_dbar <- mean(per_snapshot$dbar_pres[per_snapshot$cond == 3L], na.rm = TRUE)

cat(sprintf(
  "\nsign concordance range: %.0f%%-%.0f%%\n", 100 * sc_lo, 100 * sc_hi))
cat(sprintf(
  "signed rank rho range: %.2f-%.2f\n", rs_lo, rs_hi))
cat(sprintf(
  "Delta-bar sign preserved: %.0f%% flux-conserved, %.0f%% rate-conserved\n",
  100 * flux_dbar, 100 * rate_dbar))
