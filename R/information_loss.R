# Fig_InformationNull and Tab_InformationNull.
#
# Separates directional signal from information loss in the Redistributed
# scenario. As lineages approach fixation, Population Graphs are estimated
# from fewer informative loci, and noisier edge weights inflate |Delta|
# whatever the direction of gene flow. The null (data/information_null.csv,
# R/extract_information_null.R) rebuilds isotropic graphs -- no direction --
# from k randomly chosen polymorphic loci, giving the |Delta| and signed
# Delta-bar expected from information content alone. Redistributed snapshots
# are compared with that null at matched within-deme polymorphism (mean number
# of loci polymorphic within a deme, n_poly_pop). Matching on metapopulation-
# wide polymorphic loci instead flatters Redistributed, whose remaining
# polymorphic loci are near fixation within demes (He ~ 0.02-0.05 against
# ~ 0.23 in the null at the same count); even within-deme matching leaves the
# null's loci more informative, so late-phase excesses are upper bounds.
#
#   Fig_InformationNull  mean |Delta| and signed Delta-bar for every
#                        Redistributed snapshot (points) against within-deme
#                        polymorphic loci, over the null median (line) and
#                        90% range (band).
#   Tab_InformationNull  Redistributed vs null by polymorphism bin.
#
# The informative horizon is the within-deme polymorphism below which the
# null median of |Delta| exceeds the full-data isotropic plateau by more than
# HORIZON_TOL; it is reported as the generation (median across replicates) at
# which each Redistributed lineage first falls below it.
#
# Source data: data/forward_scenarios.csv, data/information_null.csv, and
# data/derived/isotropic_fit_params.csv (R/isotropic_baseline.R).
#
# Run from the repository root:
#   Rscript R/information_loss.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

fwd     <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
nl      <- read.csv("data/information_null.csv", stringsAsFactors = FALSE) |> filter(ok == 1)
iso_par <- read.csv("data/derived/isotropic_fit_params.csv", stringsAsFactors = FALSE)

red <- filter(fwd, Scenario == "Redistributed")
iso_abs_plateau <- iso_par$yinf[iso_par$stat == "mean_abs_delta"]

BREAKS      <- c(0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 6, 8, 10)
HORIZON_TOL <- 0.10

## ---- Null by within-deme polymorphism -------------------------------------

null_bins <- nl |>
  mutate(bin = cut(n_poly_pop, BREAKS)) |>
  filter(!is.na(bin)) |>
  group_by(bin) |>
  summarise(n_null   = n(),
            x        = mean(n_poly_pop),
            abs_med  = median(mean_abs_delta),
            abs_q05  = quantile(mean_abs_delta, 0.05),
            abs_q95  = quantile(mean_abs_delta, 0.95),
            sgn_med  = median(mean_delta),
            sgn_q05  = quantile(mean_delta, 0.05),
            sgn_q95  = quantile(mean_delta, 0.95),
            .groups  = "drop")

## ---- Tab_InformationNull --------------------------------------------------

tbl_null <- red |>
  mutate(bin = cut(n_poly_pop, BREAKS)) |>
  inner_join(null_bins, by = "bin") |>
  group_by(bin) |>
  summarise(`Snapshots`                = n(),
            `Median generation`        = median(generation),
            `Redistributed |Delta|`    = median(mean_abs_delta),
            `Null |Delta|`             = first(abs_med),
            `Ratio`                    = median(mean_abs_delta) / first(abs_med),
            `> null 95th (|Delta|)`    = mean(mean_abs_delta > abs_q95),
            `Redistributed Delta-bar`  = median(mean_delta),
            `> null 95th (Delta-bar)`  = mean(mean_delta > sgn_q95),
            `< null 5th (Delta-bar)`   = mean(mean_delta < sgn_q05),
            .groups = "drop") |>
  arrange(desc(bin)) |>
  rename(`Polymorphic loci / deme` = bin) |>
  mutate(across(c(`Redistributed |Delta|`, `Null |Delta|`, `Redistributed Delta-bar`), ~ round(.x, 4)),
         across(c(Ratio, starts_with(">"), starts_with("<")), ~ round(.x, 2)))

cat("=== Tab_InformationNull ===\n")
print(as.data.frame(tbl_null), row.names = FALSE)
cat(sprintf("Redistributed snapshots below the lowest null bin: %d\n",
            sum(red$n_poly_pop <= min(BREAKS))))

## ---- Inline text values: informative horizon ------------------------------

inflated <- null_bins |> filter(abs_med > (1 + HORIZON_TOL) * iso_abs_plateau)
rv_il_threshold <- BREAKS[max(as.integer(inflated$bin)) + 1]   # upper edge of highest inflated bin

horizon <- red |>
  group_by(replicate) |>
  summarise(gen   = min(c(generation[n_poly_pop < rv_il_threshold], NA), na.rm = TRUE),
            ae    = sumAe_total[generation == gen][1],
            .groups = "drop")

rv_il_horizon_gen <- median(horizon$gen, na.rm = TRUE)
rv_il_horizon_iqr <- quantile(horizon$gen, c(0.25, 0.75), na.rm = TRUE)
rv_il_horizon_ae  <- median(horizon$ae, na.rm = TRUE)

cat(sprintf(paste0(
  "\nrv_il_threshold   = %.2f polymorphic loci / deme (null |Delta| within %.0f%% of isotropic plateau above this)",
  "\nrv_il_horizon_gen = %.0f (median across replicates; IQR %.0f-%.0f; %d of %d lineages cross by gen 2999)",
  "\nrv_il_horizon_ae  = %.2f (median sumAe_total at the horizon)\n"),
  rv_il_threshold, 100 * HORIZON_TOL,
  rv_il_horizon_gen, rv_il_horizon_iqr[1], rv_il_horizon_iqr[2],
  sum(!is.na(horizon$gen)), nrow(horizon),
  rv_il_horizon_ae))

## ---- Fig_InformationNull --------------------------------------------------

lev <- c("mean_abs_delta", "mean_delta")
lab <- c("bar('|' * Delta * '|')", "bar(Delta)")

pts <- red |>
  pivot_longer(all_of(lev), names_to = "stat", values_to = "value") |>
  mutate(Statistic = factor(stat, lev, lab))
band <- bind_rows(
  transmute(null_bins, x, med = abs_med, lo = abs_q05, hi = abs_q95, stat = "mean_abs_delta"),
  transmute(null_bins, x, med = sgn_med, lo = sgn_q05, hi = sgn_q95, stat = "mean_delta")) |>
  mutate(Statistic = factor(stat, lev, lab))
zero <- data.frame(Statistic = factor("bar(Delta)", lab), y = 0)

fig_null <- ggplot(pts, aes(n_poly_pop, value)) +
  geom_hline(data = zero, aes(yintercept = y), colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = rv_il_threshold, linetype = "22", colour = "grey30") +
  geom_point(colour = "firebrick", alpha = 0.12, size = 0.5) +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE,
              fill = "grey30", alpha = 0.3) +
  geom_line(data = band, aes(x, med), inherit.aes = FALSE, linewidth = 0.8) +
  scale_x_reverse() +
  facet_grid(Statistic ~ ., scales = "free_y", labeller = label_parsed) +
  labs(x = "Polymorphic loci per deme (of 20)", y = NULL) +
  theme_minimal(base_size = 11)

out_png <- "data/derived/fig-information-null_reproduced.png"
ggsave(out_png, fig_null, width = 7, height = 5.5, dpi = 150)
cat(sprintf("\nSaved %s\n", out_png))
