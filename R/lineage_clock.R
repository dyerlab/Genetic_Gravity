# Fig_LineageDiversity, Fig_LineageClock, Fig_PeakAlignment,
# Tab_LineageDiversity, and Tab_LineageFixation.
#
# Each replicate is an independent evolutionary lineage: without mutation,
# drift and migration drive every lineage toward fixation, but at different
# rates, at different loci, and to different alleles. At a given generation
# the Redistributed replicates therefore occupy genetically distinct states,
# and the across-replicate band in Fig_Redistributed is variation among
# lineages, not noise around a common trajectory.
#
# The lineage clock is sumAe_total, the sum over all 20 loci of the
# metapopulation effective number of alleles, A_e = 1 / sum(p_i^2) (20 = every
# locus fixed; 40 = every locus at 50:50). It measures how much variation a
# lineage has left without reference to which loci or alleles carry it --
# appropriate because Delta_ij is built from distances and is invariant to
# allele labelling. sumAe_total - sumAe_within (the within-deme sum, averaged
# over demes) is the among-deme component: the covariance the Population Graph
# is built from.
#
#   Fig_LineageDiversity  sumAe_total and sumAe_within through time for every
#                         Redistributed replicate, over the isotropic
#                         replicates (grey).
#   Fig_LineageClock      mean |Delta| and signed Delta-bar of every
#                         Redistributed replicate against its own sumAe_total
#                         (25-generation windows per replicate, faint lines),
#                         with a penalized-spline trend across replicates;
#                         dashed line = isotropic plateau of mean |Delta|
#                         (Fig_IsotropicBaseline).
#   Fig_PeakAlignment     mean |Delta| of every Redistributed replicate
#                         aligned at its own (smoothed) peak, with the
#                         isotropic alignment control (see the section below).
#
# Source data: data/forward_scenarios.csv (R/extract_forward_scenarios.R) and
# data/derived/isotropic_fit_params.csv (R/isotropic_baseline.R).
#
# Run from the repository root:
#   Rscript R/lineage_clock.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
})

fwd     <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
iso_par <- read.csv("data/derived/isotropic_fit_params.csv", stringsAsFactors = FALSE)

red <- filter(fwd, Scenario == "Redistributed")
iso <- filter(fwd, Scenario == "Isotropic")

iso_abs_plateau <- iso_par$yinf[iso_par$stat == "mean_abs_delta"]

## ---- Tab_LineageDiversity -------------------------------------------------

tab_gens <- c(2004, 2254, 2504, 2754, 2999)

div_row <- function(d) {
  d |> filter(generation %in% tab_gens) |>
    group_by(generation) |>
    summarise(`sumAe (pooled)`       = median(sumAe_total),
              `sumAe (within)`       = median(sumAe_within),
              `Among-deme component` = median(sumAe_total - sumAe_within),
              `SD across replicates` = sd(sumAe_total),
              `Polymorphic loci / deme` = mean(n_poly_pop),
              .groups = "drop")
}

tbl_div <- bind_rows(mutate(div_row(red), Scenario = "Redistributed"),
                     mutate(div_row(iso), Scenario = "Isotropic")) |>
  relocate(Scenario) |>
  mutate(across(where(is.numeric) & !generation, ~ round(.x, 2)))

cat("=== Tab_LineageDiversity (medians across replicates) ===\n")
print(as.data.frame(tbl_div), row.names = FALSE)

## ---- Tab_LineageFixation --------------------------------------------------

# Per replicate: generation of its own mean |Delta| peak (50-generation
# windows), and how far it had progressed toward fixation by gen 2899.
per_rep <- red |>
  mutate(win = floor((generation - 2000) / 50)) |>
  group_by(replicate, win) |>
  summarise(a = mean(mean_abs_delta), .groups = "drop") |>
  group_by(replicate) |>
  summarise(peak_gen = 2000 + 50 * win[which.max(a)] + 25, .groups = "drop") |>
  inner_join(red |> filter(generation == 2899) |> select(replicate, ae_2899 = sumAe_total),
             by = "replicate") |>
  mutate(Group = cut(ae_2899, quantile(ae_2899, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                     labels = c("Closest to fixation", "Intermediate", "Least fixed")))

tbl_fix <- per_rep |>
  group_by(Group) |>
  summarise(Replicates = n(),
            `sumAe at gen 2899` = round(median(ae_2899), 2),
            `Peaked before gen 2900` = sprintf("%.0f%%", 100 * mean(peak_gen < 2900)),
            .groups = "drop")

cat("\n=== Tab_LineageFixation (terciles of sumAe_total at gen 2899) ===\n")
print(as.data.frame(tbl_fix), row.names = FALSE)

## ---- Inline text values ---------------------------------------------------

# Landmarks on the lineage clock from penalized-spline trends across replicates.
grid <- data.frame(sumAe_total = seq(max(red$sumAe_total), min(red$sumAe_total), by = -0.01))
trend_on_clock <- function(col) {
  m <- gam(reformulate("s(sumAe_total, k = 15)", col), data = red, method = "REML")
  as.numeric(predict(m, grid))
}
abs_tr <- trend_on_clock("mean_abs_delta")
sgn_tr <- trend_on_clock("mean_delta")

i_min  <- which.min(sgn_tr)
i_zero <- i_min - 1 + which(sgn_tr[i_min:length(sgn_tr)] > 0)[1]

rv_lc_sgn_min      <- sgn_tr[i_min]
rv_lc_sgn_min_ae   <- grid$sumAe_total[i_min]
rv_lc_sgn_zero_ae  <- grid$sumAe_total[i_zero]
rv_lc_abs_peak_ae  <- grid$sumAe_total[which.max(abs_tr)]
rv_lc_abs_peak     <- max(abs_tr)
rv_lc_cor_peak     <- cor(per_rep$ae_2899, per_rep$peak_gen)

# Share of between-replicate spread removed by comparing lineages at equal
# remaining diversity rather than equal generation (20 equal-size bins each).
bin_sd <- function(key, col) {
  red |> mutate(b = ntile(key, 20)) |> group_by(b) |>
    summarise(s = sd(.data[[col]]), .groups = "drop") |> pull(s) |> mean()
}
rv_lc_sd_gen_abs <- bin_sd(red$generation,   "mean_abs_delta")
rv_lc_sd_ae_abs  <- bin_sd(-red$sumAe_total, "mean_abs_delta")
rv_lc_sd_gen_sgn <- bin_sd(red$generation,   "mean_delta")
rv_lc_sd_ae_sgn  <- bin_sd(-red$sumAe_total, "mean_delta")

cat(sprintf(paste0(
  "\nrv_lc_sgn_min      = %.4f (signed Delta-bar minimum, at sumAe = %.2f)",
  "\nrv_lc_sgn_zero_ae  = %.2f (signed Delta-bar crosses zero)",
  "\nrv_lc_abs_peak     = %.4f (mean |Delta| trend peak, at sumAe = %.2f)",
  "\nrv_lc_cor_peak     = %.2f (cor of sumAe at gen 2899 with replicate peak generation)",
  "\nrv_lc_sd |Delta|   = %.4f by generation vs %.4f by sumAe (%.0f%% reduction)",
  "\nrv_lc_sd signed    = %.4f by generation vs %.4f by sumAe (%.0f%% reduction)\n"),
  rv_lc_sgn_min, rv_lc_sgn_min_ae, rv_lc_sgn_zero_ae,
  rv_lc_abs_peak, rv_lc_abs_peak_ae, rv_lc_cor_peak,
  rv_lc_sd_gen_abs, rv_lc_sd_ae_abs, 100 * (1 - rv_lc_sd_ae_abs / rv_lc_sd_gen_abs),
  rv_lc_sd_gen_sgn, rv_lc_sd_ae_sgn, 100 * (1 - rv_lc_sd_ae_sgn / rv_lc_sd_gen_sgn)))

## ---- Fig_LineageDiversity -------------------------------------------------

lev <- c("sumAe_total", "sumAe_within")
lab <- c("Metapopulation", "Within~deme~(mean)")

ae_long <- bind_rows(red, iso) |>
  select(Scenario, replicate, generation, all_of(lev)) |>
  pivot_longer(all_of(lev), names_to = "level", values_to = "sumAe") |>
  mutate(level = factor(level, lev, lab))

fig_div <- ggplot(filter(ae_long, Scenario == "Redistributed"),
                  aes(generation, sumAe, group = replicate)) +
  geom_line(data = filter(ae_long, Scenario == "Isotropic"),
            colour = "grey75", alpha = 0.4, linewidth = 0.3) +
  geom_line(colour = "firebrick", alpha = 0.45, linewidth = 0.3) +
  geom_hline(yintercept = 20, linetype = "22", colour = "grey30") +
  facet_grid(level ~ ., labeller = label_parsed) +
  labs(x = "Generation", y = expression(sum(A[e], loci))) +
  theme_minimal(base_size = 11)

## ---- Fig_LineageClock -----------------------------------------------------

win <- red |>
  mutate(w = floor((generation - 2000) / 25)) |>
  group_by(replicate, w) |>
  summarise(generation = mean(generation), sumAe = mean(sumAe_total),
            mean_abs_delta = mean(mean_abs_delta), mean_delta = mean(mean_delta),
            .groups = "drop")

clk_lev <- c("mean_abs_delta", "mean_delta")
clk_lab <- c("bar('|' * Delta * '|')", "bar(Delta)")

clk <- win |>
  pivot_longer(all_of(clk_lev), names_to = "stat", values_to = "value") |>
  mutate(Statistic = factor(stat, clk_lev, clk_lab))
trend <- bind_rows(data.frame(sumAe = grid$sumAe_total, value = abs_tr, stat = "mean_abs_delta"),
                   data.frame(sumAe = grid$sumAe_total, value = sgn_tr, stat = "mean_delta")) |>
  mutate(Statistic = factor(stat, clk_lev, clk_lab))
refs <- data.frame(stat = clk_lev, y = c(iso_abs_plateau, 0),
                   lt = c("22", "solid"), col = c("grey30", "grey60")) |>
  mutate(Statistic = factor(stat, clk_lev, clk_lab))

fig_clock <- ggplot(clk, aes(sumAe, value)) +
  geom_hline(data = refs, aes(yintercept = y, linetype = lt, colour = col), linewidth = 0.4) +
  scale_linetype_identity() +
  scale_colour_identity() +
  geom_path(aes(group = replicate), colour = "grey20", alpha = 0.3, linewidth = 0.3) +
  geom_line(data = trend, colour = "firebrick", linewidth = 0.8) +
  scale_x_reverse() +
  facet_grid(Statistic ~ ., scales = "free_y", labeller = label_parsed) +
  labs(x = expression(sum(A[e], loci) ~ "(metapopulation): remaining diversity"),
       y = NULL) +
  theme_minimal(base_size = 11)

## ---- Fig_PeakAlignment ----------------------------------------------------

# Each replicate's mean |Delta| trajectory is aligned at its own peak. The peak
# is taken from a per-replicate penalized spline rather than the raw maximum,
# since aligning noisy series at their maxima manufactures a spike and a
# post-peak decline even from pure noise. The same procedure is applied to the
# isotropic replicates, whose mean |Delta| is stationary, so their aligned
# median shows how much rise-and-fall the alignment itself produces. Aligned
# medians are drawn only where at least MIN_REPS replicates contribute.

MIN_REPS <- 5

align_at_peak <- function(d) {
  peaks <- d |>
    group_by(replicate) |>
    group_modify(function(x, key) {
      m  <- gam(mean_abs_delta ~ s(generation, k = 10), data = x, method = "REML")
      gg <- seq(min(x$generation), max(x$generation), by = 5)
      p  <- as.numeric(predict(m, data.frame(generation = gg)))
      tibble(peak_gen = gg[which.max(p)], peak_fit = max(p), end_fit = p[length(p)])
    }) |>
    ungroup()
  aligned <- d |>
    mutate(w = floor((generation - 2000) / 25)) |>
    group_by(replicate, w) |>
    summarise(generation = mean(generation), value = mean(mean_abs_delta), .groups = "drop") |>
    inner_join(peaks, by = "replicate") |>
    mutate(rel = 25 * round((generation - peak_gen) / 25))
  list(peaks = peaks, aligned = aligned,
       median = aligned |> group_by(rel) |>
         summarise(value = median(value), n = n_distinct(replicate), .groups = "drop") |>
         filter(n >= MIN_REPS))
}

pa_red <- align_at_peak(red)
pa_iso <- align_at_peak(iso)

turned <- pa_red$peaks |> filter(peak_gen <= 2950)
cat(sprintf(paste0(
  "\nrv_pa_turned      = %d of %d redistributed lineages peak (smoothed) by gen 2950",
  "\nrv_pa_peak_gen    = %.0f (median smoothed peak generation of those lineages)",
  "\nrv_pa_drop        = %.0f%% (median decline from smoothed peak to gen 2999 in those lineages)",
  "\nrv_pa_iso_turned  = %d of %d isotropic lineages (alignment control)",
  "\nrv_pa_iso_drop    = %.0f%% (median decline from smoothed peak, isotropic control)\n"),
  nrow(turned), nrow(pa_red$peaks), median(turned$peak_gen),
  100 * median(1 - turned$end_fit / turned$peak_fit),
  sum(pa_iso$peaks$peak_gen <= 2950), nrow(pa_iso$peaks),
  100 * median(1 - with(filter(pa_iso$peaks, peak_gen <= 2950), end_fit / peak_fit))))

fig_peak <- ggplot(pa_red$aligned, aes(rel, value)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey30") +
  geom_line(aes(group = replicate), colour = "grey20", alpha = 0.3, linewidth = 0.3) +
  geom_line(data = filter(pa_iso$median, rel <= max(pa_red$aligned$rel)),
            colour = "grey30", linetype = "22", linewidth = 0.7) +
  geom_line(data = pa_red$median, colour = "firebrick", linewidth = 0.9) +
  labs(x = "Generations relative to each replicate's own peak",
       y = expression(bar("|" * Delta * "|"))) +
  theme_minimal(base_size = 11)

save_fig <- function(p, name, height) {
  out_png <- sprintf("data/derived/%s_reproduced.png", name)
  ggsave(out_png, p, width = 7, height = height, dpi = 150)
  cat(sprintf("Saved %s\n", out_png))
}

save_fig(fig_div,   "fig-lineage-diversity", 5.5)
save_fig(fig_clock, "fig-lineage-clock",     5.5)
save_fig(fig_peak,  "fig-peak-alignment",    4)
