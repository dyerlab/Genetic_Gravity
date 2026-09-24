# Fig_Redistributed, Fig_Obstructed, and Fig_SignedAsymmetry.
#
# Forward-phase (gens 2004-2999, N=100) trends in graph structure and
# asymmetry under the two directional scenarios, each against the isotropic
# reference. Fig_Redistributed and Fig_Obstructed mirror Fig_IsotropicBaseline:
# C_D-bar, diameter, and mean |Delta_ij|, with the scenario's across-replicate
# +/- 1 SD as the shaded band, the across-replicate mean as a faint line, its
# fitted trend as the solid line, and the isotropic burn-in exponential fit
# (extrapolated across the forward phase, as in Fig_IsotropicBaseline) as the
# dashed reference. Fig_SignedAsymmetry stacks the signed graph mean Delta-bar
# for the two scenarios, drawn the same way, over the fitted isotropic trend.
#
# Treatment trends are penalized regression splines,
#     y ~ s(t, k = 20) + s(replicate, bs = "re"),
# fit by REML (mgcv) to every replicate snapshot, with the replicate random
# intercept excluded from the plotted prediction. A smoother is used rather
# than a parametric form because the forward trajectories are not approaches
# to a single equilibrium: Redistributed diameter is sigmoidal, its signed
# Delta-bar is biphasic (negative, then reversing to positive), and Obstructed
# diameter rises and then declines, so no single parametric family fits all
# of them without parameters running to their bounds. The isotropic signed
# Delta-bar is fit with the same smoother for Fig_SignedAsymmetry.
#
# Source data: data/forward_scenarios.csv (R/extract_forward_scenarios.R) and
# data/derived/isotropic_fit_params.csv (written by R/isotropic_baseline.R,
# which must be run first).
#
# Run from the repository root:
#   Rscript R/forward_scenarios.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
})

fwd     <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
iso_par <- read.csv("data/derived/isotropic_fit_params.csv", stringsAsFactors = FALSE)

GENS <- 2000:2999

stat_levels <- c("CD", "diameter", "mean_abs_delta")
stat_labels <- c("bar(C)[D]", "Diameter", "bar('|' * Delta * '|')")

## ---- Trend fits -----------------------------------------------------------

fit_trend <- function(df) {
  df$replicate <- factor(df$replicate)
  m <- gam(value ~ s(generation, k = 20) + s(replicate, bs = "re"),
           data = df, method = "REML")
  nd <- data.frame(generation = GENS, replicate = df$replicate[1])
  list(model = m,
       curve = data.frame(generation = GENS,
                          fit = as.numeric(predict(m, nd, exclude = "s(replicate)"))))
}

band <- function(df) {
  df |> group_by(generation) |>
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
}

# Isotropic reference: burn-in exponential fit, evaluated over the forward phase.
iso_curve <- expand_grid(stat = stat_levels, generation = GENS) |>
  left_join(iso_par, by = "stat") |>
  mutate(fit = yinf + (y0 - yinf) * exp(-generation / tau))

scenario_figure <- function(scn) {
  long <- fwd |>
    filter(Scenario == scn) |>
    pivot_longer(all_of(stat_levels), names_to = "stat", values_to = "value")

  fits <- lapply(split(long, long$stat), fit_trend)
  cat(sprintf("\n=== %s trend fits ===\n", scn))
  for (s in stat_levels)
    cat(sprintf("  %-15s edf = %5.2f   deviance explained = %.3f\n", s,
                sum(summary(fits[[s]]$model)$edf[1]),
                summary(fits[[s]]$model)$dev.expl))

  trend <- bind_rows(lapply(stat_levels, function(s)
    mutate(fits[[s]]$curve, stat = s)))
  bands <- long |> group_by(stat) |> group_modify(~ band(.x)) |> ungroup()

  lab <- function(d) mutate(d, Statistic = factor(stat, stat_levels, stat_labels))

  ggplot(lab(bands), aes(generation)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
    geom_line(aes(y = mean), colour = "grey20", alpha = 0.45, linewidth = 0.4) +
    geom_line(data = lab(iso_curve), aes(y = fit), colour = "grey30",
              linetype = "22", linewidth = 0.6) +
    geom_line(data = lab(trend), aes(y = fit), colour = "firebrick", linewidth = 0.7) +
    facet_grid(Statistic ~ ., scales = "free_y", labeller = label_parsed) +
    labs(x = "Generation", y = NULL) +
    theme_minimal(base_size = 11)
}

save_fig <- function(p, name, height) {
  out_png <- sprintf("data/derived/%s_reproduced.png", name)
  ggsave(out_png, p, width = 7, height = height, dpi = 150)
  cat(sprintf("Saved %s\n", out_png))
}

## ---- Fig_Redistributed / Fig_Obstructed -----------------------------------

save_fig(scenario_figure("Redistributed"), "fig-redistributed", 6.5)
save_fig(scenario_figure("Obstructed"),    "fig-obstructed",    6.5)

## ---- Fig_SignedAsymmetry --------------------------------------------------

signed <- fwd |> transmute(Scenario, replicate, generation, value = mean_delta)

iso_signed <- fit_trend(filter(signed, Scenario == "Isotropic"))$curve
cat(sprintf("\nIsotropic signed Delta-bar trend range: %.5f to %.5f\n",
            min(iso_signed$fit), max(iso_signed$fit)))

scn_levels <- c("Redistributed", "Obstructed")
sig_fit <- lapply(scn_levels, function(s) fit_trend(filter(signed, Scenario == s)))
names(sig_fit) <- scn_levels
cat("\n=== Signed Delta-bar trend fits ===\n")
for (s in scn_levels)
  cat(sprintf("  %-13s edf = %5.2f   deviance explained = %.3f\n", s,
              sum(summary(sig_fit[[s]]$model)$edf[1]),
              summary(sig_fit[[s]]$model)$dev.expl))

sig_trend <- bind_rows(lapply(scn_levels, function(s)
  mutate(sig_fit[[s]]$curve, Scenario = s)))
sig_band <- signed |> filter(Scenario %in% scn_levels) |>
  group_by(Scenario) |> group_modify(~ band(.x)) |> ungroup()
iso_ref <- expand_grid(Scenario = scn_levels, iso_signed)

fac <- function(d) mutate(d, Scenario = factor(Scenario, scn_levels))

fig_signed <- ggplot(fac(sig_band), aes(generation)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = mean), colour = "grey20", alpha = 0.45, linewidth = 0.4) +
  geom_line(data = fac(iso_ref), aes(y = fit), colour = "grey30",
            linetype = "22", linewidth = 0.6) +
  geom_line(data = fac(sig_trend), aes(y = fit), colour = "firebrick", linewidth = 0.7) +
  facet_grid(Scenario ~ .) +
  labs(x = "Generation", y = expression(bar(Delta))) +
  theme_minimal(base_size = 11)

save_fig(fig_signed, "fig-signed-asymmetry", 5)
