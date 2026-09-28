# exploratory/lineage_phase_windows.R
#
# Per-lineage phase windows for annotating forward-phase figures.
#
# End of reorganization: each lineage's directional-IBGD advantage,
# dR2 = R^2(pGD ~ |dx| + |dx|:reverse) - R^2(cGD ~ |dx|), rises after the
# migration matrix changes at generation 2000 and then levels off. As in
# Fig_PeakAlignment (R/lineage_clock.R), each lineage's series is smoothed with
# a penalized spline, gam(y ~ s(generation, k = 10), REML), predicted every 5
# generations. From the lineage's own burn-in baseline (mean dR2 over
# 1904-1999) to its smoothed maximum:
#   t50, t90   first generation the smoothed curve completes 50% / 90% of that rise
#   peak_gen   generation of the smoothed maximum
# A lineage counts as having reorganized only if its rise exceeds the 95th
# percentile of the rises the same procedure finds in the Isotropic lineages
# (the alignment control, as in Fig_PeakAlignment).
# The same is repeated on the within-pGD gain (direction-aware minus
# direction-blind pGD) as a cross-check.
#
# Start of falling apart: each lineage's smoothed peak in mean |Delta|, exactly
# as in Fig_PeakAlignment (turned = peak by generation 2950).
#
# Sources: data/derived/ibgd_trajectory_fit.rda (exploratory/ibgd_trajectory.R),
# data/forward_scenarios.csv.
# Run from the repository root:  Rscript exploratory/lineage_phase_windows.R
# Outputs: data/derived/lineage_phase_windows.csv; media/fig-lineage-phase-windows-v2.png
# (earlier single-row variant: fig-lineage-phase-windows.png)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv); library(ggplot2) })

load("data/derived/ibgd_trajectory_fit.rda")
fwd <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
grid <- seq(2004, 2999, 5)
scen <- c("Isotropic", "Redistributed", "Obstructed")

smooth_fit <- function(g, y) {
  m <- gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML")
  as.numeric(predict(m, data.frame(g = grid)))
}
rise_times <- function(fit, base) {
  amp <- max(fit) - base
  at <- function(f) { i <- which(fit >= base + f * amp)[1]; if (is.na(i)) NA_real_ else grid[i] }
  c(amp = amp, t50 = at(0.5), t90 = at(0.9), peak_gen = grid[which.max(fit)], peak = max(fit))
}

d <- traj |> filter(status == "ok") |>
  mutate(dR2 = r2_pgd_dir - r2_cgd, gain = r2_pgd_dir - r2_pgd_blind)
base <- d |> filter(Treatment == "Burn-in") |> group_by(Replicate) |>
  summarise(base_dR2 = mean(dR2), base_gain = mean(gain), .groups = "drop")

fits <- list(); curves <- list()
for (s in scen) for (r in 1:50) {
  x <- d |> filter(Treatment == s, Replicate == r) |> arrange(generation)
  b <- base[base$Replicate == r, ]
  f1 <- smooth_fit(x$generation, x$dR2); f2 <- smooth_fit(x$generation, x$gain)
  a <- d_abs <- fwd |> filter(Scenario == s, replicate == r) |> arrange(generation)
  f3 <- smooth_fit(a$generation, a$mean_abs_delta)
  fits[[length(fits) + 1]] <- data.frame(Scenario = s, Replicate = r,
    t(setNames(rise_times(f1, b$base_dR2), paste0("dR2_", c("amp", "t50", "t90", "peak_gen", "peak")))),
    t(setNames(rise_times(f2, b$base_gain), paste0("gain_", c("amp", "t50", "t90", "peak_gen", "peak")))),
    absdelta_peak_gen = grid[which.max(f3)])
  curves[[length(curves) + 1]] <- data.frame(Scenario = s, Replicate = r, generation = grid, dR2_fit = f1)
}
fits <- do.call(rbind, fits); curves <- do.call(rbind, curves)

thr <- fits |> filter(Scenario == "Isotropic") |>
  summarise(dR2 = quantile(dR2_amp, .95), gain = quantile(gain_amp, .95))
fits <- fits |> mutate(reorg_dR2 = dR2_amp > thr$dR2, reorg_gain = gain_amp > thr$gain,
                       turned = absdelta_peak_gen <= 2950,
                       Scenario = factor(Scenario, scen))
write.csv(fits, "data/derived/lineage_phase_windows.csv", row.names = FALSE)

q <- function(v) if (all(is.na(v))) "-" else
  sprintf("%.0f (%.0f-%.0f)", median(v, na.rm = TRUE), quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE))
cat(sprintf("Isotropic-control rise thresholds (95th pct): dR2 %.3f, within-pGD gain %.3f\n\n", thr$dR2, thr$gain))
cat("Per-lineage phase markers, median (IQR) among lineages that reorganized:\n")
print(fits |> group_by(Scenario) |>
        summarise(reorganized = sum(reorg_dR2),
                  dR2_rise = sprintf("%.3f", median(dR2_amp[reorg_dR2])),
                  t50 = q(dR2_t50[reorg_dR2]), t90 = q(dR2_t90[reorg_dR2]),
                  dR2_peak = q(dR2_peak_gen[reorg_dR2]),
                  gain_t90 = q(gain_t90[reorg_gain]),
                  absdelta_peak_turned = q(absdelta_peak_gen[turned]), n_turned = sum(turned),
                  .groups = "drop") |> as.data.frame())
cat("\nWithin lineage (reorganized): cor(dR2 t90, gain t90), share with dR2 t90 < |Delta| peak:\n")
print(fits |> filter(Scenario != "Isotropic", reorg_dR2, reorg_gain) |> group_by(Scenario) |>
        summarise(n = n(), cor_t90 = cor(dR2_t90, gain_t90, use = "complete.obs"),
                  t90_before_absdelta_peak = mean(dR2_t90 < absdelta_peak_gen),
                  median_stable_span = median(pmin(absdelta_peak_gen, 2999) - dR2_t90)) |> as.data.frame(),
      digits = 3)

## ---- figure -----------------------------------------------------------------------
## Rows: the two rise measures (each lineage's smoothed curve, its t90 as a point,
## the IQR of t90 across reorganized lineages shaded blue). Orange: IQR of the
## smoothed |Delta| peak among lineages that turned (Fig_PeakAlignment).
metrics <- c(dR2 = "\u0394R\u00b2 (pGD vs cGD)", gain = "\u0394R\u00b2 (directional pGD vs pGD)")
asym <- c("Redistributed", "Obstructed")
curves_g <- list()
for (s in asym) for (r in 1:50) {
  x <- d |> filter(Treatment == s, Replicate == r) |> arrange(generation)
  curves_g[[length(curves_g) + 1]] <- data.frame(Scenario = s, Replicate = r, generation = grid,
                                                 gain = smooth_fit(x$generation, x$gain))
}
cv <- curves |> filter(Scenario %in% asym) |> left_join(do.call(rbind, curves_g),
                                                         by = c("Scenario", "Replicate", "generation")) |>
  pivot_longer(c(dR2_fit, gain), names_to = "metric", values_to = "fit") |>
  mutate(metric = factor(ifelse(metric == "dR2_fit", "dR2", "gain"), names(metrics), metrics),
         Scenario = factor(Scenario, asym))
mk <- fits |> filter(Scenario %in% asym) |>
  select(Scenario, Replicate, dR2 = dR2_t90, gain = gain_t90, reorg_dR2, reorg_gain) |>
  pivot_longer(c(dR2, gain), names_to = "metric", values_to = "t90") |>
  filter(ifelse(metric == "dR2", reorg_dR2, reorg_gain), !is.na(t90)) |>
  mutate(metric = factor(metric, names(metrics), metrics), Scenario = factor(as.character(Scenario), asym))
pts <- inner_join(mk, cv, by = c("Scenario", "Replicate", "metric", "t90" = "generation"))
reorg <- mk |> group_by(Scenario, metric) |>
  summarise(lo = quantile(t90, .25), hi = quantile(t90, .75), .groups = "drop") |>
  mutate(band = "End of reorganization (t90, IQR)")
coll <- fits |> filter(Scenario %in% asym, turned) |> group_by(Scenario) |>
  summarise(lo = quantile(absdelta_peak_gen, .25), hi = quantile(absdelta_peak_gen, .75), .groups = "drop") |>
  mutate(Scenario = factor(as.character(Scenario), asym)) |>
  tidyr::crossing(metric = factor(metrics, metrics)) |>
  mutate(band = "Start of falling apart (|\u0394| peak, IQR)")
bands <- bind_rows(reorg, coll)
p <- ggplot() +
  geom_rect(data = bands, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = band), alpha = 0.18) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(data = cv, aes(generation, fit, group = Replicate), colour = "grey25", alpha = 0.35,
            linewidth = 0.3) +
  geom_point(data = pts, aes(t90, fit), size = 1, colour = "grey10") +
  scale_fill_manual(values = c("#2a78d6", "#eb6834")) +
  facet_grid(metric ~ Scenario, scales = "free_y", switch = "y") +
  labs(x = "Generation", y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", strip.placement = "outside", panel.grid.minor = element_blank())
ggsave("media/fig-lineage-phase-windows-v2.png", p, width = 10, height = 6.5, dpi = 150)
cat("\nFigure: media/fig-lineage-phase-windows-v2.png\n")
