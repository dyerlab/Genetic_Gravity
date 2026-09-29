# ibgd_fit_v2.R
# Phase strips use the heterozygosity-gradient reorganization boundary (2026-09-29); earlier file = gain-based boundary.
# Fig_DirectionalIBGD with phase-share strips -> media/fig-ibgd-fit-v3.png
# Top: standard IBGD R^2 (cGD ~ |dx|). Middle: the within-pGD directional gain,
# direction-aware minus direction-blind pGD adjusted R^2 (what the caption and
# text describe; the earlier media/fig-ibgd-fit.png plotted pGD - cGD Delta R^2
# instead). Median and IQR by generation, burn-in shared (as in
# exploratory/ibgd_trajectory.R). Bottom: share of Redistributed and Obstructed
# lineages in each phase (exploratory/lineage_phases.R).
# Input: data/derived/ibgd_trajectory_fit.rda. Run from the manuscript repo root.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
source("exploratory/lineage_phases.R")
load("data/derived/ibgd_trajectory_fit.rda")
scen <- c("Isotropic", "Redistributed", "Obstructed")
cols <- c(Isotropic = "#2a78d6", Redistributed = "#eb6834", Obstructed = "#1baf7a")
ink <- "#3d3d3a"
d <- traj |> filter(status == "ok") |> mutate(Treatment = as.character(Treatment), gain = r2_pgd_dir - r2_pgd_blind)
d <- bind_rows(bind_rows(lapply(scen, function(s) filter(d, Treatment == "Burn-in") |> mutate(Treatment = s))),
               filter(d, Treatment != "Burn-in")) |>
  select(Treatment, generation, r2_cgd, gain) |>
  pivot_longer(c(r2_cgd, gain), names_to = "metric") |>
  mutate(metric = factor(metric, c("r2_cgd", "gain"),
                         c("Standard IBGD R² (cGD)", "Directional gain in R² (pGD)")),
         Treatment = factor(Treatment, c("Redistributed", "Obstructed", "Isotropic"))) |>
  group_by(Treatment, metric, generation) |>
  summarise(med = median(value), lo = quantile(value, .25), hi = quantile(value, .75), .groups = "drop")

xs <- scale_x_continuous(limits = c(1900, 3000), breaks = seq(2000, 3000, 250), expand = c(0, 0))
th <- theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = ink),
                                            strip.placement = "outside", strip.text.y.left = element_text(angle = 90, colour = ink))
top <- ggplot(d, aes(generation, med, colour = Treatment, fill = Treatment)) +
  geom_vline(xintercept = 2000, colour = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols, breaks = scen) + scale_fill_manual(values = cols, breaks = scen) +
  facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left") +
  xs + labs(x = NULL, y = NULL, colour = NULL, fill = NULL) + th +
  theme(legend.position = "top", axis.text.x = element_blank())
strip <- lp_strip(c("Redistributed", "Obstructed"), xs, facet = "row", strip_labels = TRUE) +
  geom_vline(xintercept = 2000, colour = "grey60", linewidth = 0.4) +
  labs(x = "Generation") + theme(strip.text.y = element_text(colour = ink, size = 9, angle = 0, hjust = 0), panel.spacing.y = unit(0.6, "lines"))
p <- top / strip + plot_layout(heights = c(6, 1.5))
ggsave("media/fig-ibgd-fit-v3.png", p, width = 8.5, height = 7.4, dpi = 150, bg = "white")
cat("Saved media/fig-ibgd-fit-v3.png\n")
