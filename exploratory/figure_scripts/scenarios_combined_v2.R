# scenarios_combined_v2.R
# Phase strips use the heterozygosity-gradient reorganization boundary (2026-09-29); earlier file = gain-based boundary.
# Fig_ScenarioComparison with a phase-share strip: the 3x2 facet of C_D-bar,
# diameter and mean |Delta| (Redistributed | Obstructed, as in
# scenarios_combined.R), plus a bottom row showing, at each census, the share of
# the 50 lineages in reorganization, plateau and falling apart
# (exploratory/lineage_phases.R). Each lineage changes phase at its own
# generation, so the strip shows how many lineages each calendar-time median mixes.
# -> media/fig-scenarios-combined-v3.png (earlier: fig-scenarios-combined.png)
# Run from the manuscript repo root.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(mgcv); library(patchwork) })
source("exploratory/lineage_phases.R")

fwd     <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
iso_par <- read.csv("data/derived/isotropic_fit_params.csv", stringsAsFactors = FALSE)
GENS <- 2000:2999
stat_levels <- c("CD", "diameter", "mean_abs_delta")
stat_labels <- c("bar(C)[D]", "Diameter", "bar('|' * Delta * '|')")
scn_levels  <- c("Redistributed", "Obstructed")
PHASE_COL <- c(Reorganize = "#86b6ef", Plateau = "#2a78d6", `Fall apart` = "#104281")   # ordinal blue ramp (validated)
ink <- "#3d3d3a"

fit_trend <- function(df) {
  df$replicate <- factor(df$replicate)
  m <- gam(value ~ s(generation, k = 20) + s(replicate, bs = "re"), data = df, method = "REML")
  nd <- data.frame(generation = GENS, replicate = df$replicate[1])
  data.frame(generation = GENS, fit = as.numeric(predict(m, nd, exclude = "s(replicate)")))
}
long <- fwd |> filter(Scenario %in% scn_levels) |> pivot_longer(all_of(stat_levels), names_to = "stat", values_to = "value")
trend <- long |> group_by(Scenario, stat) |> group_modify(~ fit_trend(.x)) |> ungroup()
bands <- long |> group_by(Scenario, stat, generation) |> summarise(mean = mean(value), sd = sd(value), .groups = "drop")
iso_curve <- expand_grid(Scenario = scn_levels, stat = stat_levels, generation = GENS) |>
  left_join(iso_par, by = "stat") |> mutate(fit = yinf + (y0 - yinf) * exp(-generation / tau))
lab <- function(d) mutate(d, Statistic = factor(stat, stat_levels, stat_labels), Scenario = factor(Scenario, scn_levels))

xs <- scale_x_continuous(limits = c(2000, 3000), breaks = seq(2000, 2750, 250), expand = c(0, 0))
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(1.2, "lines"),
        axis.text = element_text(colour = ink), strip.text = element_text(colour = ink))

top <- ggplot(lab(bands), aes(generation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = mean), colour = "grey20", alpha = 0.45, linewidth = 0.4) +
  geom_line(data = lab(iso_curve), aes(y = fit), colour = "grey30", linetype = "22", linewidth = 0.6) +
  geom_line(data = lab(trend), aes(y = fit), colour = "firebrick", linewidth = 0.7) +
  facet_grid(Statistic ~ Scenario, scales = "free_y",
             labeller = labeller(Statistic = label_parsed, Scenario = label_value)) +
  xs + labs(x = NULL, y = NULL) + th + theme(axis.text.x = element_blank())

share <- expand.grid(Replicate = 1:50, Scenario = scn_levels, generation = seq(2004, 2999, 5), stringsAsFactors = FALSE) |>
  mutate(phase = droplevels(lp_phase(Scenario, Replicate, generation))) |>
  count(Scenario, generation, phase) |> complete(Scenario, generation, phase, fill = list(n = 0)) |>
  mutate(share = n / 50, Scenario = factor(Scenario, scn_levels), phase = factor(phase, names(PHASE_COL)))

bot <- ggplot(share, aes(generation, share, fill = phase)) +
  geom_area(position = position_stack(reverse = TRUE), colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = PHASE_COL, name = "Lineages in phase",
                    labels = c(Reorganize = "Reorganizing", Plateau = "Plateau", `Fall apart` = "Falling apart")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "25", "50"), expand = c(0, 0)) +
  facet_grid(. ~ Scenario) + xs +
  labs(x = "Generation", y = "Lineages") + th +
  theme(strip.text = element_blank(), panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_text(colour = ink, size = 10), legend.text = element_text(colour = ink, size = 10),
        legend.key.size = unit(0.4, "cm"))

p <- top / bot + plot_layout(heights = c(6.5, 1))
ggsave("media/fig-scenarios-combined-v3.png", p, width = 9, height = 7.6, dpi = 150, bg = "white")
cat("Saved media/fig-scenarios-combined-v3.png\n")
