# scenarios_combined.R
# Fig_ScenarioComparison: 3x2 facet of C_D-bar, diameter, mean |Delta| for Redistributed | Obstructed -> media/fig-scenarios-combined.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(mgcv)
})

fwd     <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
iso_par <- read.csv("data/derived/isotropic_fit_params.csv", stringsAsFactors = FALSE)

GENS <- 2000:2999
stat_levels <- c("CD", "diameter", "mean_abs_delta")
stat_labels <- c("bar(C)[D]", "Diameter", "bar('|' * Delta * '|')")
scn_levels  <- c("Redistributed", "Obstructed")

fit_trend <- function(df) {
  df$replicate <- factor(df$replicate)
  m <- gam(value ~ s(generation, k = 20) + s(replicate, bs = "re"),
           data = df, method = "REML")
  nd <- data.frame(generation = GENS, replicate = df$replicate[1])
  data.frame(generation = GENS,
             fit = as.numeric(predict(m, nd, exclude = "s(replicate)")))
}

long <- fwd |>
  filter(Scenario %in% scn_levels) |>
  pivot_longer(all_of(stat_levels), names_to = "stat", values_to = "value")

trend <- long |> group_by(Scenario, stat) |>
  group_modify(~ fit_trend(.x)) |> ungroup()
bands <- long |> group_by(Scenario, stat, generation) |>
  summarise(mean = mean(value), sd = sd(value), .groups = "drop")
iso_curve <- expand_grid(Scenario = scn_levels, stat = stat_levels, generation = GENS) |>
  left_join(iso_par, by = "stat") |>
  mutate(fit = yinf + (y0 - yinf) * exp(-generation / tau))

lab <- function(d) mutate(d,
  Statistic = factor(stat, stat_levels, stat_labels),
  Scenario  = factor(Scenario, scn_levels))

p <- ggplot(lab(bands), aes(generation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = mean), colour = "grey20", alpha = 0.45, linewidth = 0.4) +
  geom_line(data = lab(iso_curve), aes(y = fit), colour = "grey30",
            linetype = "22", linewidth = 0.6) +
  geom_line(data = lab(trend), aes(y = fit), colour = "firebrick", linewidth = 0.7) +
  facet_grid(Statistic ~ Scenario, scales = "free_y",
             labeller = labeller(Statistic = label_parsed, Scenario = label_value)) +
  labs(x = "Generation", y = NULL) +
  theme_minimal(base_size = 11)

ggsave("media/fig-scenarios-combined.png", p, width = 9, height = 6.5, dpi = 150)
cat("Saved media/fig-scenarios-combined.png\n")
