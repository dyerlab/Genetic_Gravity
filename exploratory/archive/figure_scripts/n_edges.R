# n_edges.R
# Mean edge count through time by scenario -> media/fig-n-edges.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
f <- read.csv("data/forward_scenarios.csv")
d <- f |> group_by(Scenario, generation) |>
  summarise(mean = mean(n_edges), sd = sd(n_edges), n = n(), .groups = "drop") |>
  mutate(Scenario = factor(Scenario, levels = c("Isotropic", "Redistributed", "Obstructed")))
cat("replicates per point:", range(d$n), "\n")
print(d |> filter(generation %in% c(2004, 2499, 2999)) |> select(Scenario, generation, mean, sd) |> as.data.frame(), digits = 3)
cols <- c(Isotropic = "grey40", Redistributed = "firebrick", Obstructed = "steelblue")
p <- ggplot(d, aes(generation, mean, colour = Scenario, fill = Scenario)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(x = "Generation", y = "Edges in Population Graph (mean ± 1 SD)", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-n-edges.png", p, width = 8, height = 4.5, dpi = 150)
