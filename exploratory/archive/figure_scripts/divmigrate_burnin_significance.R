# divmigrate_burnin_significance.R
# [dropped thread] rep-39 burn-in divMigrate significance profile (OUR BH z-test) -> media/fig-divmigrate-burnin-significance.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
load("~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate/firstsig/divMig.burnin_profile.39.rda")
x <- burnin_profile |> arrange(generation)
long <- x |> transmute(generation, `Significant pairs (of 300)` = n_sig, `G[ST]` = Gst, `Top pair: chain distance` = top_dist) |>
  pivot_longer(-generation, names_to = "panel", values_to = "value") |>
  mutate(panel = factor(panel, c("Significant pairs (of 300)", "G[ST]", "Top pair: chain distance")))
quiet <- x |> filter(n_sig == 0)
cat("quiet generations:", quiet$generation, "\n")
p <- ggplot(long, aes(generation, value)) +
  annotate("rect", xmin = -Inf, xmax = 16.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_line(colour = "steelblue", linewidth = 0.5, na.rm = TRUE) +
  geom_point(data = filter(long, panel == "Top pair: chain distance"), colour = "steelblue", size = 0.7, na.rm = TRUE) +
  geom_point(data = filter(long, panel == "Significant pairs (of 300)", value == 0), colour = "black", size = 1.6) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 5), breaks = c(0, 10, 25, 50, 100, 250, 500, 1000, 2000)) +
  facet_grid(panel ~ ., scales = "free_y", switch = "y", labeller = labeller(panel = c(`Significant pairs (of 300)` = "Significant pairs (of 300)",
                                                                                       `G[ST]` = "Gₛₜ", `Top pair: chain distance` = "Top pair: chain distance"))) +
  labs(x = "Burn-in generation (symmetric migration; pseudo-log scale)", y = NULL,
       caption = "Replicate 39. divMigrate Nm bootstrap (1,000 resamples of individuals), Benjamini–Hochberg across pairs at q = 0.05.\nShaded: generations 0–14, the only censuses with no significant pair (black points).") +
  theme_minimal(base_size = 11) + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))
ggsave("media/fig-divmigrate-burnin-significance.png", p, width = 8, height = 7, dpi = 150)
