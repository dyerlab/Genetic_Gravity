# divmigrate_significance_timeline.R
# [dropped thread] rep-39 full-timeline divMigrate significance (OUR BH z-test) -> media/fig-divmigrate-significance-timeline.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dir <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate/firstsig"
load(file.path(dir, "divMig.burnin_profile.39.rda")); load(file.path(dir, "divMig.forward_profile.39.rda"))
lab <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
x <- bind_rows(burnin_profile |> mutate(Series = "Burn-in (symmetric)"),
               forward_profile |> mutate(Series = unname(lab[as.character(Treatment)])) |> select(-Treatment)) |>
  mutate(Series = factor(Series, c("Burn-in (symmetric)", "Isotropic", "Redistributed", "Obstructed")))
# join each treatment to the last burn-in census so its line starts from generation 1999
start <- burnin_profile |> filter(generation == max(generation))
x <- bind_rows(x, bind_rows(lapply(lab, function(l) start |> mutate(Series = factor(l, levels(x$Series))))))
cat("== Forward phase (2004-2999), replicate 39 ==\n")
print(x |> filter(generation >= 2004) |> group_by(Series) |>
  summarise(censuses = n(), untestable = sum(is.na(n_sig)), zero_sig = sum(n_sig == 0, na.rm = TRUE),
            median_sig = median(n_sig, na.rm = TRUE), min_sig = min(n_sig, na.rm = TRUE),
            Gst_end = round(Gst[which.max(generation)], 3), median_top_dist = median(top_dist, na.rm = TRUE), .groups = "drop") |> as.data.frame())
long <- x |> transmute(generation, Series, `Significant pairs (of 300)` = n_sig, Gst = Gst, `Top pair: chain distance` = top_dist) |>
  pivot_longer(-c(generation, Series), names_to = "panel", values_to = "value") |>
  mutate(panel = factor(panel, c("Significant pairs (of 300)", "Gst", "Top pair: chain distance"),
                        c("Significant pairs (of 300)", "Gₛₜ", "Top pair: chain distance"))) |>
  arrange(Series, generation)
cols <- c("Burn-in (symmetric)" = "black", Isotropic = "grey55", Redistributed = "firebrick", Obstructed = "steelblue")
p <- ggplot(long, aes(generation, value, colour = Series)) +
  annotate("rect", xmin = -Inf, xmax = 16.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = 2000, linetype = "22", colour = "grey40") +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  scale_colour_manual(values = cols) +
  facet_grid(panel ~ ., scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = seq(0, 3000, 500)) +
  labs(x = "Generation (treatments applied at 2000, dashed line)", y = NULL, colour = NULL,
       caption = "Replicate 39. divMigrate Nm bootstrap (1,000 resamples of individuals), Benjamini–Hochberg across pairs at q = 0.05.\nShaded: generations 0–14, the only censuses with no significant pair. Gaps: censuses where the test could not be run.") +
  theme_minimal(base_size = 11) + theme(legend.position = "top", strip.placement = "outside", strip.text.y.left = element_text(angle = 90))
ggsave("media/fig-divmigrate-significance-timeline.png", p, width = 10, height = 7.5, dpi = 150)
