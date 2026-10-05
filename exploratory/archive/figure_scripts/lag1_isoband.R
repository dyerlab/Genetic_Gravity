# lag1_isoband.R
# Fig_DivMigrateLag1: lag-1 medians with isotropic mean +/- 1 SD bands (needs lag1_all.rds from lag1_plot.R) -> media/fig-divmig-lag1-isoband.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.
S_DIR <- "data/derived"; dir.create(S_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
d <- readRDS(file.path(S_DIR, "lag1_all.rds"))
med <- d |> group_by(Treatment, Matrix, Generation) |>
  summarise(Lag1 = median(Lag1, na.rm = TRUE), .groups = "drop")
# Isotropic across-replicate mean +/- 1 SD at each generation, underlaid in the
# two asymmetric panels.
iso <- d |> filter(Treatment == "Isotropic") |> group_by(Matrix, Generation) |>
  summarise(mean = mean(Lag1, na.rm = TRUE), sd = sd(Lag1, na.rm = TRUE), .groups = "drop")
band <- expand_grid(Treatment = factor(c("Redistributed", "Obstructed"), levels = levels(d$Treatment)), iso)
cat("isotropic SD, median by matrix:\n"); print(tapply(iso$sd, iso$Matrix, median))
cols  <- c("divMigrate Nm" = "steelblue", "Graph similarity S" = "firebrick")
fills <- c("divMigrate Nm" = "#cfe0ef", "Graph similarity S" = "#f3d0d0")
p <- ggplot(med, aes(Generation, Lag1)) +
  geom_ribbon(data = band, aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = Matrix),
              alpha = 0.7, colour = NA) +
  geom_line(aes(colour = Matrix), linewidth = 0.7) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = fills, guide = "none") +
  facet_grid(Treatment ~ .) +
  labs(y = expression("Median lag-1 Spearman " * rho * " (t, t + 5 generations)"), colour = NULL,
       caption = "Shaded: isotropic across-replicate mean ± 1 SD") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-divmig-lag1-isoband.png", p, width = 9, height = 7, dpi = 150)
