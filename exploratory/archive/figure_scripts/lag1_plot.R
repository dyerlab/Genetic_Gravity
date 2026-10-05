# lag1_plot.R
# Lag-1 median traces, Nm vs S (step 2; also caches lag1_all.rds used by lag1_isoband.R) -> media/fig-divmig-lag1-all.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.
S_DIR <- "data/derived"; dir.create(S_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
dm_dir <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
fs <- list.files(dm_dir, "^divMig[.][0-9]+[.]lag1[.]rda$", full.names = TRUE)
d <- bind_rows(lapply(fs, function(f) { e <- new.env(); load(f, e); e$divmig_lag1 })) |>
  mutate(Treatment = factor(as.character(Treatment), levels = c("iso", "flux", "rate"),
                            labels = c("Isotropic", "Redistributed", "Obstructed")),
         Matrix = factor(as.character(Matrix), levels = c("Nm", "S"),
                         labels = c("divMigrate Nm", "Graph similarity S")))
cat("replicates:", n_distinct(d$Replicate), " rows:", nrow(d), " NA:", sum(is.na(d$Lag1)), "\n")
med <- d |> group_by(Treatment, Matrix, Generation) |>
  summarise(Lag1 = median(Lag1, na.rm = TRUE), .groups = "drop")
p <- ggplot(med, aes(Generation, Lag1, colour = Matrix)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("divMigrate Nm" = "steelblue", "Graph similarity S" = "firebrick")) +
  facet_grid(Treatment ~ .) +
  labs(y = expression("Median lag-1 Spearman " * rho * " (t, t + 5 generations)"), colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-divmig-lag1-all.png", p, width = 9, height = 7, dpi = 150)
saveRDS(d, file.path(S_DIR, "lag1_all.rds"))
