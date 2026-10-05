# divmig_correlation_all.R
# Fig_DivMigrateCorrelation: per-census Spearman rho(S, Nm), all replicates + median -> media/fig-divmig-correlation-all.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.
S_DIR <- "data/derived"; dir.create(S_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
fs <- list.files("~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate",
                 "^divMig[.][0-9]+[.]correlation[.]rda$", full.names = TRUE)
d <- bind_rows(lapply(fs, function(f) { e <- new.env(); load(f, e); e$divmig_cor })) |>
  mutate(Treatment = factor(as.character(Treatment), levels = c("iso", "flux", "rate"),
                            labels = c("Isotropic", "Redistributed", "Obstructed")))
cat("rows:", nrow(d), " replicates:", n_distinct(d$Replicate), "\n"); print(table(d$Treatment))
cat("NA correlations:", sum(is.na(d$Correlation)), "\n")
med <- d |> group_by(Treatment, Generation) |>
  summarise(Correlation = median(Correlation, na.rm = TRUE), .groups = "drop")
p <- ggplot(d, aes(Generation, Correlation)) +
  geom_line(aes(group = Replicate), colour = "grey50", alpha = 0.2, linewidth = 0.3) +
  geom_line(data = med, colour = "firebrick", linewidth = 0.8) +
  facet_grid(Treatment ~ .) +
  labs(y = expression("Spearman " * rho * " (S, divMigrate Nm)")) +
  theme_minimal(base_size = 11)
ggsave("media/fig-divmig-correlation-all.png", p, width = 8, height = 7, dpi = 150)
saveRDS(d, file.path(S_DIR, "divmig_cor_all.rds"))
