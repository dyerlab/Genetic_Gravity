# truth_sim_overall.R
# Overall rho + RMSE, pGD similarity vs Nm -> media/fig-divmig-truth-sim-overall.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
d <- bind_rows(lapply(list.files(dm, "^divMig[.][0-9]+[.]truth_sim[.]rda$", full.names = TRUE),
                      function(f) { e <- new.env(); load(f, e); e$divmig_truth })) |>
  filter(!Excluded) |>
  mutate(Treatment = factor(as.character(Treatment), c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Method = factor(as.character(Method), c("pGD", "Nm"), c("pGD similarity 1/(1+p_ij)", "divMigrate Nm")))
med <- d |> pivot_longer(c(rho_all, rmse_all), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, c("rho_all", "rmse_all"), c("Overall Spearman rho", "Overall RMSE"))) |>
  group_by(Treatment, Method, metric, Generation) |> summarise(value = median(value), .groups = "drop")
p <- ggplot(med, aes(Generation, value, colour = Method)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("pGD similarity 1/(1+p_ij)" = "firebrick", "divMigrate Nm" = "steelblue")) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates vs. true matrix", colour = NULL,
       caption = "Censuses with any non-finite or negative Nm excluded for both methods") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-divmig-truth-sim-overall.png", p, width = 10, height = 5.5, dpi = 150)
