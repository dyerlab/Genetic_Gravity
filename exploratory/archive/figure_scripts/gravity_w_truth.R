# gravity_w_truth.R
# Gravity scores 1/(1+p), w_ij, 1/(1+w) vs truth -> media/fig-gravity-w-truth.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
d <- bind_rows(lapply(list.files(dm, "^divMig[.][0-9]+[.]truth_w[.]rda$", full.names = TRUE), function(f) { e <- new.env(); load(f, e); e$gravity_truth_w })) |>
  mutate(Treatment = factor(as.character(Treatment), c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Score = factor(as.character(Score), c("p_sim", "w_raw", "w_sim"), c("1/(1+p_ij)  (previous)", "w_ij  (raw)", "1/(1+w_ij)  (similarity)")))
mets <- c(rho_all = "Overall Spearman rho", rmse_all = "Overall RMSE", mae_all = "Overall MAE")
med <- d |> pivot_longer(names(mets), names_to = "metric", values_to = "value") |> mutate(metric = factor(metric, names(mets), mets)) |>
  group_by(Treatment, Score, metric, Generation) |> summarise(value = median(value), .groups = "drop")
cat("== Median over censuses of per-census medians (all censuses, 50 replicates) ==\n")
print(med |> group_by(metric, Score, Treatment) |> summarise(v = round(median(value), 3), .groups = "drop") |> pivot_wider(names_from = Treatment, values_from = v) |> as.data.frame())
cat("\n== Corridor metrics ==\n")
print(d |> group_by(Score, Treatment) |> summarise(rho_corr = round(median(rho_corr, na.rm = TRUE), 3), auc_corr = round(median(auc_corr), 3), .groups = "drop") |> as.data.frame())
cat("\n== Paired: w_raw rho minus previous p_sim rho ==\n")
print(d |> select(Replicate, Treatment, Generation, Score, rho_all) |> pivot_wider(names_from = Score, values_from = rho_all) |>
  mutate(gap = `w_ij  (raw)` - `1/(1+p_ij)  (previous)`) |> group_by(Treatment) |> summarise(w_raw_higher = round(mean(gap > 0), 3), median_gap = round(median(gap), 3)) |> as.data.frame())
p <- ggplot(med, aes(Generation, value, colour = Score, linetype = Score)) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("1/(1+p_ij)  (previous)" = "grey45", "w_ij  (raw)" = "firebrick", "1/(1+w_ij)  (similarity)" = "darkorange")) +
  scale_linetype_manual(values = c("1/(1+p_ij)  (previous)" = "22", "w_ij  (raw)" = "solid", "1/(1+w_ij)  (similarity)" = "solid")) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates vs. true matrix", colour = NULL, linetype = NULL,
       caption = "Gravity scores only; 0 for pairs without a graph edge; all censuses; each matrix rescaled by its maximum") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-gravity-w-truth.png", p, width = 10, height = 7.5, dpi = 150)
