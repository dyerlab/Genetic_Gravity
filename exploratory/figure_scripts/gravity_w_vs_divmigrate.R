# gravity_w_vs_divmigrate.R
# w_ij vs Nm, rho/RMSE/MAE -> media/fig-gravity-w-vs-divmigrate.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
rd <- function(pat, obj) bind_rows(lapply(list.files(dm, pat, full.names = TRUE), function(f) { e <- new.env(); load(f, e); e[[obj]] }))
w  <- rd("^divMig[.][0-9]+[.]truth_w[.]rda$", "gravity_truth_w") |> filter(Score == "w_raw")
nm <- rd("^divMig[.][0-9]+[.]truth_sim[.]rda$", "divmig_truth") |> filter(Method == "Nm")
key <- c("Replicate", "Treatment", "Generation"); mets <- c("rho_all", "rmse_all", "mae_all")
ex <- nm |> select(all_of(key), Excluded) |> mutate(Treatment = as.character(Treatment))
w  <- w |> mutate(Treatment = as.character(Treatment)) |> left_join(ex, by = key)
d <- bind_rows(
  w  |> filter(!Excluded) |> transmute(across(all_of(key)), across(all_of(mets)), Series = "Genetic gravity w_ij"),
  nm |> filter(!Excluded) |> mutate(Treatment = as.character(Treatment)) |> transmute(across(all_of(key)), across(all_of(mets)), Series = "divMigrate Nm"),
  w  |> transmute(across(all_of(key)), across(all_of(mets)), Series = "Genetic gravity w_ij, all censuses")) |>
  mutate(Treatment = factor(Treatment, c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Series = factor(Series, c("Genetic gravity w_ij", "divMigrate Nm", "Genetic gravity w_ij, all censuses")))
lab <- c(rho_all = "Overall Spearman rho", rmse_all = "Overall RMSE", mae_all = "Overall MAE")
med <- d |> pivot_longer(all_of(mets), names_to = "metric", values_to = "value") |> mutate(metric = factor(metric, names(lab), lab)) |>
  group_by(Treatment, Series, metric, Generation) |> summarise(value = median(value), .groups = "drop")
cat("== Median over censuses (paired censuses; last column = w on all censuses) ==\n")
print(med |> group_by(metric, Treatment, Series) |> summarise(v = round(median(value), 3), .groups = "drop") |> pivot_wider(names_from = Series, values_from = v) |> as.data.frame())
cat("\n== Paired per-census comparison, w_ij vs Nm ==\n")
pw <- inner_join(w |> filter(!Excluded) |> select(all_of(key), all_of(mets)), nm |> filter(!Excluded) |> mutate(Treatment = as.character(Treatment)) |> select(all_of(key), all_of(mets)), by = key, suffix = c("_w", "_nm"))
print(pw |> group_by(Treatment) |> summarise(censuses = n(), w_higher_rho = round(mean(rho_all_w > rho_all_nm), 3), w_lower_rmse = round(mean(rmse_all_w < rmse_all_nm), 3), w_lower_mae = round(mean(mae_all_w < mae_all_nm), 3)) |> as.data.frame())
p <- ggplot(med, aes(Generation, value, colour = Series, linetype = Series, linewidth = Series)) + geom_line() +
  scale_colour_manual(values = c("Genetic gravity w_ij" = "firebrick", "divMigrate Nm" = "steelblue", "Genetic gravity w_ij, all censuses" = "firebrick")) +
  scale_linetype_manual(values = c("Genetic gravity w_ij" = "solid", "divMigrate Nm" = "solid", "Genetic gravity w_ij, all censuses" = "22")) +
  scale_linewidth_manual(values = c("Genetic gravity w_ij" = 0.65, "divMigrate Nm" = 0.65, "Genetic gravity w_ij, all censuses" = 0.35)) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates vs. true matrix", colour = NULL, linetype = NULL, linewidth = NULL,
       caption = "Solid lines: censuses where Nm is intact (paired). Dashed: w_ij on all censuses. Each matrix rescaled by its maximum; diagonals excluded.") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-gravity-w-vs-divmigrate.png", p, width = 10, height = 7.5, dpi = 150)
