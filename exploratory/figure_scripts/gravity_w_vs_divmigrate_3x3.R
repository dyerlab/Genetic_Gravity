# gravity_w_vs_divmigrate_3x3.R
# Manuscript 3x3: Spearman rho, RMSE, replicates present; w_ij vs divMigrate Nm -> media/fig-gravity-w-vs-divmigrate-3x3.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
rd <- function(pat, obj) bind_rows(lapply(list.files(dm, pat, full.names = TRUE), function(f) { e <- new.env(); load(f, e); e[[obj]] }))
w  <- rd("^divMig[.][0-9]+[.]truth_w[.]rda$", "gravity_truth_w") |> filter(Score == "w_raw")
nm <- rd("^divMig[.][0-9]+[.]truth_sim[.]rda$", "divmig_truth") |> filter(Method == "Nm")
key <- c("Replicate", "Treatment", "Generation"); mets <- c("rho_all", "rmse_all")
ex <- nm |> select(all_of(key), Excluded) |> mutate(Treatment = as.character(Treatment))
w  <- w |> mutate(Treatment = as.character(Treatment)) |> left_join(ex, by = key)
d <- bind_rows(
  w  |> filter(!Excluded) |> transmute(across(all_of(key)), across(all_of(mets)), Series = "Genetic gravity w_ij"),
  nm |> filter(!Excluded) |> mutate(Treatment = as.character(Treatment)) |> transmute(across(all_of(key)), across(all_of(mets)), Series = "divMigrate Nm"),
  w  |> transmute(across(all_of(key)), across(all_of(mets)), Series = "Genetic gravity w_ij, all censuses")) |>
  mutate(Treatment = factor(Treatment, c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Series = factor(Series, c("Genetic gravity w_ij", "divMigrate Nm", "Genetic gravity w_ij, all censuses")))
lab <- c(rho_all = "Spearman\u2019s \u03c1", rmse_all = "RMSE")
med <- d |> pivot_longer(all_of(mets), names_to = "metric", values_to = "value") |> mutate(metric = factor(metric, names(lab), lab)) |>
  group_by(Treatment, Series, metric, Generation) |> summarise(value = median(value), .groups = "drop")
present <- d |> filter(Series != "Genetic gravity w_ij") |>
  group_by(Treatment, Series, Generation) |> summarise(value = 100 * n_distinct(Replicate) / 50, .groups = "drop") |>
  mutate(metric = factor("Replicates (%)", levels = c(lab, "Replicates (%)")))
med <- bind_rows(med |> mutate(metric = factor(as.character(metric), levels = levels(present$metric))), present)
cat("== Percent present, mean by window ==\n")
print(present |> mutate(w = cut(Generation, c(2000, 2250, 2500, 2750, 3000))) |> group_by(Treatment, Series, w) |> summarise(pct = round(mean(value), 1), .groups = "drop") |> pivot_wider(names_from = w, values_from = pct) |> as.data.frame())
cat("== Median over censuses (paired censuses; last column = w on all censuses) ==\n")
print(med |> group_by(metric, Treatment, Series) |> summarise(v = round(median(value), 3), .groups = "drop") |> pivot_wider(names_from = Series, values_from = v) |> as.data.frame())
cat("\n== Paired per-census comparison, w_ij vs Nm ==\n")
pw <- inner_join(w |> filter(!Excluded) |> select(all_of(key), all_of(mets)), nm |> filter(!Excluded) |> mutate(Treatment = as.character(Treatment)) |> select(all_of(key), all_of(mets)), by = key, suffix = c("_w", "_nm"))
print(pw |> group_by(Treatment) |> summarise(censuses = n(), w_higher_rho = round(mean(rho_all_w > rho_all_nm), 3), w_lower_rmse = round(mean(rmse_all_w < rmse_all_nm), 3)) |> as.data.frame())
p <- ggplot(med, aes(Generation, value, colour = Series, linetype = Series, linewidth = Series)) + geom_line() +
  scale_colour_manual(values = c("Genetic gravity w_ij" = "firebrick", "divMigrate Nm" = "steelblue", "Genetic gravity w_ij, all censuses" = "firebrick")) +
  scale_linetype_manual(values = c("Genetic gravity w_ij" = "solid", "divMigrate Nm" = "solid", "Genetic gravity w_ij, all censuses" = "22")) +
  scale_linewidth_manual(values = c("Genetic gravity w_ij" = 0.65, "divMigrate Nm" = 0.65, "Genetic gravity w_ij, all censuses" = 0.35)) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = NULL, colour = NULL, linetype = NULL, linewidth = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("media/fig-gravity-w-vs-divmigrate-3x3.png", p, width = 10, height = 7.5, dpi = 150)
