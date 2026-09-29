# gravity_w_vs_divmigrate_3x3_v3.R
# Phase strips use the heterozygosity-gradient reorganization boundary (2026-09-29); earlier file = gain-based boundary.
# v3 = v2 plus the phase-share strip (exploratory/lineage_phases.R) -> media/fig-gravity-w-vs-divmigrate-3x3-v4.png
# Fig_DivMigrateTruth, transposed-truth version (local/node_potential_reanalysis.md, Task 3a):
# w_{j|i} (W[i, j], j's weight in i's neighbourhood) is compared with m_{j->i}, the
# transposed truth, because gene flow j -> i raises it; divMigrate's relative Nm
# (row = source) is compared with m_{i->j} as before. Same layout as
# gravity_w_vs_divmigrate_3x3.R (which drew the m_{i->j} comparison, kept as
# media/fig-gravity-w-vs-divmigrate-3x3.png).
# Input: data/derived/node_potential_reanalysis.rda (rean$truth$d, from
# exploratory/node_potential_reanalysis.R). Output: media/fig-gravity-w-vs-divmigrate-3x3-v2.png
# Run from the manuscript repo root.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
source("exploratory/lineage_phases.R")
e <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e)
t3 <- e$rean$truth$d
key <- c("Replicate", "Treatment", "Generation")
d <- bind_rows(
  t3 |> filter(Exclusion == "") |> transmute(across(all_of(key)), rho_all = w_rho_tr, rmse_all = w_rmse_tr,
                                            Series = "Genetic gravity w_ij"),
  t3 |> filter(Exclusion == "") |> transmute(across(all_of(key)), rho_all = nm_rho_run, rmse_all = nm_rmse_run,
                                            Series = "divMigrate Nm"),
  t3 |> transmute(across(all_of(key)), rho_all = w_rho_tr, rmse_all = w_rmse_tr,
                  Series = "Genetic gravity w_ij, all censuses")) |>
  mutate(Treatment = factor(Treatment, c("Isotropic", "Redistributed", "Obstructed")),
         Series = factor(Series, c("Genetic gravity w_ij", "divMigrate Nm", "Genetic gravity w_ij, all censuses")))
mets <- c("rho_all", "rmse_all")
lab <- c(rho_all = "Spearman’s ρ", rmse_all = "RMSE")
med <- d |> pivot_longer(all_of(mets), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, names(lab), lab)) |>
  group_by(Treatment, Series, metric, Generation) |> summarise(value = median(value), .groups = "drop")
present <- d |> filter(Series != "Genetic gravity w_ij") |>
  group_by(Treatment, Series, Generation) |> summarise(value = 100 * n_distinct(Replicate) / 50, .groups = "drop") |>
  mutate(metric = factor("Replicates (%)", levels = c(lab, "Replicates (%)")))
med <- bind_rows(med |> mutate(metric = factor(as.character(metric), levels = levels(present$metric))), present)
p <- ggplot(med, aes(Generation, value, colour = Series, linetype = Series, linewidth = Series)) + geom_line() +
  scale_colour_manual(values = c("Genetic gravity w_ij" = "firebrick", "divMigrate Nm" = "steelblue",
                                 "Genetic gravity w_ij, all censuses" = "firebrick")) +
  scale_linetype_manual(values = c("Genetic gravity w_ij" = "solid", "divMigrate Nm" = "solid",
                                   "Genetic gravity w_ij, all censuses" = "22")) +
  scale_linewidth_manual(values = c("Genetic gravity w_ij" = 0.65, "divMigrate Nm" = 0.65,
                                    "Genetic gravity w_ij, all censuses" = 0.35)) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = NULL, colour = NULL, linetype = NULL, linewidth = NULL) +
  scale_x_continuous(limits = c(2000, 3000), breaks = seq(2000, 2750, 250), expand = c(0, 0)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top", axis.text.x = element_blank(),
                                        panel.spacing.x = unit(1, "lines")) + labs(x = NULL)
xs <- scale_x_continuous(limits = c(2000, 3000), breaks = seq(2000, 2750, 250), expand = c(0, 0))
strip <- lp_strip(c("Isotropic", "Redistributed", "Obstructed"), xs, facet = "col") +
  labs(x = "Generation") + theme(panel.spacing.x = unit(1, "lines"))
p <- p / strip + plot_layout(heights = c(7.5, 0.9))
ggsave("media/fig-gravity-w-vs-divmigrate-3x3-v4.png", p, width = 10, height = 8.6, dpi = 150, bg = "white")
cat("Saved media/fig-gravity-w-vs-divmigrate-3x3-v4.png\n")
