# truth_sim_vs_raw.R
# Truth comparison with pGD similarity; raw vs similarity vs Nm -> media/fig-divmig-truth-sim.png, fig-divmig-truth-raw-vs-sim.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
rd <- function(pat) bind_rows(lapply(list.files(dm, pat, full.names = TRUE), function(f) { e <- new.env(); load(f, e); e$divmig_truth }))
raw <- rd("^divMig[.][0-9]+[.]truth[.]rda$"); sim <- rd("^divMig[.][0-9]+[.]truth_sim[.]rda$")
stopifnot(identical(raw$Excluded, sim$Excluded))
cat("Nm rows identical across runs:", isTRUE(all.equal(raw[raw$Method == "Nm", ], sim[sim$Method == "Nm", ], check.attributes = FALSE)), "\n")
d <- bind_rows(raw |> mutate(Method = ifelse(Method == "pGD", "pGD raw p_ij", "divMigrate Nm")),
               sim |> filter(Method == "pGD") |> mutate(Method = "pGD similarity 1/(1+p_ij)")) |>
  mutate(Treatment = factor(as.character(Treatment), c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Method = factor(Method, c("pGD similarity 1/(1+p_ij)", "pGD raw p_ij", "divMigrate Nm")))
mets <- c(rho_all = "Overall rho", rmse_all = "Overall RMSE", rho_corr = "Corridor rho",
          rmse_corr = "Corridor RMSE", auc_corr = "Corridor AUC")
med <- d |> filter(!Excluded) |> pivot_longer(names(mets), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, names(mets), mets)) |>
  group_by(Treatment, Method, metric, Generation) |> summarise(value = median(value, na.rm = TRUE), .groups = "drop")
cat("\n== Median over censuses ==\n")
print(med |> group_by(metric, Method, Treatment) |> summarise(v = round(median(value, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = Treatment, values_from = v) |> arrange(metric, Method) |> as.data.frame())
cat("\n== Early (2004-2249) vs late (2754-2999), asymmetric scenarios ==\n")
print(med |> filter(Treatment != "Isotropic") |> mutate(w = ifelse(Generation < 2250, "early", ifelse(Generation >= 2754, "late", NA))) |>
  filter(!is.na(w)) |> group_by(metric, Treatment, Method, w) |> summarise(v = round(median(value, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = w, values_from = v) |> arrange(metric, Treatment, Method) |> as.data.frame())
theme_set(theme_minimal(base_size = 10) + theme(legend.position = "top"))
cap <- "Censuses with any non-finite or negative Nm excluded for all methods"
p1 <- ggplot(filter(med, Method != "pGD raw p_ij"), aes(Generation, value, colour = Method)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("pGD similarity 1/(1+p_ij)" = "firebrick", "divMigrate Nm" = "steelblue")) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates (vs. rescaled true migration matrix)", colour = NULL, caption = cap)
ggsave("media/fig-divmig-truth-sim.png", p1, width = 10, height = 10, dpi = 150)
p2 <- ggplot(med, aes(Generation, value, colour = Method, linetype = Method)) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = c("pGD similarity 1/(1+p_ij)" = "firebrick", "pGD raw p_ij" = "firebrick", "divMigrate Nm" = "steelblue")) +
  scale_linetype_manual(values = c("pGD similarity 1/(1+p_ij)" = "solid", "pGD raw p_ij" = "22", "divMigrate Nm" = "solid")) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates (vs. rescaled true migration matrix)", colour = NULL, linetype = NULL, caption = cap)
ggsave("media/fig-divmig-truth-raw-vs-sim.png", p2, width = 10, height = 10, dpi = 150)
