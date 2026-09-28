# truth_raw.R
# Truth comparison, raw p_ij vs Nm, five metrics (caches truth_all.rds) -> media/fig-divmig-truth.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.
S_DIR <- "data/derived"; dir.create(S_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
fs <- list.files(dm, "^divMig[.][0-9]+[.]truth[.]rda$", full.names = TRUE)
d <- bind_rows(lapply(fs, function(f) { e <- new.env(); load(f, e); e$divmig_truth })) |>
  mutate(Treatment = factor(as.character(Treatment), c("iso", "flux", "rate"), c("Isotropic", "Redistributed", "Obstructed")),
         Method = factor(as.character(Method), c("pGD", "Nm"), c("pGD  p_ij", "divMigrate Nm")))
cat("replicates:", n_distinct(d$Replicate), "\n\n== Censuses excluded for non-finite Nm ==\n")
ex <- d |> filter(Method == "divMigrate Nm") |> group_by(Treatment) |>
  summarise(censuses = n(), nonfinite = sum(Exclusion == "nonfinite"), negative = sum(Exclusion == "negative"),
            excluded = sum(Excluded), frac = round(mean(Excluded), 3),
            reps_affected = n_distinct(Replicate[Excluded]))
print(as.data.frame(ex))
mets <- c(rho_all = "Overall rho", rmse_all = "Overall RMSE", rho_corr = "Corridor rho",
          rmse_corr = "Corridor RMSE", auc_corr = "Corridor AUC")
long <- d |> filter(!Excluded) |> pivot_longer(names(mets), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, names(mets), mets))
med <- long |> group_by(Treatment, Method, metric, Generation) |>
  summarise(value = median(value, na.rm = TRUE), n = sum(!is.na(value)), .groups = "drop")
saveRDS(list(d = d, med = med, ex = ex), file.path(S_DIR, "truth_all.rds"))
cat("\n== Median over all censuses (of per-census medians) ==\n")
print(med |> group_by(metric, Treatment, Method) |> summarise(v = round(median(value, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = c(Treatment), values_from = v) |> as.data.frame())
cat("\n== Early vs late (2004-2249 vs 2754-2999) ==\n")
print(med |> mutate(w = ifelse(Generation < 2250, "early", ifelse(Generation >= 2754, "late", NA))) |> filter(!is.na(w)) |>
  group_by(metric, Treatment, Method, w) |> summarise(v = round(median(value, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = w, values_from = v) |> as.data.frame())
p <- ggplot(med, aes(Generation, value, colour = Method)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("pGD  p_ij" = "firebrick", "divMigrate Nm" = "steelblue")) +
  facet_grid(metric ~ Treatment, scales = "free_y") +
  labs(y = "Median across replicates (vs. rescaled true migration matrix)", colour = NULL,
       caption = "Censuses with any non-finite or negative Nm excluded for both methods") +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
ggsave("media/fig-divmig-truth.png", p, width = 10, height = 10, dpi = 150)
