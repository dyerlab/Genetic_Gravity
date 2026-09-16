# Boundary decomposition: Tab_EdgeClassDecomposition and Fig_BoundaryDiversity.
#
# Reproduces the manuscript's edge-class (leaf-incident vs. interior vs. all)
# decomposition of the asymmetry index, and the terminus-vs-interior expected
# heterozygosity figure. Source data (data/boundary_summary_edges.csv,
# boundary_summary_diversity.csv) are the per-snapshot summaries already
# computed from the individual-based simulation output by the private
# research repo's R/boundary_timecourse.R; no raw simulation data is required
# to run this script.
#
# Run from the repository root:
#   Rscript R/boundary_decomposition.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

edges <- read.csv("data/boundary_summary_edges.csv", stringsAsFactors = FALSE)
div   <- read.csv("data/boundary_summary_diversity.csv", stringsAsFactors = FALSE)

scn <- c("Isotropic", "Flux-conserved", "Rate-conserved")
edges <- edges |> mutate(scenario = factor(scenario, levels = scn))
div <- div |>
  mutate(scenario = factor(scenario, levels = scn),
         position = ifelse(is_terminus, "Terminus", "Interior"))

# Established-signal window: generation 2304 +/- 50, i.e. all sampled
# generations in [2254, 2354]. Snapshots are cached every 50 generations
# (2004-2999). A window is used rather than the single nearest generation
# because leaf-incident edges can drop out entirely at any one snapshot by
# sampling sparsity alone (a degree-1 chain terminus is not guaranteed to
# still be degree-1 in any specific estimated graph).
gen_lo <- 2304 - 50
gen_hi <- 2304 + 50

edge_at <- function(scenario_name, edge_class, col) {
  d <- edges |> filter(scenario == scenario_name, edge_class == !!edge_class,
                        generation >= gen_lo, generation <= gen_hi)
  if (nrow(d) == 0) NA_real_ else round(mean(d[[col]], na.rm = TRUE), 3)
}

## ---- Tab_EdgeClassDecomposition ----------------------------------------

tbl_boundary <- edges |>
  filter(generation >= gen_lo, generation <= gen_hi) |>
  group_by(Scenario = scenario, `Edge class` = edge_class) |>
  summarise(`Mean signed Delta` = round(mean(mean_delta), 3),
            `Mean |Delta|`      = round(mean(mean_abs), 3),
            `Edges/graph`       = round(mean(n), 1), .groups = "drop") |>
  arrange(Scenario, factor(`Edge class`, levels = c("interior", "leaf-incident", "all")))

cat("=== Tab_EdgeClassDecomposition ===\n")
print(as.data.frame(tbl_boundary), row.names = FALSE)

rv_bd_int_flux <- edge_at("Flux-conserved", "interior", "mean_delta")
rv_bd_int_rate <- edge_at("Rate-conserved", "interior", "mean_delta")
rv_bd_int_iso  <- edge_at("Isotropic",      "interior", "mean_delta")
rv_bd_leaf_abs <- round(mean(vapply(scn, edge_at, numeric(1),
                                     edge_class = "leaf-incident", col = "mean_abs"),
                              na.rm = TRUE), 3)
rv_bd_int_abs  <- round(mean(vapply(scn, edge_at, numeric(1),
                                     edge_class = "interior", col = "mean_abs"),
                              na.rm = TRUE), 3)

# NOTE: Isotropic has zero leaf-incident edges anywhere in the archived data
# (all 21 sampled generations, all replicates) -- not a windowing artifact.
# rv_bd_leaf_abs above is therefore the mean of Flux-conserved and
# Rate-conserved only; see manuscript text for how this is described.
cat(sprintf("\nrv_bd_int_flux = %s\nrv_bd_int_rate = %s\nrv_bd_int_iso  = %s\n",
            rv_bd_int_flux, rv_bd_int_rate, rv_bd_int_iso))
cat(sprintf("rv_bd_leaf_abs = %s (Flux-conserved + Rate-conserved only; Isotropic has no leaf-incident edges)\n",
            rv_bd_leaf_abs))
cat(sprintf("rv_bd_int_abs  = %s\n", rv_bd_int_abs))

## ---- Fig_BoundaryDiversity ----------------------------------------------

last_gen <- max(div$generation, na.rm = TRUE)
he_at <- function(scenario_name, position_name, gen) {
  d <- div |> filter(scenario == scenario_name, position == position_name, generation == gen)
  if (nrow(d) == 0) NA_real_ else round(mean(d$He, na.rm = TRUE), 3)
}
rv_bd_he_term <- he_at("Rate-conserved", "Terminus", last_gen)
rv_bd_he_int  <- he_at("Rate-conserved", "Interior", last_gen)
cat(sprintf("\nrv_bd_he_term (Rate-conserved, final generation) = %s\n", rv_bd_he_term))
cat(sprintf("rv_bd_he_int  (Rate-conserved, final generation) = %s\n", rv_bd_he_int))

fig_boundary_diversity <- div |>
  group_by(scenario, position, generation) |>
  summarise(He = mean(He, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(generation, He, colour = position)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scenario) +
  scale_colour_manual(values = c("Terminus" = "#d73027", "Interior" = "#4575b4")) +
  labs(x = "Generation", y = "Mean expected heterozygosity", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

if (interactive()) print(fig_boundary_diversity)
ggsave("data/derived/fig-boundary-diversity_reproduced.png", fig_boundary_diversity,
       width = 9, height = 3.6, dpi = 150)
cat("\nSaved reproduced Fig_BoundaryDiversity to data/derived/fig-boundary-diversity_reproduced.png\n")
