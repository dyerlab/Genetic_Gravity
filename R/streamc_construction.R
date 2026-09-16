# Fig_ConstructionSurvival and Fig_ConstructionErosion.
#
# Reproduces the Stream C Arm B graph-construction "survival" analysis: the
# fraction of replicates yielding a constructible Population Graph at each
# forward census generation, by migration ratio and population size, and the
# erosion of the surviving graphs' mean edge count on the approach to
# construction failure. Source data (data/graph_summary.csv) is the
# per-replicate/scenario/generation graph summary already computed from the
# individual-based simulation output by the private research repo's
# R/graph_summarizer.R; no raw simulation data is required to run this
# script. construction_survival()/survival_summary() below are ported
# verbatim from the private repo's R/construction_survival.R -- the internal
# `cond = as.integer(as.character(scenario))` coercion works unchanged
# whether `scenario` is a factor (private .rda) or a plain integer column
# (this CSV), so no adaptation was needed beyond the top-level data load.
#
# Run from the repository root:
#   Rscript R/streamc_construction.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

graph_summary <- read.csv("data/graph_summary.csv", stringsAsFactors = FALSE)

## ---- ported from R/construction_survival.R (private repo), verbatim ------

# Arm B cells (incl. the 1:1 nulls) mapping condition code -> (N, ratio).
ARM_B_CELLS <- tibble::tribble(
  ~N,  ~ratio, ~cond,
   50, 1.0,  7,  50, 1.5,  8,  50, 2.0,  9,  50, 2.5, 10,  50, 4.0, 11,
  100, 1.0,  1, 100, 1.5, 12, 100, 2.0, 13, 100, 2.5, 14, 100, 4.0,  2,
  200, 1.0, 15, 200, 1.5, 16, 200, 2.0, 17, 200, 2.5, 18, 200, 4.0, 19,
  400, 1.0, 20, 400, 1.5, 21, 400, 2.0, 22, 400, 2.5, 23, 400, 4.0, 24
)

construction_survival <- function(graph_summary, cells = ARM_B_CELLS) {
  cells <- cells |> mutate(N = as.integer(N), cond = as.integer(cond))
  gs <- graph_summary
  if (!"N" %in% names(gs)) gs$N <- 100L
  fwd <- gs |>
    filter(phase == "forward") |>
    mutate(N = as.integer(N),
           cond = as.integer(as.character(scenario))) |>
    inner_join(cells, by = c("N", "cond"))

  intended <- fwd |>
    group_by(N, ratio, cond) |>
    summarise(n_intended = n_distinct(replicate), .groups = "drop")

  per_gen <- fwd |>
    group_by(N, ratio, cond, generation) |>
    summarise(n_built      = n_distinct(replicate),
              mean_n_edges = mean(n_edges, na.rm = TRUE),
              mean_cGD     = mean(mean_cGD, na.rm = TRUE),
              .groups = "drop")

  census_gens <- sort(unique(fwd$generation))
  expand_grid(distinct(cells, N, ratio, cond), generation = census_gens) |>
    left_join(intended, by = c("N", "ratio", "cond")) |>
    left_join(per_gen,  by = c("N", "ratio", "cond", "generation")) |>
    mutate(n_built = coalesce(n_built, 0L),
           buildability = ifelse(n_intended > 0, n_built / n_intended, NA_real_)) |>
    arrange(N, ratio, generation)
}

survival_summary <- function(surv, tol = 1e-9) {
  fin <- max(surv$generation, na.rm = TRUE)
  per_cell <- surv |>
    group_by(N, ratio, cond) |>
    summarise(
      n_intended    = dplyr::first(n_intended),
      last_full_gen = { gg <- generation[buildability >= 1 - tol]
                        if (length(gg)) max(gg) else NA_real_ },
      last_valid_gen = { gg <- generation[n_built > 0]
                         if (length(gg)) max(gg) else NA_real_ },
      n_built_2999  = { v <- n_built[generation == fin]; if (length(v)) v[1] else 0L },
      build_2999    = { v <- buildability[generation == fin]; if (length(v)) v[1] else 0 },
      .groups = "drop")

  bygen <- surv |>
    group_by(generation) |>
    summarise(min_build = min(buildability, na.rm = TRUE), .groups = "drop")
  gg <- bygen$generation[bygen$min_build >= 1 - tol]
  common_gen <- if (length(gg)) max(gg) else NA_real_

  list(per_cell = per_cell, common_pre_fixation_gen = common_gen)
}

## ---- run --------------------------------------------------------------

surv <- construction_survival(graph_summary)
css  <- survival_summary(surv)
pc   <- css$per_cell

.get <- function(n, r, col) { v <- pc[[col]][pc$N == n & pc$ratio == r]; if (length(v)) v[1] else NA }
cs_common_gen   <- css$common_pre_fixation_gen
cs_n50_4_built  <- .get(50, 4.0, "n_built_2999"); cs_n50_4_int  <- .get(50, 4.0, "n_intended")
cs_n50_4_lv     <- .get(50, 4.0, "last_valid_gen")
cs_n50_25_built <- .get(50, 2.5, "n_built_2999"); cs_n50_25_int <- .get(50, 2.5, "n_intended")
e50 <- surv[surv$N == 50 & surv$ratio == 4.0, ]
cs_edges_early <- round(e50$mean_n_edges[e50$generation == min(e50$generation)][1])
cs_edges_late  <- round(e50$mean_n_edges[e50$generation == cs_n50_4_lv][1])

cat("=== Construction survival summary (per_cell) ===\n")
print(as.data.frame(pc), row.names = FALSE)

cat(sprintf("\ncs_common_gen (common pre-fixation horizon) = %s\n", cs_common_gen))
cat(sprintf("N=50, 4:1 ratio: %d of %d replicates still yield a graph at gen 2999 (last valid gen %d)\n",
            cs_n50_4_built, cs_n50_4_int, cs_n50_4_lv))
cat(sprintf("N=50, 2.5:1 ratio: %d of %d (failure count, per qmd formula: n_intended - n_built)\n",
            (cs_n50_25_int - cs_n50_25_built), cs_n50_25_int))
cat(sprintf("Mean edges at N=50, 4:1: %d near onset -> %d just before construction fails\n",
            cs_edges_early, cs_edges_late))

## ---- Fig_ConstructionSurvival ---------------------------------------------

fig_survival <- surv |>
  mutate(Npanel = factor(sprintf("N = %d", N),
                          levels = c("N = 50", "N = 100", "N = 200", "N = 400")),
         Ratio  = factor(sprintf("%.1f:1", ratio))) |>
  ggplot(aes(generation, buildability, colour = Ratio)) +
  geom_vline(xintercept = cs_common_gen, linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Npanel, nrow = 2) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(x = "Generation", y = "Replicates with a constructible graph", colour = "Ratio") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

out_png1 <- "data/derived/fig-construction-survival_reproduced.png"
ggsave(out_png1, fig_survival, width = 7, height = 6, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_ConstructionSurvival to %s\n", out_png1))

## ---- Fig_ConstructionErosion -----------------------------------------------

fig_erosion <- surv |>
  filter(is.finite(mean_n_edges)) |>
  mutate(Npanel = factor(sprintf("N = %d", N),
                          levels = c("N = 50", "N = 100", "N = 200", "N = 400")),
         Ratio  = factor(sprintf("%.1f:1", ratio))) |>
  ggplot(aes(generation, mean_n_edges, colour = Ratio)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Npanel, nrow = 2) +
  labs(x = "Generation", y = "Mean retained edges (surviving graphs)", colour = "Ratio") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

out_png2 <- "data/derived/fig-construction-erosion_reproduced.png"
ggsave(out_png2, fig_erosion, width = 7, height = 6, dpi = 150)
cat(sprintf("Saved reproduced Fig_ConstructionErosion to %s\n", out_png2))
