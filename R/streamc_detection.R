# Fig_StreamCPower, Tab_DetectionLimitsEndState, Tab_DetectionLimitsPreFixation.
#
# Reproduces the Stream C Arm B detection-limits / power-envelope analysis:
# single-snapshot power (at a fixed 5% false-positive rate, matched 1:1 null
# per Ne) across the migration-ratio x population-size grid, both at the
# final generation (2999, "end state") and at the common pre-fixation
# horizon (generation 2744, the latest generation at which every grid cell
# -- including N=50 -- still yields a constructible graph for every
# replicate). Source data (data/graph_summary.csv) is the
# per-replicate/scenario/generation graph summary already computed from the
# individual-based simulation output by the private research repo's
# R/graph_summarizer.R; no raw simulation data is required to run this
# script. armB_detection() below is ported verbatim from the private repo's
# R/detection_analysis.R (its own .auc() rank-based helper kept, matching
# the same formula reused in R/divmigrate_sensitivity.R) -- the internal
# `as.integer(as.character(scenario))` coercion works unchanged whether
# `scenario` is a factor (private .rda) or a plain integer column (this
# CSV). construction_survival()/survival_summary() are duplicated here
# (verbatim, from R/streamc_construction.R) only to compute the
# pre-fixation-horizon generation for the second table.
#
# Run from the repository root:
#   Rscript R/streamc_detection.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

graph_summary <- read.csv("data/graph_summary.csv", stringsAsFactors = FALSE)

## ---- ported from R/detection_analysis.R (private repo), verbatim ---------

ARM_B_GRID <- tibble::tribble(
  ~N,   ~ratio, ~cond,  ~null_cond,
  50,   1.5,    8,      7,
  50,   2.0,    9,      7,
  50,   2.5,    10,     7,
  50,   4.0,    11,     7,
  100,  1.5,    12,     1,
  100,  2.0,    13,     1,
  100,  2.5,    14,     1,
  100,  4.0,    2,      1,
  200,  1.5,    16,     15,
  200,  2.0,    17,     15,
  200,  2.5,    18,     15,
  200,  4.0,    19,     15,
  400,  1.5,    21,     20,
  400,  2.0,    22,     20,
  400,  2.5,    23,     20,
  400,  4.0,    24,     20
)

.auc <- function(pos, null) {
  pos  <- pos[is.finite(pos)]
  null <- null[is.finite(null)]
  n1 <- length(pos); n0 <- length(null)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(c(pos, null))
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

armB_detection <- function(graph_summary, stab_gen = 2999) {
  gs <- graph_summary
  if (!"N" %in% names(gs)) gs$N <- 100L

  code <- function(x) as.integer(as.character(x))
  gs_fwd_all <- gs |> filter(phase == "forward")
  gs_forward <- gs_fwd_all |> filter(generation == stab_gen)

  present_conds <- unique(code(gs_forward$scenario))
  grid_subset <- ARM_B_GRID |> filter(cond %in% present_conds, null_cond %in% present_conds)

  if (nrow(grid_subset) == 0L) return(list(ok = FALSE, grid = NULL, summary = NULL))

  results <- purrr::map_dfr(seq_len(nrow(grid_subset)), function(i) {
    row <- grid_subset[i, ]

    n_intended <- gs_fwd_all |>
      filter(as.integer(N) == row$N, code(scenario) == row$cond) |>
      distinct(replicate) |> nrow()

    pos_vals <- gs_forward |>
      filter(as.integer(N) == row$N, code(scenario) == row$cond) |>
      pull(mean_delta)
    null_vals <- gs_forward |>
      filter(as.integer(N) == row$N, code(scenario) == row$null_cond) |>
      pull(mean_delta)

    pos_stat  <- abs(pos_vals[is.finite(pos_vals)])
    null_stat <- abs(null_vals[is.finite(null_vals)])
    n_built   <- length(pos_stat)

    thr <- if (length(null_stat) > 0L)
      stats::quantile(null_stat, 0.95, names = FALSE, type = 7) else NA_real_

    auc_val  <- if (n_built > 0L && length(null_stat) > 0L) .auc(pos_stat, null_stat) else NA_real_
    n_detect <- if (!is.na(thr)) sum(pos_stat > thr) else NA_integer_

    power_net  <- if (n_intended > 0L && !is.na(thr)) n_detect / n_intended else NA_real_
    power_cond <- if (n_built   > 0L && !is.na(thr)) n_detect / n_built   else NA_real_

    m_total <- 0.050
    m_rev   <- m_total / (row$ratio + 1.0)
    delta_m <- (row$ratio - 1.0) * m_rev

    tibble::tibble(
      N                 = row$N,
      ratio             = row$ratio,
      four_ne_deltam    = 4.0 * row$N * delta_m,
      auc               = auc_val,
      power             = power_net,
      power_conditional = power_cond,
      buildability      = if (n_intended > 0L) n_built / n_intended else NA_real_,
      n_intended        = n_intended,
      n_built           = n_built,
      n_null            = length(null_stat)
    )
  })

  min_detectable <- results |>
    filter(power >= 0.8) |>
    group_by(N) |>
    slice_min(ratio, n = 1) |>
    select(N, min_ratio = ratio, power, four_ne_deltam) |>
    ungroup()

  list(ok = TRUE, grid = results, summary = min_detectable)
}

## ---- duplicated from R/streamc_construction.R, verbatim (pre-fixation gen) --

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
    mutate(N = as.integer(N), cond = as.integer(as.character(scenario))) |>
    inner_join(cells, by = c("N", "cond"))

  intended <- fwd |>
    group_by(N, ratio, cond) |>
    summarise(n_intended = n_distinct(replicate), .groups = "drop")

  per_gen <- fwd |>
    group_by(N, ratio, cond, generation) |>
    summarise(n_built = n_distinct(replicate), .groups = "drop")

  census_gens <- sort(unique(fwd$generation))
  expand_grid(distinct(cells, N, ratio, cond), generation = census_gens) |>
    left_join(intended, by = c("N", "ratio", "cond")) |>
    left_join(per_gen,  by = c("N", "ratio", "cond", "generation")) |>
    mutate(n_built = coalesce(n_built, 0L),
           buildability = ifelse(n_intended > 0, n_built / n_intended, NA_real_)) |>
    arrange(N, ratio, generation)
}

survival_summary <- function(surv, tol = 1e-9) {
  bygen <- surv |>
    group_by(generation) |>
    summarise(min_build = min(buildability, na.rm = TRUE), .groups = "drop")
  gg <- bygen$generation[bygen$min_build >= 1 - tol]
  list(common_pre_fixation_gen = if (length(gg)) max(gg) else NA_real_)
}

## ---- End-state grid (generation 2999) -------------------------------------

da <- armB_detection(graph_summary, stab_gen = 2999)

cat("=== Tab_DetectionLimitsEndState (min ratio reaching >=80% net power) ===\n")
tbl_endstate <- da$summary |>
  transmute(
    `Population Size (N)` = sprintf("N = %d", N),
    `Min Detectable Ratio` = sprintf("%.1f : 1", min_ratio),
    `Power at Threshold` = sprintf("%.1f%%", power * 100),
    `Composite Scale (4Ne*dm)` = sprintf("%.2f", four_ne_deltam))
print(as.data.frame(tbl_endstate), row.names = FALSE)

g <- da$grid; s <- da$summary
bN50_4 <- g$buildability[g$N == 50 & g$ratio == 4][1]
cat(sprintf(
  "\nN=50, ratio=4:1 buildability at gen 2999 = %.0f%% (caps attainable power)\n", 100 * bN50_4))
cat(sprintf(
  "Among N reaching 80%% threshold: min_ratio range = %.1f-%.1f : 1\n",
  min(s$min_ratio), max(s$min_ratio)))
cat(sprintf(
  "4Ne*dm at those thresholds: range = %.1f-%.1f (%.1f-fold)\n",
  min(s$four_ne_deltam), max(s$four_ne_deltam),
  max(s$four_ne_deltam) / min(s$four_ne_deltam)))

## ---- Fig_StreamCPower (heatmap at end state) ------------------------------

fig_power <- da$grid |>
  mutate(N_factor = factor(sprintf("N = %d", N),
                            levels = c("N = 400", "N = 200", "N = 100", "N = 50"))) |>
  ggplot(aes(x = factor(ratio), y = N_factor, fill = power)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent, option = "plasma") +
  geom_text(aes(label = sprintf("%.1f%%", power * 100)),
            colour = ifelse(da$grid$power > 0.5, "black", "white"), fontface = "bold") +
  labs(x = "Migration Asymmetry Ratio (forward : reverse)", y = "Population Size (N)",
       fill = "Power (5% FPR)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "right")

out_png <- "data/derived/fig-streamc-power_reproduced.png"
ggsave(out_png, fig_power, width = 9, height = 3.6, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_StreamCPower to %s\n", out_png))

## ---- Pre-fixation-horizon grid --------------------------------------------

gen_p <- survival_summary(construction_survival(graph_summary))$common_pre_fixation_gen
da_p  <- armB_detection(graph_summary, stab_gen = gen_p)

cat(sprintf("\n=== Tab_DetectionLimitsPreFixation (generation %s) ===\n", gen_p))
tbl_prefix <- da_p$summary |>
  transmute(
    `Population Size (N)` = sprintf("N = %d", N),
    `Min Detectable Ratio` = sprintf("%.1f : 1", min_ratio),
    `Power at Threshold` = sprintf("%.1f%%", power * 100),
    `Composite Scale (4Ne*dm)` = sprintf("%.2f", four_ne_deltam))
print(as.data.frame(tbl_prefix), row.names = FALSE)

sp <- da_p$summary
n50_p <- { v <- sp$min_ratio[sp$N == 50]; if (length(v)) v[1] else NA_real_ }
cat(sprintf(
  "\nEvaluated at generation %d (%d generations after onset). N=50 min detectable ratio at pre-fixation horizon = %.1f : 1 (end-state N=50 does not reach 80%% threshold at any tested ratio)\n",
  as.integer(gen_p), as.integer(gen_p - 2000), n50_p))
