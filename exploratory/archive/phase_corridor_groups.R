# exploratory/phase_corridor_groups.R
#
# EXPLORATORY: not (yet) behind any manuscript figure or table.
#
# How genetic gravity (pGD similarity) and divMigrate (Nm) score the true
# stepping-stone corridors versus the non-migrating pairs as each lineage moves
# through three phases after onset: reorganize -> accumulate -> fall apart.
#
# Phases, per replicate x scenario (the "hybrid, diameter-based" rule):
#   1 reorganize   generations 2004-2349 (onset is shared by all lineages);
#   3 fall apart   from the lineage's own onset: the first generation after its
#                  smoothed-diameter maximum at which smoothed diameter is 20%
#                  below that maximum (penalized spline, k = 10), never before
#                  2350; lineages where this never happens have no phase 3;
#   2 accumulate   everything in between.
# The rule fires in every redistributed lineage, but in only about half of the
# obstructed and isotropic ones, i.e. no more often than under symmetric
# migration; isotropic is the control for what the rule does on noise.
#
# Pair groups (off-diagonal only; defined by the true matrix, not by either
# method): corridor forward (i -> i+1), corridor reverse (i+1 -> i), and
# non-corridor pairs at chain distance 2, 3-5, and 6+.
# Scores: gravity = 1 / (1 + p_ij) on retained edges, 0 elsewhere; divMigrate =
# relative Nm; each rescaled by its own off-diagonal maximum, as is the truth.
# divMigrate is NA at censuses where its Nm matrix has any non-finite or
# negative entry. Gravity is reported on ALL censuses and, for comparison, on
# the censuses where Nm survives: the latter is heavily survivor-biased in the
# fall-apart phase (Nm has broken down in ~82% of redistributed fall-apart
# censuses) and hides the loss of spatial contrast visible on all censuses.
#
# Inputs: data/forward_scenarios.csv (this repo) and the saved per-census
# matrices in the private research repo's data/divMigrate/
# (divMig.{rep}.matrices.rda from R/divmigrate_matrix_census.R;
#  divMig.{rep}.pij.rda from R/divmigrate_truth_compare.R).
#
# Run from the repository root:
#   Rscript exploratory/phase_corridor_groups.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(mgcv)
})

DM_DIR <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
OUT    <- "data/derived/phase_corridor_groups.rds"
REORG_END <- 2350    # first generation of phase 2
FALL_FRAC <- 0.8     # fall apart once smoothed diameter <= 80% of its maximum
CORES <- max(1L, parallel::detectCores() - 1L)

scn   <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
rates <- list(iso = c(0.025, 0.025), flux = c(0.040, 0.010), rate = c(0.025, 0.010))

## ---- Per-lineage phase markers --------------------------------------------

fwd <- read.csv("data/forward_scenarios.csv", stringsAsFactors = FALSE)
markers <- fwd |>
  group_by(Scenario, replicate) |>
  group_modify(function(x, key) {
    pd   <- as.numeric(predict(gam(diameter ~ s(generation, k = 10), data = x, method = "REML")))
    imax <- which.max(pd)
    hit  <- which(seq_along(pd) > imax & pd <= FALL_FRAC * pd[imax])
    tibble(diam_peak_gen = x$generation[imax],
           fall_gen = if (length(hit)) max(x$generation[hit[1]], REORG_END) else NA_integer_)
  }) |>
  ungroup()

cat("=== Fall-apart onset (diameter rule) ===\n")
print(as.data.frame(markers |> group_by(Scenario) |>
  summarise(lineages_falling = sum(!is.na(fall_gen)), median_onset = median(fall_gen, na.rm = TRUE),
            q10 = quantile(fall_gen, 0.1, na.rm = TRUE), q90 = quantile(fall_gen, 0.9, na.rm = TRUE),
            .groups = "drop")), row.names = FALSE)

## ---- Pair-group scores per census ------------------------------------------

K <- 25; i <- row(diag(K)); j <- col(diag(K)); off <- i != j; d <- abs(i - j)
grp <- ifelse(j == i + 1, "corridor fwd", ifelse(i == j + 1, "corridor rev",
       ifelse(d == 2, "non-corr d=2", ifelse(d <= 5, "non-corr d=3-5", "non-corr d=6+"))))[off]

one_rep <- function(r) {
  e <- new.env(parent = emptyenv())
  load(file.path(DM_DIR, sprintf("divMig.%d.matrices.rda", r)), envir = e)
  load(file.path(DM_DIR, sprintf("divMig.%d.pij.rda", r)), envir = e)
  bind_rows(lapply(names(scn), function(s) {
    Tm <- matrix(0, K, K); Tm[j == i + 1] <- rates[[s]][1]; Tm[i == j + 1] <- rates[[s]][2]
    Tr <- (Tm / max(Tm))[off]
    fg   <- markers$fall_gen[markers$Scenario == scn[[s]] & markers$replicate == r]
    gens <- as.integer(dimnames(e$divmig_M[[s]])[[3]])
    bind_rows(lapply(seq_along(gens), function(k) {
      M <- e$divmig_M[[s]][, , k]; P <- e$divmig_P[[s]][, , k]
      broken <- any(!is.finite(M[off])) || any(M[off] < 0)
      G <- ifelse(P > 0, 1 / (1 + P), 0); G <- (G / max(G[off]))[off]
      Mv <- if (broken) rep(NA_real_, sum(off)) else (M / max(M[off]))[off]
      phase <- if (gens[k] < REORG_END) "1 reorganize"
               else if (length(fg) && !is.na(fg) && gens[k] >= fg) "3 fall apart" else "2 accumulate"
      data.frame(group = grp, gravity = G, divMigrate = Mv, truth = Tr) |>
        group_by(group) |>
        summarise(gravity = mean(gravity), divMigrate = mean(divMigrate), truth = mean(truth),
                  .groups = "drop") |>
        mutate(replicate = r, scenario = scn[[s]], generation = gens[k], phase = phase,
               nm_broken = broken)
    }))
  }))
}
res <- bind_rows(parallel::mclapply(1:50, one_rep, mc.cores = CORES))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(markers = markers, census = res), OUT)
cat(sprintf("\nSaved per-census pair-group scores to %s\n", OUT))

## ---- Summaries (lineages weighted equally) -----------------------------------

lin <- res |>
  group_by(scenario, phase, group, replicate) |>
  summarise(gravity_all = mean(gravity), gravity_nm_ok = mean(gravity[!nm_broken]),
            divMigrate = mean(divMigrate, na.rm = TRUE), truth = mean(truth),
            nm_broken = mean(nm_broken), censuses = n(), .groups = "drop") |>
  mutate(across(c(gravity_nm_ok, divMigrate), ~ ifelse(is.nan(.x), NA_real_, .x)))

cat("\n=== Share of censuses where divMigrate has broken down ===\n")
print(as.data.frame(lin |> filter(group == "corridor fwd") |>
  group_by(scenario, phase) |>
  summarise(lineages = n(), nm_broken = round(sum(nm_broken * censuses) / sum(censuses), 2),
            censuses = sum(censuses), .groups = "drop")),
  row.names = FALSE)

tab <- lin |>
  group_by(scenario, phase, group) |>
  summarise(truth = round(median(truth), 2), gravity_all = round(median(gravity_all), 3),
            gravity_nm_ok = round(median(gravity_nm_ok, na.rm = TRUE), 3),
            divMigrate_nm_ok = round(median(divMigrate, na.rm = TRUE), 3), .groups = "drop")
for (sc in scn) {
  cat(sprintf("\n=== %s: median over lineages of mean rescaled score ===\n", sc))
  print(as.data.frame(tab |> filter(scenario == sc) |> select(-scenario)), row.names = FALSE)
}

ratio <- function(v, col) {
  lin |>
    select(scenario, phase, group, replicate, value = all_of(col)) |>
    pivot_wider(names_from = group, values_from = value) |>
    mutate(fwd_rev = `corridor fwd` / `corridor rev`,
           corr_vs_d2 = `corridor fwd` / `non-corr d=2`,
           corr_vs_d6 = `corridor fwd` / `non-corr d=6+`) |>
    group_by(scenario, phase) |>
    summarise(fwd_rev = round(median(fwd_rev, na.rm = TRUE), 3),
              corr_vs_d2 = round(median(corr_vs_d2, na.rm = TRUE), 2),
              corr_vs_d6 = round(median(corr_vs_d6[is.finite(corr_vs_d6)], na.rm = TRUE), 1),
              .groups = "drop") |>
    mutate(score = v) |>
    relocate(score)
}
cat("\n=== Contrasts: forward/reverse corridor (direction) and corridor/non-corridor (spatial) ===\n")
print(as.data.frame(bind_rows(ratio("gravity, all censuses", "gravity_all"),
                              ratio("divMigrate, Nm-surviving censuses", "divMigrate")) |>
  arrange(score, scenario, phase)), row.names = FALSE)

## ---- Corridor direction ratios: p_ij : p_ji and Nm_ij : Nm_ji ----------------
#
# For each nearest-neighbour corridor (i, i+1) at each census, the forward /
# reverse ratio (i -> i+1) / (i+1 -> i) of raw p_ij (which equals the directional
# weight ratio w_{i->j} / w_{j->i}, e_ij cancelling) and of relative Nm. True
# ratios: 4 (redistributed), 2.5 (obstructed), 1 (isotropic). p ratios need the
# corridor to be a retained graph edge; Nm ratios need the census's Nm matrix
# intact (no non-finite or negative entry). Averaged on the log scale
# (geometric mean) within each lineage x phase, then across lineages.

ratio_rep <- function(r) {
  e <- new.env(parent = emptyenv())
  load(file.path(DM_DIR, sprintf("divMig.%d.matrices.rda", r)), envir = e)
  load(file.path(DM_DIR, sprintf("divMig.%d.pij.rda", r)), envir = e)
  fw <- cbind(1:(K - 1), 2:K); rv <- fw[, 2:1]
  bind_rows(lapply(names(scn), function(s) {
    fg   <- markers$fall_gen[markers$Scenario == scn[[s]] & markers$replicate == r]
    gens <- as.integer(dimnames(e$divmig_M[[s]])[[3]])
    bind_rows(lapply(seq_along(gens), function(k) {
      M <- e$divmig_M[[s]][, , k]; P <- e$divmig_P[[s]][, , k]
      broken <- any(!is.finite(M[off])) || any(M[off] < 0)
      phase <- if (gens[k] < REORG_END) "1 reorganize"
               else if (length(fg) && !is.na(fg) && gens[k] >= fg) "3 fall apart" else "2 accumulate"
      data.frame(replicate = r, scenario = scn[[s]], phase = phase, corridor = 1:(K - 1),
                 log_p  = ifelse(P[fw] > 0 & P[rv] > 0, log(P[fw] / P[rv]), NA_real_),
                 log_nm = if (broken) NA_real_ else
                          ifelse(M[fw] > 0 & M[rv] > 0, log(M[fw] / M[rv]), NA_real_))
    }))
  }))
}
rr <- bind_rows(parallel::mclapply(1:50, ratio_rep, mc.cores = CORES))

rr_lin <- rr |>
  group_by(scenario, phase, replicate) |>
  summarise(p_gm = exp(mean(log_p, na.rm = TRUE)), p_gt1 = mean(log_p > 0, na.rm = TRUE),
            p_n = sum(!is.na(log_p)),
            nm_gm = exp(mean(log_nm, na.rm = TRUE)), nm_gt1 = mean(log_nm > 0, na.rm = TRUE),
            nm_n = sum(!is.na(log_nm)), .groups = "drop") |>
  mutate(across(c(nm_gm, nm_gt1), ~ ifelse(is.nan(.x), NA_real_, .x)))

cat("\n=== Corridor forward/reverse ratios (geometric mean across lineages) ===\n")
print(as.data.frame(rr_lin |>
  group_by(scenario, phase) |>
  summarise(truth = c(Isotropic = 1, Redistributed = 4, Obstructed = 2.5)[first(scenario)],
            `p_ij:p_ji` = round(exp(mean(log(p_gm), na.rm = TRUE)), 3),
            `p frac >1` = round(mean(p_gt1, na.rm = TRUE), 3),
            `p pairs` = sum(p_n),
            `Nm_ij:Nm_ji` = round(exp(mean(log(nm_gm), na.rm = TRUE)), 3),
            `Nm frac >1` = round(mean(nm_gt1, na.rm = TRUE), 3),
            `Nm pairs` = sum(nm_n),
            .groups = "drop")), row.names = FALSE)
