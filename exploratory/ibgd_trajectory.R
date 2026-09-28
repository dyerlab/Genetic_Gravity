# exploratory/ibgd_trajectory.R
#
# Isolation by graph distance (IBGD) through the forward phase, symmetric (cGD)
# versus directional (pGD) edge weights, per lineage.
#
# Reference: each replicate's own last 100 burn-in generations (1904-1999, every
# 5 gens). Whatever directional lean a lineage carries into generation 2000 is
# its starting condition, so forward-phase values are expressed as departures
# from that lineage's baseline and, for the asymmetric treatments, from the
# same lineage's Isotropic run at the same generation.
#
# Forward phase: every census 2004-2999 (5-gen spacing), Isotropic,
# Redistributed, Obstructed, 50 replicates. Statistics from ibd_one() in
# exploratory/ibd_cgd_vs_pgd.R:
#   r_cgd    r(cGD, |dx|)                  symmetric IBGD
#   r_dir    r(pGD_ij - pGD_ji, x_j - x_i)  directional IBGD (< 0: forward cheaper)
#   rel_fwd  forward / reverse IBGD slope of pGD on |dx|
# Censuses whose graph was not saved (collapse) are counted, not dropped silently.
#
# Phase markers per asymmetric lineage (r_dir smoothed over 5 censuses):
#   onset  first generation the smoothed r_dir falls below the lineage's
#          burn-in mean - 2 SD and stays there for 5 censuses
#   peak   generation and value of the minimum smoothed r_dir
#   return last generation still below the band (NA if it never returns)
#
# Run from the repository root:  Rscript exploratory/ibgd_trajectory.R
# Outputs: data/derived/ibgd_trajectory_fit.rda; media/fig-ibgd-trajectory-fit-v2.png and
#          fig-ibgd-trajectory-fit-v2-vs-cgd.png (earlier: fig-ibgd-trajectory-fit.png)
# (earlier variant without the R^2 columns: ibgd_trajectory.rda, fig-ibgd-trajectory.png)

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(tidyr)
  library(parallel); library(ggplot2) })
source("exploratory/ibd_cgd_vs_pgd.R")          # ibd_one(); main block guarded

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
OUT     <- "data/derived/ibgd_trajectory_fit.rda"
FIG     <- "media/fig-ibgd-trajectory-fit-v2.png"
scen    <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
cols    <- c(Isotropic = "#2a78d6", Redistributed = "#eb6834", Obstructed = "#1baf7a")
stats   <- c("edges", "r_cgd", "r_pgd", "r_dir", "rel_fwd", "r2_cgd", "r2_pgd_blind", "r2_pgd_dir")

jobs <- rbind(
  expand.grid(Replicate = 1:50, Treatment = "Burn-in", generation = seq(1904L, 1999L, 5L),
              stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, Treatment = names(scen), generation = seq(2004L, 2999L, 5L),
              stringsAsFactors = FALSE))
jobs$file <- with(jobs, ifelse(Treatment == "Burn-in",
  sprintf("%s/replicate%d/rep%d-graph-%d.rda", SIM_DIR, Replicate, Replicate, generation),
  sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", SIM_DIR, Replicate, Replicate,
          scen[Treatment], generation)))

if (file.exists(OUT)) load(OUT) else {
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    st <- setNames(rep(NA_real_, length(stats)), stats); status <- "missing"
    if (file.exists(jobs$file[i])) {
      ee <- new.env(); load(jobs$file[i], envir = ee)
      out <- tryCatch(ibd_one(ee$graph), error = function(err) NULL)
      if (is.null(out)) status <- "failed" else { st <- out[stats]; status <- "ok" }
    }
    data.frame(jobs[i, 1:3], t(st), status = status)
  }, mc.cores = max(1L, detectCores() - 1L))
  traj <- do.call(rbind, res)
  dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
  save(traj, file = OUT)
}
traj$Treatment <- factor(traj$Treatment, c("Burn-in", names(scen)))

## ---- lineage baselines and departures -----------------------------------------
base <- traj |> filter(Treatment == "Burn-in", status == "ok") |> group_by(Replicate) |>
  summarise(across(c(r_cgd, r_dir, rel_fwd), list(mu = mean, sd = sd)), .groups = "drop")
iso  <- traj |> filter(Treatment == "Isotropic") |>
  select(Replicate, generation, iso_r_dir = r_dir, iso_r_cgd = r_cgd)
fwd  <- traj |> filter(Treatment != "Burn-in") |> left_join(base, by = "Replicate") |>
  left_join(iso, by = c("Replicate", "generation")) |>
  mutate(d_base = r_dir - r_dir_mu, d_iso = r_dir - iso_r_dir)

cat("Burn-in baselines (last 100 generations), across the 50 lineages:\n")
cat(sprintf("  r_cgd  median %.3f (IQR %.3f-%.3f)\n", median(base$r_cgd_mu),
            quantile(base$r_cgd_mu, .25), quantile(base$r_cgd_mu, .75)))
cat(sprintf("  r_dir  median %.3f (IQR %.3f-%.3f); lineages > 0: %d/50; within-lineage SD median %.3f\n",
            median(base$r_dir_mu), quantile(base$r_dir_mu, .25), quantile(base$r_dir_mu, .75),
            sum(base$r_dir_mu > 0), median(base$r_dir_sd)))
cat(sprintf("  rel_fwd median %.3f\n\n", median(base$rel_fwd_mu)))

fwd$block <- cut(fwd$generation, seq(2000, 3000, 100), labels = paste0(seq(2000, 2900, 100), "s"))
cat("Forward phase, medians by 100-generation block (d_base: change from own burn-in;\n",
    "d_iso: difference from own Isotropic run; missing: graph not built):\n")
print(fwd |> group_by(Treatment, block) |>
        summarise(missing = mean(status != "ok"), edges = median(edges, na.rm = TRUE),
                  r_cgd = median(r_cgd, na.rm = TRUE), r_dir = median(r_dir, na.rm = TRUE),
                  d_base = median(d_base, na.rm = TRUE), d_iso = median(d_iso, na.rm = TRUE),
                  rel_fwd = median(rel_fwd, na.rm = TRUE), .groups = "drop") |> as.data.frame(),
      digits = 3)

## ---- phase markers per lineage ---------------------------------------------------
roll <- function(x, k = 5) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
run_below <- function(b, k = 5) { r <- rle(b); ends <- cumsum(r$lengths)
  i <- which(r$values & r$lengths >= k); if (length(i)) ends[i[1]] - r$lengths[i[1]] + 1 else NA }
phase <- fwd |> filter(Treatment != "Isotropic") |> arrange(Replicate, Treatment, generation) |>
  group_by(Treatment, Replicate) |>
  summarise({
    s <- roll(r_dir); thr <- r_dir_mu[1] - 2 * r_dir_sd[1]; b <- !is.na(s) & s < thr
    on <- run_below(b); pk <- which.min(s)
    data.frame(onset = if (is.na(on)) NA else generation[on],
               peak_gen = if (length(pk)) generation[pk] else NA,
               peak = if (length(pk)) s[pk] else NA,
               last_below = if (any(b)) max(generation[b]) else NA,
               below_at_end = tail(b[!is.na(s)], 1))
  }, .groups = "drop")
cat("\nPhase markers (smoothed r_dir vs own burn-in mean - 2 SD):\n")
print(phase |> group_by(Treatment) |>
        summarise(lineages = n(), with_onset = sum(!is.na(onset)),
                  onset_med = median(onset, na.rm = TRUE), onset_iqr = paste(quantile(onset, c(.25, .75), na.rm = TRUE), collapse = "-"),
                  peak_gen_med = median(peak_gen), peak_gen_iqr = paste(quantile(peak_gen, c(.25, .75)), collapse = "-"),
                  peak_med = round(median(peak), 3), below_at_end = sum(below_at_end)) |> as.data.frame())
iso_on <- fwd |> filter(Treatment == "Isotropic") |> arrange(Replicate, generation) |> group_by(Replicate) |>
  summarise(on = !is.na(run_below({ s <- roll(r_dir); !is.na(s) & s < r_dir_mu[1] - 2 * r_dir_sd[1] })))
cat(sprintf("Isotropic lineages meeting the same onset rule (false onsets): %d/50\n", sum(iso_on$on)))
save(traj, base, phase, file = OUT)

## ---- figures -----------------------------------------------------------------------
## Top: standard IBGD fit (R^2, cGD ~ |dx|). Bottom, two variants:
##   FIG      gain from modelling direction within pGD (direction-aware minus
##            direction-blind adjusted R^2);
##   FIG_CGD  direction-aware pGD adjusted R^2 minus the cGD R^2.
## Isotropic is drawn last (on top); the legend keeps the treatment order.
FIG_CGD <- sub("\\.png$", "-vs-cgd.png", FIG)
draw_order <- c("Redistributed", "Obstructed", "Isotropic")
ibgd_fig <- function(bottom, lab, file) {
  d <- traj |> filter(status == "ok") |>
    mutate(bottom = bottom(r2_pgd_dir, r2_pgd_blind, r2_cgd), Treatment = as.character(Treatment))
  d <- bind_rows(bind_rows(lapply(names(scen), function(s) d |> filter(Treatment == "Burn-in") |>
                                    mutate(Treatment = s))),
                 d |> filter(Treatment != "Burn-in")) |>
    select(Treatment, Replicate, generation, r2_cgd, bottom) |>
    pivot_longer(c(r2_cgd, bottom), names_to = "metric") |>
    mutate(metric = factor(metric, c("r2_cgd", "bottom"), lab), Treatment = factor(Treatment, draw_order)) |>
    group_by(Treatment, metric, generation) |>
    summarise(med = median(value), lo = quantile(value, .25), hi = quantile(value, .75), .groups = "drop")
  p <- ggplot(d, aes(generation, med, colour = Treatment, fill = Treatment)) +
    geom_vline(xintercept = 2000, colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = cols, breaks = names(scen)) +
    scale_fill_manual(values = cols, breaks = names(scen)) +
    facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left") +
    labs(x = "Generation", y = NULL, colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", strip.placement = "outside", panel.grid.minor = element_blank())
  ggsave(file, p, width = 8, height = 6, dpi = 150)
  cat("\nFigure:", file, "\n")
}
ibgd_fig(function(dir, blind, cgd) dir - blind,
         c("IBGD (cGD)", "\u0394R\u00b2 (directional pGD vs pGD)"), FIG)
ibgd_fig(function(dir, blind, cgd) dir - cgd,
         c("IBGD (cGD)", "\u0394R\u00b2 (pGD vs cGD)"), FIG_CGD)

cat("\nFit by 100-generation block (medians):\n")
print(fwd |> group_by(Treatment, block) |>
        summarise(dir_beats_cgd = mean(r2_pgd_dir > r2_cgd, na.rm = TRUE),   # before r2_cgd is summarised
                  gain = median(r2_pgd_dir - r2_pgd_blind, na.rm = TRUE),
                  r2_cgd = median(r2_cgd, na.rm = TRUE), r2_blind = median(r2_pgd_blind, na.rm = TRUE),
                  r2_dir = median(r2_pgd_dir, na.rm = TRUE), .groups = "drop") |>
        as.data.frame(), digits = 3)
