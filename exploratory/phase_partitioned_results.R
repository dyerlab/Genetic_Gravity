# exploratory/phase_partitioned_results.R
#
# Every calendar-window result in the Results, re-expressed by each lineage's
# own phase (exploratory/lineage_phases.R): Reorganize / Plateau / Fall apart
# for Redistributed and Obstructed; Isotropic reported over its whole forward
# phase ("Symmetric"). Rates are pooled over censuses within a phase, with the
# number of lineages contributing; lineage-level summaries average each
# lineage's censuses in the phase first.
#
# Sections
#   desc     C_D, diameter, mean |Delta|, Delta-bar          (data/forward_scenarios.csv, every 5 gens)
#   ibgd     Tab_IBGDDetection and Tab_IBGDFitComparison     (ibgd_trajectory_fit.rda, every 5 gens)
#   sstest   Tab_SourceSinkPower                             (rean$sstest, every 50 gens)
#   edgecls  Tab_EdgeClassDecomposition                      (data/boundary_summary_edges.csv, every 50 gens)
#   truth    Nm retention, truth rho, corridor direction     (rean$truth, rean$corridor, every 5 gens)
#   bw       Tab_BandwidthSensitivity recomputed by phase    (saved graphs, every 50 gens; new computation)
#
# Run from the repository root:  Rscript exploratory/phase_partitioned_results.R
# Output: data/derived/phase_partitioned_results.rda (list `ppr`)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(parallel); library(igraph); library(gstudio) })
source("exploratory/gravity_convention.R")
source("exploratory/lineage_phases.R")
OUT <- "data/derived/phase_partitioned_results.rda"
SC  <- c("Isotropic", "Redistributed", "Obstructed")
ph  <- function(s, r, g) lp_phase(s, r, g)
e <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e); rean <- e$rean
ef <- new.env(); load("data/derived/ibgd_trajectory_fit.rda", envir = ef)
ppr <- if (file.exists(OUT)) { z <- new.env(); load(OUT, envir = z); z$ppr } else list()
set.seed(20260930)

lin_n <- function(d) n_distinct(d$Replicate)

## ---- descriptive metrics ----------------------------------------------------------------
fwd <- gc_read_forward() |> filter(Scenario %in% SC) |> rename(Replicate = replicate) |>
  mutate(phase = ph(Scenario, Replicate, generation))
ppr$desc <- fwd |> group_by(Scenario, phase, Replicate) |>
  summarise(across(c(CD, diameter, mean_abs_delta, dbar), mean), n = n(), .groups = "drop") |>
  group_by(Scenario, phase) |>
  summarise(lineages = n(), censuses = sum(n), CD = median(CD), diameter = median(diameter),
            abs_delta = median(mean_abs_delta), dbar_pos = mean(dbar > 0), dbar = median(dbar), .groups = "drop")

## ---- IBGD ------------------------------------------------------------------------------------
ib <- ef$traj |> filter(status == "ok", Treatment %in% SC) |> mutate(Treatment = as.character(Treatment)) |>
  mutate(phase = ph(Treatment, Replicate, generation), det = r2_pgd_dir > r2_cgd)
ppr$ibgd <- ib |> group_by(Treatment, phase) |>
  summarise(lineages = n_distinct(Replicate), censuses = n(), detect = mean(det), med_dR2 = median(r2_pgd_dir - r2_cgd),
            ratio_lt1 = mean(rel_fwd[det] < 1), med_ratio_det = median(rel_fwd[det]),
            r2_cgd = median(r2_cgd), r2_blind = median(r2_pgd_blind), r2_dir = median(r2_pgd_dir),
            gain = median(r2_pgd_dir - r2_pgd_blind), .groups = "drop")

## ---- source-sink test --------------------------------------------------------------------
ss <- rean$sstest$tests |> filter(null == "N1w", scenario %in% SC) |> mutate(phase = ph(scenario, Replicate, generation))
ppr$sstest <- ss |> group_by(scenario, phase) |>
  summarise(lineages = n_distinct(Replicate), censuses = n(), reject = mean(p < .05), med_r = median(r_S),
            frac_up = mean(r_S < 0), .groups = "drop") |>
  left_join(ss |> group_by(scenario, phase, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
              group_by(scenario, phase) |> summarise(lin_up = mean(m < 0), .groups = "drop"), by = c("scenario", "phase"))

## ---- edge classes ----------------------------------------------------------------------------
ec <- gc_read_stored("data/boundary_summary_edges.csv") |> rename(Replicate = replicate) |>
  mutate(phase = ph(scenario, Replicate, generation))
ppr$edgecls <- ec |> group_by(scenario, phase, edge_class) |>
  summarise(censuses = n(), lineages = n_distinct(Replicate), dbar = mean(dbar), abs = mean(mean_abs),
            edges = mean(n), .groups = "drop")

## ---- divMigrate comparison ----------------------------------------------------------------------
tr <- rean$truth$d |> mutate(phase = ph(Treatment, Replicate, Generation))
co <- rean$corridor$d |> mutate(phase = ph(Treatment, Replicate, Generation))
ppr$truth <- tr |> group_by(Treatment, phase) |>
  summarise(censuses = n(), lineages = n_distinct(Replicate), nm_defined = mean(Exclusion == ""),
            nm_defined_lineages = n_distinct(Replicate[Exclusion == ""]),
            w_rho = median(w_rho_tr), w_rho_paired = median(w_rho_tr[Exclusion == ""]), nm_rho = median(nm_rho_run, na.rm = TRUE),
            w_beats_nm = mean((w_rho_tr > nm_rho_run)[Exclusion == ""]),
            w_rmse = median(w_rmse_tr), nm_rmse = median(nm_rmse_run, na.rm = TRUE), .groups = "drop") |>
  left_join(co |> group_by(Treatment, phase) |>
              summarise(w_fwd = mean(w_fwd_gt, na.rm = TRUE), nm_fwd = mean(nm_fwd_gt, na.rm = TRUE),
                        w_fwd_paired = mean(w_fwd_gt[nm_ok], na.rm = TRUE), .groups = "drop"), by = c("Treatment", "phase"))

## ---- bandwidth sensitivity by phase (new computation) ------------------------------------
if (is.null(ppr$bw)) {
  jobs <- expand.grid(Replicate = 1:50, Scenario = c("Redistributed", "Obstructed"), generation = seq(2004L, 2954L, 50L),
                      stringsAsFactors = FALSE)
  sc <- c(Redistributed = 2L, Obstructed = 3L)
  alts <- list(global = "global", `local x1/2` = list(scale = 0.5), `local x2` = list(scale = 2), perplexity = list(perplexity = 4))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    g <- gc_load_graph(jobs$Replicate[i], jobs$generation[i], sc[[jobs$Scenario[i]]]); if (is.null(g)) return(NULL)
    can <- gc_edges(gc_weights(g), g)
    out <- list()
    for (nm in names(alts)) {
      a <- alts[[nm]]
      W <- tryCatch(if (identical(a, "global")) gc_weights(g, bandwidth = mean(E(g)$weight)) else do.call(gc_weights, c(list(g), a)),
                    error = function(e) NULL)
      if (is.null(W)) next
      alt <- gc_edges(W, g); d1 <- alt$delta; d0 <- can$delta; nz <- d0 != 0 | d1 != 0
      out[[nm]] <- data.frame(jobs[i, ], bandwidth = nm, concord = mean(sign(d1[nz]) == sign(d0[nz])),
                              rho_signed = suppressWarnings(cor(d1, d0, method = "spearman")),
                              rho_abs = suppressWarnings(cor(abs(d1), abs(d0), method = "spearman")),
                              dbar_kept = sign(mean(d1)) == sign(mean(d0)))
    }
    do.call(rbind, out)
  }, mc.cores = max(1L, detectCores() - 1L), mc.preschedule = FALSE)
  bw <- do.call(rbind, res) |> mutate(phase = ph(Scenario, Replicate, generation))
  ppr$bw_raw <- bw
  ppr$bw <- bw |> group_by(bandwidth, phase) |>
    summarise(censuses = n(), concord = mean(concord), rho_signed = mean(rho_signed, na.rm = TRUE),
              rho_abs = mean(rho_abs, na.rm = TRUE), dbar_kept = mean(dbar_kept), .groups = "drop")
  ppr$bw_wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
}
ppr$bw_raw$phase <- ph(ppr$bw_raw$Scenario, ppr$bw_raw$Replicate, ppr$bw_raw$generation)   # current phases
ppr$bw <- ppr$bw_raw |> group_by(bandwidth, phase) |>
  summarise(censuses = n(), concord = mean(concord), rho_signed = mean(rho_signed, na.rm = TRUE),
            rho_abs = mean(rho_abs, na.rm = TRUE), dbar_kept = mean(dbar_kept), .groups = "drop")
ppr$meta <- list(seed = 20260930, date = as.character(Sys.Date()), phases = "exploratory/lineage_phases.R")
save(ppr, file = OUT)

if (sys.nframe() == 0L) {
  options(width = 180)
  for (nm in c("desc", "ibgd", "sstest", "edgecls", "truth", "bw")) {
    cat("\n=====", nm, "=====\n"); print(as.data.frame(ppr[[nm]]), digits = 3)
  }
  cat("\nBandwidth wall time (min):", round(ppr$bw_wall_min, 1), "\n")
  cat("\nExisting Tab_BandwidthSensitivity is gen 2999 only; phase of those censuses:\n")
  sc2 <- rep(c("Redistributed", "Obstructed"), each = 50)
  print(table(sc2, lp_phase(sc2, rep(1:50, 2), rep(2999, 100))))
}
